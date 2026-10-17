# Map reviewed seeds to SIDER without collapsing flat groups into exact entities.
library(readr);library(dplyr);library(tidyr);library(stringr);library(tibble);library(httr2)
if (!all(dir.exists(c("scripts", "config", "data"))))stop('Run from project root.')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
tsv<-function(f,names)read_tsv(file.path('data/raw/sider',f),col_names=names,
 col_types=cols(.default=col_character()),quote='',na=c('','NA'),progress=FALSE)
x<-r('data/curated/seed_identity_crosscheck_reviewed.csv')
aliases<-r('data/curated/seed_aliases_reviewed.csv')
se<-tsv('meddra_all_se.tsv.gz',c('stitch_flat_id','stitch_stereo_id','label_umls_cui','meddra_type','meddra_umls_cui','side_effect_name'))
names<-tsv('drug_names.tsv',c('stitch_flat_id','sider_drug_name'))
freq<-tsv('meddra_freq.tsv.gz',c('stitch_flat_id','stitch_stereo_id','label_umls_cui','placebo_raw',
 'frequency_description','frequency_lower_raw','frequency_upper_raw','meddra_type','meddra_umls_cui','side_effect_name'))
ind<-tsv('meddra_all_indications.tsv.gz',c('stitch_flat_id','label_umls_cui','detection_method',
 'concept_name','meddra_type','meddra_umls_cui','meddra_concept_name'))
stopifnot(nrow(problems(se))==0,nrow(problems(freq))==0,nrow(problems(ind))==0,
 all(str_detect(se$stitch_stereo_id,'^CID0[0-9]+$')),!anyDuplicated(names$stitch_flat_id))
pairs<-se |> distinct(stitch_flat_id,stitch_stereo_id) |>
 mutate(sider_pubchem_cid=as.character(as.numeric(str_remove(stitch_stereo_id,'^CID')))) |>
 left_join(names,by='stitch_flat_id')
# Decode stereo CID only. Flat identifiers are grouped SIDER entities, never exact CID joins.
cid_leads<-bind_rows(
 x |> transmute(drugbank_id,cid=reviewed_pubchem_cid,lead_method='reviewed_CID'),
 x |> transmute(drugbank_id,cid=pubchem_cid,lead_method='original_source_CID_lead')) |>
 filter(!is.na(cid)) |> inner_join(pairs,by=c('cid'='sider_pubchem_cid'),relationship='many-to-many') |>
 mutate(sider_pubchem_cid=cid) |> select(-cid)
name_leads<-aliases |> transmute(drugbank_id=seed_drugbank_id,
 name_key=str_to_lower(str_squish(alias))) |> distinct() |>
 inner_join(names |> mutate(name_key=str_to_lower(str_squish(sider_drug_name))),
 by='name_key',relationship='many-to-many') |>
 select(drugbank_id,stitch_flat_id) |> distinct() |>
 inner_join(pairs,by='stitch_flat_id',relationship='many-to-many') |> mutate(lead_method='exact_normalised_alias_lead')
candidates<-bind_rows(cid_leads,name_leads) |>
 group_by(drugbank_id,stitch_flat_id,stitch_stereo_id,sider_pubchem_cid,sider_drug_name) |>
 summarise(lead_methods=paste(sort(unique(lead_method)),collapse='|'),.groups='drop')
write_csv(candidates,'outputs/11b_sider_mapping_candidates.csv',na='')
message('Candidate seed/SIDER pairs: ',nrow(candidates),' for ',n_distinct(candidates$drugbank_id),' seed drugs')

# Current PubChem structures independently check the historical SIDER CID leads.
pc<-r('data/curated/pubchem_seed_identity_records.csv') |>
 transmute(cid,pc_inchikey,pc_inchi,pc_formula,pc_title)
needed<-sort(setdiff(unique(candidates$sider_pubchem_cid),pc$cid))
cache_dir<-'data/raw/sider/identity_api';dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
log<-list();scalar<-function(z)if(is.null(z)||!length(z))NA_character_ else as.character(z[[1]])
fetch<-function(ids,label){
 url<-paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',paste(ids,collapse=','),
 '/property/InChIKey,InChI,MolecularFormula,Title/JSON')
 path<-file.path(cache_dir,paste0(label,'.rds'))
 if(file.exists(path)){z<-readRDS(path);if(!identical(z$url,url))stop('Cache URL changed: ',label)}else{
 Sys.sleep(.25);resp<-request(url)|>req_timeout(40)|>req_retry(max_tries=3,retry_on_failure=TRUE)|>
 req_error(is_error=function(r)FALSE)|>req_perform()
 status<-resp_status(resp);if(!status %in% c(200L,404L))stop('PubChem failed: ',status)
 z<-list(url=url,status=status,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),body=resp_body_json(resp));saveRDS(z,path)}
 log[[length(log)+1L]]<<-tibble(url=z$url,status=z$status,retrieved_utc=z$retrieved_utc,cache_file=path)
 z
}
parse<-function(z){
 if(z$status!=200L)return(tibble(cid=character(),pc_inchikey=character(),pc_inchi=character(),pc_formula=character(),pc_title=character()))
 bind_rows(lapply(z$body$PropertyTable$Properties,function(p)tibble(cid=scalar(p$CID),
 pc_inchikey=scalar(p$InChIKey),pc_inchi=scalar(p$InChI),pc_formula=scalar(p$MolecularFormula),pc_title=scalar(p$Title))))}
groups<-split(needed,ceiling(seq_along(needed)/20))
for(i in seq_along(groups)){
 z<-fetch(groups[[i]],paste0('batch_',i))
 if(z$status==404L){for(id in groups[[i]])pc<-bind_rows(pc,parse(fetch(id,paste0('cid_',id))))}
 else pc<-bind_rows(pc,parse(z))
 message('SIDER structure batch ',i,' / ',length(groups),' checked')
}
pc<-pc |> distinct(cid,.keep_all=TRUE)
mapping<-candidates |> left_join(pc,by=c('sider_pubchem_cid'='cid')) |>
 left_join(x |> select(drugbank_id,drug_name,inchikey,drug_type,seed_tiers),by='drugbank_id') |>
 mutate(mapping_status=case_when(
 is.na(inchikey)~'review_name_lead_no_seed_structure',
 is.na(pc_inchikey)~'review_CID_structure_not_returned',
 inchikey==pc_inchikey~'accepted_full_InChIKey_match',
 substr(inchikey,1,14)==substr(pc_inchikey,1,14)~'review_connectivity_only_DO_NOT_pool',
 TRUE~'review_different_structure_or_representation_DO_NOT_pool'),
 evidence_note='SIDER 4.1 label associations; flat-group origin and historical identifiers retained; full key alone is not coordination-geometry proof')
accepted<-mapping |> filter(mapping_status=='accepted_full_InChIKey_match')
review<-mapping |> filter(mapping_status!='accepted_full_InChIKey_match')
# PT only for the main binary feature table; raw LLT and frequencies remain separate.
pt<-se |> filter(meddra_type=='PT') |>
 distinct(stitch_flat_id,stitch_stereo_id,meddra_umls_cui,side_effect_name)
features<-accepted |> select(drugbank_id,drug_name,stitch_flat_id,stitch_stereo_id,mapping_status) |>
 inner_join(pt,by=c('stitch_flat_id','stitch_stereo_id'),relationship='many-to-many')
binary<-features |> distinct(drugbank_id,meddra_umls_cui) |> mutate(present=1L)
# No zero-filled matrix: lack of SIDER coverage is not an observed negative.
counts<-binary |> count(drugbank_id,name='unique_PT_side_effects')
candidate_counts<-mapping |> count(drugbank_id,name='candidate_SIDER_records')
accepted_counts<-accepted |> count(drugbank_id,name='accepted_SIDER_records')
coverage<-x |> select(drugbank_id,drug_name,drug_type,seed_tiers,cancers) |>
 left_join(candidate_counts,by='drugbank_id') |> left_join(accepted_counts,by='drugbank_id') |>
 left_join(counts,by='drugbank_id') |>
 mutate(candidate_SIDER_records=replace_na(candidate_SIDER_records,0L),
 accepted_SIDER_records=replace_na(accepted_SIDER_records,0L),
 coverage_status=case_when(!is.na(unique_PT_side_effects)~'mapped_with_PT_side_effect_data',
 candidate_SIDER_records>0~'candidate_requires_identity_review',TRUE~'not_mapped_in_this_release'))
freq_seed<-accepted |> select(drugbank_id,drug_name,stitch_flat_id,stitch_stereo_id) |>
 inner_join(freq,by=c('stitch_flat_id','stitch_stereo_id'),relationship='many-to-many') |>
 mutate(placebo=is.na(placebo_raw)==FALSE & placebo_raw=='placebo',
 frequency_lower=suppressWarnings(as.numeric(frequency_lower_raw)),
 frequency_upper=suppressWarnings(as.numeric(frequency_upper_raw)))
# Flat-only indication/ATC records are contextual leads, not exact-structure evidence.
ind_seed<-accepted |> distinct(drugbank_id,drug_name,stitch_flat_id) |>
 inner_join(ind,by='stitch_flat_id',relationship='many-to-many') |>
 mutate(mapping_scope='flat_group_context_only_not_seed_indication_validation')
primary_ids<-r('data/curated/primary_compound_mapping_candidates.csv')$drugbank_id
summary<-bind_rows(coverage |> count(coverage_status,name='drugs') |> mutate(scope='all_170_seeds'),
 coverage |> filter(drugbank_id %in% primary_ids) |> count(coverage_status,name='drugs') |> mutate(scope='primary_55_structure_queue'))
stopifnot(nrow(coverage)==170,!anyDuplicated(coverage$drugbank_id),
 all(accepted$inchikey==accepted$pc_inchikey),!anyDuplicated(binary[c('drugbank_id','meddra_umls_cui')]),
 all(!is.na(binary$meddra_umls_cui)),all(binary$present==1L))
write_csv(pairs,'data/curated/sider_compound_index.csv',na='')
write_csv(mapping,'data/curated/seed_sider_mapping_all.csv',na='')
write_csv(accepted,'data/curated/seed_sider_mapping_accepted.csv',na='')
write_csv(review,'outputs/11b_sider_identity_review_queue.csv',na='')
write_csv(pc,'data/curated/sider_candidate_PubChem_structures.csv',na='')
write_csv(features,'data/curated/seed_sider_PT_associations.csv',na='')
write_csv(binary,'data/curated/seed_sider_PT_binary_long.csv',na='')
write_csv(freq_seed,'data/curated/seed_sider_frequencies.csv',na='')
write_csv(ind_seed,'data/curated/seed_sider_indication_context.csv',na='')
write_csv(coverage,'outputs/11b_sider_seed_coverage.csv',na='')
write_csv(summary,'outputs/11b_sider_coverage_summary.csv',na='')
if(length(log))write_csv(bind_rows(log),'outputs/11b_sider_identity_API_provenance.csv')
report<-c('SCRIPT 11b - SIDER MAPPING AND COVERAGE',capture.output(print(summary)),
 paste('Accepted seed/SIDER pairs:',nrow(accepted)),paste('Distinct seed/PT associations:',nrow(binary)),
 paste('Candidate pairs retained for review:',nrow(review)),
 'Inputs: reviewed 10d seeds/aliases, original SIDER tables and cached PubChem structures.',
 'Outputs: data/curated/seed_sider_*; sider_compound_index.csv; candidate structures.',
 'Coverage and identity-review queue: outputs/11b_*. Primary queue is 55, not all 75 primary drugs.',
 'Only full-key stereo-CID mappings enter the current main side-effect feature table.',
 'Other names/structures and biologics remain review candidates, not negative results.',
 'PT CUIs are deduplicated per seed. LLT rows do not become additional independent events.',
 'Frequency rows retain placebo status, raw descriptions and ranges; no pooling across trials.',
 'SIDER is a 2015 label-derived snapshot. Absence is missing coverage, not evidence of safety.',
 'Flat-group indication mentions (including NLP preconditions/text mentions) do not alter seed approval decisions.',
 'No similarity scores, p-values, BH correction or eEF1A binding claims are generated here.',
 'Next: review identity exceptions before computing side-effect similarity on mapped seeds.')
writeLines(report,'outputs/11b_sider_mapping_summary.txt')
cat(paste0('\nSCRIPT 11b RUN - ',format(Sys.time(),tz='Europe/Rome',usetz=TRUE),'\n',
 'Script: scripts/11b_map_sider_seeds.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
