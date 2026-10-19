# Follow stereo CID leads; preserve historical groups and related structures separately.
library(readr);library(dplyr);library(stringr);library(tibble);library(xml2)
if(!file.exists('outputs/13_identity_review_queue_summary.csv'))stop('Run from project root after script 13.')
source('scripts/helpers/cached_public_get.R')
source('scripts/helpers/stream_drugbank_fields.R')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
q<-r('outputs/13_identity_review_queue_summary.csv')
index<-r('data/curated/sider_compound_index.csv')
leads<-q |> inner_join(index,by='stitch_flat_id')
stopifnot(n_distinct(leads$stitch_flat_id)==nrow(q),!anyDuplicated(q$stitch_flat_id),
 all(grepl('^CID0[0-9]+$',leads$stitch_stereo_id)))
cids<-sort(unique(leads$sider_pubchem_cid));properties<-list();synonyms<-list();logs<-list()
scalar<-function(z)if(length(z))as.character(z[[1]])else NA_character_
for(i in seq_along(cids)){
 cid<-cids[i];base<-paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',cid)
 z<-public_cache_get(paste0(base,'/property/InChIKey,InChI,MolecularFormula,MolecularWeight,IUPACName,Title/JSON'),
  paste0('data/raw/candidate_identity/pubchem_properties_',cid,'.rds'))
 logs[[length(logs)+1L]]<-tibble(url=z$url,status=z$status,retrieved_utc=z$retrieved_utc)
 if(z$status==200L)for(p in z$body$PropertyTable$Properties){
 properties[[length(properties)+1L]]<-tibble(requested_cid=cid,returned_cid=scalar(p$CID),
  current_inchikey=scalar(p$InChIKey),current_inchi=scalar(p$InChI),molecular_formula=scalar(p$MolecularFormula),
  molecular_weight=scalar(p$MolecularWeight),pubchem_title=scalar(p$Title),iupac_name=scalar(p$IUPACName))}
 z<-public_cache_get(paste0(base,'/synonyms/JSON'),paste0('data/raw/candidate_identity/pubchem_synonyms_',cid,'.rds'))
 logs[[length(logs)+1L]]<-tibble(url=z$url,status=z$status,retrieved_utc=z$retrieved_utc)
 if(z$status==200L)for(p in z$body$InformationList$Information){
 v<-as.character(unlist(p$Synonym,use.names=FALSE))
 synonyms[[length(synonyms)+1L]]<-tibble(requested_cid=rep(cid,length(v)),returned_cid=rep(scalar(p$CID),length(v)),alias=v)}
 message('Candidate CID checked: ',i,' / ',length(cids))
}
pc<-bind_rows(properties);syn<-bind_rows(synonyms)
# No CID returned by the API is silently substituted for a historical requested CID.
identity<-leads |> left_join(pc,by=c('sider_pubchem_cid'='requested_cid'),relationship='many-to-many')
local_index<-r('data/curated/drugbank_structure_identity_index.csv')
links<-identity |> filter(!is.na(current_inchikey)) |>
 inner_join(local_index |> filter(!is.na(record_inchikey)),by=c('current_inchikey'='record_inchikey'),relationship='many-to-many') |>
 distinct(stitch_flat_id,sider_pubchem_cid,returned_cid,current_inchikey,record_drugbank_id,record_name)
source_manifest<-r('config/database_manifest.csv') |> filter(source_id=='drugbank')
xml_path<-source_manifest$local_file[[1]]
ids<-sort(unique(links$record_drugbank_id));cache<-'data/raw/candidate_identity/drugbank_record_enrichment.rds'
signature<-paste(normalizePath(xml_path),file.info(xml_path)$size,as.numeric(file.info(xml_path)$mtime),paste(ids,collapse='|'),sep='|')
if(file.exists(cache)){
 local<-readRDS(cache);if(!identical(local$signature,signature))stop('Local enrichment cache source/IDs changed.')
}else{
 records<-list();aliases<-list();xrefs<-list()
 process<-function(node){
  get<-function(p)xml_text(xml_find_first(node,p))
  id<-get('./drugbank-id[@primary="true"]');if(!id %in% ids)return(invisible(NULL))
  records[[length(records)+1L]]<<-tibble(record_drugbank_id=id,record_name=get('./name'),
   unii=get('./unii'),cas_number=get('./cas-number'),drug_type=xml_attr(node,'type'),
   description=get('./description'),indication=get('./indication'),mechanism_of_action=get('./mechanism-of-action'))
  v<-unique(c(get('./name'),xml_text(xml_find_all(node,'./synonyms/synonym'))))
  aliases[[length(aliases)+1L]]<<-tibble(record_drugbank_id=id,alias=v)
  e<-xml_find_all(node,'./external-identifiers/external-identifier')
  xrefs[[length(xrefs)+1L]]<<-tibble(record_drugbank_id=rep(id,length(e)),
   resource=xml_text(xml_find_first(e,'./resource')),identifier=xml_text(xml_find_first(e,'./identifier')))
 }
 stream_drugbank_fields(xml_path,c('drugbank-id','name','unii','cas-number','description','indication','mechanism-of-action','synonyms','external-identifiers'),process)
 local<-list(signature=signature,records=bind_rows(records),aliases=bind_rows(aliases),xrefs=bind_rows(xrefs));saveRDS(local,cache)
}
stopifnot(all(ids %in% local$records$record_drugbank_id))
local_links<-links |> left_join(local$records,by=c('record_drugbank_id','record_name'))
counts<-links |> group_by(stitch_flat_id) |> summarise(exact_key_DrugBank_records=n_distinct(record_drugbank_id),
 matched_names=paste(sort(unique(record_name)),collapse='|'),.groups='drop')
seed<-r('data/curated/seed_identity_crosscheck_reviewed.csv')
identity<-identity |> left_join(counts,by='stitch_flat_id') |>
 mutate(identity_status=case_when(is.na(current_inchikey)~'current_CID_without_structure_key_review',
 !is.na(exact_key_DrugBank_records)~'current_PubChem_DrugBank_full_key_agreement',TRUE~'current_PubChem_record_only_review'),
  already_seed_full_key=!is.na(current_inchikey)&current_inchikey %in% na.omit(seed$inchikey),
  stereo_or_charge_key_note='Full standard key agreement is representation-level; preserve stereochemistry, charge, metal/isotope and formulation context')
aliases<-bind_rows(identity |> select(stitch_flat_id,sider_pubchem_cid) |>
 inner_join(syn,by=c('sider_pubchem_cid'='requested_cid'),relationship='many-to-many') |>
 transmute(stitch_flat_id,alias,source='PubChem',source_record=returned_cid,alias_scope='record_label_not_individual_equivalence_proof'),
 links |> inner_join(local$aliases,by='record_drugbank_id',relationship='many-to-many') |>
 transmute(stitch_flat_id,alias,source='DrugBank supplied XML',source_record=record_drugbank_id,
 alias_scope='full_key_matched_record_label_not_individual_equivalence_proof')) |> distinct()
label_support<-identity |> select(stitch_flat_id,sider_drug_name) |>
 left_join(aliases |> mutate(label_norm=str_to_lower(str_squish(alias))) |>
 select(stitch_flat_id,label_norm) |> distinct(),by='stitch_flat_id',relationship='many-to-many') |>
 group_by(stitch_flat_id) |> summarise(SIDER_name_exact_alias=any(str_to_lower(str_squish(sider_drug_name))==label_norm,na.rm=TRUE),.groups='drop')
exceptions<-r('config/candidate_identity_exceptions.csv')
stopifnot(!anyDuplicated(exceptions$stitch_flat_id),all(exceptions$stitch_flat_id %in% identity$stitch_flat_id))
identity<-identity |> left_join(label_support,by='stitch_flat_id') |>
 left_join(exceptions,by='stitch_flat_id') |>
 mutate(historical_identity_status=case_when(!is.na(identity_exception)~'conflicting_historical_label_and_current_record',
 SIDER_name_exact_alias~'alias_supported_historical_formulation_still_unverified',TRUE~'historical_label_to_current_record_review_required'),
 automatic_label_mapping_blocked=!is.na(identity_exception))
identifiers<-links |> inner_join(local$xrefs,by='record_drugbank_id',relationship='many-to-many')
# Group shared current full keys for navigation, without pooling SIDER profiles or related salts.
structure_groups<-identity |> filter(!is.na(current_inchikey)) |> group_by(current_inchikey) |>
 summarise(SIDER_groups=n_distinct(stitch_flat_id),group_ids=paste(sort(unique(stitch_flat_id)),collapse='|'),.groups='drop')
write_csv(identity,'data/curated/candidate_identity_registry.csv',na='')
write_csv(local_links,'data/curated/candidate_DrugBank_record_enrichment.csv',na='')
write_csv(identifiers,'data/curated/candidate_external_identifiers.csv',na='')
write_csv(aliases,'data/curated/candidate_aliases.csv',na='')
write_csv(identity |> select(stitch_flat_id,sider_drug_name,stitch_stereo_id,sider_pubchem_cid,returned_cid,
 pubchem_title,current_inchikey,molecular_formula,identity_status,exact_key_DrugBank_records,matched_names,
 SIDER_name_exact_alias,already_seed_full_key,historical_identity_status,automatic_label_mapping_blocked,identity_exception,overlap_context), 'outputs/14a_candidate_identity_summary.csv',na='')
write_csv(structure_groups,'outputs/14a_structure_group_summary.csv',na='')
write_csv(bind_rows(logs),'data/curated/candidate_identity_API_provenance.csv')
stopifnot(setequal(identity$stitch_flat_id,q$stitch_flat_id),!anyDuplicated(links[c('stitch_flat_id','record_drugbank_id')]))
report<-c('SCRIPT 14a - CANDIDATE IDENTITY CROSS-CHECK',paste('SIDER groups retained:',n_distinct(identity$stitch_flat_id)),
 paste('Stereo CID leads requested:',length(cids)),paste('Groups with current structure keys:',n_distinct(identity$stitch_flat_id[!is.na(identity$current_inchikey)])),
 paste('Groups linked to DrugBank by full key:',n_distinct(links$stitch_flat_id)),
 paste('DrugBank source records enriched:',length(ids)),paste('Provenance-bearing alias rows:',nrow(aliases)),
 'Outputs: candidate_identity_registry.csv; candidate_DrugBank_record_enrichment.csv; candidate_aliases.csv;',
 'candidate_external_identifiers.csv and API provenance (data/curated/); compact identity/structure summaries (outputs/14a_*).',
 'Historical SIDER IDs and requested/returned CIDs remain separate. No flat CID arithmetic or name-only merging.',
 'Full-key matches link current representations, not historical label/formulation equivalence or independent experimental identity.',
 'Matched source records retain names, UNIIs, CAS numbers, descriptions, indications, mechanisms and cross-references locally.',
 'Alias labels remain search leads. Same-key groups are navigational only; original SIDER side-effect profiles are not pooled.',
 paste('Groups with exact historical-name alias support:',sum(identity$SIDER_name_exact_alias)),
 paste('Groups already represented by a seed full key:',sum(identity$already_seed_full_key)),
 paste('Conflicting historical/current identities blocked from automatic label mapping:',sum(identity$automatic_label_mapping_blocked)),
 'Sacrosidase CID lead 135113 resolves to 2-Methyl-5-tert-butyl-1,3-benzodioxole (C12H16O2); preserve as a conflict, not an enzyme structure.',
 'Current-record structure resolution is not confirmation of the historical SIDER drug identity.',
 'Next: verify current label ingredients/active moieties and preserve scope before constructing comparable PT profiles.')
writeLines(report,'outputs/14a_identity_resolution_summary.txt')
cat(paste0('\nSCRIPT 14a RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/14a_resolve_candidate_identities.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
