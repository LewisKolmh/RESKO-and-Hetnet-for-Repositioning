# Complete human STITCH 5 detailed-file extraction for checked target protein IDs.
# Association/channel scores are not measured binding or statistical p-values.
library(readr);library(dplyr);library(tidyr);library(stringr);library(tibble);library(httr2)
if(!file.exists('config/eef1a_STITCH_protein_crosswalk.csv'))stop('Run from project root after 17b.')
source('scripts/helpers/stream_stitch_targets.R')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
root<-'data/raw/stitch';dir.create(root,recursive=TRUE,showWarnings=FALSE)
url<-'http://stitch.embl.de/download/protein_chemical.links.detailed.v5.0/9606.protein_chemical.links.detailed.v5.0.tsv.gz'
path<-file.path(root,'9606.protein_chemical.links.detailed.v5.0.tsv.gz');meta_path<-paste0(path,'.metadata.rds')
if(!file.exists(path)){
 resp<-request(url) |> req_timeout(300) |> req_retry(max_tries=2,retry_on_failure=TRUE) |> req_perform(path=paste0(path,'.partial'))
 stopifnot(resp_status(resp)==200L)
 con<-file(paste0(path,'.partial'),'rb');magic<-readBin(con,'raw',2L);close(con);stopifnot(identical(magic,as.raw(c(31,139))))
 stopifnot(file.rename(paste0(path,'.partial'),path))
 saveRDS(list(url=url,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),last_modified=resp_header(resp,'last-modified'),
 file_bytes=file.info(path)$size,md5=unname(tools::md5sum(path))),meta_path)
}
meta<-readRDS(meta_path);stopifnot(identical(meta$url,url),identical(meta$md5,unname(tools::md5sum(path))),file.info(path)$size==meta$file_bytes)
# Cache the publisher's namespace and score definitions, not an inferred CID conversion.
readme_path<-file.path(root,'README.v5.rds')
if(file.exists(readme_path))readme<-readRDS(readme_path)else{
 ru<-'http://stitch.embl.de/download/README';resp<-request(ru) |> req_timeout(30) |> req_perform()
 readme<-list(url=ru,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),body=resp_body_string(resp));saveRDS(readme,readme_path)
}
stopifnot(str_detect(readme$body,fixed('CIDs / CID0')),str_detect(readme$body,fixed('CIDm / CID1')))
page_path<-file.path(root,'download_page.v5.rds')
if(!file.exists(page_path)){
 pu<-'http://stitch.embl.de/cgi/download.pl';resp<-request(pu) |> req_timeout(30) |> req_retry(max_tries=3,retry_on_failure=TRUE,
 is_transient=function(r)resp_status(r) %in% c(408L,429L,500L,502L,503L,504L,522L)) |> req_perform()
 saveRDS(list(url=pu,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),body=resp_body_string(resp)),page_path)
}
con<-file(path,'rb');magic<-readBin(con,'raw',2L);close(con);stopifnot(identical(magic,as.raw(c(31,139))))
inputs<-c('config/eef1a_STITCH_protein_crosswalk.csv','outputs/16f_target_evidence_handoff_summary.csv',
 'data/curated/seed_identity_crosscheck_reviewed.csv','data/curated/candidate_identity_registry.csv',
 'data/curated/seed_sider_mapping_accepted.csv',path,meta_path,readme_path,page_path,'scripts/helpers/stream_stitch_targets.R')
before<-tools::md5sum(inputs);cw<-r(inputs[1]);inventory<-r(inputs[2]);seed<-r(inputs[3]);candidate<-r(inputs[4]);sider<-r(inputs[5])
stopifnot(nrow(cw)==2L,all(cw$crosswalk_supported=='TRUE'),all(grepl('^9606[.]ENSP[0-9]+$',cw$stitch_protein_id)))
sig<-paste(meta$md5,paste(sort(cw$stitch_protein_id),collapse='|'),sep='|');cache<-file.path(root,'eef1a_filtered_v5.rds')
if(file.exists(cache)){cached<-readRDS(cache);stopifnot(identical(cached$signature,sig));z<-cached$data}else{
 z<-stream_stitch_targets(path,cw$stitch_protein_id);saveRDS(list(signature=sig,data=z),cache)
}
edges<-read_tsv(I(paste(c(z$header,z$lines),collapse='\n')),col_types=cols(chemical=col_character(),protein=col_character(),.default=col_integer()))
score_cols<-c('experimental','prediction','database','textmining','combined_score')
stopifnot(all(edges$protein %in% cw$stitch_protein_id),!anyDuplicated(edges[c('chemical','protein')]),
 all(grepl('^CID[ms][0-9]+$',edges$chemical)),all(vapply(edges[score_cols],function(v)all(!is.na(v)&v>=0L&v<=1000L),logical(1))))
edges<-edges |> left_join(cw |> select(gene,uniprot,protein=stitch_protein_id,canonical_STRING_crossreference_match,
 current_UniProt_Ensembl_protein_match,sequence_variant_equivalence_established),by='protein') |>
 mutate(chemical_namespace=if_else(startsWith(chemical,'CIDs'),'stereo_specific','merged_salt_isomer_group'),
  textmining_only=experimental==0L & prediction==0L & database==0L & textmining>0L,
  experimental_channel_present=experimental>0L,database_channel_present=database>0L,
  evidence_status='STITCH_association_not_confirmed_direct_binding',source_release='5.0',
  source_url=url,retrieved_utc=meta$retrieved_utc)
# Only stereo-specific documented suffixes become PubChem-ID leads, supported by existing key crosschecks.
make_stereo<-function(cid){stopifnot(all(grepl('^[0-9]+$',cid)));paste0('CIDs',sprintf('%08.0f',as.numeric(cid)))}
map<-list()
for(i in seq_len(nrow(seed))){
 a<-seed[i,];if(is.na(a$inchikey)||is.na(a$pc_inchikey)||a$inchikey!=a$pc_inchikey)next
 ids<-unique(na.omit(c(a$resolved_pubchem_cid,a$reviewed_pubchem_cid)));ids<-ids[nzchar(ids)]
 if(length(ids))map[[length(map)+1L]]<-tibble(entity_id=paste0('seed:',a$drugbank_id),chemical=make_stereo(ids),
 mapping_relation='stereo_ID_with_existing_current_full_key_support_historical_structure_not_revalidated',mapping_priority=3L,
 mapping_source='reviewed_seed_PubChem_key_crosscheck')
}
for(i in seq_len(nrow(candidate))){
 a<-candidate[i,];if(a$automatic_label_mapping_blocked=='TRUE'||is.na(a$returned_cid))next
 verified<-a$identity_status=='current_PubChem_DrugBank_full_key_agreement'
 map[[length(map)+1L]]<-tibble(entity_id=paste0('sider:',a$stitch_flat_id),chemical=make_stereo(a$returned_cid),
 mapping_relation=if(verified)'stereo_ID_with_existing_current_full_key_support_historical_structure_not_revalidated'else 'stereo_ID_current_PubChem_record_only_identity_review',
 mapping_priority=if(verified)3L else 1L,mapping_source='candidate_current_PubChem_record')
}
# Exact published namespace translation links GROUP CONTEXT only. No flat numeric ID becomes a molecular structure.
flat<-bind_rows(sider |> transmute(entity_id=paste0('seed:',drugbank_id),stitch_flat_id),
 candidate |> filter(automatic_label_mapping_blocked!='TRUE') |> transmute(entity_id=paste0('sider:',stitch_flat_id),stitch_flat_id)) |>
 distinct() |> mutate(chemical=sub('^CID1','CIDm',stitch_flat_id),
 mapping_relation='merged_group_cross_release_context_only_NOT_structure_equivalence',mapping_priority=0L,
 mapping_source='SIDER_group_ID_publisher_documented_namespace') |>
 select(entity_id,chemical,mapping_relation,mapping_priority,mapping_source)
map<-bind_rows(bind_rows(map),flat) |> distinct() |> arrange(desc(mapping_priority)) |>
 distinct(entity_id,chemical,.keep_all=TRUE)
stopifnot(!any(map$entity_id=='sider:CID100135113'),all(map$entity_id %in% inventory$entity_id))
links<-map |> inner_join(edges,by='chemical',relationship='many-to-many') |>
 left_join(inventory |> select(entity_id,origin,source_name,current_record_name,identity_status,structure_mapping_blocked),by='entity_id') |>
 mutate(identity_link_requires_historical_review=TRUE,direct_binding_established=NA,
 isoform_selectivity_established=NA,RT_PPI_disruption_established=NA,allostery_established=NA)
stopifnot(!anyDuplicated(links[c('entity_id','chemical','protein')]),!any(links$structure_mapping_blocked=='TRUE'))
source_coverage<-edges |> group_by(gene,uniprot,chemical_namespace) |>
 summarise(association_rows=n(),distinct_chemical_IDs=n_distinct(chemical),
 experimental_channel_rows=sum(experimental_channel_present),database_channel_rows=sum(database_channel_present),
 prediction_channel_rows=sum(prediction>0),textmining_only_rows=sum(textmining_only),.groups='drop')
counts<-links |> group_by(entity_id,gene) |> summarise(associated_chemical_IDs=n_distinct(chemical),
 stereo_IDs_with_current_key_support=sum(chemical_namespace=='stereo_specific' & mapping_priority==3L),
 stereo_IDs_needing_identity_review=sum(chemical_namespace=='stereo_specific' & mapping_priority<3L),
 merged_group_context_IDs=sum(chemical_namespace=='merged_salt_isomer_group'),
 experimental_channel_IDs=sum(experimental_channel_present),database_channel_IDs=sum(database_channel_present),
 textmining_only_IDs=sum(textmining_only),max_combined_association_score=max(combined_score),.groups='drop')
coverage<-crossing(entity_id=inventory$entity_id,gene=cw$gene) |>
 left_join(inventory |> select(entity_id,origin,source_name,current_record_name,structure_mapping_blocked),by='entity_id') |>
 left_join(map |> count(entity_id,name='usable_source_ID_leads'),by='entity_id') |> left_join(counts,by=c('entity_id','gene')) |>
 mutate(across(c(usable_source_ID_leads,associated_chemical_IDs,stereo_IDs_with_current_key_support,stereo_IDs_needing_identity_review,
 merged_group_context_IDs,experimental_channel_IDs,database_channel_IDs,textmining_only_IDs),~replace_na(.x,0L)),
 source_coverage_status=case_when(structure_mapping_blocked=='TRUE'~'quarantined_structure_named_biologic_context_only',
 usable_source_ID_leads==0L~'no_usable_source_identifier_crosswalk',associated_chemical_IDs==0L~'no_matching_ID_in_this_release_NOT_negative_evidence',
 TRUE~'source_association_for_review_NOT_binding_or_selectivity'))
stopifnot(nrow(coverage)==2L*nrow(inventory),!anyDuplicated(coverage[c('entity_id','gene')]),identical(before,tools::md5sum(inputs)))
write_csv(edges,'data/curated/eef1a_STITCH_association_channels.csv',na='')
write_csv(map,'data/curated/STITCH_inventory_identifier_leads.csv',na='')
write_csv(links,'data/curated/eef1a_STITCH_inventory_association_evidence.csv',na='')
write_csv(source_coverage,'outputs/17c_STITCH_source_coverage.csv',na='')
write_csv(links |> group_by(entity_id,origin,source_name,current_record_name) |>
 summarise(target_genes=paste(sort(unique(gene)),collapse='|'),
 stereo_ID_with_current_key_support_any=any(mapping_priority==3L),
 experimental_channel_present_any=any(experimental_channel_present),database_channel_present_any=any(database_channel_present),
 all_linked_rows_textmining_only=all(textmining_only),max_combined_association_score=max(combined_score),.groups='drop'),
 'outputs/17c_linked_inventory_summary.csv',na='')
write_csv(coverage,'outputs/17c_STITCH_inventory_coverage.csv',na='')
write_csv(links |> count(origin,gene,mapping_relation,experimental_channel_present,database_channel_present,textmining_only,name='source_ID_associations'),
 'outputs/17c_inventory_evidence_summary.csv',na='')
write_csv(tibble(file=inputs,md5=unname(before)),'outputs/17c_input_summary.csv')
write_csv(tibble(source_id='stitch',release='5.0',url=url,transport='HTTP_public_download',file_bytes=meta$file_bytes,
 retrieved_utc=meta$retrieved_utc,last_modified=meta$last_modified,md5=meta$md5,rows_scanned=z$rows_scanned,
 local_file=path),'data/raw/stitch/acquisition_manifest.csv')
manifest<-r('config/database_manifest.csv');k<-manifest$source_id=='stitch'
manifest$status[k]<-'human_v5_detailed_acquired_eef1a_subset_extracted'
manifest$release_version[k]<-'5.0 (legacy/unsupported; source file modified 2016-08-16)'
manifest$retrieved_date[k]<-format(Sys.Date());manifest$download_url[k]<-url;manifest$local_file[k]<-'data/raw/stitch/acquisition_manifest.csv'
manifest$licence_note[k]<-'Academic research; publisher states commercial users require a licence for detailed scores. Raw files and source-rich rows remain local; preserve STITCH attribution.'
write_csv(manifest,'config/database_manifest.csv',na='')
report<-c('SCRIPT 17c - HUMAN STITCH DETAILED ASSOCIATION ACQUISITION',
 paste('Human source rows scanned:',z$rows_scanned),paste('Exact target-protein association rows:',nrow(edges)),
 paste('Distinct STITCH chemical IDs across both targets:',n_distinct(edges$chemical)),
 paste('Inventory entities with source-ID association links:',n_distinct(links$entity_id)),
 paste('Inventory entities with stereo-ID links and existing current full-key support:',n_distinct(links$entity_id[links$mapping_priority==3L])),
 paste('Inventory source-ID/protein association rows:',nrow(links)),
 paste('Target rows supported only by textmining channel:',sum(edges$textmining_only)),
 paste('Target rows with experimental-channel support:',sum(edges$experimental_channel_present)),
 paste('Inventory-linked rows with experimental-channel support:',sum(links$experimental_channel_present)),
 paste('Inventory-linked drug names:',paste(sort(unique(links$current_record_name)),collapse='; ')),
 'Complete human detailed source file scanned; no API neighbour cap or additional score cutoff was applied.',
 'Source is legacy STITCH 5.0, not a current release; file modification date and retrieval date are distinct.',
 'Experimental, prediction, database, textmining and combined scores stay separate in their original 0..1000 scale.',
 'STITCH scores are association confidence, not affinity, p-values, direct binding or inhibition/selectivity measurements.',
 'The detailed file does not distinguish direct versus transferred evidence; experimental-channel support alone cannot establish direct human isoform binding.',
 'Stereo-ID inventory links use documented PubChem suffixes plus existing identity checks; original STITCH 5 InChIKeys were not independently revalidated.',
 'Merged salt/isomer IDs remain cross-release group-context leads only; merged/stereo rows are not independent biological evidence.',
 'No sacrosidase wrong-structure link accepted. Missing records or identity coverage is not negative activity evidence.',
 'Public HTTP transport is explicit; HTTPS certificate checks are not bypassed. Raw download and source-rich scores remain local.',
 'Outputs locally: eef1a_STITCH_association_channels.csv, STITCH_inventory_identifier_leads.csv, eef1a_STITCH_inventory_association_evidence.csv.',
 'Compact source/inventory coverage, evidence-layer counts and hashes in outputs/17c_*; source acquisition manifest is under data/raw/stitch/.',
 'Next: primary-source review of experimental/database-supported association leads, resolve historical chemical forms and provenance overlap with ChEMBL/DrugBank before binding prioritisation.')
writeLines(report,'outputs/17c_STITCH_acquisition_summary.txt')
writeLines(c('PASS: source gzip/MD5 and exact target protein fields; unique chemical/protein IDs, score ranges and source schema.',
 'PASS: unique inventory chemical/protein links, documented namespace mapping, quarantined structure excluded and all 400 entity/target coverage rows retained.',
 'Scope: ID links are not proof of historical structure equivalence or direct binding; source channel provenance remains explicit.'),
 'outputs/17c_STITCH_validation.txt')
cat(paste0('\nSCRIPT 17c RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/optional/eef1a_evidence/17c_acquire_stitch_associations.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
