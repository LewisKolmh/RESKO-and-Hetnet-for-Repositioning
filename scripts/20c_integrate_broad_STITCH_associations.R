# Core RESKO: full chemical-neighbour extraction, legacy protein aliases and human pathways.
library(readr);library(dplyr);library(tidyr);library(tibble);library(httr2)
source('scripts/helpers/stream_stitch_chemicals.R')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
inputs<-c('outputs/16f_target_evidence_handoff_summary.csv',
 'data/curated/seed_identity_crosscheck_reviewed.csv','data/curated/candidate_identity_registry.csv',
 'data/curated/seed_sider_mapping_accepted.csv',
 'data/curated/Reactome_human_lowest_level_annotations.csv',
 'data/curated/Reactome_human_all_level_annotations.csv')
before<-tools::md5sum(inputs)
inventory<-r(inputs[1]);seed<-r(inputs[2]);candidate<-r(inputs[3]);sider<-r(inputs[4])
root<-'data/raw/stitch';dir.create(root,recursive=TRUE,showWarnings=FALSE)
acquire<-function(url,path){
 metadata<-paste0(path,'.metadata.rds')
 if(!file.exists(path)){
  z<-request(url)|>req_timeout(300)|>req_retry(max_tries=3,retry_on_failure=TRUE)|>req_perform(path=paste0(path,'.partial'))
  stopifnot(resp_status(z)==200L,file.rename(paste0(path,'.partial'),path))
  saveRDS(list(url=url,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
   last_modified=resp_header(z,'last-modified'),file_bytes=file.info(path)$size,
   md5=unname(tools::md5sum(path))),metadata)
 }
 stopifnot(file.exists(metadata));m<-readRDS(metadata)
 stopifnot(identical(m$url,url),identical(m$md5,unname(tools::md5sum(path))))
 con<-file(path,'rb');magic<-readBin(con,'raw',n=2L);close(con)
 stopifnot(identical(magic,as.raw(c(31,139))));m
}
source_url<-'http://stitch.embl.de/download/protein_chemical.links.detailed.v5.0/9606.protein_chemical.links.detailed.v5.0.tsv.gz'
source_path<-file.path(root,'9606.protein_chemical.links.detailed.v5.0.tsv.gz')
source_meta<-acquire(source_url,source_path)
alias_url<-'https://version-10.string-db.org/download/protein.aliases.v10/9606.protein.aliases.v10.txt.gz'
alias_path<-file.path(root,'9606.protein.aliases.v10.txt.gz');alias_meta<-acquire(alias_url,alias_path)
# Independently derive compound leads from earlier core identity tables, not optional target studies.
make_stereo<-function(cid){stopifnot(all(grepl('^[0-9]+$',cid)));paste0('CIDs',sprintf('%08.0f',as.numeric(cid)))}
maps<-list()
for(i in seq_len(nrow(seed))){
 a<-seed[i,];if(is.na(a$inchikey)||is.na(a$pc_inchikey)||a$inchikey!=a$pc_inchikey)next
 ids<-unique(na.omit(c(a$resolved_pubchem_cid,a$reviewed_pubchem_cid)));ids<-ids[nzchar(ids)]
 if(length(ids))maps[[length(maps)+1L]]<-tibble(entity_id=paste0('seed:',a$drugbank_id),chemical=make_stereo(ids),
 mapping_relation='stereo_ID_with_existing_current_full_key_support_historical_structure_not_revalidated',mapping_priority=3L,
 mapping_source='reviewed_seed_PubChem_key_crosscheck')
}
for(i in seq_len(nrow(candidate))){
 a<-candidate[i,];if(a$automatic_label_mapping_blocked=='TRUE'||is.na(a$returned_cid))next
 verified<-a$identity_status=='current_PubChem_DrugBank_full_key_agreement'
 maps[[length(maps)+1L]]<-tibble(entity_id=paste0('sider:',a$stitch_flat_id),chemical=make_stereo(a$returned_cid),
 mapping_relation=if(verified)'stereo_ID_with_existing_current_full_key_support_historical_structure_not_revalidated'else 'stereo_ID_current_PubChem_record_only_identity_review',
 mapping_priority=if(verified)3L else 1L,mapping_source='candidate_current_PubChem_record')
}
flat<-bind_rows(sider |> transmute(entity_id=paste0('seed:',drugbank_id),stitch_flat_id),
 candidate |> filter(automatic_label_mapping_blocked!='TRUE') |> transmute(entity_id=paste0('sider:',stitch_flat_id),stitch_flat_id)) |>
 distinct() |> mutate(chemical=sub('^CID1','CIDm',stitch_flat_id),
 mapping_relation='merged_group_cross_release_context_only_NOT_structure_equivalence',mapping_priority=0L,
 mapping_source='SIDER_group_ID_publisher_documented_namespace') |>
 select(entity_id,chemical,mapping_relation,mapping_priority,mapping_source)
map<-bind_rows(bind_rows(maps),flat) |> distinct() |> arrange(desc(mapping_priority)) |>
 distinct(entity_id,chemical,.keep_all=TRUE)
stopifnot(!anyNA(map$chemical),!any(map$entity_id=='sider:CID100135113'),all(map$entity_id %in% inventory$entity_id))
write_csv(map,'data/curated/inventory_broad_STITCH_chemical_leads.csv')
cache<-file.path(root,'inventory_chemical_neighbours_v5.rds')
signature<-paste(source_meta$md5,paste(sort(unique(map$chemical)),collapse='|'),
 unname(tools::md5sum('scripts/helpers/stream_stitch_chemicals.R')),sep=':')
if(file.exists(cache)){
 scan<-readRDS(cache);stopifnot(identical(scan$signature,signature))
}else{
 scan<-stream_stitch_chemicals(source_path,unique(map$chemical));scan$signature<-signature;saveRDS(scan,cache)
}
stopifnot(length(scan$lines)>0L)
edges<-read_tsv(I(paste(c(scan$header,scan$lines),collapse='\n')),
 col_types=cols(chemical=col_character(),protein=col_character(),.default=col_integer()),quote='',progress=FALSE)
score_cols<-c('experimental','prediction','database','textmining','combined_score')
stopifnot(nrow(problems(edges))==0L,!anyDuplicated(edges[c('chemical','protein')]),
 all(edges$chemical %in% map$chemical),all(grepl('^9606[.]ENSP[0-9]+$',edges$protein)),
 all(vapply(edges[score_cols],function(x)all(!is.na(x)&x>=0L&x<=1000L),logical(1))))
edges<-edges |> mutate(chemical_namespace=if_else(startsWith(chemical,'CIDs'),'stereo_specific','merged_salt_isomer_group'),
 textmining_only=experimental==0L & prediction==0L & database==0L & textmining>0L,
 evidence_status='STITCH_association_not_measured_binding',source_release='STITCH 5.0 legacy',
 source_url=source_url,source_md5=source_meta$md5,retrieved_utc=source_meta$retrieved_utc)
links<-map |> inner_join(edges,by='chemical',relationship='many-to-many') |>
 left_join(inventory |> select(entity_id,origin,source_name,seed_tiers,cancers,structure_mapping_blocked),
 by='entity_id',relationship='many-to-one') |>
 mutate(compound_identity_layer=case_when(mapping_priority==3L~'current_key_supported_stereo_ID_historical_identity_not_revalidated',
 mapping_priority==1L~'current_PubChem_only_identity_review',TRUE~'merged_group_context_only'),
 experimental_channel_alone_establishes_direct_binding=FALSE,affinity_required_for_inclusion=FALSE)
stopifnot(!any(links$structure_mapping_blocked=='TRUE'),!anyDuplicated(links[c('entity_id','chemical','protein')]))
# Only direct Ensembl/UniProt accession aliases; BLAST-only matches are never accepted.
aliases<-read_tsv(alias_path,skip=1,col_names=c('protein','alias','alias_sources'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
stopifnot(nrow(problems(aliases))==0L,all(grepl('^9606[.]',aliases$protein)))
uniprot_pattern<-'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$'
accession_aliases<-aliases |> filter(protein %in% edges$protein,grepl(uniprot_pattern,alias),
 grepl('(^| )Ensembl_UniProt_AC( |$)|(^| )Ensembl_HGNC_UniProt_ID[(]supplied_by_UniProt[)]_AC( |$)',alias_sources)) |>
 transmute(protein,uniprot=alias,alias_sources) |> distinct()
lowest<-r(inputs[5]);all_levels<-r(inputs[6])
# Legacy-to-current projections retain every candidate. Ambiguous aliases do not enter automatic pathway joins.
review<-accession_aliases |> mutate(present_in_current_human_Reactome=uniprot %in% all_levels$uniprot,
 alias_release='STRING 10.0',alias_source_url=alias_url,alias_file_md5=alias_meta$md5,
 sequence_equivalence_established=FALSE)
counts<-review |> group_by(protein) |> summarise(direct_accession_aliases=n_distinct(uniprot),
 current_Reactome_supported_aliases=n_distinct(uniprot[present_in_current_human_Reactome]),.groups='drop')
protein_coverage<-tibble(protein=sort(unique(edges$protein))) |>
 left_join(counts,by='protein',relationship='one-to-one') |>
 mutate(across(c(direct_accession_aliases,current_Reactome_supported_aliases),~coalesce(.x,0L)),
 pathway_projection_status=case_when(direct_accession_aliases==0L~'no_direct_UniProt_alias',
 current_Reactome_supported_aliases==0L~'no_exact_current_Reactome_accession',
 current_Reactome_supported_aliases>1L~'ambiguous_current_accessions_hold',
 TRUE~'unique_Reactome_supported_historical_alias_projection'))
accepted<-review |> filter(present_in_current_human_Reactome) |>
 inner_join(protein_coverage |> filter(current_Reactome_supported_aliases==1L),by='protein',relationship='many-to-one')
stopifnot(!anyDuplicated(accepted$protein))
# Projections remain association-context tables, separate from DrugBank curated target roles.
projected<-links |> inner_join(accepted,by='protein',relationship='many-to-one')
write_csv(edges,'data/curated/STITCH_inventory_chemical_neighbour_source_rows.csv')
write_csv(links,'data/curated/inventory_broad_STITCH_association_evidence.csv',na='')
write_csv(review,'data/curated/STITCH_historical_UniProt_alias_review.csv')
write_csv(protein_coverage,'data/curated/STITCH_historical_protein_mapping_coverage.csv')
write_csv(projected,'data/curated/inventory_STITCH_unique_Reactome_alias_projections.csv',na='')
# Deduplicated entity/protein features prevent stereo/merged rows from becoming independent observations.
# Chemical identity strata are retained in pathway features and never collapsed across these layers.
protein_features<-projected |> distinct(entity_id,protein,uniprot,compound_identity_layer)
lowest_features<-protein_features |> inner_join(lowest,by='uniprot',relationship='many-to-many') |>
 distinct(entity_id,protein,uniprot,compound_identity_layer,pathway_id,pathway_name)
all_features<-protein_features |> inner_join(all_levels,by='uniprot',relationship='many-to-many') |>
 distinct(entity_id,protein,uniprot,compound_identity_layer,pathway_id,pathway_name)
write_csv(lowest_features,'data/curated/inventory_STITCH_Reactome_lowest_level_context.csv')
write_csv(all_features,'data/curated/inventory_STITCH_Reactome_all_level_context.csv')
edge_counts<-links |> group_by(entity_id) |> summarise(source_ID_protein_links=n(),
 distinct_legacy_proteins=n_distinct(protein),experimental_channel_links=sum(experimental>0L),
 database_channel_links=sum(database>0L),textmining_only_links=sum(textmining_only),
 current_key_supported_stereo_links=sum(mapping_priority==3L),merged_group_context_links=sum(mapping_priority==0L),.groups='drop')
coverage<-inventory |> select(entity_id,origin,source_name,seed_tiers,cancers) |>
 left_join(map |> count(entity_id,name='chemical_ID_leads'),by='entity_id',relationship='one-to-one') |>
 left_join(edge_counts,by='entity_id',relationship='one-to-one') |>
 mutate(across(where(is.numeric),~coalesce(.x,0L)),
 coverage_status=case_when(chemical_ID_leads==0L~'no_accepted_chemical_identifier_lead',
 source_ID_protein_links==0L~'identifier_lead_no_source_rows',TRUE~'source_associations_available'),
 missing_coverage_is_not_negative=TRUE)
stopifnot(nrow(coverage)==200L,!anyDuplicated(coverage$entity_id),
 all(coverage$entity_id %in% inventory$entity_id),identical(before,tools::md5sum(inputs)),
 sum(coverage$source_ID_protein_links)==nrow(links),
 !any(all_features$entity_id=='sider:CID100135113'))
write_csv(coverage,'data/curated/inventory_broad_STITCH_entity_coverage.csv')
write_csv(coverage |> count(origin,coverage_status,name='inventory_entities'),'outputs/20c_entity_coverage.csv')
write_csv(protein_coverage |> count(pathway_projection_status,name='legacy_proteins'),'outputs/20c_protein_mapping_coverage.csv')
write_csv(links |> group_by(compound_identity_layer) |> summarise(inventory_entities=n_distinct(entity_id),
 source_ID_protein_links=n(),distinct_proteins=n_distinct(protein),experimental_channel_links=sum(experimental>0L),
 database_channel_links=sum(database>0L),textmining_only_links=sum(textmining_only),.groups='drop'),
 'outputs/20c_channel_summary.csv')
write_csv(all_features |> group_by(compound_identity_layer) |> summarise(inventory_entities=n_distinct(entity_id),
 legacy_proteins=n_distinct(protein),current_UniProt_accessions=n_distinct(uniprot),
 pathway_IDs=n_distinct(pathway_id),.groups='drop'),'outputs/20c_pathway_coverage.csv')
write_csv(tibble(source=c('human_STITCH_detailed','human_STRING_10_aliases'),
 url=c(source_url,alias_url),md5=c(source_meta$md5,alias_meta$md5),
 retrieved_utc=c(source_meta$retrieved_utc,alias_meta$retrieved_utc),
 file_bytes=c(source_meta$file_bytes,alias_meta$file_bytes)),'outputs/20c_source_summary.csv')
write_csv(tibble(file=inputs,md5=unname(before)),'outputs/20c_input_summary.csv')
manifest<-r('config/database_manifest.csv');manifest$status[manifest$source_id=='stitch']<-
 'human_v5_detailed_inventory_neighbours_acquired_historical_alias_projection_reviewed'
write_csv(manifest,'config/database_manifest.csv',na='')
checklist<-r('config/RESKO_method_completion_checklist.csv')
checklist$project_status[checklist$step=='06']<-'DrugBank_and_broad_STITCH_inventory_links_available_STRING_network_pending'
checklist$next_work[checklist$step=='06']<-'Acquire human STRING PPI network; retain historical/current identifier and evidence layers separately'
write_csv(checklist,'config/RESKO_method_completion_checklist.csv')
report<-c('SCRIPT 20c - BROAD STITCH CHEMICAL NEIGHBOURS AND HISTORICAL PROTEIN ALIASES',
 paste('Inventory entities retained:',nrow(coverage)),paste('Entities with chemical identifier leads:',n_distinct(map$entity_id)),
 paste('Complete human source rows scanned:',scan$rows_scanned),paste('Distinct selected chemical/protein source rows:',nrow(edges)),
 paste('Inventory entities with associations:',n_distinct(links$entity_id)),paste('Inventory source-ID/protein links:',nrow(links)),
 paste('Distinct legacy human protein IDs:',n_distinct(edges$protein)),
 paste('Legacy proteins with unique current Reactome-supported accession projection:',nrow(accepted)),
 'All source score channels retained separately on original 0..1000 scale; no API neighbour cap or new score threshold.',
 'Merged chemical groups, current-PubChem-only identities and current-key-supported stereo leads remain separate.',
 'Experimental channel is association support, not proof of direct human binding. Source overlap with DrugBank/ChEMBL is not independent confirmation.',
 'Only direct Ensembl/UniProt accession aliases accepted; BLAST-only mappings excluded; multiple current accessions held.',
 'Historical alias to current Reactome is a labelled projection, not proof of sequence/isoform equivalence. Missing/ambiguous mappings remain explicit.',
 'Source-rich rows/projections/pathway features stay local; compact source, channel, mapping and coverage reports in outputs/20c_*.',
 'No drug excluded for missing affinity, no enrichment/BH test and no allostery claim.',
 'Core stage derives chemical leads from identity tables independently of optional eEF1A binding/literature scripts.',
 'Next: human STRING network acquisition and network-ready identifier harmonisation, followed by GO/DO integration.')
writeLines(report,'outputs/20c_integration_summary.txt')
writeLines(c('PASS: complete source scan, exact chemical/protein fields and bounded score validation; unique entity/source-ID/protein links.',
 'PASS: full 200-entity coverage and protected inputs unchanged; blocked sacrosidase structure excluded.',
 'PASS: chemical identity strata, direct-versus-BLAST alias provenance and ambiguous historical mappings retained.'),
 'outputs/20c_integration_validation.txt')
notes<-readLines('Notes.txt',warn=FALSE);start<-grep('^CURRENT CHECKPOINT',notes)[1];end<-grep('^PROJECT AIM',notes)[1]
stopifnot(start<end);writeLines(c(notes[seq_len(start-1)],'CURRENT CHECKPOINT - SCRIPT 20c BROAD STITCH INTEGRATION COMPLETE',report,'',notes[end:length(notes)]),'Notes.txt')
cat(paste0('\nSCRIPT 20c RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/20c_integrate_broad_STITCH_associations.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
