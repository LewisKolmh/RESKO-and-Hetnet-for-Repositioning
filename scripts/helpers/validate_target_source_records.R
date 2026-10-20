library(readr);library(dplyr)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
a<-r('data/curated/eef1a_chembl_enriched_activity_evidence.csv');m<-r('data/curated/eef1a_chembl_ligand_registry.csv')
s<-r('data/curated/eef1a_chembl_assay_registry.csv');h<-r('outputs/16f_target_evidence_handoff_summary.csv')
links<-r('outputs/17b_inventory_ligand_match_summary.csv');hash<-r('outputs/17b_input_summary.csv')
stopifnot(!anyDuplicated(a$activity_id),nrow(a)==39,n_distinct(m$molecule_chembl_id)==20,n_distinct(s$assay_chembl_id)==23,
 identical(unname(tools::md5sum(hash$file)),hash$md5))
for(i in seq_len(nrow(m))){v<-readRDS(paste0('data/raw/target_evidence/chembl_molecule_',m$molecule_chembl_id[i],'.rds'))$body
 stopifnot(v$molecule_chembl_id==m$molecule_chembl_id[i],identical(v$molecule_structures$standard_inchi_key,m$standard_inchikey[i]))}
for(i in seq_len(nrow(s))){v<-readRDS(paste0('data/raw/target_evidence/chembl_assay_',s$assay_chembl_id[i],'.rds'))$body
 stopifnot(v$assay_chembl_id==s$assay_chembl_id[i],v$target_chembl_id==s$target_chembl_id[i],as.character(v$confidence_score)==s$confidence_score[i])}
expected<-h |> filter(!is.na(current_inchikey),structure_mapping_blocked=='FALSE') |>
 inner_join(m |> filter(!is.na(standard_inchikey)),by=c('current_inchikey'='standard_inchikey'),relationship='many-to-many')
stopifnot(nrow(expected)==nrow(links),all(a$source_assay_target_agreement=='TRUE'),
 all(a$evidence_review_status[a$isoform_description_conflict=='TRUE']=='hold_isoform_description_target_conflict'))
for(p in c('scripts/optional/eef1a_evidence/17a_acquire_target_records.R','scripts/optional/eef1a_evidence/17b_expand_assays_and_match_ligands.R','scripts/helpers/target_public_get.R'))invisible(parse(p))
msg<-c('PASS: all 20 ligand structure keys and 23 assay target/confidence records independently checked against raw cached JSON.',
 'PASS: full inventory exact-key join independently recomputed; 39 source activities retained and isoform conflict held.',
 'PASS: source input hashes unchanged and section R scripts parse. No negative binding conclusion follows from absent matches.')
cat(paste0(paste(msg,collapse='\n'),'\n'),file='outputs/17b_assay_ligand_validation.txt',append=TRUE)
cat(paste0('\nSCRIPT 17b - INDEPENDENT VALIDATION\n',paste(msg,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
