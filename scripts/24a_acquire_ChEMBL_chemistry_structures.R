# Exact-key ChEMBL chemistry source acquisition. No fingerprints or ranking yet.
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr)})
source('scripts/helpers/target_public_get.R')
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
inputs <- c('outputs/16f_target_evidence_handoff_summary.csv',
 'data/curated/RESKO23_completed_chemistry_proposals.csv',
 'data/curated/RESKO23_candidate_identity_registry.csv',
 'config/RESKO23_supplemental_known_treatments.csv',
 'data/curated/RESKO23_completed_SIDER_handoff.csv')
stopifnot(all(file.exists(inputs)));hashes <- tools::md5sum(inputs)
policy <- tibble(parameter=c('structure_match','entity_scope','missingness','forms',
 'parent_hierarchy','fingerprint_plan','encoding_limit','stereochemistry_limit','source','affinity'),
 value=c('exact full standard InChIKey only; retain all returned exact-key records',
 'original 200 entries plus unrepresented 23f proposals; supplemental thalidomide as separate comparator',
 'all entries retained in coverage; missing structures are not inactive or zero similarity',
 'no fragment stripping, neutralisation, stereochemistry removal or name-based merge',
 'record parent ID only; never replace queried structure with its parent',
 'planned ChemmineR sdf2ap -> desc2fp(descnames=2048) -> binary Tanimoto',
 'paper specifies atom-pair 2048 bits but not encoding/vocabulary; this implementation is an explicit choice',
 '2D atom-pair similarities cannot establish stereo/form equivalence; full keys remain separate',
 'live ChEMBL snapshot cached with status release before/after and hashes; not a 2018 snapshot',
 'no binding-affinity filter'))
write_csv(policy,'config/RESKO_chemistry_structure_policy.csv')
inv <- r(inputs[1]);proposal <- r(inputs[2]);id <- r(inputs[3]);supp <- r(inputs[4])
stopifnot(nrow(inv)==200L,!anyDuplicated(inv$entity_id))
entities <- inv |> transmute(entity_id,name=source_name,origin,
 expected_inchikey=current_inchikey,structure_mapping_blocked=structure_mapping_blocked=='TRUE',
 seed_tiers,cancers,identity_context=identity_status,
 historical_equivalence_proven=FALSE,scope='original_inventory')
# Proposals matching an existing full key remain navigational links rather than duplicate entities.
new <- proposal |> filter(is.na(existing_inventory_entities)|!nzchar(existing_inventory_entities)) |>
 transmute(entity_id=paste0('proposal23:',flat),name=historical_name,origin='23f_chemistry_proposal',
 expected_inchikey=current_inchikeys,structure_mapping_blocked=FALSE,
 seed_tiers=NA_character_,cancers=NA_character_,identity_context="current_representation_historical_form_unverified",
 historical_equivalence_proven=FALSE,scope='proposed_current_representation')
control <- id |> filter(flat %in% supp$flat,historical_name_alias_supported=='TRUE',
 requested_returned_CID_equal=='TRUE') |>
 transmute(entity_id=paste0('comparator23:',flat),name=historical_name,
 origin='label_verified_known_treatment_comparator',expected_inchikey=current_inchikey,
 structure_mapping_blocked=FALSE,seed_tiers='supplemental_known_labelled_treatment',
 cancers='Multiple myeloma',identity_context='current_representation_historical_form_unverified',
 historical_equivalence_proven=FALSE,scope='supplemental_known_treatment_comparator')
entities <- bind_rows(entities,new,control)
stopifnot(!anyDuplicated(entities$entity_id))
entities <- entities |> mutate(key_format_valid=!is.na(expected_inchikey)&
 grepl('^[A-Z]{14}-[A-Z]{10}-[A-Z]$',expected_inchikey),
 query_eligible=key_format_valid & !structure_mapping_blocked)
keys <- sort(unique(entities$expected_inchikey[entities$query_eligible]))
dir.create('data/raw/chemistry24',recursive=TRUE,showWarnings=FALSE)
status_url <- 'https://www.ebi.ac.uk/chembl/api/data/status.json'
before <- target_public_get(status_url,'data/raw/chemistry24/status_before.rds')
writeLines(capture.output(str(before$body)),'data/raw/chemistry24/status_before.txt')
scalar <- function(z) if(length(z)) as.character(z[[1]]) else NA_character_
records <- logs <- list()
batches <- split(keys,ceiling(seq_along(keys)/10L))
for(i in seq_along(batches)) {
 batch <- batches[[i]]
 url <- paste0('https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__standard_inchi_key__in=',
 paste(batch,collapse=','),'&limit=1000')
 # Cache filenames contain requested key signature, avoiding silent positional-batch reuse.
 sigfile <- tempfile();writeLines(batch,sigfile);sig <- unname(tools::md5sum(sigfile));unlink(sigfile)
 path <- paste0('data/raw/chemistry24/exact_keys_',sig,'.rds')
 z <- target_public_get(url,path)
 stopifnot(is.null(z$body$page_meta[["next"]]),
 as.integer(z$body$page_meta$total_count)==length(z$body$molecules))
 logs[[i]] <- tibble(url=z$url,status=z$status,retrieved_utc=z$retrieved_utc,cache_path=path,
 requested_keys=length(batch),returned_records=length(z$body$molecules))
 for(m in z$body$molecules) {
  s <- m$molecule_structures
  stopifnot(s$standard_inchi_key %in% batch)
  records[[length(records)+1L]] <- tibble(chembl_id=scalar(m$molecule_chembl_id),
   pref_name=scalar(m$pref_name),molecule_type=scalar(m$molecule_type),
   standard_inchikey=scalar(s$standard_inchi_key),standard_inchi=scalar(s$standard_inchi),
   canonical_smiles=scalar(s$canonical_smiles),molfile=scalar(s$molfile),
   parent_chembl_id=scalar(m$molecule_hierarchy$parent_chembl_id),
   source_url=z$url,retrieved_utc=z$retrieved_utc,
   chembl_release=scalar(before$body$chembl_db_version))
 }
 message('24a ChEMBL exact-key batch ',i,'/',length(batches))
}
after <- target_public_get(status_url,'data/raw/chemistry24/status_after.rds')
stopifnot(identical(before$body$chembl_db_version,after$body$chembl_db_version))
rec <- bind_rows(records)
stopifnot(!anyDuplicated(rec$chembl_id),all(rec$standard_inchikey %in% keys))
links <- entities |> filter(query_eligible) |>
 inner_join(rec,by=c('expected_inchikey'='standard_inchikey'),relationship='many-to-many')
counts <- links |> group_by(entity_id) |> summarise(ChEMBL_records=n_distinct(chembl_id),
 SMILES_records=sum(!is.na(canonical_smiles)&nzchar(canonical_smiles)),
 molfile_records=sum(!is.na(molfile)&nzchar(molfile)),.groups='drop')
coverage <- entities |> left_join(counts,by='entity_id') |>
 mutate(across(c(ChEMBL_records,SMILES_records,molfile_records),~replace_na(.x,0L)),
 structure_status=case_when(structure_mapping_blocked~'blocked_identity_hold',
 !key_format_valid~'no_valid_current_structure_key',ChEMBL_records==0L~'no_exact_ChEMBL_match',
 molfile_records>0L~'exact_ChEMBL_molfile_available',
 SMILES_records>0L~'exact_ChEMBL_SMILES_only',TRUE~'record_without_usable_structure'))
write_csv(coverage,'data/curated/RESKO24_chemistry_entity_coverage.csv')
# Preserve candidate cancer/cohort context separately from original seed indications.
context <- r(inputs[5]) |>
 filter(final_section23_action=='advance_current_representation_to_chemistry_annotation') |>
 select(flat,historical_name,cancer_id,cohort_id,current_inchikeys) |>
 inner_join(entities |> select(entity_id,expected_inchikey),
 by=c('current_inchikeys'='expected_inchikey'),relationship='many-to-many')
write_csv(context,'data/curated/RESKO24_proposal_cancer_cohort_context.csv')
saveRDS(rec,'data/curated/RESKO24_ChEMBL_structure_records.rds')
write_csv(links |> select(-molfile),'data/curated/RESKO24_entity_ChEMBL_structure_links.csv')
write_csv(bind_rows(logs),'data/curated/RESKO24_ChEMBL_API_provenance.csv')
# Preserve raw supplied mol blocks, one SDF per record. Parsing/chemistry checks are 24b.
dir.create('data/raw/chemistry24/sdf',showWarnings=FALSE)
for(i in which(!is.na(rec$molfile)&nzchar(rec$molfile))) {
 mol <- sub('[\r\n]+$','',rec$molfile[i])
 writeLines(c(mol,'$$$$'),paste0('data/raw/chemistry24/sdf/',rec$chembl_id[i],'.sdf'),useBytes=TRUE)
}
packages <- c('ChemmineR','ChemmineOB','rcdk','fingerprint','BiocManager')
write_csv(tibble(package=packages,available=vapply(packages,requireNamespace,quietly=TRUE,FUN.VALUE=logical(1))),
 'outputs/24a_chemistry_tool_coverage.csv')
write_csv(coverage |> count(scope,structure_status),'outputs/24a_structure_coverage_summary.csv')
write_csv(tibble(chembl_release=scalar(before$body$chembl_db_version),
 requested_full_keys=length(keys),returned_records=nrow(rec),
 matched_full_keys=n_distinct(rec$standard_inchikey),entities=nrow(entities)),
 'outputs/24a_ChEMBL_acquisition_summary.csv')
stopifnot(identical(hashes,tools::md5sum(inputs)))
protected <- c(inputs,'config/RESKO_chemistry_structure_policy.csv',
 list.files('data/raw/chemistry24',pattern='\\.rds$',full.names=TRUE))
write_csv(tibble(path=protected,md5=unname(tools::md5sum(protected))),
 'outputs/24a_input_hashes_summary.csv')
print(coverage |> count(scope,structure_status),width=Inf)
message('24a complete: structures acquired; parsing, fingerprints and Tanimoto remain 24b work.')
