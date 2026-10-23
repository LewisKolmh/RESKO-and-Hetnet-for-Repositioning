# Reconcile qualifying mechanism queues; weak broad patterns remain context only.
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(stringr)})
source('scripts/helpers/cached_public_get.R')
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
inputs <- c('data/curated/RESKO_mechanism_SIDER_primary_group_review_queue.csv',
 'data/curated/RESKO_mechanism_SIDER_exploratory_group_review_queue.csv',
 'data/curated/sider_compound_index.csv','data/curated/drugbank_structure_identity_index.csv',
 'outputs/16f_target_evidence_handoff_summary.csv','data/curated/cancer_primary_seeds_expanded.csv',
 'data/curated/cancer_investigational_seeds_expanded.csv','data/curated/cancer_documented_seeds_expanded.csv',
 'data/curated/drugbank_drugs.csv','config/candidate_identity_exceptions.csv',
 'data/curated/RESKO_SIDER_new_candidate_group_review_queue.csv')
stopifnot(all(file.exists(inputs))); hashes <- tools::md5sum(inputs)
q <- bind_rows(r(inputs[1]),r(inputs[2])) |> distinct(flat,historical_name)
augmented_queue <- 'data/curated/RESKO23_augmented_IMiD_group_review_queue.csv'
if(file.exists(augmented_queue)) {
 inputs <- c(inputs,augmented_queue); hashes <- tools::md5sum(inputs)
 extra <- r(augmented_queue) |> distinct(flat) |> left_join(
 r(inputs[3]) |> distinct(stitch_flat_id,sider_drug_name),
 by=c('flat'='stitch_flat_id')) |> transmute(flat,historical_name=sider_drug_name)
 q <- bind_rows(q,extra) |> distinct(flat,historical_name)
}
index <- r(inputs[3]) |> distinct()
leads <- q |> left_join(index,by=c('flat'='stitch_flat_id'),relationship='one-to-many') |>
 mutate(lead_id=paste(flat,sider_pubchem_cid,sep=':'))
stopifnot(n_distinct(leads$flat)==nrow(q),!anyDuplicated(leads$lead_id),!anyNA(leads$sider_pubchem_cid))
scalar <- function(x) if(length(x)) as.character(x[[1]]) else NA_character_
normal <- function(x) str_to_lower(str_squish(x))
properties <- aliases <- logs <- list()
for(i in seq_len(nrow(leads))) {
 cid <- leads$sider_pubchem_cid[i]; base <- paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',cid)
 url <- paste0(base,'/property/InChIKey,InChI,MolecularFormula,MolecularWeight,IUPACName,Title/JSON')
 z <- public_cache_get(url,paste0('data/raw/section23_identity/properties_',cid,'.rds'))
 logs[[length(logs)+1L]] <- tibble(lead_id=leads$lead_id[i],flat=leads$flat[i],requested_cid=cid,
 url=z$url,status=z$status,retrieved_utc=z$retrieved_utc)
 if(z$status==200L) for(p in z$body$PropertyTable$Properties)
  properties[[length(properties)+1L]] <- tibble(lead_id=leads$lead_id[i],flat=leads$flat[i],requested_cid=cid,
   returned_cid=scalar(p$CID),current_inchikey=scalar(p$InChIKey),
   current_inchi=scalar(p$InChI),pubchem_title=scalar(p$Title),
   molecular_formula=scalar(p$MolecularFormula),molecular_weight=scalar(p$MolecularWeight))
 z <- public_cache_get(paste0(base,'/synonyms/JSON'),paste0('data/raw/section23_identity/synonyms_',cid,'.rds'))
 logs[[length(logs)+1L]] <- tibble(lead_id=leads$lead_id[i],flat=leads$flat[i],requested_cid=cid,
 url=z$url,status=z$status,retrieved_utc=z$retrieved_utc)
 if(z$status==200L) for(p in z$body$InformationList$Information)
  aliases[[length(aliases)+1L]] <- tibble(lead_id=leads$lead_id[i],flat=leads$flat[i],returned_cid=scalar(p$CID),
    alias=as.character(unlist(p$Synonym,use.names=FALSE)))
 message('23c identity ',i,'/',nrow(leads),': ',leads$historical_name[i])
}
pc <- bind_rows(properties); syn <- bind_rows(aliases)
stopifnot(!anyDuplicated(pc$lead_id))
registry <- leads |> left_join(pc |> select(-flat),by='lead_id')
support <- registry |> select(lead_id,flat,historical_name,returned_cid) |>
 left_join(syn |> select(-flat) |> rename(synonym_returned_cid=returned_cid),by='lead_id',relationship='many-to-many') |>
 group_by(lead_id) |> summarise(historical_name_alias_supported=any(
  normal(historical_name)==normal(alias) & returned_cid==synonym_returned_cid,na.rm=TRUE),.groups='drop')
db_index <- r(inputs[4]) |> filter(!is.na(record_inchikey),nzchar(record_inchikey))
links <- registry |> filter(!is.na(current_inchikey)) |>
 inner_join(db_index,by=c('current_inchikey'='record_inchikey'),relationship='many-to-many') |>
 select(lead_id,flat,current_inchikey,record_drugbank_id,record_name) |> distinct()
collapse <- function(x) paste(sort(unique(na.omit(x[nzchar(x)]))),collapse='|')
db_summary <- links |> group_by(lead_id) |> summarise(
 exact_key_DrugBank_ids=collapse(record_drugbank_id),exact_key_DrugBank_names=collapse(record_name),
 exact_key_DrugBank_records=n_distinct(record_drugbank_id),.groups='drop')
inventory <- r(inputs[5]); stopifnot(nrow(inventory)==200L,!anyDuplicated(inventory$entity_id))
inv_links <- registry |> filter(!is.na(current_inchikey)) |>
 inner_join(inventory |> filter(!is.na(current_inchikey),nzchar(current_inchikey),
 structure_mapping_blocked=='FALSE'),by='current_inchikey',relationship='many-to-many') |>
 select(lead_id,flat,entity_id,current_inchikey,origin,seed_tiers,cancers,identity_status)
inv_summary <- inv_links |> group_by(lead_id) |> summarise(
 existing_inventory_entities=collapse(entity_id),existing_inventory_cancers=collapse(cancers),
 existing_inventory_tiers=collapse(seed_tiers),.groups='drop')
tiers <- bind_rows(lapply(inputs[6:8],r)) |>
 distinct(cancer_id,cancer_name,drugbank_id,drug_name,decision,evidence_url,evidence_excerpt,
 regulatory_status,source_export_date)
known <- links |> inner_join(tiers,by=c('record_drugbank_id'='drugbank_id'),relationship='many-to-many') |>
 mutate(link_basis='exact current full InChIKey; historical group equivalence remains unverified')
db <- r(inputs[9])
indications <- links |> left_join(db,by=c('record_drugbank_id'='drugbank_id'),relationship='many-to-one')
exceptions <- r(inputs[10])
registry <- registry |> left_join(support,by='lead_id') |> left_join(db_summary,by='lead_id') |>
 left_join(inv_summary,by='lead_id') |>
 left_join(exceptions,by=c('flat'='stitch_flat_id')) |>
 mutate(requested_returned_CID_equal=!is.na(returned_cid)&sider_pubchem_cid==returned_cid,
  structure_key_available=!is.na(current_inchikey)&grepl('^[A-Z]{14}-[A-Z]{10}-[A-Z]$',current_inchikey),
  exact_key_DrugBank_records=replace_na(exact_key_DrugBank_records,0L),
  historical_formulation_equivalence_proven=FALSE,
  identity_disposition=case_when(!is.na(identity_exception)~'conflict_hold',
   !structure_key_available~'structure_unavailable_hold',
   !requested_returned_CID_equal~'CID_redirect_review_hold',
   !historical_name_alias_supported~'historical_name_support_review_hold',
   exact_key_DrugBank_records==0L~'current_structure_only_no_DrugBank_match',
   TRUE~'alias_supported_current_full_key_match_historical_form_unverified'),
  inventory_navigation=if_else(is.na(existing_inventory_entities),
   'not_in_original_inventory_by_accepted_full_key','existing_inventory_current_structure_link'),
  evidence_limit='current representation and alias evidence, not historical formulation proof or target binding')
# Include weak broad-pattern groups in coverage accounting, without advancing them.
weak <- r(inputs[11]) |> distinct(stitch_flat_id) |>
 mutate(mechanism_identity_reviewed=stitch_flat_id %in% registry$flat,
 status=if_else(mechanism_identity_reviewed,'also_in_qualifying_mechanism_queue',
 'weak_two_term_pattern_only_not_advanced'))
write_csv(registry,'data/curated/RESKO23_candidate_identity_registry.csv')
write_csv(links,'data/curated/RESKO23_candidate_DrugBank_structure_links.csv')
write_csv(inv_links,'data/curated/RESKO23_existing_inventory_structure_links.csv')
write_csv(known,'data/curated/RESKO23_known_seed_treatment_evidence.csv')
write_csv(indications,'data/curated/RESKO23_DrugBank_indication_context.csv')
write_csv(syn,'data/curated/RESKO23_candidate_alias_evidence.csv')
write_csv(bind_rows(logs),'data/curated/RESKO23_identity_API_provenance.csv')
write_csv(weak,'data/curated/RESKO23_weak_broad_pattern_dispositions.csv')
write_csv(registry |> count(identity_disposition,inventory_navigation),
 'outputs/23c_identity_coverage_summary.csv')
write_csv(weak |> count(status),'outputs/23c_weak_pattern_coverage_summary.csv')
stopifnot(n_distinct(registry$flat)==nrow(q),!anyDuplicated(registry$lead_id),identical(hashes,tools::md5sum(inputs)))
cache <- list.files('data/raw/section23_identity',pattern='\\.rds$',full.names=TRUE)
write_csv(tibble(path=c(inputs,cache),md5=unname(tools::md5sum(c(inputs,cache)))),
 'outputs/23c_input_hashes_summary.csv')
print(registry |> count(identity_disposition,inventory_navigation),width=Inf)
message('23c complete: all qualifying groups accounted for; weak patterns retained separately.')
