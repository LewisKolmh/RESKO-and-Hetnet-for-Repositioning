# Core RESKO: DrugBank protein roles for the full existing inventory, then human pathways.
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(xml2)
source('scripts/helpers/stream_selected_drugbank_fields.R')
source('scripts/helpers/validate_selected_drugbank_stream.R')
r <- function(path) read_csv(path, col_types = cols(.default = col_character()))
stopifnot(validate_selected_drugbank_stream())
inputs <- c('outputs/16f_target_evidence_handoff_summary.csv',
 'data/curated/drugbank_structure_identity_index.csv',
 'data/curated/Reactome_human_lowest_level_annotations.csv',
 'data/curated/Reactome_human_all_level_annotations.csv',
 'config/RESKO_evidence_scope.csv')
before <- tools::md5sum(inputs)
inventory <- r(inputs[1]); index <- r(inputs[2])
manifest <- r('config/database_manifest.csv'); db <- manifest |> filter(source_id == 'drugbank')
stopifnot(nrow(db) == 1L, file.exists(db$local_file), nrow(inventory) == 200L,
 !anyDuplicated(inventory$entity_id))
xml_path <- db$local_file[[1]]
# Explicit DB identifiers already reviewed in earlier identity stages; no name-only merge.
entity_map <- inventory |> select(entity_id, origin, source_name, drugbank_ids, current_inchikey,
 identity_status, structure_mapping_blocked, seed_tiers, cancers) |>
 filter(!is.na(drugbank_ids), nzchar(drugbank_ids)) |>
 separate_longer_delim(drugbank_ids, delim = '|') |> rename(drugbank_id = drugbank_ids) |>
 left_join(index, by = c('drugbank_id' = 'record_drugbank_id'), relationship = 'many-to-one') |>
 mutate(mapping_scope = case_when(origin == 'seed' ~ 'reviewed_seed_primary_DrugBank_record',
 identity_status == 'named_biologic_UNII_verified_historical_CID_conflict' ~ 'named_biologic_context_only',
 !is.na(current_inchikey) & current_inchikey == record_inchikey & structure_mapping_blocked != 'TRUE' ~
 'candidate_current_full_InChIKey_agreement', TRUE ~ 'identity_hold'))
stopifnot(!anyNA(entity_map$record_name), !anyDuplicated(entity_map[c('entity_id','drugbank_id')]),
 !any(entity_map$mapping_scope == 'identity_hold'), all(grepl('^DB[0-9]+$', entity_map$drugbank_id)))
ids <- sort(unique(entity_map$drugbank_id))
message('DrugBank records selected: ', length(ids), ' for ', n_distinct(entity_map$entity_id), ' entities')
root <- 'data/raw/drugbank/protein_inventory'; dir.create(root, recursive = TRUE, showWarnings = FALSE)
message('Hashing source XML for reproducible extraction/cache validation')
source_md5 <- unname(tools::md5sum(xml_path))
source_stat <- file.info(xml_path)
signature <- paste(source_md5, paste(ids, collapse = '|'),
 unname(tools::md5sum('scripts/helpers/stream_selected_drugbank_fields.R')), 'extractor-v1', sep = ':')
cache <- file.path(root, 'selected_protein_roles.rds')
if (file.exists(cache)) {
 extracted <- readRDS(cache)
 if (!identical(extracted$signature, signature)) stop('Protein extraction cache/input changed; inspect before refreshing')
} else {
 record_rows <- list(); role_rows <- list(); protein_rows <- list(); article_rows <- list(); link_rows <- list()
 scalar <- function(node, path) {
  value <- xml_text(xml_find_first(node, path))
  if (length(value) && !is.na(value) && nzchar(trimws(value))) value else NA_character_
 }
 process <- function(node) {
  id <- scalar(node, './drugbank-id[@primary="true"]')
  record_rows[[length(record_rows) + 1L]] <<- tibble(drugbank_id = id,
   record_name = scalar(node, './name'), record_type = xml_attr(node, 'type'),
   record_updated = xml_attr(node, 'updated'))
  for (role in c('target', 'enzyme', 'carrier', 'transporter')) {
   nodes <- xml_find_all(node, paste0('./', role, 's/', role))
   for (j in seq_along(nodes)) {
    interaction <- nodes[[j]]; interaction_id <- scalar(interaction, './id')
    row_id <- paste(id, role, j, sep = ':')
    actions <- xml_text(xml_find_all(interaction, './actions/action'))
    articles <- xml_find_all(interaction, './references/articles/article')
    urls <- xml_find_all(interaction, './references/links/link')
    role_rows[[length(role_rows) + 1L]] <<- tibble(source_role_row_id = row_id,
     drugbank_id = id, protein_role = role, role_position = as.character(j),
     source_interaction_id = interaction_id, source_interaction_name = scalar(interaction, './name'),
     entry_organism = scalar(interaction, './organism'), source_known_action = scalar(interaction, './known-action'),
     source_actions = if (length(actions)) paste(actions, collapse = '|') else NA_character_,
     source_article_count = length(articles), source_link_count = length(urls))
    for (article in articles) article_rows[[length(article_rows) + 1L]] <<- tibble(source_role_row_id = row_id,
     reference_id = scalar(article, './ref-id'), pubmed_id = scalar(article, './pubmed-id'),
     citation = scalar(article, './citation'))
    for (link in urls) link_rows[[length(link_rows) + 1L]] <<- tibble(source_role_row_id = row_id,
     reference_id = scalar(link, './ref-id'), title = scalar(link, './title'), url = scalar(link, './url'))
    polypeptides <- xml_find_all(interaction, './polypeptide')
    for (k in seq_along(polypeptides)) {
     protein <- polypeptides[[k]]
     protein_rows[[length(protein_rows) + 1L]] <<- tibble(source_role_row_id = row_id,
      source_protein_row_id = paste(row_id, k, sep = ':'), polypeptide_id = xml_attr(protein, 'id'),
      polypeptide_id_source = xml_attr(protein, 'source'),
      protein_name = scalar(protein, './name'), gene_name = scalar(protein, './gene-name'),
      protein_organism = scalar(protein, './organism'),
      protein_taxon = xml_attr(xml_find_first(protein, './organism'), 'ncbi-taxonomy-id'))
    }
   }
  }
 }
 scan <- stream_selected_drugbank_fields(xml_path, ids,
  c('drugbank-id','name','targets','enzymes','carriers','transporters'), process)
 extracted <- list(signature = signature, scan = scan, records = bind_rows(record_rows),
  roles = bind_rows(role_rows), proteins = bind_rows(protein_rows),
  articles = bind_rows(article_rows), links = bind_rows(link_rows),
  source_md5 = source_md5, source_file_bytes = as.character(source_stat$size),
  extracted_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE))
 saveRDS(extracted, cache)
}
stopifnot(nrow(extracted$records) == length(ids), setequal(extracted$records$drugbank_id, ids),
 !anyDuplicated(extracted$roles$source_role_row_id), !anyDuplicated(extracted$proteins$source_protein_row_id),
 all(extracted$proteins$source_role_row_id %in% extracted$roles$source_role_row_id))
# Accept exact UniProt accessions only; do not strip isoform suffixes or infer human from gene symbols.
uniprot_pattern <- '^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$'
proteins <- extracted$roles |> left_join(extracted$proteins, by = 'source_role_row_id', relationship = 'one-to-many') |>
 mutate(accession_is_UniProt = !is.na(polypeptide_id_source) &
 polypeptide_id_source %in% c('Swiss-Prot','TrEMBL','UniProtKB') &
 !is.na(polypeptide_id) & grepl(uniprot_pattern, polypeptide_id),
 taxon_label_conflict = coalesce(protein_taxon == '9606' &
 (!protein_organism %in% c('Humans','Homo sapiens') | !entry_organism %in% c('Humans','Homo sapiens')), FALSE),
 human_mapping_status = case_when(is.na(source_protein_row_id) ~ 'no_polypeptide_identifier',
 taxon_label_conflict ~ 'organism_conflict_hold',
 is.na(protein_taxon) ~ 'taxon_missing_hold',
 protein_taxon != '9606' ~ 'nonhuman_context',
 !accession_is_UniProt ~ 'human_accession_unresolved', TRUE ~ 'human_exact_UniProt_accession'),
 uniprot = if_else(human_mapping_status == 'human_exact_UniProt_accession', polypeptide_id, NA_character_),
 evidence_layer = 'DrugBank_curated_role_annotation_not_new_binding_measurement',
 role_interpretation = if_else(protein_role == 'target',
 'Source-labelled target; not necessarily selective or primary therapeutic target',
 'Metabolism_transport_or_carrier_context_not_automatically_therapeutic_target'),
 source_export_date = db$export_date[[1]], source_xml_md5 = source_md5,
 source_drug_url = paste0('https://go.drugbank.com/drugs/', drugbank_id),
 direct_human_binding_affinity = NA_real_, affinity_missing_is_inactive = FALSE)
entity_proteins <- entity_map |> inner_join(proteins, by = 'drugbank_id', relationship = 'many-to-many')
human <- entity_proteins |> filter(human_mapping_status == 'human_exact_UniProt_accession')
# Named-biologic/historical identity context is retained separately from strict candidate links.
strict_human <- human |> filter(mapping_scope != 'named_biologic_context_only')
context_human <- human |> filter(mapping_scope == 'named_biologic_context_only')
lowest <- r(inputs[3]); all_levels <- r(inputs[4])
stopifnot(!anyDuplicated(lowest), !anyDuplicated(all_levels), all(lowest$species == 'Homo sapiens'))
join_pathways <- function(annotation, level) {
 strict_human |> inner_join(annotation, by = 'uniprot', relationship = 'many-to-many') |>
  mutate(pathway_level_scope = level, pathway_interpretation =
   'Protein membership only; no inference of drug pathway inhibition, selectivity or allostery')
}
lowest_links <- join_pathways(lowest, 'lowest_level')
all_links <- join_pathways(all_levels, 'all_levels_including_ancestors')
# Distinct features remain per role: pharmacokinetic proteins are never silently mixed into target sets.
features <- all_links |> distinct(entity_id, mapping_scope, protein_role, pathway_id, pathway_name)
counts <- entity_proteins |> group_by(entity_id) |> summarise(
 DrugBank_source_role_rows = n_distinct(source_role_row_id),
 human_protein_accessions = n_distinct(uniprot[!is.na(uniprot)]),
 human_target_accessions = n_distinct(uniprot[protein_role == 'target' & !is.na(uniprot)]),
 organism_conflict_rows = sum(human_mapping_status == 'organism_conflict_hold'), .groups = 'drop')
pathway_counts <- features |> group_by(entity_id) |> summarise(
 all_role_pathway_features = n_distinct(pathway_id),
 target_role_pathway_features = n_distinct(pathway_id[protein_role == 'target']), .groups = 'drop')
coverage <- inventory |> select(entity_id,origin,source_name,seed_tiers,cancers,identity_status,structure_mapping_blocked) |>
 left_join(entity_map |> count(entity_id,name='DrugBank_record_links'),by='entity_id',relationship='one-to-one') |>
 left_join(counts,by='entity_id',relationship='one-to-one') |>
 left_join(pathway_counts,by='entity_id',relationship='one-to-one') |>
 mutate(across(c(DrugBank_record_links,DrugBank_source_role_rows,human_protein_accessions,human_target_accessions,
 organism_conflict_rows,all_role_pathway_features,target_role_pathway_features),~coalesce(.x,0L)),
 protein_coverage_status=case_when(DrugBank_record_links==0L~'no_reviewed_DrugBank_identifier',
 DrugBank_source_role_rows==0L~'record_present_no_protein_roles',
 human_protein_accessions==0L~'roles_present_no_accepted_human_UniProt',TRUE~'human_protein_annotations_available'),
 affinity_required_for_RESKO=FALSE,missing_coverage_is_not_negative=TRUE)
stopifnot(nrow(coverage)==200L,!anyDuplicated(coverage$entity_id),
 all(lowest_links$protein_taxon=='9606'),all(lowest_links$species=='Homo sapiens'),
 !any(lowest_links$mapping_scope=='named_biologic_context_only'),
 !any(lowest_links$entity_id=='sider:CID100135113'),
 all(lowest_links$source_role_row_id %in% extracted$roles$source_role_row_id),
 identical(before,tools::md5sum(inputs)),
 file.info(xml_path)$size==source_stat$size,as.numeric(file.info(xml_path)$mtime)==as.numeric(source_stat$mtime))
write_csv(entity_map,'data/curated/inventory_DrugBank_record_links.csv',na='')
write_csv(extracted$records,'data/curated/inventory_DrugBank_record_metadata.csv',na='')
write_csv(proteins,'data/curated/DrugBank_selected_protein_role_annotations.csv',na='')
write_csv(extracted$articles,'data/curated/DrugBank_selected_protein_role_references.csv',na='')
write_csv(extracted$links,'data/curated/DrugBank_selected_protein_role_reference_links.csv',na='')
write_csv(entity_proteins,'data/curated/inventory_DrugBank_protein_role_evidence.csv',na='')
write_csv(strict_human,'data/curated/inventory_DrugBank_human_protein_links.csv',na='')
write_csv(context_human,'data/curated/inventory_DrugBank_named_biologic_context.csv',na='')
write_csv(lowest_links,'data/curated/inventory_DrugBank_Reactome_lowest_level_evidence.csv',na='')
write_csv(all_links,'data/curated/inventory_DrugBank_Reactome_all_level_evidence.csv',na='')
write_csv(features,'data/curated/inventory_DrugBank_Reactome_role_features.csv',na='')
write_csv(coverage,'data/curated/inventory_DrugBank_Reactome_coverage.csv',na='')
write_csv(coverage |> count(origin,protein_coverage_status,name='inventory_entities'),'outputs/20b_entity_coverage.csv')
write_csv(proteins |> count(protein_role,human_mapping_status,name='source_annotation_rows'),'outputs/20b_role_coverage.csv')
write_csv(strict_human |> group_by(protein_role) |> summarise(inventory_entities=n_distinct(entity_id),
 distinct_human_proteins=n_distinct(uniprot),entity_annotation_rows=n(),.groups='drop'),
 'outputs/20b_human_role_summary.csv')
write_csv(features |> group_by(protein_role) |> summarise(inventory_entities=n_distinct(entity_id),
 distinct_pathways=n_distinct(pathway_id),entity_pathway_features=n(),.groups='drop'),
 'outputs/20b_pathway_coverage.csv')
write_csv(tibble(source_xml_md5=source_md5,source_file_bytes=extracted$source_file_bytes,
 export_date=db$export_date[[1]],extracted_utc=extracted$extracted_utc,
 records_scanned=extracted$scan$records_scanned,selected_DrugBank_records=nrow(extracted$records),
 source_role_rows=nrow(extracted$roles),source_reference_rows=nrow(extracted$articles),
 source_reference_link_rows=nrow(extracted$links)),'outputs/20b_source_summary.csv')
write_csv(tibble(file=inputs,md5=unname(before)),'outputs/20b_input_summary.csv')
manifest$status[manifest$source_id=='drugbank']<-'registered_export_inventory_protein_roles_extracted'
write_csv(manifest,'config/database_manifest.csv',na='')
checklist<-r('config/RESKO_method_completion_checklist.csv')
checklist$project_status[checklist$step=='06']<-'DrugBank_inventory_roles_integrated_broad_STITCH_STRING_pending'
checklist$next_work[checklist$step=='06']<-'Expand STITCH chemical associations with historical protein crosswalk; acquire human STRING network'
checklist$project_status[checklist$step=='07']<-'DrugBank_human_protein_pathway_links_integrated_remaining_sources_pending'
checklist$next_work[checklist$step=='07']<-'Integrate remaining source protein IDs, construct source-aware matrices and specify enrichment universe'
write_csv(checklist,'config/RESKO_method_completion_checklist.csv')
report<-c('SCRIPT 20b - BROAD DRUGBANK PROTEIN ROLES AND HUMAN REACTOME JOINS',
 paste('Inventory entities retained:',nrow(coverage)),paste('Entities with reviewed DrugBank record links:',n_distinct(entity_map$entity_id)),
 paste('Selected DrugBank records:',length(ids)),paste('Source protein-role entries:',nrow(extracted$roles)),
 paste('Source article references preserved:',nrow(extracted$articles)),
 paste('Entities with strict human protein links:',n_distinct(strict_human$entity_id)),
 paste('Distinct strict human UniProt accessions across all roles:',n_distinct(strict_human$uniprot)),
 paste('Entities with target-role human protein links:',n_distinct(strict_human$entity_id[strict_human$protein_role=='target'])),
 paste('Entities with all-level Reactome features:',n_distinct(features$entity_id)),
 paste('Distinct all-level pathways across all roles:',n_distinct(features$pathway_id)),
 'Targets, enzymes, carriers and transporters remain separate; pharmacokinetic proteins are not automatically therapeutic targets.',
 'Human mapping requires explicit taxon 9606, consistent organism labels and an exact UniProt accession; isoform suffixes are not stripped.',
 'Candidate current-key links, named-biologic context and identity/missing-coverage holds stay separate; sacrosidase wrong structure excluded.',
 'Article references are source citations, not independently adjudicated affinity measurements. Known-action/action annotations do not prove selective binding.',
 'Lowest-level and ancestor-expanded pathways stay separate. Pathway membership does not show drug inhibition or allostery.',
 'No inventory filtering by affinity, no enrichment p-values and no BH correction at this acquisition/integration stage.',
 'Outputs: source-rich protein/reference/pathway joins locally; compact counts, source hashes and validation in outputs/20b_*.',
 'Next: broad STITCH chemical/protein links and historical protein crosswalks, followed by human STRING and GO/DO integration.')
writeLines(report,'outputs/20b_integration_summary.txt')
writeLines(c('PASS: selected-record streaming fixture crosses chunk boundaries, skips unselected records and drops unwanted fields.',
 'PASS: requested records all present, unique role/protein row IDs, source-reference and human pathway joins validated.',
 'PASS: all 200 identities/tiers retained; original protected input hashes unchanged; source size/mtime unchanged.',
 'PASS: source roles, identity contexts, taxon restrictions and exact accession matching retained; no affinity exclusion.'),
 'outputs/20b_integration_validation.txt')
notes<-readLines('Notes.txt',warn=FALSE);start<-grep('^CURRENT CHECKPOINT',notes)[1];end<-grep('^PROJECT AIM',notes)[1]
stopifnot(start<end);writeLines(c(notes[seq_len(start-1)],'CURRENT CHECKPOINT - SCRIPT 20b BROAD DRUGBANK/HUMAN PATHWAY INTEGRATION COMPLETE',report,'',notes[end:length(notes)]),'Notes.txt')
cat(paste0('\nSCRIPT 20b RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/20b_integrate_drugbank_proteins_pathways.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
