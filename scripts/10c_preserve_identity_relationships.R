# Preserve source-versus-resolved identity relationships without merging variants.
library(readr); library(dplyr); library(tibble)
if (!all(dir.exists(c("scripts", "config", "data")))) stop('Run from the project root.')
r <- function(p) read_csv(p, col_types=cols(.default=col_character()))
x <- r('data/curated/seed_identity_crosscheck.csv')
pc <- r('data/curated/pubchem_seed_identity_records.csv')
local_ids <- r('data/curated/seed_equivalent_record_identifiers.csv')
source_comparison <- x |>
  select(drugbank_id,drug_name,inchikey,pubchem_cid,resolved_pubchem_cid) |>
  left_join(pc,by=c('pubchem_cid'='cid')) |>
  mutate(source_cid_structure_status=case_when(
    is.na(pubchem_cid) ~ 'no_source_CID',
    is.na(pc_inchikey) ~ 'source_CID_not_returned',
    is.na(inchikey) ~ 'no_DrugBank_key_for_comparison',
    inchikey==pc_inchikey ~ 'full_key_match',
    substr(inchikey,1,14)==substr(pc_inchikey,1,14) ~ 'connectivity_only_DO_NOT_merge',
    TRUE ~ 'different_representation_or_structure_DO_NOT_merge'
  ))
roles <- bind_rows(
  x |> transmute(drugbank_id,drug_name,source_inchikey=inchikey,cid=pubchem_cid,record_role='source_CID_lead'),
  x |> transmute(drugbank_id,drug_name,source_inchikey=inchikey,cid=resolved_pubchem_cid,record_role='resolved_full_key_CID')
) |> filter(!is.na(cid)) |> left_join(pc,by='cid') |>
  mutate(full_key_match=!is.na(source_inchikey)&!is.na(pc_inchikey)&source_inchikey==pc_inchikey,
    interpretation=if_else(full_key_match,'full_structure_key_agreement',
      'retained_for_identity_review_NOT_equivalent'))
# One-to-many identifiers remain long, with resource namespaces and record provenance.
identifiers <- bind_rows(
  local_ids |> transmute(seed_drugbank_id,resource,identifier,source='DrugBank supplied XML',
    source_record=record_drugbank_id,relationship),
  x |> filter(chembl_key_status=='full_key_match') |>
    transmute(seed_drugbank_id=drugbank_id,resource='ChEMBL',identifier=chembl_id,
      source='ChEMBL API',source_record=chembl_id,relationship='external_full_key_match'),
  x |> filter(pubchem_key_status=='full_key_match') |>
    transmute(seed_drugbank_id=drugbank_id,resource='PubChem Compound',identifier=resolved_pubchem_cid,
      source='PubChem API',source_record=resolved_pubchem_cid,relationship='external_full_key_match')
) |> distinct()
stopifnot(nrow(x)==170,!anyDuplicated(x$drugbank_id),nrow(source_comparison)==170,
  nrow(r('data/curated/cancer_drug_review_all.csv'))==238)
write_csv(source_comparison,'outputs/10c_source_CID_structure_comparison.csv',na='')
write_csv(roles,'data/curated/seed_pubchem_record_relationships.csv',na='')
write_csv(identifiers,'data/curated/seed_external_identifiers_enriched.csv',na='')
write_csv(source_comparison |> filter(source_cid_structure_status %in%
  c('connectivity_only_DO_NOT_merge','different_representation_or_structure_DO_NOT_merge')),
  'outputs/10c_source_CID_disagreements.csv',na='')
report <- c('SCRIPT 10c - RETAIN SOURCE AND RESOLVED STRUCTURE RELATIONSHIPS',
  capture.output(print(count(source_comparison,source_cid_structure_status))),
  paste('Provenance-bearing identifier rows:',nrow(identifiers)),
  'Different source/resolved CIDs are NOT assumed to be retired-ID redirects.',
  'Source records include stereochemical, protonation/coordination or other differences.',
  'These require chemical review; related forms are retained separately, not pooled.',
  'Cisplatin, motexafin lutetium and DOTMP HO-166 have source keys without unique',
  'exact-key PubChem resolution. The latter two are outside the primary queue.',
  'Carboplatin has an exact-key PubChem record but no ChEMBL structure returned.',
  'No missing structures or unverified forms are counted as negative binding results.',
  'Outputs: seed_pubchem_record_relationships.csv; seed_external_identifiers_enriched.csv',
  '(data/curated/); source CID comparison and disagreement reports (outputs/10c_*).',
  'Result meaning: names, cross-references and alternative source representations are',
  'available for future searches, with entity equivalence kept separate from relatedness.',
  'No biological activity, approval or eEF1A selectivity is inferred from these mappings.')
writeLines(report,'outputs/10c_identity_relationship_summary.txt')
cat(paste0('\nSCRIPT 10c RUN - ',format(Sys.time(),tz='Europe/Rome',usetz=TRUE),'\n',
  'Script: scripts/10c_preserve_identity_relationships.R\n',
  'Inputs: cached 10a/10b curated identity and cross-reference outputs. No new API calls.\n',
  paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
