# Bounded, explicit assay review; source endpoints remain unchanged.
library(readr)
library(dplyr)
library(tibble)
r <- function(p) read_csv(p, col_types = cols(.default = col_character()))
inputs <- c('data/curated/eef1a_chembl_enriched_activity_evidence.csv',
 'data/curated/eef1a_chembl_document_registry.csv',
 'data/curated/target_primary_literature_registry.csv')
before <- tools::md5sum(inputs)
a <- r(inputs[1]); docs <- r(inputs[2]); lit <- r(inputs[3])
# Decisions apply to exact assays, not keyword hits or ChEMBL assay-type labels.
decisions <- tribble(
 ~assay_chembl_id, ~review_layer, ~review_reason,
 'CHEMBL1228364','functional_partner_readout','tRNA ternary-complex readout does not measure compound binding to eEF1A',
 'CHEMBL1228365','translation_function','Peptidyl-transfer readout does not measure compound binding to eEF1A',
 'CHEMBL1805138','affinity_capture','Biotinylated-probe pull-down; direct affinity and probe/parent identity need full methods',
 'CHEMBL1805291','reported_SPR_pending_methods','Reported EEF1A1 SPR; full primary methods and probe/parent identities not verified',
 'CHEMBL1805281','reported_SPR_pending_methods','Reported EEF1A2 SPR; full primary methods and probe/parent identities not verified',
 'CHEMBL1805284','translation_function','Cellular methionine incorporation is not physical binding',
 'CHEMBL1805285','translation_function','Cellular methionine incorporation is not physical binding',
 'CHEMBL1805286','isoform_conflict_hold','EEF1A2 target annotation conflicts with EEF1A1 assay description',
 'CHEMBL1805287','translation_function','Cellular methionine incorporation is not physical binding',
 'CHEMBL1805288','downstream_signalling','Akt phosphorylation is not physical binding',
 'CHEMBL1805289','downstream_signalling','Akt phosphorylation is not physical binding',
 'CHEMBL3803738','chemoproteomic_enrichment','C646-yne enrichment; reactive-probe and abundant-protein liabilities, no affinity measurement',
 'CHEMBL4024027','source_integrity_hold','Expression of Concern PMID 31496240 linked to source PMID 29498519',
 'CHEMBL4024030','source_integrity_hold','Expression of Concern PMID 31496240 linked to source PMID 29498519',
 'CHEMBL4024033','source_integrity_hold','ITC does not override the linked Expression of Concern; quarantine derived activities',
 'CHEMBL4270827','chemoproteomic_enrichment','Gymnemic acid I ribosome-complex proteomics; enrichment is not direct eEF1A affinity',
 'CHEMBL5697473','competition_profile_not_positive','Molibresib competition profile; near-unit fold changes do not establish positive binding',
 'CHEMBL5698173','competition_profile_not_positive','Molibresib competition profile; preserve source concentrations and ratios',
 'CHEMBL5698873','competition_profile_not_positive','Molibresib competition profile; preserve source concentrations and ratios',
 'CHEMBL5699573','competition_profile_not_positive','Molibresib competition profile; preserve source concentrations and ratios',
 'CHEMBL5700273','competition_profile_not_positive','Molibresib competition profile; preserve source concentrations and ratios',
 'CHEMBL5700973','competition_profile_not_positive','Censored IC50 >=10000 nM is not a positive inhibition hit or proof of absent binding',
 'CHEMBL5650978','kinobead_competition','Same molecule/assay/value with Kd and ED50 labels: one competition-derived observation, not purified affinity'
)
stopifnot(nrow(decisions)==23L, !anyDuplicated(decisions$assay_chembl_id),
 setequal(a$assay_chembl_id, decisions$assay_chembl_id))
write_csv(decisions, 'config/17e_assay_review_decisions.csv')
review <- a |> left_join(decisions, by='assay_chembl_id',relationship='many-to-one') |>
 left_join(docs |> select(document_chembl_id,pubmed_id,doi),by='document_chembl_id',relationship='many-to-one') |>
 left_join(lit |> select(pmid,full_text_status,source_integrity_status,correction_PMIDs,pubmed_url),
 by=c('pubmed_id'='pmid'),relationship='many-to-one') |>
 mutate(review_date='2026-10-20',review_basis='ChEMBL assay record and primary publication metadata/abstract where available; full-text access status explicit',
 binding_conclusion=case_when(review_layer=='source_integrity_hold' ~ 'quarantined_source_integrity',
 review_layer=='reported_SPR_pending_methods' ~ 'reported_physical_binding_pending_methods',
 TRUE ~ 'not_established_by_this_record'),
 selectivity_conclusion='not_established_in_this_review',allostery_conclusion='not_established_in_this_review',
 RT_PPI_conclusion='not_established_by_this_record',translation_sparing_conclusion='not_established_in_this_review',
 observation_group=if_else(assay_chembl_id=='CHEMBL5650978',
 paste(assay_chembl_id,molecule_chembl_id,standard_value,sep=':'),activity_id))
stopifnot(nrow(review)==39L,!anyDuplicated(review$activity_id),!anyNA(review$review_layer),
 sum(review$review_layer=='source_integrity_hold')==11L,
 all(review$review_layer[review$pubmed_id=='29498519']=='source_integrity_hold'),
 n_distinct(review$observation_group[review$assay_chembl_id=='CHEMBL5650978'])==1L)
write_csv(review,'data/curated/eef1a_reviewed_activity_evidence.csv',na='')
write_csv(review |> count(review_layer,name='source_records'),'outputs/17e_assay_review_summary.csv')
write_csv(review |> distinct(assay_chembl_id,review_layer,review_reason,pubmed_id,full_text_status),
 'outputs/17e_assay_coverage.csv',na='')
# Positive control stays outside the original seed/candidate inventory.
index <- r('data/curated/drugbank_structure_identity_index.csv')
pc <- index |> filter(record_drugbank_id=='DB04977')
stopifnot(nrow(pc)==1L,pc$record_inchikey=='UUSZLLQJYRSZIS-LXNNNBEUSA-N')
z <- readRDS('data/raw/target_evidence/literature/fulltext_PMC5054363.rds')
stopifnot(z$status==200L)
p <- xml2::xml_text(xml2::xml_find_all(xml2::read_xml(z$body),'//body//p'))
stopifnot(any(grepl('eEF1A2 purified from rabbit muscle',p,fixed=TRUE)))
control <- tibble(control_id='control:DB04977',compound='Plitidepsin',drugbank_id='DB04977',
 inchikey=pc$record_inchikey,identity_basis='Local DrugBank structure index; independent cross-source identity check still pending',
 role='separate_positive_binding_reference_not_inventory_seed',source_pmcid='PMC5054363',doi='10.1038/srep35100',
 assay='Radioligand saturation binding',protein='Rabbit muscle eEF1A2',reported_Kd_nM=79,reported_uncertainty_nM=30,
 uncertainty_note='Source reports 79 +/- 30 nM; retain original uncertainty definition in paper',
 human_context='Human tumour-cell target engagement and lysate evidence; not a purified-human affinity measurement',
 selectivity='A2-exclusive binding not established; both isoforms detected in binding lysate fractions',
 allostery='Domain 1/2 site proposed by modelling; experimental allostery not established',
 RT_PPI_disruption='Not established by this control paper',translation_sparing='Not established by this review',
 source_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC5054363/')
write_csv(control,'data/curated/eef1a_positive_binding_controls.csv')
write_csv(control,'outputs/17e_control_summary.csv')
# Source-level association evidence cannot be upgraded without upstream studies.
stitch <- r('data/curated/eef1a_STITCH_inventory_association_evidence.csv')
write_csv(stitch |> mutate(review_status='association_only_not_direct_binding',
 upstream_study_resolution='unresolved_detailed_export_has_no_original_study_identifiers',
 independence='Potential source overlap with ChEMBL/DrugBank; not independent confirmation'),
 'data/curated/eef1a_STITCH_reviewed_inventory_evidence.csv',na='')
stopifnot(identical(before,tools::md5sum(inputs)))
write_csv(tibble(file=inputs,md5=unname(before)),'outputs/17e_input_summary.csv')
report <- c('SCRIPT 17e - BOUNDED ASSAY ADJUDICATION AND CONTROL',
 '39 source activity records accounted for across 23 assays and 20 ligands.',
 '11 records from PMID 29498519 quarantined because of linked Expression of Concern; this is not described as a retraction.',
 '6 SPR records retained as reported binding pending full methods; no inventory identity match.',
 'One isoform contradiction held. Two kinobead endpoint labels represent one observation.',
 'Plitidepsin is a separate positive reference: reported rabbit eEF1A2 Kd 79 +/- 30 nM, not proven A2-exclusive or allosteric.',
 '13 inventory STITCH links retain association-only status and unresolved upstream study provenance.',
 'Primary full-text retrieval remains incomplete; this is a bounded review, not exhaustive adjudication.',
 'Outputs: reviewed activity evidence, STITCH review and control registry locally; compact review/control reports in outputs/17e_*.',
 'No BH adjustment: no valid statistical test family or null/background is defined here.')
writeLines(report,'outputs/17e_review_summary.txt')
writeLines('PASS: exact assay coverage, all source rows accounted for, integrity quarantine, duplicate grouping, source XML control check and protected input hashes.',
 'outputs/17e_review_validation.txt')
cat(paste0('\nSCRIPT 17e RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
