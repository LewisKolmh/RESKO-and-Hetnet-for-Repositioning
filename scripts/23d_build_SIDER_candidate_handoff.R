# Source-aware robustness and known-treatment handoff; no integrated RESKO ranking.
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr)})
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
inputs <- c('data/curated/RESKO23_candidate_identity_registry.csv',
 'data/curated/RESKO23_known_seed_treatment_evidence.csv',
 'data/curated/RESKO23_DrugBank_indication_context.csv',
 'data/curated/RESKO_mechanism_SIDER_primary_group_review_queue.csv',
 'data/curated/RESKO_mechanism_SIDER_exploratory_group_review_queue.csv',
 'data/curated/RESKO_mechanism_SIDER_all_scores.rds',
 'outputs/23b_mechanism_SIDER_cohort_summary.csv',
 'config/RESKO_mechanism_cohort_membership.csv',
 'config/RESKO23_supplemental_known_treatments.csv')
stopifnot(all(file.exists(inputs))); hashes <- tools::md5sum(inputs)
collapse <- function(x) paste(sort(unique(na.omit(x[nzchar(x)]))),collapse='|')
id <- r(inputs[1]); known <- r(inputs[2]); indication <- r(inputs[3])
sc <- readRDS(inputs[6]); cohort <- read_csv(inputs[7],show_col_types=FALSE)
member <- r(inputs[8]); supplemental <- r(inputs[9]); q <- bind_rows(r(inputs[4]),r(inputs[5])) |> distinct(flat,cohort_id)
# Count independent settings transparently; minimum-term sensitivity is not replication.
robust <- sc |> semi_join(q,by=c('flat','cohort_id')) |>
 filter(primary_gate | exploratory_covered_subset_gate) |>
 group_by(flat,cohort_id) |> summarise(historical_name=first(historical_name),
  full_coverage_cohort=all(primary_gate),
  settings_evaluated=n(),pruning_interpretations_evaluated=n_distinct(branch),
  settings_pass_ge_50=sum(pass_ge_50),settings_pass_gt_50=sum(pass_gt_50),
  passes_all_pruning_and_seed_boundaries=all(pass_ge_50)&n_distinct(branch)==2L&n()==4L,
  passes_all_strict_boundaries=all(pass_gt_50)&n_distinct(branch)==2L&n()==4L,
  min_pattern_coverage=min(coverage_fraction),max_pattern_coverage=max(coverage_fraction),
  joint_PT_min=min(joint_PT_count),joint_PT_max=max(joint_PT_count),
  .groups='drop') |>
 left_join(member |> distinct(cohort_id,cancer_id,mechanism_class),by='cohort_id')
# Retain every lead in registry. This summary holds unresolved multi-structure groups.
groups <- id |> group_by(flat) |> summarise(
 historical_CID_leads=n_distinct(lead_id),
 current_structure_keys=n_distinct(current_inchikey[!is.na(current_inchikey)]),
 current_inchikeys=collapse(current_inchikey),
 current_DrugBank_ids=collapse(exact_key_DrugBank_ids),
 current_DrugBank_names=collapse(exact_key_DrugBank_names),
 existing_inventory_entities=collapse(existing_inventory_entities),
 all_aliases_supported=all(historical_name_alias_supported=='TRUE'),
 all_CIDs_unchanged=all(requested_returned_CID_equal=='TRUE'),
 current_structure_available=all(structure_key_available=='TRUE'),
 any_exact_DrugBank_match=any(as.integer(exact_key_DrugBank_records)>0L),
 identity_dispositions=collapse(identity_disposition),
 .groups='drop') |>
 mutate(identity_scope=case_when(current_structure_keys>1L~'multiple_current_forms_do_not_pool',
 !current_structure_available~'structure_missing_hold',!all_CIDs_unchanged~'CID_redirect_hold',
 !all_aliases_supported~'historical_name_support_hold',
 TRUE~'single_current_representation_historical_formulation_unverified'))
write_csv(groups,'data/curated/RESKO23_group_identity_summary.csv')
known_context <- known |> group_by(flat,cancer_id) |> summarise(
 known_seed_tiers=collapse(decision),known_seed_DrugBank_ids=collapse(record_drugbank_id),
 known_seed_evidence_urls=collapse(evidence_url),.groups='drop')
# Historical/name-only leads are annotations, never promoted to structural equivalence.
# Do not claim absence of a clinical indication from absence in the seed review.
indication_context <- indication |> group_by(flat) |> summarise(
 DrugBank_indication_available=any(!is.na(indication)&nzchar(indication)),.groups='drop')
handoff <- robust |> left_join(groups,by='flat') |>
 left_join(known_context,by=c('flat','cancer_id')) |>
 left_join(indication_context,by='flat') |>
 left_join(supplemental |> select(flat,cancer_id,supplemental_known_status=decision,
 supplemental_known_evidence_url=evidence_url),by=c('flat','cancer_id')) |>
 mutate(known_treatment_status=case_when(
 !is.na(supplemental_known_status)~'known_labelled_treatment_supplemental_source',
 grepl('include_primary',known_seed_tiers)~'known_primary_current_structure_link',
 grepl('include_documented',known_seed_tiers)~'known_documented_use_current_structure_link',
 grepl('include_investigational',known_seed_tiers)~'known_investigational_current_structure_link',
 TRUE~'not_in_reviewed_seed_tiers_not_proof_of_novelty'),
 evidence_lane=if_else(full_coverage_cohort,'fully_covered_declared_mechanism_cohort',
 'exploratory_incomplete_mechanism_cohort'),
 next_action=case_when(identity_scope!='single_current_representation_historical_formulation_unverified'~
 'retain_form_specific_identity_hold',
 known_treatment_status!='not_in_reviewed_seed_tiers_not_proof_of_novelty'~
 'retain_known_treatment_comparator_not_novel_repurposing',
 !full_coverage_cohort~'exploratory_context_only',
 !passes_all_strict_boundaries~'threshold_sensitive_lead_preserve_for_sensitivity',
 TRUE~'advance_current_representation_to_chemistry_annotation'),
 binding_evidence_required=FALSE,affinity_missing_excludes_candidate=FALSE,
 historical_formulation_equivalence_proven=FALSE,
 interpretation='side-effect candidate context; not human eEF1A binding, allostery, RT-PPI or efficacy evidence')
write_csv(handoff,'data/curated/RESKO23_source_aware_candidate_handoff.csv')
write_csv(handoff |> count(evidence_lane,next_action),
 'outputs/23d_candidate_disposition_summary.csv')
write_csv(handoff |> group_by(evidence_lane) |> summarise(
 unique_historical_groups=n_distinct(flat),candidate_cohort_pairs=n(),
 strict_all_settings_groups=n_distinct(flat[passes_all_strict_boundaries]),.groups='drop'),
 'outputs/23d_candidate_robustness_summary.csv')
# Export current structures only as proposals. Never mutate the 200-entry inventory.
proposal <- handoff |> filter(next_action=='advance_current_representation_to_chemistry_annotation') |>
 select(flat,historical_name,current_inchikeys,current_DrugBank_ids,existing_inventory_entities) |>
 distinct() |> mutate(proposal_only=TRUE,identity_scope='current_representation_historical_form_unverified')
write_csv(proposal,'data/curated/RESKO23_chemistry_annotation_proposals.csv')
write_csv(handoff |> count(cancer_id,mechanism_class,evidence_lane,next_action),
 'outputs/23d_cohort_disposition_summary.csv')
stopifnot(!anyDuplicated(handoff[c('flat','cohort_id')]),
 setequal(paste(handoff$flat,handoff$cohort_id),paste(q$flat,q$cohort_id)),
 identical(hashes,tools::md5sum(inputs)))
write_csv(tibble(path=inputs,md5=unname(hashes)),'outputs/23d_input_hashes_summary.csv')
print(handoff |> count(evidence_lane,next_action),width=Inf)
message('23d complete. Novelty is unproven; current structures are proposals, not inventory additions.')
