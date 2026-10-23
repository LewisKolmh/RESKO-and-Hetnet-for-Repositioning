# Final bounded SIDER handoff, source-corrected sensitivity and section audit.
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr)})
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
inputs <- c('data/curated/RESKO23_final_candidate_handoff.csv',
 'data/curated/RESKO23_augmented_IMiD_scores.rds','data/curated/RESKO23_group_identity_summary.csv',
 'data/curated/RESKO23_known_seed_treatment_evidence.csv',
 'config/RESKO23_supplemental_known_treatments.csv',
 'data/curated/RESKO23_candidate_identity_registry.csv',
 'outputs/16f_target_evidence_handoff_summary.csv')
stopifnot(all(file.exists(inputs))); hashes <- tools::md5sum(inputs)
original <- read_csv(inputs[1],show_col_types=FALSE)
sc <- readRDS(inputs[2]); groups <- read_csv(inputs[3],show_col_types=FALSE)
known <- r(inputs[4]); supplemental <- r(inputs[5]); id <- r(inputs[6])
collapse <- function(x) paste(sort(unique(na.omit(x[nzchar(x)]))),collapse='|')
q <- sc |> filter(augmented_cohort_gate,pass_ge_50,!known_or_held_treatment,!identity_conflict) |>
 distinct(flat)
aug <- sc |> semi_join(q,by='flat') |> group_by(flat) |> summarise(
 settings_evaluated=n(),settings_pass_ge_50=sum(pass_ge_50),settings_pass_gt_50=sum(pass_gt_50),
 passes_all_pruning_and_seed_boundaries=all(pass_ge_50)&n()==4L&n_distinct(branch)==2L,
 passes_all_strict_boundaries=all(pass_gt_50)&n()==4L&n_distinct(branch)==2L,
 min_pattern_coverage=min(coverage_fraction),max_pattern_coverage=max(coverage_fraction),
 joint_PT_min=min(joint_PT_count),joint_PT_max=max(joint_PT_count),.groups='drop') |>
 left_join(id |> distinct(flat,historical_name),by='flat') |>
 left_join(groups,by='flat') |>
 mutate(cohort_id='multiple_myeloma__cereblon_IMiD_augmented',cancer_id='multiple_myeloma',
 mechanism_class='cereblon_IMiD',evidence_lane='label_corrected_declared_mechanism_cohort')
k <- known |> filter(cancer_id=='multiple_myeloma') |> group_by(flat) |>
 summarise(known_seed_tiers=collapse(decision),.groups='drop')
aug <- aug |> left_join(k,by='flat') |> mutate(
 known_treatment_status=case_when(flat %in% supplemental$flat~'known_labelled_treatment_supplemental_source',
 !is.na(known_seed_tiers)~'known_reviewed_seed_use_current_structure_link',
 TRUE~'not_in_reviewed_seed_tiers_not_proof_of_novelty'),
 final_section23_action=case_when(identity_scope!='single_current_representation_historical_formulation_unverified'~
 'retain_form_specific_identity_hold',
 known_treatment_status!='not_in_reviewed_seed_tiers_not_proof_of_novelty'~
 'retain_known_treatment_comparator_not_novel_repurposing',
 !passes_all_strict_boundaries~'threshold_sensitive_lead_preserve_for_sensitivity',
 TRUE~'advance_current_representation_to_chemistry_annotation'),
 seed_completeness_note='thalidomide_added_using_external_label_not_outcome_selected_subset',
 historical_formulation_equivalence_proven=FALSE,binding_evidence_required=FALSE,
 affinity_missing_excludes_candidate=FALSE,
 interpretation='source-corrected side-effect pattern, not efficacy or human eEF1A binding evidence')
# Original IMiD entries remain visible, explicitly sensitivity-only. No favourable union hidden.
final <- bind_rows(original,aug)
write_csv(final,'data/curated/RESKO23_completed_SIDER_handoff.csv')
write_csv(final |> count(evidence_lane,final_section23_action),
 'outputs/23f_completed_SIDER_disposition_summary.csv')
write_csv(final |> group_by(evidence_lane) |> summarise(candidate_cohort_pairs=n(),
 historical_groups=n_distinct(flat),.groups='drop'),
 'outputs/23f_completed_SIDER_coverage_summary.csv')
proposal <- final |> filter(final_section23_action=='advance_current_representation_to_chemistry_annotation') |>
 select(flat,historical_name,current_inchikeys,current_DrugBank_ids,existing_inventory_entities) |>
 distinct() |> mutate(proposal_only=TRUE,historical_formulation_equivalence_proven=FALSE,
 source_stage='23f; current representations only; prior 200-entry inventory not mutated')
write_csv(proposal,'data/curated/RESKO23_completed_chemistry_proposals.csv')
# Compact approved candidate names/counts, no raw side-effect or licensed indication text.
write_csv(proposal |> select(historical_name,flat),
 'outputs/23f_chemistry_proposal_summary.csv')
# Validate original handoff setting summaries independently, using all source score rows.
base_scores <- readRDS('data/curated/RESKO_mechanism_SIDER_all_scores.rds')
for(i in seq_len(nrow(original))) {
 z <- original[i,]; rows <- base_scores[base_scores$flat==z$flat & base_scores$cohort_id==z$cohort_id &
 (base_scores$primary_gate | base_scores$exploratory_covered_subset_gate),]
 stopifnot(nrow(rows)==z$settings_evaluated,sum(rows$pass_ge_50)==z$settings_pass_ge_50,
 sum(rows$pass_gt_50)==z$settings_pass_gt_50,
 z$passes_all_strict_boundaries==(all(rows$pass_gt_50)&nrow(rows)==4L&length(unique(rows$branch))==2L))
}
# Independently reconstruct the corrected three-drug intersections from raw PT source.
raw <- read_tsv('data/raw/sider/meddra_all_se.tsv.gz',
 col_names=c('flat','stereo','label','type','pt','term'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
x <- unique(as.data.frame(raw[raw$type=='PT',c('flat','pt')]))
f <- table(x$pt);N <- length(unique(x$flat))
member <- r('data/curated/RESKO23_augmented_IMiD_membership.csv')
mp <- r('data/curated/seed_sider_mapping_accepted.csv')
mp <- rbind(as.data.frame(mp[c('drugbank_id','stitch_flat_id')]),
 data.frame(drugbank_id=supplemental$drugbank_id,stitch_flat_id=supplemental$flat))
sm <- read_csv('outputs/23e_augmented_IMiD_summary.csv',show_col_types=FALSE)
for(i in seq_len(nrow(sm))) {
 z <- sm[i,]; removed <- if(z$branch=='top_decile_types') names(f)[f>=sort(f,decreasing=TRUE)[ceiling(length(f)*.1)]] else names(f)[f/N>.1]
 keep <- x[!x$pt %in% removed,];sets <- lapply(member$drugbank_id,function(d)
 unique(keep$pt[keep$flat %in% mp$stitch_flat_id[mp$drugbank_id==d]]))
 stopifnot(all(lengths(sets)>=z$minimum_seed_terms))
 joint <- Reduce(intersect,sets);ss <- sc[sc$branch==z$branch & sc$minimum_seed_terms==z$minimum_seed_terms,]
 expected <- vapply(ss$flat,function(g)sum(unique(keep$pt[keep$flat==g]) %in% joint),integer(1))
 stopifnot(length(joint)==z$joint_PT_count,all(ss$shared_PT==expected),
 all(ss$pass_ge_50==(2L*expected>=length(joint))),all(ss$pass_gt_50==(2L*expected>length(joint))))
}
stopifnot(!anyDuplicated(final[c('flat','cohort_id')]),
 all(proposal$historical_formulation_equivalence_proven==FALSE),nrow(r(inputs[7]))==200L,
 !any(final$binding_evidence_required),!any(final$affinity_missing_excludes_candidate),
 !any(proposal$flat=='CID100135113'),
 !any(final$flat=='CID100005426' & final$final_section23_action=='advance_current_representation_to_chemistry_annotation'),
 identical(hashes,tools::md5sum(inputs)))
for(stage in c('23a','23b','23c','23d','23e')) {
 h <- r(paste0('outputs/',stage,'_input_hashes_summary.csv'))
 stopifnot(all(unname(tools::md5sum(h$path))==h$md5))
}
invisible(lapply(list.files('scripts',pattern='\\.R$',recursive=TRUE,full.names=TRUE),parse))
write_csv(tibble(path=inputs,md5=unname(hashes)),'outputs/23f_input_hashes_summary.csv')
report <- c('SECTION 23 COMPLETE: bounded SIDER analysis and source-aware handoff, not full RESKO completion.',
 '23a broad cancer baselines: coverage limited, no five-term pattern.',
 '23b all 75 original primary drugs / 85 associations retained in 53 predeclared mechanism cohorts.',
 paste('23c identity dispositions cover',n_distinct(id$flat),'historical groups /',nrow(id),'stereo-CID leads.'),
 '23d setting robustness, known-treatment comparators and historical form holds retained.',
 '23e label-verified thalidomide seed gap: added sensitivity gives 6 / 12 joint PT; original IMiD cohort sensitivity-only.',
 paste('23f completed handoff candidate/cohort contexts:',nrow(final)),
 paste('Distinct current-representation chemistry proposals:',nrow(proposal)),
 'PASS: all original setting summaries and corrected IMiD raw-source overlaps independently verified.',
 'PASS: all protected stage hashes unchanged; all R scripts parse; original inventory remains 200.',
 'Unresolved historical forms are explicit holds. Known seed gaps do not establish exhaustive clinical cohort coverage.',
 'No p-values/BH, affinity filter, clinical recommendation, target-binding or allostery claim.',
 'NEXT: paper-matched chemical fingerprints and similarity matrices using source-aware proposals; then ontology/network integration and valid enrichment families.')
writeLines(report,'outputs/23f_section_completion_summary.txt')
writeLines(report[grepl('PASS:',report)],'outputs/23f_section_validation.txt')
cat(report,sep='\n')
