# Mandatory response to a source-discovered missing known treatment; not hit optimisation.
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr)})
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
inputs <- c('config/RESKO23_supplemental_known_treatments.csv',
 'config/RESKO_mechanism_cohort_membership.csv',
 'data/curated/seed_sider_mapping_accepted.csv','data/curated/seed_sider_mapping_all.csv',
 'data/curated/RESKO23_candidate_identity_registry.csv','data/raw/sider/meddra_all_se.tsv.gz',
 'data/curated/RESKO23_source_aware_candidate_handoff.csv',
 'data/raw/section23_identity/thalidomide_label.rds')
stopifnot(all(file.exists(inputs)));hashes <- tools::md5sum(inputs)
supplement <- r(inputs[1]); member <- r(inputs[2]); map <- r(inputs[3]); all_map <- r(inputs[4])
id <- r(inputs[5]); raw <- read_tsv(inputs[6],
 col_names=c('flat','stereo','label','type','pt','term'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
profile <- raw |> filter(type=='PT',!is.na(pt),nzchar(pt)) |> distinct(flat,pt)
f <- profile |> count(pt,name='n');cut <- sort(f$n,decreasing=TRUE)[ceiling(nrow(f)*.1)]
# Identity support is required before using an additional known drug's flat profile.
extra <- id |> filter(flat==supplement$flat[1],historical_name_alias_supported=='TRUE',
 requested_returned_CID_equal=='TRUE',exact_key_DrugBank_ids==supplement$drugbank_id[1])
stopifnot(nrow(extra)==1L)
label <- readRDS(inputs[8]);stopifnot(label$status==200L,
 grepl('newly diagnosed multiple myeloma',label$body,fixed=TRUE))
cohort <- 'multiple_myeloma__cereblon_IMiD'
original <- member |> filter(cohort_id==cohort) |> select(drugbank_id,drug_name)
augmented <- bind_rows(original,supplement |> select(drugbank_id,drug_name))
aug_map <- bind_rows(map |> filter(drugbank_id %in% original$drugbank_id) |> select(drugbank_id,flat=stitch_flat_id),
 supplement |> select(drugbank_id,flat)) |> distinct()
write_csv(augmented |> mutate(cohort_id=cohort,source='23e_label_corrected_sensitivity_original_membership_preserved'),
 'data/curated/RESKO23_augmented_IMiD_membership.csv')
known_ids <- member$drugbank_id[member$cancer_id=='multiple_myeloma']
exclude <- union(all_map$stitch_flat_id[all_map$drugbank_id %in% known_ids],supplement$flat)
summary <- patterns <- scores <- list()
for(branch in c('top_decile_types','prevalence_gt_0.10')) {
 removed <- if(branch=='top_decile_types') f$pt[f$n>=cut] else f$pt[f$n/n_distinct(profile$flat)>.1]
 kept <- profile |> filter(!pt %in% removed)
 for(minimum in c(3L,4L)) {
  sets <- lapply(augmented$drugbank_id,function(id) kept$pt[kept$flat %in% aug_map$flat[aug_map$drugbank_id==id]])
  eligible <- lengths(sets)>=minimum
  joint <- if(all(eligible)) Reduce(intersect,sets) else character()
  J <- length(joint);gate <- all(eligible)&n_distinct(aug_map$flat)>=2L&J>=5L
  counts <- kept |> filter(pt %in% joint) |> count(flat,name='shared_PT')
  sc <- tibble(flat=sort(unique(profile$flat))) |> left_join(counts,by='flat') |>
   mutate(shared_PT=replace_na(shared_PT,0L),joint_PT_count=J,
    coverage_fraction=if(J>0)shared_PT/J else NA_real_,
    pass_ge_50=J>0 & 2L*shared_PT>=J,pass_gt_50=J>0 & 2L*shared_PT>J,
    known_or_held_treatment=flat %in% exclude,identity_conflict=flat=='CID100135113',
    branch=branch,minimum_seed_terms=minimum,augmented_cohort_gate=gate)
  scores[[length(scores)+1L]] <- sc
  patterns[[length(patterns)+1L]] <- tibble(branch=branch,minimum_seed_terms=minimum,pt=joint)
  summary[[length(summary)+1L]] <- tibble(branch=branch,minimum_seed_terms=minimum,
   intended_seeds=nrow(augmented),eligible_seeds=sum(eligible),joint_PT_count=J,
   augmented_cohort_gate=gate,
   nonknown_ge_50=sum(sc$pass_ge_50 & !sc$known_or_held_treatment & !sc$identity_conflict),
   nonknown_gt_50=sum(sc$pass_gt_50 & !sc$known_or_held_treatment & !sc$identity_conflict))
 }
}
sm <- bind_rows(summary); sc <- bind_rows(scores)
write_csv(sm,'outputs/23e_augmented_IMiD_summary.csv')
write_csv(bind_rows(patterns),'data/curated/RESKO23_augmented_IMiD_joint_patterns.csv')
saveRDS(sc,'data/curated/RESKO23_augmented_IMiD_scores.rds')
write_csv(sc |> filter(augmented_cohort_gate,pass_ge_50,!known_or_held_treatment,!identity_conflict),
 'data/curated/RESKO23_augmented_IMiD_group_review_queue.csv')
# Original handoff is immutable here; annotate explicitly in a separate final handoff.
h <- read_csv(inputs[7],show_col_types=FALSE) |>
 mutate(seed_completeness_note=if_else(cohort_id==cohort,
 'original_two_seed_IMiD_cohort_omits_label_verified_thalidomide_see_augmented_sensitivity',
 'coverage_of_declared_seeds_not_exhaustive_clinical_cohort'),
 final_section23_action=if_else(cohort_id==cohort &
 next_action!='retain_known_treatment_comparator_not_novel_repurposing',
 'retain_original_IMiD_sensitivity_only_use_augmented_analysis_for_followup',next_action))
write_csv(h,'data/curated/RESKO23_final_candidate_handoff.csv')
write_csv(h |> count(evidence_lane,final_section23_action),'outputs/23e_final_disposition_summary.csv')
stopifnot(identical(hashes,tools::md5sum(inputs)))
write_csv(tibble(path=inputs,md5=unname(hashes)),'outputs/23e_input_hashes_summary.csv')
print(sm,width=Inf)
message('23e complete. Missing known treatment recorded; original cohort not silently rewritten.')
