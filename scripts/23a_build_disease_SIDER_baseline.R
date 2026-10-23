# Disease-forward SIDER-only baseline. Run from project root.
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr)})
r <- function(p) read_csv(p, col_types=cols(.default=col_character()))
inputs <- c('data/raw/sider/meddra_all_se.tsv.gz',
 'data/curated/cancer_primary_seeds_expanded.csv',
 'data/curated/seed_sider_mapping_accepted.csv', 'config/cancers.csv')
stopifnot(all(file.exists(inputs)))
hashes <- tools::md5sum(inputs)
# Declare before calculating results; no outcome-dependent fallback.
settings <- tibble(parameter=c('source','cohort','pruning_branches','top_decile_ties',
 'minimum_seed_terms','candidate_boundary','coverage_denominator','missing_seed_policy',
 'fallback','small_pattern_policy','identity_unit'),
 value=c('SIDER 4.1 PT only','reviewed primary conventional seeds, separately per cancer',
 'top_decile_types|prevalence_gt_0.10','remove all tied at ceiling(0.10 * number of PT types)',
 '>=3 and >3 after pruning','>=0.50 and >0.50',
 'distinct joint PT count; overlap / joint PT count',
 'report missing; covered-subset adaptation does not establish full-cohort pattern',
 'none; no seed-halving selected after results',
 'score nonempty patterns; flag joint count <5 as insufficient for Figure 3',
 'historical SIDER flat group, not an exact modern structure'))
write_csv(settings,'config/RESKO_SIDER_baseline_settings.csv')
message('Reading SIDER PT associations...')
raw <- read_tsv(inputs[1],col_names=c('stitch_flat_id','stereo','label_cui','type','pt_cui','term'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
profiles <- raw |> filter(type=='PT',!is.na(pt_cui),nzchar(pt_cui)) |>
 distinct(stitch_flat_id,pt_cui)
N <- n_distinct(profiles$stitch_flat_id)
freq <- profiles |> count(pt_cui,name='profiles_with_term') |>
 mutate(profile_fraction=profiles_with_term/N)
cut <- sort(freq$profiles_with_term,decreasing=TRUE)[ceiling(nrow(freq)*0.10)]
freq <- freq |> mutate(remove_top_decile_types=profiles_with_term>=cut,
 remove_prevalence_gt_0.10=profile_fraction>0.10)
labels <- raw |> filter(type=='PT') |> distinct(pt_cui,term)
seeds <- r(inputs[2]) |> semi_join(r(inputs[4]) |> filter(include=='TRUE'),by='cancer_id') |>
 distinct(cancer_id,cancer_name,drugbank_id,drug_name)
map <- r(inputs[3]) |> distinct(drugbank_id,stitch_flat_id,mapping_status,evidence_note)
stopifnot(!anyDuplicated(seeds[c('cancer_id','drugbank_id')]))
# Preserve distinct mapped flat groups but count each seed/PT only once.
seed_map <- seeds |> left_join(map,by='drugbank_id',relationship='many-to-many')
coverage <- list(); summaries <- list(); patterns <- list(); scores <- list()
for(branch in c('top_decile_types','prevalence_gt_0.10')) {
 removed <- freq$pt_cui[freq[[paste0('remove_',branch)]]]
 kept <- profiles |> filter(!pt_cui %in% removed)
 seed_pt <- seed_map |> filter(!is.na(stitch_flat_id)) |>
 inner_join(kept,by='stitch_flat_id',relationship='many-to-many') |>
 distinct(cancer_id,drugbank_id,pt_cui)
 counts <- seed_pt |> count(cancer_id,drugbank_id,name='retained_PT_count')
 for(minimum in c(3L,4L)) {
  cv <- seeds |> left_join(counts,by=c('cancer_id','drugbank_id')) |>
   mutate(retained_PT_count=replace_na(retained_PT_count,0L),
    mapping_present=drugbank_id %in% map$drugbank_id,
    eligible=mapping_present & retained_PT_count>=minimum,
    status=case_when(!mapping_present~'no_accepted_SIDER_mapping',
      !eligible~'too_few_retained_PT',TRUE~'eligible_covered_seed'),
    branch=branch,minimum_seed_terms=minimum)
  coverage[[length(coverage)+1L]] <- cv
  for(cancer in unique(seeds$cancer_id)) {
   cc <- cv |> filter(cancer_id==cancer)
   ids <- cc$drugbank_id[cc$eligible]
   joint <- seed_pt |> filter(cancer_id==cancer,drugbank_id %in% ids) |>
    count(pt_cui) |> filter(n==length(ids),length(ids)>0) |> pull(pt_cui)
   J <- length(joint)
   full <- all(cc$eligible)
   context <- tibble(cancer_id=cancer,branch=branch,minimum_seed_terms=minimum,
    total_primary_seeds=nrow(cc),eligible_seeds=length(ids),
    missing_mapping=sum(!cc$mapping_present),joint_PT_count=J,
    full_cohort_assessable=full,figure3_minimum_met=J>=5L,
    scope=if(full)'full_primary_cohort' else 'eligible_covered_subset_adaptation')
   patterns[[length(patterns)+1L]] <- tibble(cancer_id=cancer,branch=branch,
    minimum_seed_terms=minimum,pt_cui=joint)
   known_groups <- seed_map$stitch_flat_id[seed_map$cancer_id==cancer]
   hits <- kept |> filter(pt_cui %in% joint) |> count(stitch_flat_id,name='shared_PT')
   sc <- tibble(stitch_flat_id=sort(unique(profiles$stitch_flat_id))) |>
    left_join(hits,by='stitch_flat_id') |>
    mutate(shared_PT=replace_na(shared_PT,0L),joint_PT_count=J,
     coverage_fraction=if(J>0)shared_PT/J else NA_real_,
     pass_ge_50=J>0 & 2L*shared_PT>=J,pass_gt_50=J>0 & 2L*shared_PT>J,
     known_primary_group=stitch_flat_id %in% known_groups,
     known_identity_conflict=stitch_flat_id=='CID100135113',
     cancer_id=cancer,branch=branch,minimum_seed_terms=minimum,
     full_cohort_assessable=full,figure3_minimum_met=J>=5L,
     identity_scope='historical_SIDER_flat_group_requires_identity_review')
   scores[[length(scores)+1L]] <- sc
   summaries[[length(summaries)+1L]] <- context |>
    mutate(new_groups_ge_50=sum(sc$pass_ge_50 & !sc$known_primary_group & !sc$known_identity_conflict),
     new_groups_gt_50=sum(sc$pass_gt_50 & !sc$known_primary_group & !sc$known_identity_conflict))
  }
 }
}
cv <- bind_rows(coverage); sm <- bind_rows(summaries); sc <- bind_rows(scores)
pat <- bind_rows(patterns) |> left_join(labels,by='pt_cui',relationship='many-to-many')
write_csv(freq,'data/curated/RESKO_SIDER_global_PT_prevalence.csv')
write_csv(seed_map,'data/curated/RESKO_SIDER_primary_seed_mapping.csv')
write_csv(cv,'data/curated/RESKO_SIDER_disease_seed_coverage.csv')
write_csv(pat,'data/curated/RESKO_SIDER_disease_joint_patterns.csv')
write_csv(sc,'data/curated/RESKO_SIDER_disease_candidate_scores.csv')
write_csv(sc |> filter(pass_ge_50,!known_primary_group,!known_identity_conflict),
 'data/curated/RESKO_SIDER_new_candidate_group_review_queue.csv')
write_csv(sm,'outputs/23a_SIDER_disease_baseline_summary.csv')
write_csv(cv |> count(cancer_id,branch,minimum_seed_terms,status),
 'outputs/23a_SIDER_seed_coverage_summary.csv')
write_csv(tibble(source_profiles=N,PT_types=nrow(freq),
 top_decile_rank=ceiling(nrow(freq)*0.10),tie_count_cutoff=cut,
 top_decile_removed=sum(freq$remove_top_decile_types),
 prevalence_gt_10pct_removed=sum(freq$remove_prevalence_gt_0.10)),
 'outputs/23a_SIDER_pruning_summary.csv')
stopifnot(identical(hashes,tools::md5sum(inputs)))
write_csv(tibble(path=inputs,md5=unname(hashes)),'outputs/23a_input_hashes_summary.csv')
print(sm,width=Inf)
message('23a complete. Scores are pattern coverage, not p-values or evidence of eEF1A binding.')
