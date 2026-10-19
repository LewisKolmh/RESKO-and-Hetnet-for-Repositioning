# Reassessment with reviewed context. Descriptive set overlap is not a validated prediction.
library(readr);library(dplyr);library(tidyr);library(stringr);library(tibble)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
paths<-c('outputs/16c_pair_compatibility_summary.csv','outputs/16c_profile_compatibility_summary.csv',
 'data/curated/reviewed_candidate_full_context.csv','data/curated/profile_comparison_source_context.csv',
 'outputs/16a_candidate_PT_features.csv','outputs/16b_seed_PT_features.csv','data/curated/side_effect_term_prevalence.csv')
before<-tools::md5sum(paths);pairs<-r(paths[1]);profiles<-r(paths[2]);context<-r(paths[3]);prefix<-r(paths[4])
f<-bind_rows(r(paths[5]),r(paths[6]));prevalence<-r(paths[7])
seed_context<-tribble(~entity_id,~threshold_operator,~threshold_percent,~marker,~population_scope,
 'seed:DB01101','>',30,'metastatic breast cancer who received capecitabine as a single agent','metastatic_breast_cancer',
 'seed:DB01181','>=',10,'ifosfamide monotherapy','ifosfamide_trials',
 'seed:DB01229','>=',20,'single-agent use of paclitaxel protein-bound','albumin_bound_metastatic_breast',
 'seed:DB09074','>=',10,'pooled single agent safety population','pooled_olaparib_safety')
for(i in seq_len(nrow(seed_context))){s<-prefix |> filter(entity_id==seed_context$entity_id[i]);stopifnot(nrow(s)==1L,
 str_detect(s$source_list_prefix,fixed(seed_context$marker[i])))}
write_csv(seed_context,'config/seed_single_agent_frequency_review.csv')
settings<-tibble(setting=c('descriptive_minimum_PT_terms','descriptive_require_same_route','descriptive_require_same_nominal_cutoff',
 'candidate_descriptive_treatment_scope'),value=c('5','TRUE','TRUE','single_investigational_treatment_arm_supported_by_source'),
 interpretation=c('Sparsity screen retained from earlier pilot','One shared label route only','Same cutoff does not imply equal populations or ascertainment',
 'Permits reviewed study-arm context; does not assert absence of supportive care or other background treatment'))
write_csv(settings,'config/context_overlap_settings.csv')
# Explicit reviewed source decisions, not a name/keyword classifier.
context<-context |> mutate(candidate_study_arm_supported=entity_id %in% c('sider:CID100040632','sider:CID100148121'))
write_csv(context,'data/curated/candidate_context_reassessment.csv',na='')
pairs<-pairs |> left_join(context |> select(candidate_profile_id=profile_id,candidate_regimen_review=regimen_review,
 candidate_cohort=cohort_review,candidate_operator=threshold_operator,candidate_threshold=threshold_percent,
 comparator_rule,frequency_limitation,candidate_study_arm_supported),by='candidate_profile_id') |>
 left_join(seed_context |> select(seed_entity_id=entity_id,seed_operator=threshold_operator,seed_threshold=threshold_percent),by='seed_entity_id') |>
 mutate(same_nominal_cutoff=!is.na(seed_threshold)&seed_threshold==candidate_threshold&seed_operator==candidate_operator,
 seed_explicit_scope_reviewed=seed_entity_id %in% seed_context$entity_id,
 descriptive_aligned_summary_eligible=seed_explicit_scope_reviewed & candidate_study_arm_supported &
 exact_same_single_route=='TRUE' & both_meet_term_minimum=='TRUE' & same_nominal_cutoff & same_current_full_key=='FALSE',
 strict_pair_ready=FALSE,
 assessment_scope=if_else(descriptive_aligned_summary_eligible,'descriptive_aligned_cutoff_summary_only','cross_context_description_only'))
# Counts and Jaccard are descriptive only. No ranking or FDR adjustment is generated.
sets<-split(f$pt_cui,f$entity_id);common<-prevalence$pt_cui[prevalence$common_term=='TRUE']
counts<-lapply(seq_len(nrow(pairs)),function(i){
 a<-unique(sets[[pairs$seed_entity_id[i]]]);b<-unique(sets[[pairs$candidate_entity_id[i]]]);shared<-intersect(a,b)
 tibble(shared_PT_terms=length(shared),union_PT_terms=length(union(a,b)),
 descriptive_jaccard=length(shared)/length(union(a,b)),shared_common_terms=sum(shared %in% common),
 shared_below_common_cutoff=sum(!shared %in% common))})
pairs<-bind_cols(pairs,bind_rows(counts));eligible<-pairs |> filter(descriptive_aligned_summary_eligible)
shared_terms<-bind_rows(lapply(seq_len(nrow(eligible)),function(i){
 ids<-intersect(sets[[eligible$seed_entity_id[i]]],sets[[eligible$candidate_entity_id[i]]])
 tibble(seed_entity_id=eligible$seed_entity_id[i],candidate_entity_id=eligible$candidate_entity_id[i],pt_cui=ids)
})) |> left_join(prevalence |> select(pt_cui,pt_name,SIDER_profile_fraction,common_term),by='pt_cui')
stopifnot(nrow(pairs)==384L,!anyDuplicated(pairs[c('seed_profile_id','candidate_profile_id')]),
 all(pairs$descriptive_jaccard>=0 & pairs$descriptive_jaccard<=1),all(!pairs$strict_pair_ready),
 all(!pairs$descriptive_aligned_summary_eligible | pairs$same_nominal_cutoff),identical(before,tools::md5sum(paths)))
write_csv(pairs,'data/curated/context_annotated_profile_pairs.csv',na='')
write_csv(pairs,'outputs/16e_context_pair_summary.csv',na='')
write_csv(eligible,'outputs/16e_aligned_summary_overlap_summary.csv',na='')
write_csv(shared_terms,'outputs/16e_aligned_shared_terms_summary.csv',na='')
write_csv(context |> select(-source_url,-retrieved_utc),'outputs/16e_candidate_reassessment_summary.csv',na='')
write_csv(tibble(file=paths,md5=unname(before)),'outputs/16e_input_summary.csv')
report<-c('SCRIPT 16e - CONTEXT REASSESSMENT AND DESCRIPTIVE OVERLAP',
 paste('Full reviewed profile pairs retained:',nrow(pairs)),paste('Strict-ready pairs:',sum(pairs$strict_pair_ready)),
 paste('Pairs in declared descriptive aligned-cutoff stratum:',nrow(eligible)),
 'Four explicitly single-agent seed lists have reviewed nominal frequency cutoffs; all six candidate contexts are retained.',
 'Pirfenidone and pralatrexate have source-supported single investigational treatment arms, without asserting no background/supportive treatment.',
 'Descriptive stratum requires reviewed seed single-agent scope, candidate study-arm support, same label route, >=5 terms on both sides, same nominal cutoff and no self full-key comparison.',
 'This relaxed descriptive stratum does not supersede the original stricter audit or prove cross-disease cohort comparability.',
 'All 384 raw set overlaps are labelled cross-context or aligned-summary descriptions; no rank, p-value or BH-adjusted claim is made.',
 'Shared common terms use the previously declared >=20% prevalence threshold across original SIDER profiles; lower prevalence is not target specificity.',
 'Source frequency definitions, comparator filtering, era, indication and population remain confounders even when nominal cutoffs match.',
 'Outputs: context_annotated_profile_pairs.csv (local); outputs/16e_* pair, eligible-pair and shared-term summaries.',
 'Next: final section-16 source/evidence ledger, consolidated feature export and target-evidence handoff without hiding pending review work.')
if(nrow(eligible))report<-c(report,vapply(seq_len(nrow(eligible)),function(i)paste0(
 'Aligned descriptive pair: ',eligible$seed_name[i],' / ',eligible$candidate_current_name[i],
 '; shared=',eligible$shared_PT_terms[i],'; union=',eligible$union_PT_terms[i],
 '; Jaccard=',format(round(eligible$descriptive_jaccard[i],4),nsmall=4),
 '; shared common terms=',eligible$shared_common_terms[i],
 '; shared below common cutoff=',eligible$shared_below_common_cutoff[i],'.'),character(1)))
writeLines(report,'outputs/16e_reassessment_summary.txt')
writeLines('PASS: four seed source-scope markers, all 384 unique pairs, bounded set scores, explicit analysis strata and unchanged protected inputs.',
 'outputs/16e_reassessment_validation.txt')
cat(paste0('\nSCRIPT 16e RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/16e_reassess_context_and_overlap.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
