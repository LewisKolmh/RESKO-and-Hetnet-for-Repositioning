# Audit descriptive overlap contributions; no target-binding inference or significance testing.
library(readr);library(dplyr);library(tidyr);library(stringr);library(tibble)
if(!file.exists('data/curated/current_label_SIDER_similarity_pilot.csv'))stop('Run from project root after script 12.')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
settings<-tibble(parameter=c('common_term_profile_fraction','minimum_retained_terms','top_rank_review_window','minimum_shared_PT_terms'),
 value=c('0.20','3','10','3'),purpose=c('Declared sensitivity cutoff; not optimised on ranks',
 'Both reduced profiles must retain this many terms','Review window includes score ties',
 'Descriptive corroboration threshold; not statistical significance'))
config_path<-'config/side_effect_audit_settings.csv'
if(file.exists(config_path))settings<-r(config_path)
get_setting<-function(k)as.numeric(settings$value[match(k,settings$parameter)])
cutoff<-get_setting('common_term_profile_fraction');minimum<-get_setting('minimum_retained_terms')
window<-get_setting('top_rank_review_window');shared_min<-get_setting('minimum_shared_PT_terms')
stopifnot(cutoff>0,cutoff<1,minimum>=1,window>=1,shared_min>=1)
write_csv(settings,config_path)
pilot<-read_csv('data/curated/current_label_SIDER_similarity_pilot.csv',show_col_types=FALSE)
f<-r('data/curated/current_label_PT_supplement.csv')
seed<-f |> filter(comparison_stratum=='explicit_monotherapy_or_single_agent',drugbank_id %in% pilot$drugbank_id) |>
 distinct(drugbank_id,pt_cui)
se<-read_tsv('data/raw/sider/meddra_all_se.tsv.gz',col_names=c('stitch_flat_id','stereo','label_cui','type','pt_cui','term'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
ref<-se |> filter(type=='PT') |> distinct(stitch_flat_id,pt_cui)
N<-n_distinct(ref$stitch_flat_id)
labels<-se |> filter(type=='PT') |> group_by(pt_cui) |> summarise(pt_name=first(term),.groups='drop')
weights<-ref |> count(pt_cui,name='SIDER_profiles_with_term') |>
 mutate(SIDER_profile_fraction=SIDER_profiles_with_term/N,
 idf_weight=1+log((N+1)/(SIDER_profiles_with_term+1)),
 common_term=SIDER_profile_fraction>=cutoff) |> left_join(labels,by='pt_cui')
stopifnot(all(seed$pt_cui %in% weights$pt_cui))
seed_w<-seed |> left_join(weights,by='pt_cui')
ref_w<-ref |> left_join(weights,by='pt_cui')
seed_mass<-seed_w |> group_by(drugbank_id) |> summarise(seed_IDF_mass=sum(idf_weight),
 seed_retained_terms=sum(!common_term),.groups='drop')
ref_mass<-ref_w |> group_by(stitch_flat_id) |> summarise(reference_IDF_mass=sum(idf_weight),
 reference_retained_terms=sum(!common_term),.groups='drop')
shared<-seed |> inner_join(ref_w,by='pt_cui',relationship='many-to-many')
shared_mass<-shared |> group_by(drugbank_id,stitch_flat_id) |>
 summarise(shared_IDF_mass=sum(idf_weight),shared_common_terms=sum(common_term),
 shared_retained_terms=sum(!common_term),.groups='drop')
audit<-pilot |> left_join(seed_mass,by='drugbank_id') |> left_join(ref_mass,by='stitch_flat_id') |>
 left_join(shared_mass,by=c('drugbank_id','stitch_flat_id')) |>
 mutate(across(c(shared_IDF_mass,shared_common_terms,shared_retained_terms),~replace_na(.x,0)),
 weighted_jaccard=shared_IDF_mass/(seed_IDF_mass+reference_IDF_mass-shared_IDF_mass),
 retained_profile_status=case_when(seed_retained_terms<minimum~'seed_retains_too_few_terms',
 reference_retained_terms<minimum~'reference_retains_too_few_terms',TRUE~'eligible'),
 reduced_jaccard=if_else(retained_profile_status=='eligible',
 shared_retained_terms/(seed_retained_terms+reference_retained_terms-shared_retained_terms),NA_real_),
 common_term_fraction_of_overlap=if_else(shared_PT_terms>0,shared_common_terms/shared_PT_terms,NA_real_)) |>
 group_by(drugbank_id) |> mutate(baseline_rank_with_ties=min_rank(desc(jaccard)),
 weighted_rank_with_ties=min_rank(desc(weighted_jaccard)),reduced_rank_with_ties=min_rank(desc(reduced_jaccard)),
 weighted_rank_change=baseline_rank_with_ties-weighted_rank_with_ties,
 review_union=baseline_rank_with_ties<=window|weighted_rank_with_ties<=window|
 coalesce(reduced_rank_with_ties<=window,FALSE),
 retained_in_both_main_top_windows=baseline_rank_with_ties<=window&weighted_rank_with_ties<=window,
 priority_for_identity_review=retained_in_both_main_top_windows&shared_PT_terms>=shared_min) |> ungroup()
review<-audit |> filter(review_union)
contributions<-shared |> inner_join(review |> select(drugbank_id,stitch_flat_id,drug_name,SIDER_drug_name,
 shared_PT_terms,shared_IDF_mass,jaccard,weighted_jaccard),by=c('drugbank_id','stitch_flat_id')) |>
 mutate(unweighted_overlap_share=1/shared_PT_terms,
 IDF_overlap_share=idf_weight/shared_IDF_mass,
 weighted_jaccard_contribution=idf_weight/(seed_mass$seed_IDF_mass[match(drugbank_id,seed_mass$drugbank_id)]+
 ref_mass$reference_IDF_mass[match(stitch_flat_id,ref_mass$stitch_flat_id)]-shared_IDF_mass)) |>
 arrange(drugbank_id,stitch_flat_id,desc(IDF_overlap_share),pt_cui)
# Rank consistency is descriptive and includes ties; tiny numerical differences are not significant.
stability<-audit |> group_by(drugbank_id,drug_name) |> summarise(
 nonself_comparisons=n(),baseline_top_window=sum(baseline_rank_with_ties<=window),
 weighted_top_window=sum(weighted_rank_with_ties<=window),
 both_top_windows=sum(retained_in_both_main_top_windows),
 overlap_fraction_of_baseline_window=both_top_windows/baseline_top_window,
 baseline_top_window_median_common_overlap=median(common_term_fraction_of_overlap[baseline_rank_with_ties<=window],na.rm=TRUE),
 baseline_top_window_median_reference_PT=median(SIDER_PT_terms[baseline_rank_with_ties<=window]),
 weighted_top_window_median_reference_PT=median(SIDER_PT_terms[weighted_rank_with_ties<=window]),
 reference_size_vs_baseline_spearman=cor(SIDER_PT_terms,jaccard,method='spearman'),
 reference_size_vs_weighted_spearman=cor(SIDER_PT_terms,weighted_jaccard,method='spearman'),
 reduced_seed_terms=first(seed_retained_terms),reduced_eligible_comparisons=sum(retained_profile_status=='eligible'),.groups='drop')
queue<-review |> filter(priority_for_identity_review) |>
 group_by(stitch_flat_id,SIDER_drug_name) |>
 summarise(seed_profiles_for_descriptive_review=n_distinct(drugbank_id),
 seed_names=paste(sort(unique(drug_name)),collapse='|'),best_baseline_rank=min(baseline_rank_with_ties),
 best_weighted_rank=min(weighted_rank_with_ties),max_shared_PT_terms=max(shared_PT_terms),
 max_lower_frequency_shared_terms=max(shared_retained_terms),
 max_weighted_jaccard=max(weighted_jaccard),.groups='drop') |>
 mutate(name_label_flag=if_else(is.na(SIDER_drug_name)|str_detect(str_to_lower(SIDER_drug_name),'^(nsc|unknown|compound|drug)$'),
 'uninformative_label_requires_identifier_resolution','label_present_identity_still_unverified'),
 overlap_context=if_else(max_lower_frequency_shared_terms>0,
 'some_overlap_below_common_term_cutoff_NOT_target_specific',
 'all_overlap_common_terms_only'),
 disposition='identity_and_comparable_label_review_only_NOT_eEF1A_evidence') |>
 arrange(best_weighted_rank,best_baseline_rank,stitch_flat_id)
# Independent set calculations validate all three scoring methods for one nonzero pair.
t<-review |> arrange(desc(shared_PT_terms)) |> slice(1)
A<-unique(seed$pt_cui[seed$drugbank_id==t$drugbank_id]);B<-unique(ref$pt_cui[ref$stitch_flat_id==t$stitch_flat_id])
w<-setNames(weights$idf_weight,weights$pt_cui);common<-weights$pt_cui[weights$common_term]
stopifnot(abs(t$jaccard-length(intersect(A,B))/length(union(A,B)))<1e-12,
 abs(t$weighted_jaccard-sum(w[intersect(A,B)])/sum(w[union(A,B)]))<1e-12,
 !anyDuplicated(audit[c('drugbank_id','stitch_flat_id')]),nrow(audit)==nrow(pilot),
 all(audit$weighted_jaccard>=0 & audit$weighted_jaccard<=1),
 all(is.na(audit$reduced_jaccard[audit$retained_profile_status!='eligible'])))
if(t$retained_profile_status=='eligible')stopifnot(abs(t$reduced_jaccard-
 length(intersect(setdiff(A,common),setdiff(B,common)))/length(union(setdiff(A,common),setdiff(B,common))))<1e-12)
checks<-contributions |> group_by(drugbank_id,stitch_flat_id) |>
 summarise(total_fraction=sum(IDF_overlap_share),total_score=sum(weighted_jaccard_contribution),
 score=first(weighted_jaccard),.groups='drop')
stopifnot(all(abs(checks$total_fraction-1)<1e-12),all(abs(checks$total_score-checks$score)<1e-12))
# Check that the cached pilot still describes the current source sets for every pair.
source_seed_count<-seed |> count(drugbank_id,name='source_seed_terms')
source_ref_count<-ref |> count(stitch_flat_id,name='source_reference_terms')
source_shared_count<-shared |> count(drugbank_id,stitch_flat_id,name='source_shared_terms')
verify<-audit |> left_join(source_seed_count,by='drugbank_id') |>
 left_join(source_ref_count,by='stitch_flat_id') |>
 left_join(source_shared_count,by=c('drugbank_id','stitch_flat_id')) |>
 mutate(source_shared_terms=replace_na(source_shared_terms,0L))
stopifnot(all(verify$seed_PT_terms==verify$source_seed_terms),
 all(verify$SIDER_PT_terms==verify$source_reference_terms),
 all(verify$shared_PT_terms==verify$source_shared_terms))
inputs<-c('scripts/13_audit_side_effect_similarity.R',config_path,
 'data/curated/current_label_SIDER_similarity_pilot.csv','data/curated/current_label_PT_supplement.csv',
 'data/raw/sider/meddra_all_se.tsv.gz')
write_csv(tibble(input_file=inputs,md5=unname(tools::md5sum(inputs))),
 'outputs/13_input_provenance_summary.csv')
writeLines(trimws(capture.output(sessionInfo()),which='right'),'outputs/13_environment_summary.txt')
leading<-audit |> filter(baseline_rank_with_ties<=window)
write_csv(weights,'data/curated/side_effect_term_prevalence.csv')
write_csv(audit,'data/curated/side_effect_similarity_audit.csv',na='')
write_csv(contributions,'data/curated/side_effect_overlap_contributions.csv',na='')
write_csv(review,'outputs/13_rank_comparison_summary.csv',na='')
write_csv(contributions,'outputs/13_shared_term_contributions_summary.csv',na='')
write_csv(stability,'outputs/13_rank_stability_summary.csv',na='')
write_csv(queue,'outputs/13_identity_review_queue_summary.csv',na='')
write_csv(weights |> semi_join(seed,by='pt_cui') |> arrange(desc(SIDER_profile_fraction)),
 'outputs/13_seed_term_prevalence_summary.csv')
report<-c('SCRIPT 13 - SIDE-EFFECT OVERLAP AND RANK SENSITIVITY AUDIT',
 paste('Seed profiles:',n_distinct(audit$drugbank_id)),paste('Non-self comparisons:',nrow(audit)),
 paste('Common terms at declared cutoff:',sum(weights$common_term),'of',nrow(weights)),
 paste('Review-union seed/reference pairs (including ties):',nrow(review)),
 paste('Shared-term contribution rows:',nrow(contributions)),
 paste('Unique SIDER groups queued for identity/comparable-label review:',nrow(queue)),
 paste('Queued groups with any overlap below common-term cutoff:',sum(queue$max_lower_frequency_shared_terms>0)),
 paste('Leading baseline pairs with common-term-only overlap:',sum(leading$shared_retained_terms==0),'of',nrow(leading)),
 paste('Seeds with >=',minimum,'terms after common-term removal:',sum(stability$reduced_seed_terms>=minimum)),
 'Inputs: original script-12 comparisons, reviewed single-agent supplementary profiles, raw SIDER PT profiles.',
 'IDF weight = 1 + log((N+1)/(df+1)); N is all SIDER flat-group profiles, df is profile frequency.',
 'Weighted Jaccard sums term weights in intersection divided by union; it is descriptive, not a probability.',
 paste('Sensitivity removes terms found in >=',100*cutoff,'% of profiles; retain >=',minimum,'terms in BOTH profiles.'),
 'Excluded reduced profiles receive NA, not a zero score. No ranking cutoff was optimised after inspecting results.',
 paste('Top-window review uses tied ranks <=',window,'and >=',shared_min,'shared PTs in both raw/IDF windows; workload rule, not significance.'),
 'Shared-term contributions sum to weighted scores. Profile-size correlations are descriptive, with no causal interpretation.',
 'Outputs: full audit/term prevalence/contributions locally; compact rank, contribution, stability and review-queue CSVs in outputs/13_*.',
 'The review queue consists of unverified SIDER flat-group identities. Multiple seed appearances are not independent target evidence.',
 'Common toxicities, sparse common-reaction lists, historical labels, route/formulation and source-era differences remain confounders.',
 'No p-values, BH correction, binding, inhibition, allostery or PPI-disruption claims are generated.',
 'Next: verify queued chemical identities and obtain comparably processed current labels before target-evidence prioritisation.')
writeLines(report,'outputs/13_similarity_audit_summary.txt')
writeLines(c('Script 13 checks passed: independent set scores; IDF contribution sums;',
 'unchanged comparison universe; unique pair keys; bounded weights/scores;',
 'tied-rank review windows; too-small reduced profiles remain NA.'),'outputs/13_similarity_audit_validation.txt')
cat(paste0('\nSCRIPT 13 RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/13_audit_side_effect_similarity.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
