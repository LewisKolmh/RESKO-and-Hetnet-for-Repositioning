# Audit reviewed profile compatibility; no similarity scores or activity claims.
library(readr);library(dplyr);library(tidyr);library(stringr);library(tibble)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
paths<-c('outputs/16a_candidate_PT_features.csv','outputs/16b_seed_PT_features.csv',
 'data/curated/candidate_label_identity_review.csv','data/curated/candidate_identity_registry.csv',
 'data/curated/seed_identity_crosscheck_reviewed.csv','outputs/15_full_inventory_coverage.csv',
 'data/curated/candidate_common_summary_segment_review.csv','data/curated/seed_representative_summary_review.csv',
 'data/curated/comparable_label_evidence_segments.csv')
if(!all(file.exists(paths)))stop('Run from project root after 16b.')
before<-tools::md5sum(paths)
ca<-r(paths[1]);se<-r(paths[2]);labels<-r(paths[3]);ci<-r(paths[4]);si<-r(paths[5]);inventory<-r(paths[6])
settings<-tibble(setting=c('minimum_reviewed_PT_terms','strict_comparison_scope','frequency_alignment_required'),
 value=c('5','current_label_explicit_single_agent_same_route','TRUE'),
 rationale=c('Retains the earlier script-12 minimum as a declared sparsity screen, not statistical significance',
 'A deliberately narrow sensitivity stratum; different routes are not universally scientifically incomparable',
 'Common-reaction thresholds and clinical cohorts differ; alignment must be reviewed before strict comparison'))
write_csv(settings,'config/profile_compatibility_settings.csv')
minimum<-as.integer(settings$value[settings$setting=='minimum_reviewed_PT_terms'])
# Re-evaluate from source qualifiers rather than simply copying prior comparison_ready flags.
seed_profiles<-se |> group_by(entity_id,entity_name,drugbank_id,set_id,label_id,effective_time,label_version,
 regimen_scope,population_scope,routes,seed_tiers,cancers,identity_scope) |>
 summarise(reviewed_PT_terms=n_distinct(pt_cui),.groups='drop') |> mutate(cohort='seed') |>
 left_join(si |> select(drugbank_id,current_inchikey=inchikey,current_record_name=drug_name,source_drug_type=drug_type),by='drugbank_id')
candidate_profiles<-ca |> group_by(entity_id,entity_name,stitch_flat_id,set_id,label_id,effective_time,label_version,
 regimen_scope,population_scope,identity_scope) |> summarise(reviewed_PT_terms=n_distinct(pt_cui),.groups='drop') |>
 left_join(labels |> select(stitch_flat_id,label_id,routes,supported_identity),by=c('stitch_flat_id','label_id')) |>
 left_join(ci |> select(stitch_flat_id,current_inchikey,current_record_name=matched_names,historical_identity_status,
 already_seed_full_key),by='stitch_flat_id') |> mutate(cohort='candidate',drugbank_id=NA_character_)
stopifnot(nrow(candidate_profiles)==n_distinct(ca$entity_id),all(candidate_profiles$supported_identity=='TRUE'))
profiles<-bind_rows(seed_profiles,candidate_profiles) |>
 mutate(profile_id=paste(entity_id,set_id,population_scope,regimen_scope,sep=':'),
 route_status=case_when(is.na(routes)|!nzchar(routes)~'unknown',str_detect(routes,fixed('|'))~'multiple_label_routes',TRUE~'single_label_route'),
 regimen_class=case_when(regimen_scope %in% c('explicit_single_agent','explicit_monotherapy','explicit_pooled_single_agent')~'explicit_single_agent',
 str_detect(regimen_scope,'combination')~'explicit_combination',TRUE~'regimen_unresolved_or_pooled'),
 frequency_alignment_status='not_yet_reviewed',identity_layer='current_label_attribution_supported',
 sparse_profile=reviewed_PT_terms<minimum,
 non_cancer_population=str_detect(population_scope,'NON_CANCER'),
 local_product_context=str_detect(population_scope,'intravesical|ocular|local_product'),
 formulation_specific=str_detect(population_scope,'albumin_bound|product_only|selected_doxorubicin'),
 strict_profile_ready=regimen_class=='explicit_single_agent' & route_status=='single_label_route' & !sparse_profile &
 !non_cancer_population & frequency_alignment_status=='reviewed_compatible')
stopifnot(!anyDuplicated(profiles$profile_id))
# Keep full source prefix metadata locally to support the next frequency/cohort review; do not infer thresholds from PT counts.
segments<-r(paths[9]);csp<-r(paths[7]);ssp<-r(paths[8])
spans<-bind_rows(csp |> select(entity_id,label_id,source_segment_id,accepted_span_start,accepted_span_end),
 ssp |> select(entity_id,label_id,source_segment_id,accepted_span_start,accepted_span_end)) |>
 filter(!is.na(accepted_span_start)) |>
 left_join(segments |> select(source_segment_id,source_text),by='source_segment_id') |>
 mutate(source_list_prefix=str_sub(source_text,1L,as.integer(accepted_span_start)-1L))
# Offsets refer to normalised text; recompute prefix on that representation.
spans<-spans |> mutate(source_list_prefix=str_sub(str_to_lower(str_squish(str_replace_all(source_text,'[\u2010-\u2015]','-'))),
 1L,as.integer(accepted_span_start)-1L))
write_csv(spans,'data/curated/profile_comparison_source_context.csv',na='')
# Cartesian pairs retain all reviewed profiles, including those failing the declared strict stratum.
s<-profiles |> filter(cohort=='seed');c<-profiles |> filter(cohort=='candidate')
pairs<-crossing(seed_profile_id=s$profile_id,candidate_profile_id=c$profile_id) |>
 left_join(s |> select(seed_profile_id=profile_id,seed_entity_id=entity_id,seed_name=entity_name,
 seed_key=current_inchikey,seed_route=routes,seed_regimen=regimen_class,seed_terms=reviewed_PT_terms,
 seed_sparse=sparse_profile,seed_non_cancer=non_cancer_population,seed_route_status=route_status),by='seed_profile_id') |>
 left_join(c |> select(candidate_profile_id=profile_id,candidate_entity_id=entity_id,candidate_name=entity_name,
 candidate_current_name=current_record_name,candidate_key=current_inchikey,candidate_route=routes,
 candidate_regimen=regimen_class,candidate_terms=reviewed_PT_terms,candidate_sparse=sparse_profile,
 candidate_route_status=route_status),by='candidate_profile_id') |>
 mutate(same_current_full_key=!is.na(seed_key)&!is.na(candidate_key)&seed_key==candidate_key,
 both_explicit_single_agent=seed_regimen=='explicit_single_agent'&candidate_regimen=='explicit_single_agent',
 exact_same_single_route=seed_route_status=='single_label_route'&candidate_route_status=='single_label_route'&seed_route==candidate_route,
 both_meet_term_minimum=!seed_sparse&!candidate_sparse,
 frequency_and_cohort_alignment_reviewed=FALSE,
 strict_pair_ready=!same_current_full_key&both_explicit_single_agent&exact_same_single_route&both_meet_term_minimum&
 !seed_non_cancer&frequency_and_cohort_alignment_reviewed)
barriers<-function(self,mono,route,terms,noncancer){paste(c(if(self)'same_current_full_key_self_comparison',
 if(!mono)'regimen_not_explicit_single_agent_on_both_sides',if(!route)'outside_same_single_route_stratum',
 if(!terms)'one_or_both_profiles_below_declared_term_minimum',if(noncancer)'non_cancer_seed_population',
 'frequency_and_cohort_alignment_pending'),collapse='|')}
pairs$strict_stratum_barriers<-mapply(barriers,pairs$same_current_full_key,pairs$both_explicit_single_agent,
 pairs$exact_same_single_route,pairs$both_meet_term_minimum,pairs$seed_non_cancer,USE.NAMES=FALSE)
# Actionable profile queue: one row per profile, not hundreds of duplicated pair tasks.
actions<-profiles |> mutate(next_action=case_when(sparse_profile~'Review additional same-context label evidence to assess sparse coverage',
 regimen_class!='explicit_single_agent'~'Inspect full source cohort and regimen before single-agent attribution',
 route_status!='single_label_route'~'Resolve route/formulation for the reviewed source cohort',
 TRUE~'Review frequency cutoff, source cohort and population alignment'),
 action_priority=if_else(cohort=='candidate','candidate_context_first','seed_context_followup'))
# Coverage remains explicit for all 200 inventory entities, not just the 70 reviewed profiles.
full<-inventory |> left_join(profiles |> group_by(entity_id) |> summarise(reviewed_profiles=n(),
 reviewed_profile_terms=sum(reviewed_PT_terms),strict_profiles_ready=sum(strict_profile_ready),.groups='drop'),by='entity_id') |>
 mutate(across(c(reviewed_profiles,reviewed_profile_terms,strict_profiles_ready),~replace_na(.x,0L)),
 reviewed_profile_status=case_when(reviewed_profiles>0~'reviewed_source_context_profile',
 label_coverage_status=='no_collected_label'~'no_collected_label',TRUE~'collected_label_without_reviewed_profile'))
stopifnot(nrow(pairs)==nrow(s)*nrow(c),!anyDuplicated(pairs[c('seed_profile_id','candidate_profile_id')]),
 nrow(full)==nrow(inventory),!anyDuplicated(full$entity_id),
 all(!pairs$strict_pair_ready|(!pairs$same_current_full_key & pairs$both_explicit_single_agent & pairs$exact_same_single_route &
 pairs$both_meet_term_minimum & pairs$frequency_and_cohort_alignment_reviewed)),identical(before,tools::md5sum(paths)))
write_csv(profiles,'data/curated/reviewed_profile_compatibility_registry.csv',na='')
write_csv(pairs,'data/curated/reviewed_profile_pair_compatibility.csv',na='')
write_csv(profiles,'outputs/16c_profile_compatibility_summary.csv',na='')
write_csv(pairs,'outputs/16c_pair_compatibility_summary.csv',na='')
write_csv(actions |> select(profile_id,entity_id,entity_name,cohort,current_record_name,reviewed_PT_terms,routes,
 regimen_scope,population_scope,sparse_profile,next_action,action_priority),
 'outputs/16c_profile_followup_summary.csv',na='')
write_csv(full,'outputs/16c_full_inventory_coverage.csv',na='')
write_csv(profiles |> count(cohort,regimen_class,route_status,sparse_profile,name='profiles'),
 'outputs/16c_readiness_summary.csv',na='')
write_csv(tibble(file=paths,md5=unname(before)),'outputs/16c_input_summary.csv')
report<-c('SCRIPT 16c - REVIEWED PROFILE COMPATIBILITY AUDIT',
 paste('Reviewed seed profiles:',nrow(s)),paste('Reviewed candidate profiles:',nrow(c)),
 paste('Seed/candidate profile pairs audited:',nrow(pairs)),
 paste('Seed profiles explicitly described as single-agent:',sum(s$regimen_class=='explicit_single_agent')),
 paste('Candidate profiles explicitly described as single-agent:',sum(c$regimen_class=='explicit_single_agent')),
 paste('Profiles below declared minimum of',minimum,'PTs:',sum(profiles$sparse_profile)),
 paste('Pairs satisfying the same-single-route condition:',sum(pairs$exact_same_single_route)),
 paste('Pairs with frequency/cohort alignment already reviewed:',sum(pairs$frequency_and_cohort_alignment_reviewed)),
 paste('Pairs ready for the declared strict stratum:',sum(pairs$strict_pair_ready)),
 paste('Full inventory entities retained:',nrow(full)),
 'Inputs: reviewed 16a/16b features, label routes, current structure registry and original inventory coverage.',
 'Outputs: reviewed_profile_compatibility_registry.csv and reviewed_profile_pair_compatibility.csv locally;',
 'outputs/16c_* profile/pair audit, follow-up queue, full inventory coverage, hashes and validation.',
 'The strict stratum requires explicit single-agent context, one identical label route, at least five reviewed PTs, no self full-key comparison and reviewed frequency/cohort alignment.',
 'These are declared analysis constraints, not universal biological incompatibility rules or significance tests.',
 'Single active ingredient does not establish a monotherapy trial. Source label routes may be broader than trial routes and require context review.',
 'Current named DrugBank identities accompany abbreviated SIDER labels; historical SIDER identity/formulation uncertainty is not repaired by this audit.',
 'Current-label attribution is separate from strict structural equivalence; related salts, active moieties and nonidentical stereochemistry are not silently merged.',
 'All barriers are retained per pair; incomplete coverage is not scored as zero similarity or negative target evidence.',
 'Regimen and frequency/cohort gaps explain zero strict-ready comparisons; prior FALSE readiness flags were not used as the decision rule.',
 'No similarity ranks, p-values, BH correction or eEF1A binding/inhibition/PPI claims are generated.',
 'Next: inspect full candidate source cohorts/regimens and threshold definitions using the profile-level follow-up queue; then reassess the declared stratum.')
writeLines(report,'outputs/16c_profile_compatibility_summary.txt')
writeLines(c('PASS: complete reviewed-profile Cartesian product, unique pair/profile keys, all inventory entities retained.',
 'PASS: candidate routes joined by exact group/label key; source identity support checked; raw input hashes unchanged.',
 'PASS: strict readiness derived from explicit evidence fields; all failed pair constraints recorded rather than hidden.'),
 'outputs/16c_profile_compatibility_validation.txt')
cat(paste0('\nSCRIPT 16c RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/16c_audit_profile_compatibility.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
