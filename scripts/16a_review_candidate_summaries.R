# Assistant-reviewed bounded reaction lists; explicit context, never whole-paragraph acceptance.
library(readr);library(dplyr);library(stringr);library(tibble)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
paths<-c('data/curated/comparable_label_common_summary_review_queue.csv',
 'data/curated/comparable_label_evidence_segments.csv','data/curated/comparable_label_PT_candidate_evidence.csv')
if(!all(file.exists(paths)))stop('Run from project root after script 15.')
before<-tools::md5sum(paths);q<-r(paths[1]);sg<-r(paths[2]);all_hits<-r(paths[3])
norm<-function(z)str_to_lower(str_squish(str_replace_all(z,'[\u2010-\u2015]','-')))
# These exact source clauses were inspected before writing the rules. Matching guards against source drift.
rules<-tribble(~entity_id,~label_id,~segment_index,~reaction_list,~regimen_scope,~population_scope,
 'sider:CID100006116','5d3fdea6-fca8-4dba-e063-6394a90a4265',2L,
 'hypercalcemia, nausea and vomiting','label_summary_regimen_unspecified','calcium_acetate_label_population',
 'sider:CID100040632','24b60d05-c1f7-4fba-ba4f-009eceb0c409',2L,
 'nausea, rash, abdominal pain, upper respiratory tract infection, diarrhea, fatigue, headache, decreased appetite, dyspepsia, dizziness, vomiting, gastro-esophageal reflux disease, sinusitis, insomnia, weight decreased, and arthralgia',
 'label_summary_regimen_unspecified','pirfenidone_label_population',
 'sider:CID100148121','a10f6330-f643-4d2b-b4b0-541d5e8ec387',2L,
 'mucositis, thrombocytopenia, nausea, and fatigue','label_summary_regimen_unspecified','pralatrexate_label_population',
 'sider:CID109912092','8d1d49ca-6dc1-4580-b2d5-f69ccb4f1c8e',2L,
 'diarrhea, nausea, and headache','label_summary_regimen_unspecified','psoriatic_arthritis_only',
 'sider:CID110482134','8e05336e-2d56-44e5-a76b-5b39f2f26a0f',2L,
 'diarrhea, flatulence, and headache','label_summary_regimen_unspecified','adult_glycerol_phenylbutyrate_only',
 'sider:CID116126651','43280dcb-52ae-40e4-bd52-723f71f74ab7',2L,
 'thrombocytopenia, diarrhea, fatigue, musculoskeletal pain, neutropenia, rash, anemia, bruising, and nausea',
 'pooled_B_cell_label_summary_regimen_unspecified','B_cell_malignancies_only')
rules<-rules |> mutate(review_method='assistant_reviewed_exact_positive_list',review_batch='16a_candidate_common_summaries')
write_csv(rules,'config/candidate_positive_list_review_rules.csv',na='')
candidate_segments<-q |> filter(cohort=='candidate') |> distinct(entity_id,entity_name,label_id,source_segment_id) |>
 left_join(sg,by=c('source_segment_id','label_id'))
# Exact PT-name mappings reviewed against the supplied dictionary; do not choose broader LLT-derived alternatives.
resolutions<-tribble(~alias_norm,~reviewed_pt_cui,
 'abdominal pain','C0000737','rash','C0015230','fatigue','C0015672','arthralgia','C0003862') |>
 mutate(reason='Literal source term exactly names the supplied SIDER PT; alternate LLT-derived PT is not asserted')
dictionary<-r('data/curated/current_label_PT_dictionary.csv')
stopifnot(nrow(inner_join(resolutions,dictionary |> filter(alias_source=='SIDER_PT_name'),
 by=c('alias_norm','reviewed_pt_cui'='pt_cui')))==nrow(resolutions))
write_csv(resolutions,'config/candidate_PT_mapping_review.csv',na='')
accepted<-list();clauses<-list()
for(i in seq_len(nrow(rules))){
 a<-rules[i,];s<-candidate_segments |> filter(entity_id==a$entity_id,label_id==a$label_id,
 source_field=='adverse_reactions',as.integer(segment_index)==a$segment_index)
 stopifnot(nrow(s)==1L)
 text<-norm(s$source_text);literal<-norm(a$reaction_list);loc<-str_locate(text,fixed(literal))
 stopifnot(!is.na(loc[1]),str_count(text,fixed(literal))==1L)
 h<-all_hits |> filter(entity_id==a$entity_id,source_segment_id==s$source_segment_id,
 as.integer(match_start)>=loc[1],as.integer(match_end)<=loc[2])
 stopifnot(nrow(h)>0L,all(h$identity_eligible=='TRUE'))
 h<-h |> left_join(resolutions |> select(alias_norm,reviewed_pt_cui),by='alias_norm') |>
  mutate(review_decision=case_when(broad_term=='TRUE'~'hold_dictionary_scope',
   ambiguous_dictionary_mapping=='TRUE' & !is.na(reviewed_pt_cui) & pt_cui==reviewed_pt_cui~'accept_positive_list_context',
   ambiguous_dictionary_mapping=='TRUE' & !is.na(reviewed_pt_cui)~'reject_alternate_dictionary_mapping',
   ambiguous_dictionary_mapping=='TRUE'~'hold_dictionary_scope',TRUE~'accept_positive_list_context'),
  review_reason=if_else(!is.na(reviewed_pt_cui),
   'Reviewed exact PT name selected; alternate LLT-derived mapping not asserted',
   'Literal mention wholly within assistant-reviewed positive reaction list'),
 regimen_scope=a$regimen_scope,population_scope=a$population_scope,review_method=a$review_method,
 review_batch=a$review_batch,comparison_ready=FALSE,evidence_group=paste0('SPL:',set_id))
 accepted[[i]]<-h
 clauses[[i]]<-tibble(entity_id=a$entity_id,source_segment_id=s$source_segment_id,
  accepted_span_start=loc[1],accepted_span_end=loc[2],regimen_scope=a$regimen_scope,population_scope=a$population_scope)
}
reviewed_hits<-bind_rows(accepted);spans<-bind_rows(clauses)
segment_decisions<-candidate_segments |> select(entity_id,entity_name,label_id,source_segment_id) |>
 left_join(spans,by=c('entity_id','source_segment_id')) |>
 mutate(segment_review_decision=if_else(!is.na(accepted_span_start),'accept_bounded_list_only','hold_outside_selected_profile'),
  reason=case_when(!is.na(accepted_span_start)~'Only mentions inside the exact reviewed reaction list are accepted',
   entity_id=='sider:CID100003053'~'Technetium Tc-99m pentetate formulation; do not attribute complex-specific reactions to free pentetic acid',
   entity_id=='sider:CID111556711'~'Isatuximab/carfilzomib/dexamethasone combination context; no single-drug attribution',
   TRUE~'Duplicate label support, distinct cohort, discontinuation/dose-change endpoint or causality-unspecified summary; not pooled into chosen profile'))
# Exhaustive decisions for this bounded batch; the much larger seed/table queue remains pending.
queue_decisions<-q |> filter(cohort=='candidate') |>
 left_join(segment_decisions |> select(entity_id,source_segment_id,segment_review_decision,reason),
 by=c('entity_id','source_segment_id')) |>
 left_join(reviewed_hits |> select(evidence_key,bounded_decision=review_decision),by='evidence_key') |>
 mutate(batch_review_decision=coalesce(bounded_decision,'hold_outside_accepted_span'))
features<-reviewed_hits |> filter(review_decision=='accept_positive_list_context') |>
 distinct(entity_id,entity_name,drugbank_id,stitch_flat_id,pt_cui,pt_name,set_id,label_id,effective_time,label_version,
 source_url,identity_scope,comparison_layer,regimen_scope,population_scope,evidence_group,review_method,review_batch,comparison_ready) |>
 mutate(present=1L,feature_layer='reviewed_candidate_positive_list_context')
stopifnot(!anyNA(queue_decisions$segment_review_decision),nrow(queue_decisions)==sum(q$cohort=='candidate'),
 !anyDuplicated(features[c('entity_id','pt_cui','regimen_scope','population_scope')]),
 !any(features$entity_id %in% c('sider:CID100003053','sider:CID111556711','sider:CID100135113')),
 all(!features$comparison_ready),identical(before,tools::md5sum(paths)))
write_csv(segment_decisions,'data/curated/candidate_common_summary_segment_review.csv',na='')
write_csv(queue_decisions,'data/curated/candidate_common_summary_batch_review.csv',na='')
write_csv(reviewed_hits,'data/curated/candidate_positive_list_reviewed_evidence.csv',na='')
write_csv(features,'data/curated/candidate_positive_list_PT_features.csv',na='')
write_csv(features,'outputs/16a_candidate_PT_features.csv',na='')
write_csv(features |> count(entity_id,entity_name,regimen_scope,population_scope,name='reviewed_PT_terms'),
 'outputs/16a_candidate_profile_coverage.csv',na='')
write_csv(queue_decisions |> count(batch_review_decision,reason,name='evidence_rows'),
 'outputs/16a_candidate_review_summary.csv',na='')
write_csv(tibble(file=paths,md5=unname(before)),'outputs/16a_input_summary.csv')
report<-c('SCRIPT 16a - CANDIDATE COMMON-SUMMARY CONTEXT REVIEW',
 paste('Candidate summary segments reviewed:',nrow(segment_decisions)),paste('Candidate queue evidence rows adjudicated:',nrow(queue_decisions)),
 paste('Bounded positive lists retained:',nrow(rules)),paste('Reviewed distinct entity/PT/context features:',nrow(features)),
 paste('Candidate profiles represented:',n_distinct(features$entity_id)),
 paste('Seed common-summary evidence rows still pending in this batch:',sum(q$cohort=='seed')),
 'Inputs: script-15 candidate common-summary queue, full source segments and full term evidence.',
 'Outputs: candidate_common_summary_* review tables and candidate_positive_list_* (local); compact features/reports in outputs/16a_*.',
 'Exact reviewed reaction-list spans prevent indication prefixes, adjacent tables and unreviewed surrounding sentences becoming accepted features.',
 'Six representative profiles use one selected source list each; repeated label/summary mentions are not independent evidence.',
 'Four ambiguous aliases reviewed explicitly: abdominal pain, rash, fatigue and arthralgia use exact source PT names, not alternate LLT-derived concepts.',
 'Accepted lists retain population/regimen qualifiers. None is promoted to confirmed single-agent, translation-sparing or target-binding evidence.',
 'Technetium-pentetate formulation and carfilzomib combination summaries are held; historical SIDER structure/formulation uncertainty persists.',
 'Discontinuation, dose modification, other cohorts and causality-unspecified lists remain outside the selected profiles.',
 'All comparison_ready flags are FALSE: seed/candidate sampling, regimen and route compatibility must be resolved before scoring.',
 'Only this candidate common-summary batch is adjudicated; seed summaries, dictionary ambiguities, narrative/table evidence and sacrosidase context remain pending.',
 'Original extracted evidence and prior reviewed features are unchanged. No ranks, p-values, BH adjustment or eEF1A activity claims generated.',
 'Next: review seed common-summary batches using the same bounded-list approach, then check comparable identity/regimen/route strata.')
writeLines(report,'outputs/16a_candidate_context_review_summary.txt')
writeLines(c('PASS: all candidate common-summary rows adjudicated; exact reviewed source spans present once and accepted mentions wholly contained.',
 'PASS: accepted identity support, unique entity/PT/context features, explicit pending comparison status and unchanged source hashes.'),
 'outputs/16a_candidate_context_validation.txt')
cat(paste0('\nSCRIPT 16a RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/16a_review_candidate_summaries.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
