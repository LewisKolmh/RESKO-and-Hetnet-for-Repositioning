# Review one pinned representative seed passage per common-summary entity.
# Exact bounded lists are accepted as source-context features, never inferred single-agent effects.
library(readr);library(dplyr);library(stringr);library(tibble);library(tidyr)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
paths<-c('data/curated/comparable_label_common_summary_review_queue.csv',
 'data/curated/comparable_label_evidence_segments.csv','data/curated/comparable_label_PT_candidate_evidence.csv',
 'data/curated/current_label_selection_supplemented.csv','data/curated/seed_identity_crosscheck_reviewed.csv',
 'data/curated/current_label_PT_dictionary.csv','config/candidate_PT_mapping_review.csv')
if(!all(file.exists(paths)))stop('Run from project root after 16a.')
before<-tools::md5sum(paths);q<-r(paths[1]);sg<-r(paths[2]);all_hits<-r(paths[3]);labels<-r(paths[4]);seeds<-r(paths[5]);dictionary<-r(paths[6])
norm<-function(z)str_to_lower(str_squish(str_replace_all(z,'[\u2010-\u2015]','-')))
# Source boundaries inspected against the representative passages; exact markers guard against source drift.
spec_text<-'drugbank_id|list_start|list_end|population_scope|regimen_scope
DB00014|hot flushes, headache|and peripheral edema|female_breast_uterine_bleeding_endometriosis|label_summary_regimen_unspecified
DB00072|headache, diarrhea, nausea|and chills|adjuvant_breast_cancer|label_summary_regimen_unspecified
DB00112|epistaxis, headache|and exfoliative dermatitis|pooled_label_indications|label_summary_regimen_unspecified
DB00188|nausea, diarrhea, thrombocytopenia|and anorexia|pooled_clinical_studies|label_summary_regimen_unspecified
DB00262|cerebral edema, asthenia|and depression|newly_diagnosed_high_grade_glioma_local_product|label_summary_regimen_unspecified
DB00351|nausea, diarrhea, impotence|and asthenia|oral_suspension_800mg_per_20mL_trials|label_summary_regimen_unspecified
DB00361|leukopenia, neutropenia, anemia|and peripheral neuropathy|vinorelbine_label_population|label_summary_regimen_unspecified
DB00441|urinary frequency, urinary tract infection|and increased ALT|intravesical_product_only|label_summary_regimen_unspecified
DB00445|leukopenia, neutropenia, anemia|and rash/itch|epirubicin_label_population|label_summary_regimen_unspecified
DB00480|diarrhea, fatigue, anemia|and tremor|multiple_myeloma_label_summary|label_summary_regimen_unspecified
DB00515|nephrotoxicity, peripheral neuropathy|and ototoxicity|cisplatin_label_population|label_summary_regimen_unspecified
DB00531|neutropenia, febrile neutropenia|and diarrhea|cyclophosphamide_trial_context|label_summary_regimen_unspecified
DB00539|hot flashes, sweating|and vaginal discharge|toremifene_label_population|label_summary_regimen_unspecified
DB00563|ulcerative stomatitis, leukopenia|and abdominal distress|methotrexate_label_population|label_summary_regimen_unspecified
DB00783|headache, breast tenderness|and intermenstrual bleeding|MINIVELLE_product_only|label_summary_regimen_unspecified
DB00947|injection site pain (11.6%|bone pain (9.4% of patients)|fulvestrant_500mg_group_only|label_summary_regimen_unspecified
DB00997|asthenia, fatigue, fever|and anemia|selected_doxorubicin_product_only|label_summary_regimen_unspecified
DB01006|hot flashes, arthralgia|and musculoskeletal|letrozole_label_population|label_summary_regimen_unspecified
DB01101|lymphopenia, anemia, diarrhea|and dermatitis|metastatic_breast_cancer|explicit_single_agent
DB01181|alopecia, nausea/vomiting|and infection|ifosfamide_trials|explicit_monotherapy
DB01217|hot flashes, asthenia, arthritis|and lymphedema|women_taking_anastrozole_tablets|label_summary_regimen_unspecified
DB01229|alopecia, neutropenia, sensory neuropathy|and diarrhea|metastatic_breast_albumin_bound_paclitaxel_only|explicit_single_agent
DB01234|cataract, increased intraocular pressure|and conjunctival hemorrhage|selected_ocular_product_only|label_summary_regimen_unspecified
DB01248|infections, neutropenia, anemia|and myalgia|pooled_docetaxel_indications|mixed_or_unspecified_regimens
DB01590|peripheral edema, constipation|and hyperlipidemia|kidney_transplantation_NON_CANCER_context|label_summary_regimen_unspecified
DB04845|peripheral sensory neuropathy|and musculoskeletal pain|ixabepilone_label_population|label_summary_regimen_unspecified
DB05109|nausea, fatigue, vomiting|and headache|trabectedin_label_population|label_summary_regimen_unspecified
DB05773|fatigue, nausea, musculoskeletal pain|and epistaxis|first_KADCYLA_indication_summary_only|label_summary_regimen_unspecified
DB06287|rash, asthenia, mucositis|and anorexia|TORISEL_label_population|label_summary_regimen_unspecified
DB06317|fatigue, diarrhea, pyrexia|pneumonia|elotuzumab_lenalidomide_dexamethasone|explicit_combination_lenalidomide_dexamethasone
DB06374|musculoskeletal pain, nausea|and dyspepsia|ORSERDU_label_population|label_summary_regimen_unspecified
DB06603|diarrhea, fatigue, nausea|and vomiting|panobinostat_clinical_studies|label_summary_regimen_unspecified
DB08871|neutropenia, anemia, asthenia/fatigue|and constipation|metastatic_breast_cancer|label_summary_regimen_unspecified
DB08889|upper respiratory tract infection, fatigue|and edema peripheral|carfilzomib_selected_trial_context|regimen_requires_wider_source_review
DB08910|fatigue and asthenia, neutropenia|and pyrexia|pomalidomide_multiple_myeloma_summary|label_summary_regimen_unspecified
DB09073|decreased neutrophils, decreased hemoglobin|and headache|palbociclib_selected_trial_context|regimen_requires_wider_source_review
DB09074|nausea (60%), fatigue (55%)|thrombocytopenia (10%)|pooled_olaparib_safety_population|explicit_pooled_single_agent
DB09331|upper respiratory infection, neutropenia|and asthenia|daratumumab_pooled_label_summary|label_summary_regimen_unspecified
DB09570|thrombocytopenia, neutropenia, diarrhea|and bronchitis|ixazomib_label_summary|label_summary_regimen_unspecified
DB11652|diarrhea, fatigue, rash|and pyrexia|TUKYSA_trastuzumab_trial_context|explicit_combination_trastuzumab
DB11730|lymphocytes decreased, leukocytes decreased|and fatigue|early_breast_cancer_only|regimen_requires_wider_source_review
DB11760|hemoglobin decreased, neutrophils decreased|and decreased appetite|EMBRACA_label_context|regimen_requires_wider_source_review
DB11793|nausea, thrombocytopenia, anemia|and hypomagnesemia|ZEJULA_label_population|label_summary_regimen_unspecified
DB11828|diarrhea, nausea, abdominal pain|and urinary tract infection|NERLYNX_selected_trial_context|regimen_requires_wider_source_review
DB11942|fatigue, nausea, decreased appetite|and vomiting|multiple_myeloma_XVd_only|explicit_XVd_combination
DB12001|diarrhea, infections, neutropenia|and vomiting|abemaciclib_selected_label_summary|label_summary_regimen_unspecified
DB12015|glucose increased, creatinine increased|and alopecia|PIQRAY_label_summary|regimen_requires_wider_source_review
DB12218|diarrhea, cutaneous adverse reactions|vomiting and stomatitis|HR_positive_HER2_negative_breast_only|regimen_requires_wider_source_review
DB12332|nausea, fatigue (including asthenia)|and leukopenia|ovarian_cancer_only|label_summary_regimen_unspecified
DB12489|increased aspartate aminotransferase, fatigue|and decreased neutrophils|mirvetuximab_label_summary|label_summary_regimen_unspecified
DB12893|decreased leukocyte count (83%)|decreased potassium (26%)|pooled_sacituzumab_safety_population|label_summary_regimen_unspecified
DB13944|hematocrit increased, hypertension|and headache|XYOSTED_product_only|label_summary_regimen_unspecified
DB14811|upper respiratory tract infection, fatigue|and diarrhea|isatuximab_selected_trial_context|regimen_requires_wider_source_review
DB14962|decreased white blood cell count (73%)|musculoskeletal pain (24%)|pooled_ENHERTU_safety_population|label_summary_regimen_unspecified
DB14976|decreased hemoglobin, decreased neutrophils|and decreased appetite|relacorilant_selected_trial_context|regimen_requires_wider_source_review
DB15275|decreased neutrophils, decreased hemoglobin|and headache|ITOVEBI_selected_label_summary|regimen_requires_wider_source_review
DB15395|CRS, fatigue, injection site reaction|and pyrexia|elranatamab_label_summary|label_summary_regimen_unspecified
DB15719|reduction in BCVA, corneal exam findings|and COVID-19|BLENREP_selected_trial_context|regimen_requires_wider_source_review
DB16410|stomatitis, nausea, alopecia|and rash|EGFR_mutated_NSCLC_only|label_summary_regimen_unspecified
DB16655|hypogammaglobulinemia, upper respiratory tract infection|and weight decreased|teclistamab_selected_trial_context|label_summary_regimen_unspecified
DB16678|pyrexia, CRS, dysgeusia|and headache|talquetamab_label_summary|label_summary_regimen_unspecified
DB17447|musculoskeletal pain, cytokine release syndrome|and dyspnea|linvoseltamab_label_summary|label_summary_regimen_unspecified
DB19006|decreased white blood cells, increased AST|and constipation|VEPPANU_label_summary|label_summary_regimen_unspecified
DB19043|decreased hemoglobin, musculoskeletal pain|and abdominal pain|imlunestrant_selected_trial_context|regimen_requires_wider_source_review'
spec<-read_delim(I(spec_text),delim='|',col_types=cols(.default=col_character()),trim_ws=TRUE)
holds<-tribble(~drugbank_id,~hold_reason,
 'DB00286','Summary points to a table without listing reactions; table review required',
 'DB00570','Underlining determines common reactions; plain text loses this formatting',
 'DB00990','Sentence segmentation truncates the comparator expression; wider source review required',
 'DB01259','Selected passage concerns discontinuation and Grade 3/4 effects; not a general common-reaction profile',
 'DB06366','Selected passage is Grade 3/4 only; not a general all-grade profile',
 'DB09037','Selected passage concerns permanent discontinuation only; endpoint-specific context held')
representatives<-q |> filter(cohort=='seed') |> distinct(entity_id,entity_name,label_id,source_segment_id) |>
 left_join(sg,by=c('source_segment_id','label_id')) |> filter(source_field=='adverse_reactions') |>
 group_by(entity_id) |> slice_min(as.integer(segment_index),n=1,with_ties=FALSE) |> ungroup() |>
 mutate(drugbank_id=sub('^seed:','',entity_id))
# Pin source keys on first execution; later runs must use exactly the reviewed representative keys.
config_path<-'config/seed_positive_list_review_rules.csv'
rules<-spec |> left_join(representatives |> select(drugbank_id,entity_id,label_id,source_segment_id),by='drugbank_id')
stopifnot(!anyNA(rules$source_segment_id),setequal(representatives$drugbank_id,c(spec$drugbank_id,holds$drugbank_id)),
 !anyDuplicated(rules$drugbank_id))
if(file.exists(config_path)){
 previous<-r(config_path);stopifnot(identical(previous$source_segment_id,rules$source_segment_id))
}
write_csv(rules,config_path,na='');write_csv(holds,'config/seed_summary_context_holds.csv',na='')
resolution<-r(paths[7]) |> select(alias_norm,reviewed_pt_cui)
# Hot-flush variant was previously adjudicated against supplied term/CUI rows in 11g; revalidate current mappings.
hot<-tibble(alias_norm=c('hot flashes','hot flushes','hot flush'),reviewed_pt_cui='C0600142')
stopifnot(all(hot$alias_norm %in% dictionary$alias_norm[dictionary$pt_cui=='C0600142']))
resolution<-bind_rows(resolution,hot) |> distinct()
write_csv(resolution,'config/seed_PT_mapping_review.csv')
reviewed<-list();span_rows<-list()
for(i in seq_len(nrow(rules))){
 a<-rules[i,];s<-representatives |> filter(source_segment_id==a$source_segment_id,entity_id==a$entity_id)
 stopifnot(nrow(s)==1L)
 text<-norm(s$source_text);start<-str_locate(text,fixed(norm(a$list_start)))[1]
 stopifnot(!is.na(start))
 tail<-str_sub(text,start);ending<-str_locate(tail,fixed(norm(a$list_end)))[2]
 stopifnot(!is.na(ending));end<-start+ending-1L
 h<-all_hits |> filter(entity_id==a$entity_id,source_segment_id==a$source_segment_id,
 as.integer(match_start)>=start,as.integer(match_end)<=end)
 stopifnot(nrow(h)>0L,all(h$identity_eligible=='TRUE'))
 source_label<-labels |> filter(drugbank_id==a$drugbank_id,label_id==a$label_id)
 stopifnot(nrow(source_label)==1L)
 tier<-seeds |> filter(drugbank_id==a$drugbank_id)
 h<-h |> left_join(resolution,by='alias_norm') |>
 mutate(review_decision=case_when(broad_term=='TRUE'~'hold_dictionary_scope',
  ambiguous_dictionary_mapping=='TRUE' & !is.na(reviewed_pt_cui) & pt_cui==reviewed_pt_cui~'accept_positive_list_context',
  ambiguous_dictionary_mapping=='TRUE' & !is.na(reviewed_pt_cui)~'reject_alternate_dictionary_mapping',
  ambiguous_dictionary_mapping=='TRUE'~'hold_dictionary_scope',TRUE~'accept_positive_list_context'),
  review_reason='Mention wholly within inspected representative positive list; unresolved dictionary scope remains held',
  regimen_scope=a$regimen_scope,population_scope=a$population_scope,routes=source_label$routes,
  seed_tiers=tier$seed_tiers,cancers=tier$cancers,
  review_method='assistant_reviewed_exact_source_list_boundaries',review_batch='16b_seed_representatives',
  comparison_ready=FALSE,evidence_group=paste0('SPL:',set_id))
 reviewed[[i]]<-h
 span_rows[[i]]<-tibble(entity_id=a$entity_id,source_segment_id=a$source_segment_id,
 accepted_span_start=start,accepted_span_end=end,population_scope=a$population_scope,regimen_scope=a$regimen_scope)
}
e<-bind_rows(reviewed);spans<-bind_rows(span_rows)
features<-e |> filter(review_decision=='accept_positive_list_context') |>
 distinct(entity_id,entity_name,drugbank_id,pt_cui,pt_name,set_id,label_id,effective_time,label_version,source_url,
 identity_scope,comparison_layer,regimen_scope,population_scope,routes,seed_tiers,cancers,evidence_group,
 review_method,review_batch,comparison_ready) |> mutate(present=1L,feature_layer='reviewed_seed_positive_list_context')
representative_review<-representatives |> select(entity_id,entity_name,drugbank_id,source_segment_id,label_id) |>
 left_join(spans,by=c('entity_id','source_segment_id')) |> left_join(holds,by='drugbank_id') |>
 mutate(segment_review_decision=if_else(!is.na(accepted_span_start),'accept_bounded_list_only','hold_source_context'),
 reason=coalesce(hold_reason,'Only mentions inside the inspected positive-list span accepted; other text is not adjudicated'))
# Every queued row has an explicit workflow status, but only inspected representatives are called reviewed.
queue<-q |> filter(cohort=='seed') |>
 left_join(representative_review |> select(entity_id,source_segment_id,segment_review_decision,reason),
 by=c('entity_id','source_segment_id')) |>
 left_join(e |> select(evidence_key,bounded_decision=review_decision),by='evidence_key') |>
 mutate(batch_review_decision=case_when(!is.na(bounded_decision)~bounded_decision,
  segment_review_decision=='hold_source_context'~'hold_source_context',
  segment_review_decision=='accept_bounded_list_only'~'not_reviewed_outside_accepted_span',
  TRUE~'pending_nonrepresentative_source_review'))
stopifnot(!anyDuplicated(features[c('entity_id','pt_cui','regimen_scope','population_scope')]),
 all(!features$comparison_ready),!any(features$drugbank_id %in% holds$drugbank_id),
 nrow(queue)==sum(q$cohort=='seed'),!anyNA(queue$batch_review_decision),identical(before,tools::md5sum(paths)))
write_csv(e,'data/curated/seed_positive_list_reviewed_evidence.csv',na='')
write_csv(features,'data/curated/seed_positive_list_PT_features.csv',na='')
write_csv(representative_review,'data/curated/seed_representative_summary_review.csv',na='')
write_csv(queue,'data/curated/seed_common_summary_batch_review.csv',na='')
write_csv(features,'outputs/16b_seed_PT_features.csv',na='')
write_csv(features |> count(entity_id,entity_name,regimen_scope,population_scope,routes,seed_tiers,name='reviewed_PT_terms'),
 'outputs/16b_seed_profile_coverage.csv',na='')
write_csv(queue |> count(batch_review_decision,name='evidence_rows'),'outputs/16b_seed_review_workload_summary.csv')
write_csv(representative_review |> select(entity_id,entity_name,segment_review_decision,reason),
 'outputs/16b_representative_review_summary.csv',na='')
write_csv(tibble(file=paths,md5=unname(before)),'outputs/16b_input_summary.csv')
report<-c('SCRIPT 16b - SEED REPRESENTATIVE COMMON-SUMMARY REVIEW',
 paste('Representative seed passages inspected:',nrow(representative_review)),
 paste('Positive source-list profiles retained:',n_distinct(features$entity_id)),
 paste('Representative profiles held for source context:',nrow(holds)),
 paste('Reviewed entity/PT/context features:',nrow(features)),
 paste('Seed common-summary evidence rows tracked:',nrow(queue)),
 paste('Nonrepresentative evidence rows still pending:',sum(queue$batch_review_decision=='pending_nonrepresentative_source_review')),
 'Inputs: script-15 seed summary queue/segments/evidence; label route metadata, seed tiers and supplied dictionary.',
 'Outputs locally: seed_positive_list_* tables, seed_representative_summary_review.csv, seed_common_summary_batch_review.csv.',
 'Public compact feature, coverage, representative decisions, workload, hashes and validation: outputs/16b_*.',
 'One representative passage per seed with a queued narrative summary; no uninspected duplicate is labelled reviewed.',
 'Exact list boundaries exclude indication prefixes, warning lists, other cohorts, adjacent tables and laboratory summaries beyond the chosen list.',
 'Formulation/population labels and original routes are retained: intravesical gemcitabine, ocular dexamethasone, albumin-bound paclitaxel and non-cancer transplant everolimus remain distinct.',
 'Six source-context holds require tables, lost underlining, comparator repair or endpoint/grade-specific review.',
 'Four exact PT-name decisions plus previously reviewed hot-flush variant mapping are explicit; other ambiguous/broad dictionary hits stay held.',
 'All comparison_ready flags are FALSE. Regimen-unspecified lists are not inferred single-agent; no ranking or eEF1A binding claim is produced.',
 'Existing seed/candidate identity, original SIDER data and prior reviewed feature layers are unchanged.',
 'Next: audit seed/candidate identity, route, regimen and list coverage to define comparison strata; retain pending narrative/table and named-biologic evidence separately.')
writeLines(report,'outputs/16b_seed_context_review_summary.txt')
writeLines(c('PASS: exact representative source keys and boundary markers; accepted hits contained in reviewed spans.',
 'PASS: unique feature keys, valid source routes and seed tiers, explicit holds/pending statuses, unchanged source input hashes.'),
 'outputs/16b_seed_context_validation.txt')
cat(paste0('\nSCRIPT 16b RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/16b_review_seed_summaries.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
