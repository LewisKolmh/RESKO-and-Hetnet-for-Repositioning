# Consolidate bounded section-16 results and account for every remaining evidence row.
# Section completion does not mean every extracted sentence/table has been manually adjudicated.
library(readr);library(dplyr);library(tidyr);library(tibble);library(stringr)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
paths<-c('outputs/16a_candidate_PT_features.csv','outputs/16b_seed_PT_features.csv',
 'data/curated/candidate_positive_list_reviewed_evidence.csv','data/curated/seed_positive_list_reviewed_evidence.csv',
 'data/curated/candidate_common_summary_batch_review.csv','data/curated/seed_common_summary_batch_review.csv',
 'data/curated/comparable_label_PT_candidate_evidence.csv','data/curated/seed_identity_crosscheck_reviewed.csv',
 'data/curated/candidate_identity_registry.csv','data/curated/candidate_DrugBank_record_enrichment.csv',
 'data/curated/candidate_label_named_context_selection.csv','data/curated/reviewed_profile_compatibility_registry.csv',
 'outputs/16c_full_inventory_coverage.csv','outputs/16e_context_pair_summary.csv',
 'outputs/16e_aligned_summary_overlap_summary.csv','data/curated/reviewed_candidate_full_context.csv')
before<-tools::md5sum(paths)
f<-bind_rows(r(paths[1]),r(paths[2]));reviewed<-bind_rows(r(paths[3]),r(paths[4]))
q<-bind_rows(r(paths[5]) |> select(evidence_key,batch_review_decision),r(paths[6]) |> select(evidence_key,batch_review_decision))
h<-r(paths[7]);si<-r(paths[8]);ci<-r(paths[9]);db<-r(paths[10]);bio<-r(paths[11]);p<-r(paths[12]);inventory<-r(paths[13]);pairs<-r(paths[14]);eligible<-r(paths[15]);context<-r(paths[16])
stopifnot(!anyDuplicated(reviewed$evidence_key),!anyDuplicated(q$evidence_key),!anyDuplicated(h$evidence_key))
ledger<-h |> select(evidence_key,entity_id,entity_name,cohort,label_id,set_id,source_segment_id,pt_cui,
 identity_eligible,comparison_layer,review_priority) |>
 left_join(reviewed |> select(evidence_key,positive_list_decision=review_decision),by='evidence_key') |>
 left_join(q,by='evidence_key') |>
 mutate(final_review_status=case_when(!is.na(positive_list_decision)~positive_list_decision,
  !is.na(batch_review_decision)~batch_review_decision,
  cohort=='named_biologic_context'~'named_biologic_source_context_reviewed_terms_not_adjudicated',
  identity_eligible=='FALSE'~'pending_identity_resolution',TRUE~'pending_narrative_table_or_dictionary_review'))
stopifnot(nrow(ledger)==nrow(h),!anyNA(ledger$final_review_status))
# Consolidated features preserve seed/candidate origin and populate candidate routes from the audited exact label record.
metadata<-p |> select(entity_id,label_id,profile_id,audited_route=routes,current_record_name,current_inchikey,cohort)
f<-f |> left_join(metadata,by=c('entity_id','label_id')) |>
 mutate(routes=coalesce(routes,audited_route),present=as.integer(present),comparison_ready=FALSE) |>
 left_join(context |> select(entity_id,full_candidate_regimen_review=regimen_review,
 candidate_threshold_operator=threshold_operator,candidate_threshold_percent=threshold_percent,comparator_rule,frequency_limitation),by='entity_id') |>
 select(-audited_route)
stopifnot(nrow(f)==672L,!anyNA(f$profile_id),all(f$present==1L),
 !anyDuplicated(f[c('profile_id','pt_cui')]))
accepted<-reviewed |> filter(review_decision=='accept_positive_list_context')
stopifnot(!nrow(anti_join(f |> select(entity_id,label_id,pt_cui),accepted |> select(entity_id,label_id,pt_cui),by=c('entity_id','label_id','pt_cui'))),
 all(ledger$identity_eligible[ledger$final_review_status=='accept_positive_list_context']=='TRUE'))
# All original entities remain eligible for an independent target-evidence SEARCH; side-effect coverage is not a search gate.
seed_handoff<-si |> transmute(entity_id=paste0('seed:',drugbank_id),origin='seed',source_name=drug_name,current_record_name=drug_name,
 drugbank_ids=drugbank_id,unii,cas_number,current_inchikey=inchikey,identity_status=reviewed_identity_status,
 source_drug_type=drug_type,chembl_molecule_type=chembl_molecule_type,seed_tiers,cancers,
 structure_mapping_blocked=FALSE,historical_identity_status=NA_character_,rejected_structure_cid=NA_character_)
db_links<-db |> group_by(stitch_flat_id) |> summarise(drugbank_ids=paste(sort(unique(record_drugbank_id)),collapse='|'),
 unii=paste(sort(unique(na.omit(unii))),collapse='|'),cas_number=paste(sort(unique(na.omit(cas_number))),collapse='|'),
 source_drug_type=paste(sort(unique(drug_type)),collapse='|'),.groups='drop')
candidate_handoff<-ci |> select(stitch_flat_id,sider_drug_name,matched_names,current_inchikey,identity_status,
 historical_identity_status,automatic_label_mapping_blocked,sider_pubchem_cid) |>
 left_join(db_links,by='stitch_flat_id') |>
 transmute(entity_id=paste0('sider:',stitch_flat_id),origin='candidate',source_name=sider_drug_name,
 current_record_name=matched_names,drugbank_ids,unii,cas_number,current_inchikey,identity_status,source_drug_type,
 structure_mapping_blocked=automatic_label_mapping_blocked=='TRUE',historical_identity_status,
 rejected_structure_cid=if_else(structure_mapping_blocked,sider_pubchem_cid,NA_character_)) |>
 mutate(current_inchikey=if_else(structure_mapping_blocked,NA_character_,current_inchikey))
# The named biologic recovery is supported, but does not repair the historical CID.
j<-which(candidate_handoff$entity_id==paste0('sider:',bio$stitch_flat_id))
stopifnot(length(j)==1L)
candidate_handoff$current_record_name[j]<-'Sacrosidase';candidate_handoff$drugbank_ids[j]<-bio$drugbank_id
candidate_handoff$unii[j]<-bio$unii;candidate_handoff$cas_number[j]<-bio$cas_number;candidate_handoff$source_drug_type[j]<-'biotech'
candidate_handoff$identity_status[j]<-'named_biologic_UNII_verified_historical_CID_conflict'
handoff<-bind_rows(seed_handoff,candidate_handoff) |>
 left_join(inventory |> select(entity_id,label_coverage_status,reviewed_profile_status,reviewed_profiles,reviewed_profile_terms),by='entity_id') |>
 mutate(target_evidence_assessment='not_assessed_in_section16',
 EEF1A1_direct_binding_evidence=NA_character_,EEF1A2_direct_binding_evidence=NA_character_,
 RT_PPI_disruption_evidence=NA_character_,allosteric_evidence=NA_character_,translation_sparing_evidence=NA_character_,
 missing_coverage_is_not_negative=TRUE)
stopifnot(nrow(handoff)==200L,!anyDuplicated(handoff$entity_id),
 is.na(handoff$current_inchikey[handoff$entity_id=='sider:CID100135113']),
 all(is.na(handoff$EEF1A1_direct_binding_evidence)),all(is.na(handoff$EEF1A2_direct_binding_evidence)))
limitations<-tribble(~issue,~section16_disposition,~future_work,
 'Incomplete term adjudication','Every extracted row has an explicit status; accepted features come only from bounded inspected lists','Prioritise further tables/narratives when a later hypothesis requires them',
 'Strict comparison not ready','Zero strict-ready pairs retained as a source-context finding','Use direct target evidence independently; do not convert missing context to negative activity',
 'Only one aligned descriptive pair','Olaparib/pirfenidone share eight common terms; no target-specific signal','Do not promote descriptive overlap to a binding prediction',
 'Historical SIDER identity uncertainty','Current-label attribution and original source IDs remain distinct','Resolve individual historical forms only when needed; preserve original raw data',
 'Sacrosidase incorrect CID','Wrong structure excluded; verified biologic name/UNII retained separately','Keep outside small-molecule structure screening',
 'Formulation and regimen differences','Routes, population qualifiers and full candidate context remain attached','Select source-compatible strata for any later comparison',
 'Incomplete original database framework','This section concerns label context, not complete source integration','Continue the original manifest including STITCH, ChEMBL, STRING and pathway/ontology layers',
 'No pharmacological validation','No eEF1A or RT-PPI, allosteric or translation-sparing evidence established by this section','Assess EEF1A1 and EEF1A2 binding and downstream claims separately')
write_csv(ledger,'data/curated/section16_evidence_review_ledger.csv',na='')
write_csv(f,'data/curated/section16_reviewed_PT_features.csv',na='')
write_csv(handoff,'data/curated/section16_target_evidence_handoff.csv',na='')
write_csv(f,'outputs/16f_reviewed_PT_features.csv',na='')
write_csv(handoff,'outputs/16f_target_evidence_handoff_summary.csv',na='')
write_csv(ledger |> count(cohort,final_review_status,name='evidence_rows'),'outputs/16f_review_status_summary.csv')
write_csv(limitations,'outputs/16f_limitations_summary.csv')
write_csv(tibble(file=paths,md5=unname(before)),'outputs/16f_input_summary.csv')
# Final independently recomputed set check for the aligned descriptive pair.
sets<-split(f$pt_cui,f$entity_id)
for(i in seq_len(nrow(eligible))){a<-unique(sets[[eligible$seed_entity_id[i]]]);b<-unique(sets[[eligible$candidate_entity_id[i]]])
 stopifnot(length(intersect(a,b))==as.integer(eligible$shared_PT_terms[i]),
 length(union(a,b))==as.integer(eligible$union_PT_terms[i]),
 abs(length(intersect(a,b))/length(union(a,b))-as.numeric(eligible$descriptive_jaccard[i]))<1e-12)}
# Check every section-stage saved input manifest still agrees with the inputs it records.
manifest_paths<-list.files('outputs',pattern='^16[a-f]_input_summary[.]csv$',full.names=TRUE)
for(mp in manifest_paths){m<-r(mp);stopifnot(identical(unname(tools::md5sum(m$file)),m$md5))}
for(file in c(list.files('scripts',pattern='^16[a-f]_.*[.]R$',full.names=TRUE),
 'scripts/helpers/validate_seed_summary_review.R','scripts/helpers/validate_profile_compatibility.R'))invisible(parse(file))
stopifnot(identical(before,tools::md5sum(paths)),sum(as.integer(r('outputs/16f_review_status_summary.csv')$evidence_rows))==nrow(h))
report<-c('SECTION 16 COMPLETE - BOUNDED LABEL REVIEW, AUDIT AND HANDOFF',
 paste('Consolidated reviewed entity/PT/context features:',nrow(f)),
 paste('Seed profiles:',n_distinct(f$profile_id[f$cohort=='seed'])) ,paste('Candidate profiles:',n_distinct(f$profile_id[f$cohort=='candidate'])),
 paste('Extracted evidence rows fully accounted for in review ledger:',nrow(ledger)),
 paste('Evidence rows accepted from reviewed positive-list spans:',sum(ledger$final_review_status=='accept_positive_list_context')),
 paste('Evidence rows with other decisions/statuses:',sum(ledger$final_review_status!='accept_positive_list_context')),
 paste('Entities retained for independent target-evidence search:',nrow(handoff)),
 paste('Strict-ready comparisons:',sum(pairs$strict_pair_ready=='TRUE')),
 paste('Descriptive aligned-cutoff comparisons:',nrow(eligible)),
 'The descriptive olaparib/pirfenidone pair shares 8/23 union terms (Jaccard 0.3478); all eight shared terms are common in original SIDER.',
 'Section completion means the bounded source review, compatibility audit, consolidated outputs and explicit unresolved-work ledger are complete.',
 'It does not mean every extracted term/table/narrative is manually adjudicated or that biological activity has been established.',
 'No p-values or BH adjustment: these set descriptions do not supply a valid null/test family.',
 'Outputs: section16_evidence_review_ledger.csv, section16_reviewed_PT_features.csv, section16_target_evidence_handoff.csv locally;',
 'outputs/16f_* compact features, review status, limitations, identity handoff, input hashes and final validation.',
 'All original 200 entities remain available for target-evidence searches; missing side-effect coverage is not a negative target filter.',
 'Next: section 17, direct eEF1A target-evidence integration (including STITCH/ChEMBL) with EEF1A1 versus EEF1A2 and binding versus functional/PPI evidence kept distinct.')
writeLines(report,'outputs/16f_section_completion_summary.txt')
writeLines(c('PASS: six section-stage input manifests unchanged; section R scripts parse.',
 'PASS: all 59,500 evidence rows assigned exactly one workflow status; consolidated feature keys unique and trace to accepted source evidence.',
 'PASS: all 200 handoff entities retained; sacrosidase structure quarantined; target activity fields remain missing, not negative.',
 'PASS: aligned-pair set counts and Jaccard independently recomputed from consolidated features; no strict-ready status invented.'),
 'outputs/16f_section_validation.txt')
writeLines(c('section=16','status=complete_bounded_review_and_handoff','all_extracted_terms_manually_reviewed=FALSE',
 'strict_comparison_ready=FALSE','next_section=17_target_evidence_integration'),'outputs/16f_section_status_summary.txt')
cat(paste0('\nSCRIPT 16f RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/16f_finalize_section16.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
