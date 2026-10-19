# Independent cross-table checks and full-inventory coverage; run after script 15.
library(readr);library(dplyr);library(tidyr);library(stringr);library(tibble)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
registry<-r('data/curated/comparable_label_registry.csv')
h<-r('data/curated/comparable_label_PT_candidate_evidence.csv')
segments<-r('data/curated/comparable_label_evidence_segments.csv')
dict<-r('data/curated/current_label_PT_dictionary.csv')
cov<-r('outputs/15_comparable_label_coverage.csv');hash<-r('outputs/15_input_summary.csv')
stopifnot(identical(unname(tools::md5sum(hash$file)),hash$md5),!anyDuplicated(h$evidence_key),
 !anyDuplicated(segments$source_segment_id),all(h$source_segment_id %in% segments$source_segment_id),
 !nrow(anti_join(h |> distinct(alias_norm,pt_cui),dict |> distinct(alias_norm,pt_cui),by=c('alias_norm','pt_cui'))))
check<-h |> select(source_segment_id,match_start,match_end,alias_norm) |>
 left_join(segments |> select(source_segment_id,source_text),by='source_segment_id',relationship='many-to-one')
text<-str_to_lower(str_squish(str_replace_all(check$source_text,'[\u2010-\u2015]','-')))
stopifnot(all(str_sub(text,as.integer(check$match_start),as.integer(check$match_end))==check$alias_norm),
 all(h$identity_eligible[h$comparison_layer %in% c('candidate_identity_unresolved','named_biologic_context_only')]=='FALSE'),
 all(h$extraction_status=='candidate_requires_context_review'),all(is.na(h$review_decision)),
 sum(as.integer(cov$evidence_rows))==nrow(h))
seeds<-r('data/curated/seed_identity_crosscheck_reviewed.csv')
candidates<-r('data/curated/candidate_identity_registry.csv')
all_entities<-bind_rows(seeds |> transmute(entity_id=paste0('seed:',drugbank_id),entity_name=drug_name,inventory='seed'),
 candidates |> transmute(entity_id=paste0('sider:',stitch_flat_id),entity_name=sider_drug_name,inventory='candidate'))
entity_cov<-cov |> group_by(entity_id) |>
 summarise(comparison_layers=paste(sort(unique(comparison_layer)),collapse='|'),
  across(c(label_records,source_segments,literal_mentions,evidence_rows),~sum(as.integer(.x))),.groups='drop') |>
 left_join(h |> group_by(entity_id) |> summarise(candidate_PT_terms=n_distinct(pt_cui),.groups='drop'),by='entity_id')
all_cov<-all_entities |> left_join(entity_cov,by='entity_id') |>
 mutate(label_coverage_status=case_when(is.na(label_records)~'no_collected_label',as.integer(evidence_rows)==0L~'collected_label_no_dictionary_hit',TRUE~'dictionary_candidates_require_review'),
  across(c(label_records,source_segments,literal_mentions,candidate_PT_terms,evidence_rows),~replace_na(as.integer(.x),0L)))
stopifnot(nrow(all_cov)==nrow(all_entities),!anyDuplicated(all_cov$entity_id))
write_csv(all_cov,'outputs/15_full_inventory_coverage.csv',na='')
msg<-c('PASS: every evidence span rechecked against its full source segment and dictionary/PT link.',
 'PASS: all unresolved and named-biologic contexts remain identity-ineligible; all review decisions remain unassigned.',
 'PASS: input checksums and cross-table evidence counts agree; no unmatched segment keys.',
 paste('Full inventory coverage rows:',nrow(all_cov)),paste('Entities without collected labels:',sum(all_cov$label_coverage_status=='no_collected_label')),
 'Zero counts describe extraction/collection coverage only, never an absence of adverse effects or target activity.')
cat(paste0(paste(msg,collapse='\n'),'\n'),file='outputs/15_comparable_label_validation.txt',append=TRUE)
cat(paste0('\nSCRIPT 15 - INDEPENDENT VALIDATION AND INVENTORY COVERAGE\n',paste(msg,collapse='\n'),
 '\nOutput: outputs/15_full_inventory_coverage.csv (all seed and queued SIDER entities, including missing-label coverage).\n'),file='Notes.txt',append=TRUE)
cat(paste(msg,collapse='\n'),'\n')
