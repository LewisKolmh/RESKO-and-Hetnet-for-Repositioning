library(readr);library(dplyr);library(stringr)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
f<-r('outputs/16b_seed_PT_features.csv');e<-r('data/curated/seed_positive_list_reviewed_evidence.csv')
s<-r('data/curated/seed_representative_summary_review.csv');g<-r('data/curated/comparable_label_evidence_segments.csv')
labels<-r('data/curated/current_label_selection_supplemented.csv');hash<-r('outputs/16b_input_summary.csv')
stopifnot(identical(unname(tools::md5sum(hash$file)),hash$md5),!anyDuplicated(f[c('entity_id','pt_cui','regimen_scope','population_scope')]),
 all(f$comparison_ready=='FALSE'),all(f$review_batch=='16b_seed_representatives'))
accepted<-e |> filter(review_decision=='accept_positive_list_context') |>
 left_join(s |> select(entity_id,source_segment_id,accepted_span_start,accepted_span_end),by=c('entity_id','source_segment_id')) |>
 left_join(g |> select(source_segment_id,source_text),by='source_segment_id')
text<-str_to_lower(str_squish(str_replace_all(accepted$source_text,'[\u2010-\u2015]','-')))
stopifnot(all(as.integer(accepted$match_start)>=as.integer(accepted$accepted_span_start)),
 all(as.integer(accepted$match_end)<=as.integer(accepted$accepted_span_end)),
 all(str_sub(text,as.integer(accepted$match_start),as.integer(accepted$match_end))==accepted$alias_norm),
 !nrow(anti_join(f |> select(entity_id,pt_cui),accepted |> select(entity_id,pt_cui),by=c('entity_id','pt_cui'))))
route_check<-f |> distinct(drugbank_id,label_id,routes) |> left_join(labels |> select(drugbank_id,label_id,source_routes=routes),by=c('drugbank_id','label_id'))
stopifnot(all((is.na(route_check$routes)&is.na(route_check$source_routes))|coalesce(route_check$routes==route_check$source_routes,FALSE)))
cat('PASS: every accepted mention independently checked against full source text and approved span; exported PTs trace to accepted evidence.\nPASS: route metadata agrees with original selected labels; protected input hashes unchanged and all comparison-ready flags remain false.\n',file='outputs/16b_seed_context_validation.txt',append=TRUE)
cat('\nSCRIPT 16b - INDEPENDENT VALIDATION\nEvery accepted mention was independently checked against its full source text and approved boundaries. Exported PT features trace to accepted evidence; original route metadata and input hashes agree. No profile is comparison-ready.\nOutput: outputs/16b_seed_context_validation.txt.\n',file='Notes.txt',append=TRUE)
