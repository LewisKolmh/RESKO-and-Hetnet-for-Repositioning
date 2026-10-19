# Independent raw-feature and source-route checks for section 16c.
library(readr);library(dplyr);library(tidyr)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
p<-r('outputs/16c_profile_compatibility_summary.csv');pairs<-r('outputs/16c_pair_compatibility_summary.csv')
f<-bind_rows(r('outputs/16a_candidate_PT_features.csv'),r('outputs/16b_seed_PT_features.csv'))
hash<-r('outputs/16c_input_summary.csv');coverage<-r('outputs/16c_full_inventory_coverage.csv')
stopifnot(identical(unname(tools::md5sum(hash$file)),hash$md5),!anyDuplicated(p$profile_id),nrow(coverage)==200,
 !anyDuplicated(coverage$entity_id))
counts<-f |> group_by(entity_id,set_id,label_id,population_scope,regimen_scope) |>
 summarise(independent_terms=n_distinct(pt_cui),.groups='drop')
checked<-p |> left_join(counts,by=c('entity_id','set_id','label_id','population_scope','regimen_scope'))
stopifnot(all(as.integer(checked$reviewed_PT_terms)==checked$independent_terms))
expected<-expand_grid(seed_profile_id=p$profile_id[p$cohort=='seed'],candidate_profile_id=p$profile_id[p$cohort=='candidate'])
stopifnot(!nrow(anti_join(expected,pairs,by=c('seed_profile_id','candidate_profile_id'))),
 !nrow(anti_join(pairs,expected,by=c('seed_profile_id','candidate_profile_id'))),
 !anyDuplicated(pairs[c('seed_profile_id','candidate_profile_id')]))
labels<-r('data/curated/candidate_label_identity_review.csv')
routes<-p |> filter(cohort=='candidate') |> left_join(labels |> select(stitch_flat_id,label_id,original_routes=routes,raw_supported_identity=supported_identity),
 by=c('stitch_flat_id','label_id'))
stopifnot(all(routes$routes==routes$original_routes),length(routes$raw_supported_identity)==nrow(routes),all(routes$raw_supported_identity=='TRUE'),
 all(pairs$strict_pair_ready=='FALSE'),all(nzchar(pairs$strict_stratum_barriers)))
msg<-c('PASS: independent feature counts agree for every profile; complete seed/candidate pair coverage checked both directions.',
 'PASS: candidate route and identity support agree with original label records; unique full-inventory rows and input hashes confirmed.',
 'Scope: no ready strict pair is a coverage/context finding, not negative pharmacological evidence.')
cat(paste0(paste(msg,collapse='\n'),'\n'),file='outputs/16c_profile_compatibility_validation.txt',append=TRUE)
cat(paste0('\nSCRIPT 16c - INDEPENDENT VALIDATION\n',paste(msg,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
