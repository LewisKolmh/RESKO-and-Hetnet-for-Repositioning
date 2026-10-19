# Shared extraction for all collected seed/candidate label associations.
# Identity/context eligibility remains separate from dictionary matching.
library(readr);library(dplyr);library(tidyr);library(stringr);library(tibble);library(xml2)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
inputs<-c('data/curated/current_label_selection_supplemented.csv','data/curated/candidate_label_identity_review.csv',
 'data/curated/candidate_label_named_context_selection.csv','data/curated/current_label_PT_dictionary.csv')
if(!all(file.exists(inputs)))stop('Run from project root after 14c and the existing dictionary extraction.')
before<-tools::md5sum(inputs)
a<-r(inputs[1]);b<-r(inputs[2]);c<-r(inputs[3]);dict<-r(inputs[4])
selection<-bind_rows(
 a |> transmute(entity_id=paste0('seed:',drugbank_id),cohort='seed',entity_name=drug_name,
  drugbank_id,stitch_flat_id=NA_character_,label_id,set_id,effective_time,label_version,cache_file,
  source_url=label_url,identity_scope=label_identity_status,identity_eligible=TRUE,comparison_layer='seed_selected_label'),
 b |> transmute(entity_id=paste0('sider:',stitch_flat_id),cohort='candidate',entity_name=sider_drug_name,
  drugbank_id=NA_character_,stitch_flat_id,label_id,set_id,effective_time,label_version,cache_file,
  source_url=url,identity_scope=label_identity_status,identity_eligible=supported_identity=='TRUE',
  comparison_layer=if_else(supported_identity=='TRUE','candidate_identity_supported','candidate_identity_unresolved')),
 c |> transmute(entity_id=paste0('sider:',stitch_flat_id),cohort='named_biologic_context',entity_name=sider_drug_name,
  drugbank_id,stitch_flat_id,label_id,set_id,effective_time,label_version,cache_file,
  source_url=label_url,identity_scope=mapping_scope,identity_eligible=FALSE,comparison_layer='named_biologic_context_only')
) |> distinct(entity_id,label_id,.keep_all=TRUE)
stopifnot(!anyDuplicated(selection[c('entity_id','label_id')]),all(file.exists(selection$cache_file)))
selection$cache_md5<-unname(tools::md5sum(selection$cache_file))
norm<-function(z)str_to_lower(str_squish(str_replace_all(z,'[\u2010-\u2015]','-')))
terms<-unique(dict$alias_norm);terms<-terms[order(nchar(terms),decreasing=TRUE)]
pattern<-regex(paste0('(?<![\\p{L}\\p{N}])(?:',paste(str_escape(terms),collapse='|'),')(?![\\p{L}\\p{N}])'))
flags<-function(text){
 t<-norm(text)
 tibble(negation_flag=str_detect(t,'\\b(no|not|without|absence of|neither|none|did not)\\b'),
 comparator_flag=str_detect(t,'\\b(placebo|comparator|control arm|control group)\\b'),
 combination_flag=str_detect(t,'\\b(combination|combined|plus|coadministered|co-administered)\\b'),
 grade_frequency_flag=str_detect(t,'grade [0-9]|%|percent|incidence'),
 common_summary_flag=str_detect(t,'most common.{0,80}adverse (reactions|events)|common adverse reactions|most frequently reported adverse'),
 indication_context_flag=str_detect(t,'\\b(indicated|indication|treatment of|patients with|disease progression)\\b'))
}
stopifnot(!str_detect('painting',pattern),flags('No nausea in placebo patients')$negation_flag,
 flags('No nausea in placebo patients')$comparator_flag,flags('Nausea in combination therapy')$combination_flag)
flat<-function(z)as.character(unlist(z,use.names=FALSE))
rows<-list();segments_out<-list();counts<-list()
for(i in seq_len(nrow(selection))){
 s<-selection[i,];z<-readRDS(s$cache_file)
 matches<-Filter(function(v)identical(v$id,s$label_id),z$body$results)
 stopifnot(length(matches)==1L);rec<-matches[[1]];stopifnot(identical(rec$set_id,s$set_id))
 blocks<-list();j<-0L
 for(text in flat(rec$adverse_reactions)){
  for(sentence in str_split(str_squish(text),'(?<=[.!?])\\s+(?=[A-Z0-9\u2022])')[[1]]){
   j<-j+1L;blocks[[j]]<-tibble(source_field='adverse_reactions',segment_index=j,text=sentence,header=NA_character_)
  }
 }
 tn<-0L
 for(html in flat(rec$adverse_reactions_table)){
  tn<-tn+1L;doc<-read_html(html);tr<-xml_find_all(doc,'//tr')
  header<-str_squish(paste(xml_text(xml_find_all(doc,'//th')),collapse=' | '))
  if(!nzchar(header))header<-str_squish(paste(head(xml_text(tr),3),collapse=' | '))
  for(trrow in tr){
   j<-j+1L;blocks[[j]]<-tibble(source_field=paste0('adverse_reactions_table_',tn),segment_index=j,
    text=paste(str_squish(xml_text(xml_find_all(trrow,'./td|./th'))),collapse=' | '),header=header)
  }
 }
 sg<-bind_rows(blocks);matched<-0L
 if(nrow(sg)){
  locations<-str_locate_all(norm(sg$text),pattern)
  for(k in seq_len(nrow(sg))){
   loc<-locations[[k]];if(!nrow(loc))next
   v<-sg[k,];text<-norm(v$text);key<-paste(s$label_id,s$cache_md5,v$source_field,v$segment_index,sep=':')
   segments_out[[length(segments_out)+1L]]<-tibble(source_segment_id=key,label_id=s$label_id,
    source_field=v$source_field,segment_index=v$segment_index,source_text=v$text,table_header_context=v$header)
   found<-tibble(alias_norm=str_sub(text,loc[,1],loc[,2]),match_start=loc[,1],match_end=loc[,2],
    evidence_excerpt=str_sub(text,pmax(1,loc[,1]-180),pmin(nchar(text),loc[,2]+220)))
   rows[[length(rows)+1L]]<-bind_cols(s[rep(1,nrow(found)),] |> select(-cache_file),
    tibble(source_segment_id=rep(key,nrow(found)),source_field=rep(v$source_field,nrow(found))),found,
    flags(paste(text,ifelse(is.na(v$header),'',v$header)))[rep(1,nrow(found)),])
   matched<-matched+nrow(found)
  }
 }
 counts[[i]]<-s |> mutate(source_segments=nrow(sg),literal_mentions=matched)
 if(i%%10L==0L||i==nrow(selection))message('Shared extraction: ',i,' / ',nrow(selection),' label associations')
}
hits<-bind_rows(rows) |> inner_join(dict,by='alias_norm',relationship='many-to-many') |>
 distinct(entity_id,label_id,source_segment_id,match_start,pt_cui,.keep_all=TRUE) |>
 mutate(evidence_key=paste(entity_id,source_segment_id,match_start,pt_cui,sep=':'),
  extraction_status='candidate_requires_context_review',review_decision=NA_character_,review_reason=NA_character_,
  review_priority=case_when(!identity_eligible~'identity_or_named_context_review',
   ambiguous_dictionary_mapping=='TRUE'|broad_term=='TRUE'~'dictionary_scope_review',
   common_summary_flag & !negation_flag & !comparator_flag & !combination_flag~'explicit_common_summary_review',
   TRUE~'narrative_or_table_context_review'))
stopifnot(!anyDuplicated(hits$evidence_key),identical(before,tools::md5sum(inputs)))
write_csv(selection,'data/curated/comparable_label_registry.csv',na='')
write_csv(hits,'data/curated/comparable_label_PT_candidate_evidence.csv',na='')
write_csv(bind_rows(segments_out) |> distinct(source_segment_id,.keep_all=TRUE),
 'data/curated/comparable_label_evidence_segments.csv',na='')
write_csv(hits |> filter(review_priority=='explicit_common_summary_review'),
 'data/curated/comparable_label_common_summary_review_queue.csv',na='')
# Keep entities with no matches in coverage; absence of a dictionary hit is not absence of effects.
cov<-bind_rows(counts) |> group_by(entity_id,cohort,entity_name,comparison_layer) |>
 summarise(label_records=n(),source_segments=sum(source_segments),literal_mentions=sum(literal_mentions),.groups='drop') |>
 left_join(hits |> group_by(entity_id,comparison_layer) |> summarise(candidate_PT_terms=n_distinct(pt_cui),evidence_rows=n(),.groups='drop'),by=c('entity_id','comparison_layer')) |>
 mutate(candidate_PT_terms=replace_na(candidate_PT_terms,0L),evidence_rows=replace_na(evidence_rows,0L))
write_csv(cov,'outputs/15_comparable_label_coverage.csv',na='')
write_csv(hits |> count(cohort,comparison_layer,review_priority,name='evidence_rows'),
 'outputs/15_review_workload_summary.csv',na='')
write_csv(tibble(file=inputs,md5=unname(before)),'outputs/15_input_summary.csv')
report<-c('SCRIPT 15 - SHARED CURRENT-LABEL PT CANDIDATE EXTRACTION',
 paste('Entity/label associations processed:',nrow(selection)),paste('Distinct label documents:',n_distinct(selection$label_id)),
 paste('Entities with collected labels:',n_distinct(selection$entity_id)),paste('Candidate evidence rows:',nrow(hits)),
 paste('Entities with at least one dictionary match:',n_distinct(hits$entity_id)),
 paste('Explicit common-summary review rows:',sum(hits$review_priority=='explicit_common_summary_review')),
 'Inputs: all 155 seed selections, all 54 candidate group/label review records and the recovered sacrosidase context label; existing SIDER PT dictionary.',
 'Outputs locally: comparable_label_registry.csv, comparable_label_PT_candidate_evidence.csv, comparable_label_evidence_segments.csv, comparable_label_common_summary_review_queue.csv.',
 'Public reports: outputs/15_* coverage, review workload, input hashes and validation.',
 'Same normalisation, exact dictionary boundaries, sentence/table segmentation and context flags applied across cohorts.',
 'Unresolved candidate identities are processed for audit but remain outside identity-eligible comparisons; sacrosidase stays named-biologic context only.',
 'Shared label/segment IDs make repeated source documents visible; multiple labels/mentions are not independent experiments.',
 'All hits require context review: flags are screening aids, not proof of positive adverse-event attribution.',
 'Tables may contain placebo columns, mixed regimens, grades and frequencies; narrative may mention indications, negation and outcomes.',
 'Dictionary is limited to SIDER mappings; newer/unmapped terms may be missed. Longest non-overlapping literal matches are extraction candidates.',
 'Prior reviewed features and historical SIDER tables are unchanged. No accepted features, new ranks, p-values or BH correction generated.',
 'Next: review bounded source clauses/table contexts to build comparable PT profiles, retaining identity, regimen and route strata before scoring.')
writeLines(report,'outputs/15_comparable_label_extraction_summary.txt')
writeLines(c('PASS: unique entity/label and evidence keys; each cache resolves exactly one selected label and its set ID.',
 'PASS: word boundaries and distinct negation/comparator/combination flag controls.',
 'PASS: original input checksums unchanged; zero-hit selected entities retained in coverage.'),
 'outputs/15_comparable_label_validation.txt')
cat(paste0('\nSCRIPT 15 RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/15_extract_comparable_label_terms.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
