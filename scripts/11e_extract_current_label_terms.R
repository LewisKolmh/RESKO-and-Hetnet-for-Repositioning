# Map current label text to SIDER PT CUIs with an auditable dictionary and evidence.
# Extraction produces candidates, not automatically validated adverse reactions.
library(readr);library(dplyr);library(tidyr);library(stringr);library(tibble);library(xml2)
selection_path<-if(file.exists('data/curated/current_label_selection_supplemented.csv'))
 'data/curated/current_label_selection_supplemented.csv' else 'data/curated/current_label_selection.csv'
required<-c(selection_path,'data/raw/sider/meddra.tsv.gz')
if(!all(file.exists(required)))stop('Run from project root after script 11d.')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
se<-read_tsv('data/raw/sider/meddra_all_se.tsv.gz',col_names=c('flat','stereo','label_cui','type','pt_cui','term'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
med<-read_tsv('data/raw/sider/meddra.tsv.gz',col_names=c('cui','kind','meddra_id','term'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
norm<-function(z)str_to_lower(str_squish(str_replace_all(z,'[\u2010-\u2015]','-')))
pt<-se |> filter(type=='PT') |> distinct(pt_cui,pt_name=term)
# Existing SIDER label-CUI --> PT-CUI mapping supplies LLT synonyms, without invented synonyms.
links<-se |> filter(type=='PT') |> distinct(label_cui,pt_cui)
dictionary<-bind_rows(
 pt |> transmute(alias=pt_name,pt_cui,alias_source='SIDER_PT_name'),
 med |> filter(kind=='PT') |> inner_join(pt |> distinct(pt_cui),by=c('cui'='pt_cui')) |>
 transmute(alias=term,pt_cui=cui,alias_source='MedDRA_PT_same_CUI'),
 med |> filter(kind %in% c('LT','LLT')) |> inner_join(links,by=c('cui'='label_cui'),relationship='many-to-many') |>
 transmute(alias=term,pt_cui,alias_source='SIDER_label_CUI_to_PT_CUI')
) |> filter(!is.na(alias),nchar(alias)>=3L) |> mutate(alias_norm=norm(alias)) |>
 distinct(alias_norm,pt_cui,.keep_all=TRUE) |>
 left_join(pt |> group_by(pt_cui) |> summarise(pt_name=first(pt_name),.groups='drop'),by='pt_cui') |>
 group_by(alias_norm) |> mutate(ambiguous_dictionary_mapping=n_distinct(pt_cui)>1L) |> ungroup()
# Short broad terms remain dictionary entries, but are explicitly flagged for review.
dictionary<-dictionary |> mutate(broad_term=alias_norm %in% c('pain','infection','injury','disorder','death','cancer','neoplasm','lesion','toxicity','disease'))
terms<-unique(dictionary$alias_norm);terms<-terms[order(nchar(terms),decreasing=TRUE)]
pattern<-paste0('(?<![\\p{L}\\p{N}])(?:',paste(str_escape(terms),collapse='|'),')(?![\\p{L}\\p{N}])')
selected<-r(required[1]);rows<-list();segment_rows<-list();n<-0L
flatten<-function(z)as.character(unlist(z,use.names=FALSE))
context_flags<-function(text){
 t<-str_to_lower(text)
 tibble(negation_flag=str_detect(t,'\\b(no|not|without|absence of|neither|none|did not)\\b'),
 comparator_flag=str_detect(t,'\\b(placebo|comparator|control arm|control group)\\b'),
 combination_flag=str_detect(t,'\\b(combination|combined|plus|coadministered|co-administered)\\b'),
 grade_or_frequency_context=str_detect(t,'grade [0-9]|%|percent|incidence'),
 common_summary_flag=str_detect(t,'most common.{0,80}adverse (reactions|events)|common adverse reactions|most frequently reported adverse'))
}
# Boundary and context controls exercise distinct failure modes of text extraction.
stopifnot(!str_detect('painting',regex('(?<![\\p{L}\\p{N}])pain(?![\\p{L}\\p{N}])')),
 context_flags('No nausea was observed')$negation_flag,
 context_flags('Placebo patients reported nausea')$comparator_flag,
 context_flags('Nausea occurred in combination therapy')$combination_flag)
for(i in seq_len(nrow(selected))){
 s<-selected[i,];z<-readRDS(s$cache_file);rec<-Filter(function(a)identical(a$id,s$label_id),z$body$results)[[1]]
 sections<-list();j<-0L
 for(text in flatten(rec$adverse_reactions)){
   # Sentence segmentation is only an evidence window; full source sections are cached.
   sentences<-str_split(str_squish(text),'(?<=[.!?])\\s+(?=[A-Z0-9\u2022])')[[1]]
   for(t in sentences){j<-j+1L;sections[[j]]<-tibble(source_field='adverse_reactions',segment_index=j,text=t,
     table_header_context=NA_character_)}
 }
 table_no<-0L
 for(html in flatten(rec$adverse_reactions_table)){
   table_no<-table_no+1L;doc<-read_html(html)
   table_rows<-xml_find_all(doc,'//tr');header<-str_squish(paste(xml_text(xml_find_all(doc,'//th')),collapse=' | '))
   if(!nzchar(header))header<-str_squish(paste(head(xml_text(table_rows),3),collapse=' | '))
   for(row in table_rows){j<-j+1L;sections[[j]]<-tibble(source_field=paste0('adverse_reactions_table_',table_no),
      segment_index=j,text=paste(str_squish(xml_text(xml_find_all(row,'./td|./th'))),collapse=' | '),table_header_context=header)}
 }
 segments<-bind_rows(sections)
 locations<-str_locate_all(norm(segments$text),regex(pattern))
 for(k in seq_len(nrow(segments))){
   sg<-segments[k,];text<-norm(sg$text);loc<-locations[[k]]
   if(!nrow(loc)||is.na(loc[1,1]))next
   segment_id<-paste(s$drugbank_id,s$label_id,sg$source_field,sg$segment_index,sep=':')
   segment_rows[[length(segment_rows)+1L]]<-tibble(segment_id,drugbank_id=s$drugbank_id,label_id=s$label_id,
     source_field=sg$source_field,segment_index=sg$segment_index,source_text=sg$text,table_header_context=sg$table_header_context)
   for(h in seq_len(nrow(loc))){
     n<-n+1L;start<-loc[h,1];end<-loc[h,2]
     excerpt<-str_sub(text,max(1,start-180),min(nchar(text),end+220))
     flags<-context_flags(paste(text,ifelse(is.na(sg$table_header_context),'',sg$table_header_context)))
     rows[[n]]<-bind_cols(tibble(drugbank_id=s$drugbank_id,drug_name=s$drug_name,label_id=s$label_id,
       set_id=s$set_id,effective_time=s$effective_time,label_url=s$label_url,
       label_identity_status=s$label_identity_status,segment_id,source_field=sg$source_field,
       alias_norm=str_sub(text,start,end),evidence_excerpt=excerpt,match_start=start,match_end=end),flags)
   }
 }
 if(i%%10L==0L||i==nrow(selected))message('Labels extracted: ',i,' / ',nrow(selected))
}
hits<-bind_rows(rows) |> inner_join(dictionary,by='alias_norm',relationship='many-to-many') |>
 mutate(extraction_status='candidate_requires_context_review') |>
 distinct(drugbank_id,label_id,segment_id,alias_norm,pt_cui,.keep_all=TRUE)
# Stable IDs ensure a later review cannot drift silently after dictionary/label changes.
hits<-hits |> arrange(drugbank_id,label_id,segment_id,match_start,pt_cui) |> mutate(evidence_id=sprintf('LABELPT%07d',row_number()))
write_csv(dictionary,'data/curated/current_label_PT_dictionary.csv',na='')
write_csv(hits,'data/curated/current_label_PT_candidate_evidence.csv',na='')
write_csv(bind_rows(segment_rows) |> distinct(segment_id,.keep_all=TRUE),
 'data/curated/current_label_evidence_segments.csv',na='')
summary<-hits |> group_by(drugbank_id,drug_name) |> summarise(candidate_PT_terms=n_distinct(pt_cui),
 evidence_rows=n(),selected_label_records=n_distinct(label_id),.groups='drop')
write_csv(summary,'outputs/11e_current_label_term_coverage.csv')
write_csv(hits |> filter(common_summary_flag,!ambiguous_dictionary_mapping,!broad_term) |>
 mutate(review_decision=NA_character_,review_reason=NA_character_),
 'outputs/11e_common_reaction_review_queue.csv',na='')
report<-c('SCRIPT 11e - CURRENT-LABEL TERM EXTRACTION',paste('Selection input:',selection_path),paste('Dictionary aliases/PT links:',nrow(dictionary)),
 paste('Evidence rows:',nrow(hits)),paste('Seeds with candidate PT matches:',nrow(summary)),
 'Outputs: current_label_PT_dictionary.csv, current_label_PT_candidate_evidence.csv,',
 'current_label_evidence_segments.csv (data/curated/); coverage and summary review queue (outputs/11e_*).',
 'Map PT names and existing LLT-to-PT CUI links only. No embedding/fuzzy matches or invented synonym equivalences.',
 'Every hit retains label ID/date, source field, evidence segment and context flags.',
 'Narrative sentences and HTML table rows are retained separately; duplicate labels do not become independent evidence.',
 'All matches are candidates until reviewed; placebo, negation, combinations, broad terms and ambiguous dictionary links are flagged.',
 'Labels may contain outcomes/diseases unrelated to adverse reactions even within this section.',
 'Dictionary coverage is limited to PTs/mappings represented in SIDER; novel terms can be missed.',
 'No current-label candidates are inserted into the original SIDER feature table.',
 'Next: review explicit common-reaction statements and create a separately sourced supplementary PT feature table.')
writeLines(report,'outputs/11e_current_label_extraction_summary.txt')
cat(paste0('\nSCRIPT 11e RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/11e_extract_current_label_terms.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
