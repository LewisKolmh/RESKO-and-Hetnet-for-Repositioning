# Collect a bounded, provenance-bearing current-label panel; do not infer PTs or binding.
library(readr);library(dplyr);library(tidyr);library(stringr);library(tibble);library(xml2);library(httr2)
if(!file.exists('data/curated/candidate_identity_registry.csv'))stop('Run from project root after 14a.')
source('scripts/helpers/cached_public_get.R')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
x<-r('data/curated/candidate_identity_registry.csv')
db<-r('data/curated/candidate_DrugBank_record_enrichment.csv')
root<-'data/raw/candidate_labels';dir.create(root,recursive=TRUE,showWarnings=FALSE)
logs<-list();queries<-list();decisions<-list();coverage<-list()
scalar<-function(z)if(length(z))as.character(z[[1]])else NA_character_
join<-function(z)paste(sort(unique(na.omit(z[nzchar(z)]))),collapse='|')
norm<-function(z)str_to_lower(str_squish(z))
fetch<-function(url,path,json=TRUE){
 z<-public_cache_get(url,path,json)
 logs[[length(logs)+1L]]<<-tibble(url=z$url,status=z$status,retrieved_utc=z$retrieved_utc,cache_file=path)
 z
}
for(i in seq_len(nrow(x))){
 a<-x[i,];id<-a$stitch_flat_id;records<-db |> filter(stitch_flat_id==id)
 if(a$automatic_label_mapping_blocked=='TRUE'){
 coverage[[i]]<-tibble(stitch_flat_id=id,sider_drug_name=a$sider_drug_name,search_status='blocked_identity_conflict',
 queries=0L,query_result_rows=0L,sets_examined=0L,supported_AR_sets=0L,review_AR_sets=0L,search_truncated=FALSE)
 next
 }
 # Source-record UNII is evidence; names are search leads and never sufficient for acceptance.
 uniis<-sort(unique(na.omit(records$unii)));uniis<-uniis[nzchar(uniis)]
 names<-unique(na.omit(c(records$record_name,a$pubchem_title,a$sider_drug_name)))
 names<-names[nzchar(names)];names<-head(names,3L)
 tasks<-bind_rows(tibble(kind=rep('unii_code',length(uniis)),value=uniis),
  tibble(kind=rep('drug_name',length(head(names,2L))),value=head(names,2L)))
 hits<-list();truncated<-FALSE;query_n<-0L
 for(j in seq_len(nrow(tasks))){
  task<-tasks[j,]
  if(task$kind=='drug_name'&&length(uniis)&&length(hits))next
  req<-request('https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json')
  args<-setNames(list(task$value),task$kind)
  req<-do.call(req_url_query,c(list(req),args,list(pagesize=100,page=1)))
  z<-fetch(req$url,file.path(root,paste0('search_',id,'_',j,'.rds')))
  query_n<-query_n+1L
  rows<-if(z$status==200L)bind_rows(lapply(z$body$data,function(v)tibble(set_id=scalar(v$setid),title=scalar(v$title),published_date=scalar(v$published_date))))else tibble()
  total<-suppressWarnings(as.integer(scalar(z$body$metadata$total_elements)))
  if(!is.na(total)&&total>nrow(rows))truncated<-TRUE
  queries[[length(queries)+1L]]<-tibble(stitch_flat_id=id,query_kind=task$kind,query_value=task$value,
   result_rows=nrow(rows),total_elements=total,url=z$url,retrieved_utc=z$retrieved_utc)
  if(nrow(rows))hits[[length(hits)+1L]]<-rows
 }
 # At most three most recently published candidate sets from the first API page(s).
 # This is a representative panel, not exhaustive label coverage or proof of approval.
 all_hits<-bind_rows(hits)
 if(nrow(all_hits)){
  all_hits<-all_hits |> mutate(published_day=parse_date(published_date,format='%b %d, %Y',locale=locale('en')))
  stopifnot(!anyNA(all_hits$published_day))
 }
 sets<-if(nrow(all_hits))all_hits |> arrange(desc(published_day),set_id) |> distinct(set_id,.keep_all=TRUE) |> slice_head(n=3)else tibble()
 if(nrow(all_hits)&&n_distinct(all_hits$set_id)>3L)truncated<-TRUE
 group_decisions<-list()
 for(k in seq_len(nrow(sets))){
  setid<-sets$set_id[k];url<-paste0('https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/',setid,'.xml')
  path<-file.path(root,paste0('spl_',setid,'.rds'));z<-fetch(url,path,FALSE)
  if(z$status!=200L)next
  doc<-read_xml(z$body);xml_ns_strip(doc)
  actual_set<-xml_attr(xml_find_first(doc,'/document/setId'),'root')
  stopifnot(!is.na(actual_set),tolower(actual_set)==tolower(setid))
  ing<-xml_find_all(doc,'//ingredient[starts-with(@classCode,"ACTI")]/ingredientSubstance')
  ing_ids<-unique(na.omit(xml_attr(xml_find_all(ing,'./code'),'code')))
  moi_ids<-unique(na.omit(xml_attr(xml_find_all(ing,'.//activeMoiety/code'),'code')))
  ing_names<-unique(xml_text(xml_find_all(ing,'./name')));moi_names<-unique(xml_text(xml_find_all(ing,'.//activeMoiety/name')))
  source_unii<-intersect(uniis,c(ing_ids,moi_ids))
  name_support<-any(norm(names) %in% norm(c(ing_names,moi_names)))
  single<-length(ing_ids)==1L&&length(moi_ids)<=1L
  doctype<-xml_attr(xml_find_first(doc,'/document/code'),'code')
  human<-!is.na(doctype)&&doctype %in% c('34391-3','34390-5')
  supported<-length(source_unii)>0L&&name_support&&single&&human
  status<-if(supported)'source_UNII_and_name_single_ingredient_or_explicit_moiety'else
   if(!human)'review_non_human_drug_document'else if(!single)'review_combination_or_unresolved_ingredients'else
   if(!length(source_unii))'review_no_source_record_UNII_agreement'else 'review_name_agreement_required'
  ar<-xml_find_all(doc,'//section[code[@code="34084-4"]]')
  texts<-unique(xml_text(xml_find_all(ar,'.//text')))
  tables<-unique(vapply(xml_find_all(ar,'.//table'),as.character,character(1)))
  text_chars<-nchar(paste(texts,collapse=' '));labelid<-xml_attr(xml_find_first(doc,'/document/id'),'root')
  effective<-xml_attr(xml_find_first(doc,'/document/effectiveTime'),'value')
  version<-xml_attr(xml_find_first(doc,'/document/versionNumber'),'value')
  normal_path<-file.path(root,paste0('normalised_',setid,'.rds'))
  # Store even unresolved source text locally; acceptance is explicit in the selection table.
  saveRDS(list(url=z$url,retrieved_utc=z$retrieved_utc,status=200L,body=list(results=list(
   list(id=labelid,set_id=setid,effective_time=effective,version=version,
    adverse_reactions=as.list(texts),adverse_reactions_table=as.list(tables))))),normal_path)
  row<-tibble(stitch_flat_id=id,sider_drug_name=a$sider_drug_name,record_drugbank_ids=join(records$record_drugbank_id),
   set_id=setid,label_id=labelid,effective_time=effective,label_version=version,document_type=doctype,
   label_title=sets$title[k],ingredient_UNIIs=join(ing_ids),moiety_UNIIs=join(moi_ids),
   ingredient_names=join(ing_names),moiety_names=join(moi_names),matched_source_UNIIs=join(source_unii),
   name_agreement=name_support,single_ingredient=single,label_identity_status=status,
   supported_identity=supported,adverse_text_chars=text_chars,adverse_tables=length(tables),
   selected_for_context_review=supported&&text_chars>0L,
   routes=join(xml_attr(xml_find_all(doc,'//routeCode'),'displayName')),
   cache_file=normal_path,raw_cache_file=path,url=z$url,retrieved_utc=z$retrieved_utc)
  decisions[[length(decisions)+1L]]<-row;group_decisions[[length(group_decisions)+1L]]<-row
 }
 gd<-bind_rows(group_decisions)
 coverage[[i]]<-tibble(stitch_flat_id=id,sider_drug_name=a$sider_drug_name,
  search_status=if(!nrow(all_hits))'no_results_in_bounded_search'else 'bounded_candidate_sets_checked',
  queries=query_n,query_result_rows=nrow(all_hits),sets_examined=nrow(gd),
  supported_AR_sets=if(nrow(gd))sum(gd$selected_for_context_review)else 0L,
  review_AR_sets=if(nrow(gd))sum(!gd$supported_identity & gd$adverse_text_chars>0L)else 0L,search_truncated=truncated)
 message('Candidate label groups checked: ',i,' / ',nrow(x),' (',a$sider_drug_name,')')
}
d<-bind_rows(decisions);cov<-bind_rows(coverage);selection<-d |> filter(selected_for_context_review)
stopifnot(nrow(cov)==nrow(x),!anyDuplicated(cov$stitch_flat_id),setequal(cov$stitch_flat_id,x$stitch_flat_id),
 !anyDuplicated(d[c('stitch_flat_id','set_id')]),all(selection$supported_identity),all(selection$adverse_text_chars>0),
 !any(selection$stitch_flat_id %in% x$stitch_flat_id[x$automatic_label_mapping_blocked=='TRUE']),
 all(file.exists(d$cache_file)))
write_csv(d,'data/curated/candidate_label_identity_review.csv',na='')
write_csv(selection,'data/curated/candidate_label_selection.csv',na='')
write_csv(bind_rows(queries),'data/curated/candidate_label_query_provenance.csv',na='')
write_csv(bind_rows(logs),'data/curated/candidate_label_download_provenance.csv',na='')
write_csv(cov,'outputs/14b_candidate_label_coverage.csv',na='')
write_csv(d |> select(stitch_flat_id,set_id,label_id,effective_time,label_version,document_type,label_identity_status,
 supported_identity,adverse_text_chars,adverse_tables,selected_for_context_review,url,retrieved_utc),
 'outputs/14b_candidate_label_crosswalk_summary.csv',na='')
report<-c('SCRIPT 14b - CURRENT CANDIDATE LABEL COLLECTION AND INGREDIENT CHECKS',
 paste('Historical SIDER groups retained:',nrow(cov)),paste('Identity-conflict groups blocked:',sum(cov$search_status=='blocked_identity_conflict')),
 paste('Groups with any downloaded SPL:',n_distinct(d$stitch_flat_id)),paste('Group/label records checked:',nrow(d)),
 paste('Distinct SPL sets downloaded:',n_distinct(d$set_id)),
 paste('Groups with UNII/name-supported single-ingredient adverse-reaction text:',n_distinct(selection$stitch_flat_id)),
 paste('Selected group/label records:',nrow(selection)),paste('Groups with truncated representative searches:',sum(cov$search_truncated)),
 'Inputs: 14a identity registry and exact-full-key DrugBank enrichment; DailyMed search and current SPL XML.',
 'Outputs: candidate_label_identity_review.csv, candidate_label_selection.csv, query/download provenance (data/curated/).',
 'Coverage and compact source crosswalk: outputs/14b_*; raw XML and normalised AR text/tables remain local.',
 'UNII search first; name fallback if no UNII result. Without a matched source UNII, name hits remain review-only.',
 'One page (up to 100 rows) per query; at most three most recently published sets checked per group. Coverage is bounded.',
 'Exact source-record UNII plus literal ingredient/moiety name and single-ingredient human Rx/OTC scope are required.',
 'Explicit active-moiety mapping permits clinical label attribution, not salt/free-base structural equivalence.',
 'SPL presence is not regulatory approval or marketing-status verification. Combination/nonhuman/unmatched records remain review-only.',
 'Historical SIDER identity is not established by current-label agreement. Conflicting sacrosidase CID is excluded automatically.',
 'No adverse-reaction section or missing labels is missing coverage, not evidence of safety or absence of effects.',
 'This script preserves AR sections, not accepted PT features; regimen, route, placebo, negation and frequency need context review.',
 'No similarity ranks, p-values, BH correction, direct binding, inhibition, allostery or PPI claims are generated.',
 'Next: script 15, comparably extract PT candidates from collected seed/candidate labels and review context before scoring.')
writeLines(report,'outputs/14b_candidate_label_summary.txt')
writeLines(c('PASS: all 30 original groups retained once in coverage; identity conflict excluded from selection.',
 'PASS: group/set keys unique; selected records have supported identity and nonempty AR text; normalised caches exist.',
 'PASS: downloaded SPL set IDs equal requested IDs; source versions, retrieval times and query truncation recorded.'),
 'outputs/14b_candidate_label_validation.txt')
cat(paste0('\nSCRIPT 14b RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/14b_collect_candidate_labels.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
