# Recover current labels with explicit active-moiety evidence, retaining salt/combination distinctions.
library(readr);library(dplyr);library(tidyr);library(stringr);library(tibble);library(httr2);library(xml2)
if(!file.exists('data/curated/current_label_selection.csv'))stop('Run from project root after 11d.')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
x<-r('data/curated/seed_identity_crosscheck_reviewed.csv');base<-r('data/curated/current_label_selection.csv')
candidates<-r('data/curated/current_label_candidates.csv');# Gemcitabine needs the intravenous cancer-label route as well as the selected intravesical formulation.
targets<-x |> filter(!drugbank_id %in% base$drugbank_id | drugbank_id=='DB00441')
root<-'data/raw/current_labels/dailymed';dir.create(root,recursive=TRUE,showWarnings=FALSE)
logs<-list();decisions<-list();new<-list()
one<-function(z)if(length(z))as.character(z[[1]])else NA_character_
fetch<-function(url,label,json=FALSE){
 path<-file.path(root,paste0(label,'.rds'))
 if(file.exists(path)){z<-readRDS(path);stopifnot(identical(z$url,url))}else{
 Sys.sleep(.25);resp<-request(url)|>req_timeout(60)|>req_retry(max_tries=3,retry_on_failure=TRUE)|>
 req_error(is_error=function(r)FALSE)|>req_perform();status<-resp_status(resp)
 if(!status %in% c(200L,404L))stop('DailyMed status ',status,' at ',url)
 z<-list(url=url,status=status,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
 body=if(json)resp_body_json(resp)else resp_body_string(resp));saveRDS(z,path)}
 logs[[length(logs)+1L]]<<-tibble(url=z$url,status=z$status,retrieved_utc=z$retrieved_utc,cache_file=path)
 z
}
for(i in seq_len(nrow(targets))){
 seed<-targets[i,]
 choices<-candidates |> filter(drugbank_id==seed$drugbank_id) |>
 arrange(desc(original_packager),desc(effective_time)) |> distinct(set_id,.keep_all=TRUE) |> slice_head(n=3)
 if(seed$drugbank_id=='DB00441') choices<-candidates |> filter(drugbank_id=='DB00441',str_detect(routes,'INTRAVENOUS')) |>
   arrange(desc(original_packager),desc(effective_time)) |> distinct(set_id,.keep_all=TRUE) |> slice_head(n=3)
 sets<-unique(na.omit(choices$set_id))
 if(!length(sets)){
   req<-request('https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json')|>
    req_url_query(drug_name=seed$drug_name,pagesize=100)
   z<-fetch(req$url,paste0('search_',seed$drugbank_id),TRUE)
   sets<-if(z$status==200L)head(unique(vapply(z$body$data,function(a)a$setid,character(1))),3L)else character()
 }
 for(setid in sets){
   z<-fetch(paste0('https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/',setid,'.xml'),paste0('spl_',setid))
   if(z$status!=200L)next
   doc<-read_xml(z$body);xml_ns_strip(doc)
   ing<-xml_find_all(doc,'//ingredient[starts-with(@classCode,"ACTI")]/ingredientSubstance')
   ingredient_ids<-unique(na.omit(xml_attr(xml_find_all(ing,'./code'),'code')))
   moiety_ids<-unique(na.omit(xml_attr(xml_find_all(ing,'.//activeMoiety/code'),'code')))
   ingredient_names<-unique(xml_text(xml_find_all(ing,'./name')))
   moiety_names<-unique(xml_text(xml_find_all(ing,'.//activeMoiety/name')))
   supported<-!is.na(seed$unii)&&seed$unii %in% c(ingredient_ids,moiety_ids)
   named<-str_to_lower(seed$drug_name) %in% str_to_lower(c(ingredient_names,moiety_names))
   single<-length(ingredient_ids)==1L && length(moiety_ids)<=1L
   rx<-identical(xml_attr(xml_find_first(doc,'/document/code'),'code'),'34391-3')
   accepted<-supported&&named&&single&&rx
   status<-if(accepted)'explicit_single_ingredient_or_active_moiety_UNII_and_name'else
    if(!single)'review_combination_or_ingredient_structure'else 'review_name_UNII_or_product_type'
   adverse<-xml_find_all(doc,'//section[code[@code="34084-4"]]')
   texts<-unique(xml_text(xml_find_all(adverse,'.//text')))
   tables<-unique(vapply(xml_find_all(adverse,'.//table'),as.character,character(1)))
   decisions[[length(decisions)+1L]]<-tibble(drugbank_id=seed$drugbank_id,drug_name=seed$drug_name,set_id=setid,
    ingredient_UNIIs=paste(ingredient_ids,collapse='|'),moiety_UNIIs=paste(moiety_ids,collapse='|'),
    ingredient_names=paste(ingredient_names,collapse='|'),moiety_names=paste(moiety_names,collapse='|'),
    label_identity_status=status,adverse_text_chars=nchar(paste(texts,collapse=' ')),
    url=z$url,retrieved_utc=z$retrieved_utc)
   if(!accepted||!length(texts))next
   labelid<-xml_attr(xml_find_first(doc,'/document/id'),'root')
   date<-xml_attr(xml_find_first(doc,'/document/effectiveTime'),'value')
   version<-xml_attr(xml_find_first(doc,'/document/versionNumber'),'value')
   # Normalise the SPL into the same evidence schema as openFDA; keep original XML cached.
   normal_path<-file.path(root,paste0('normalised_',setid,'.rds'))
   rec<-list(id=labelid,set_id=setid,effective_time=date,version=version,
    adverse_reactions=as.list(texts),adverse_reactions_table=as.list(tables))
   saveRDS(list(url=z$url,retrieved_utc=z$retrieved_utc,status=200L,body=list(results=list(rec))),normal_path)
   new[[length(new)+1L]]<-tibble(drugbank_id=seed$drugbank_id,drug_name=seed$drug_name,label_id=labelid,set_id=setid,
    effective_time=date,label_version=version,brand_names=NA_character_,
    generic_names=paste(moiety_names,collapse='|'),substance_names=paste(ingredient_names,collapse='|'),
    label_uniis=paste(ingredient_ids,collapse='|'),routes=paste(unique(na.omit(xml_attr(xml_find_all(doc,'//routeCode'),'displayName'))),collapse='|'),application_numbers=NA_character_,
    original_packager=NA_character_,label_identity_status=status,adverse_text_chars=as.character(nchar(paste(texts,collapse=' '))),
    adverse_tables=as.character(length(tables)),cache_file=normal_path,retrieved_utc=z$retrieved_utc,
    label_url=paste0('https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=',setid))
 }
 if(i%%10L==0L||i==nrow(targets))message('DailyMed fallback seeds checked: ',i,' / ',nrow(targets))
}
extra<-bind_rows(new)
selection<-bind_rows(base,extra) |> distinct(drugbank_id,label_id,.keep_all=TRUE)
write_csv(selection,'data/curated/current_label_selection_supplemented.csv',na='')
write_csv(bind_rows(decisions),'data/curated/current_label_moiety_review.csv',na='')
write_csv(bind_rows(logs),'outputs/11f_DailyMed_query_provenance.csv',na='')
report<-c('SCRIPT 11f - DAILYMED ACTIVE-MOIETY LABEL SUPPLEMENT',
 paste('Seeds checked:',nrow(targets)),paste('Additional seeds with supported labels:',n_distinct(extra$drugbank_id[!extra$drugbank_id %in% base$drugbank_id])),
 paste('Total seeds with selected labels:',n_distinct(selection$drugbank_id)),
 'Output: current_label_selection_supplemented.csv and current_label_moiety_review.csv (data/curated/).',
 'Queries/versions and raw XML are retained; up to three candidate sets examined per seed.',
 'Accept only literal seed-name plus UNII agreement with a single ingredient or its explicit active moiety.',
 'An active-moiety relation supports clinical label attribution, not salt/free-base structural equivalence.',
 'Combination products and unmatched/ambiguous ingredients remain outside the accepted extraction selection.',
 'Gemcitabine intravenous labels are explicitly recovered; intravesical-device data do not substitute for systemic cancer-label data.',
 'DailyMed and openFDA describe overlapping SPL labels, not independent experimental evidence.',
 'Next: rerun 11e on the supplemented selection, then context-review term candidates.')
writeLines(report,'outputs/11f_DailyMed_supplement_summary.txt')
cat(paste0('\nSCRIPT 11f RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/11f_resolve_label_moieties.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
