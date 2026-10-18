# Retrieve current openFDA label candidates for all seeds, including SIDER gaps.
# These are supplementary label data, not an updated release of SIDER.
library(readr);library(dplyr);library(stringr);library(tibble);library(httr2);library(tidyr)
inputs<-c('data/curated/seed_identity_crosscheck_reviewed.csv','outputs/11c_sider_reviewed_coverage.csv')
if(!all(file.exists(inputs)))stop('Run from the project root containing config/, scripts/ and data/.')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
x<-r(inputs[1]);coverage<-r(inputs[2]);stopifnot(!anyDuplicated(x$drugbank_id))
root<-'data/raw/current_labels/openfda';dir.create(root,recursive=TRUE,showWarnings=FALSE)
flat<-function(z)unique(as.character(unlist(z,use.names=FALSE)))
one<-function(z)if(length(z))as.character(z[[1]])else NA_character_
join<-function(z)paste(flat(z),collapse='|')
norm<-function(z)str_to_lower(str_squish(z))
meta<-list();records<-list();calls<-list()
for(i in seq_len(nrow(x))){
 seed<-x[i,];name<-gsub('["\\\\]',' ',seed$drug_name)
 q<-paste0('(openfda.generic_name:"',name,'" OR openfda.substance_name:"',name,'"',
 if(!is.na(seed$unii))paste0(' OR openfda.unii:"',seed$unii,'"')else '',
 ') AND openfda.product_type:"HUMAN PRESCRIPTION DRUG"')
 path<-file.path(root,paste0(seed$drugbank_id,'.rds'))
 if(file.exists(path)){
   z<-readRDS(path);if(!identical(z$query,q))stop('Cached query changed: ',seed$drugbank_id)
 }else{
   Sys.sleep(.3)
   req<-request('https://api.fda.gov/drug/label.json')|>req_url_query(search=q,limit=50,sort='effective_time:desc')
   resp<-req|>req_timeout(60)|>req_retry(max_tries=3,retry_on_failure=TRUE)|>
     req_error(is_error=function(resp)FALSE)|>req_perform()
   status<-resp_status(resp)
   if(!status %in% c(200L,404L))stop('openFDA status ',status,' for ',seed$drugbank_id)
   z<-list(query=q,url=resp_url(resp),retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
     status=status,body=resp_body_json(resp));saveRDS(z,path)
 }
 total<-if(z$status==200L)as.integer(z$body$meta$results$total)else 0L
 calls[[i]]<-tibble(drugbank_id=seed$drugbank_id,query=q,url=z$url,status=z$status,
  retrieved_utc=z$retrieved_utc,api_last_updated=one(z$body$meta$last_updated),
  total_matching_records=total,returned_records=length(z$body$results),search_truncated=total>50L,cache_file=path)
 for(rec in z$body$results){
   f<-rec$openfda;uniis<-flat(f$unii);generic<-flat(f$generic_name);substance<-flat(f$substance_name)
   uniimatch<-!is.na(seed$unii)&&seed$unii %in% uniis
   exact<-norm(seed$drug_name) %in% norm(c(generic,substance))
   single<-length(uniis)==1L&&length(substance)<=1L&&length(generic)<=1L
   identity<-if(uniimatch&&exact&&single)'UNII_and_exact_primary_name_supported_single_ingredient' else
     if(!single)'review_multiple_or_unharmonised_ingredients' else
     if(uniimatch)'review_UNII_match_name_or_form_differs' else 'review_name_only_or_UNII_disagreement'
   id<-one(rec$id);setid<-one(rec$set_id)
   metadata<-tibble(drugbank_id=seed$drugbank_id,drug_name=seed$drug_name,label_id=id,set_id=setid,
    effective_time=one(rec$effective_time),label_version=one(rec$version),brand_names=join(f$brand_name),
    generic_names=join(f$generic_name),substance_names=join(f$substance_name),label_uniis=join(f$unii),
    routes=join(f$route),application_numbers=join(f$application_number),
    original_packager=identical(f$is_original_packager[[1]],TRUE),
    label_identity_status=identity,adverse_text_chars=nchar(paste(flat(rec$adverse_reactions),collapse=' ')),
    adverse_tables=length(rec$adverse_reactions_table),cache_file=path,retrieved_utc=z$retrieved_utc,
    label_url=paste0('https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=',setid))
   meta[[length(meta)+1L]]<-metadata
   records[[length(records)+1L]]<-list(drugbank_id=seed$drugbank_id,label_id=id,record=rec)
 }
 if(i%%10L==0L||i==nrow(x))message('Label searches: ',i,' / ',nrow(x))
}
m<-bind_rows(meta) |> distinct(drugbank_id,label_id,.keep_all=TRUE)
# Recent representative labels, favouring original packagers. Full query results remain cached.
selected<-m |> filter(label_identity_status=='UNII_and_exact_primary_name_supported_single_ingredient',
 adverse_text_chars>0L|adverse_tables>0L) |>
 arrange(drugbank_id,desc(original_packager),desc(effective_time),desc(adverse_text_chars),label_id) |>
 distinct(drugbank_id,set_id,.keep_all=TRUE) |> group_by(drugbank_id) |> slice_head(n=3L) |> ungroup()
m<-m |> mutate(selected_for_extraction=paste(drugbank_id,label_id) %in% paste(selected$drugbank_id,selected$label_id))
call<-bind_rows(calls)
cover<-coverage |> left_join(m |> group_by(drugbank_id) |>
 summarise(candidate_labels=n(),supported_labels=sum(label_identity_status=='UNII_and_exact_primary_name_supported_single_ingredient'),
 selected_labels=sum(selected_for_extraction),.groups='drop'),by='drugbank_id') |>
 mutate(across(c(candidate_labels,supported_labels,selected_labels),~replace_na(as.integer(.x),0L)))
# Explicit namespace avoids requiring tidyr to be attached for the coverage conversion.
write_csv(m,'data/curated/current_label_candidates.csv',na='')
write_csv(selected,'data/curated/current_label_selection.csv',na='')
write_csv(call,'outputs/11d_current_label_query_provenance.csv',na='')
write_csv(cover,'outputs/11d_current_label_coverage.csv',na='')
saveRDS(records,'data/raw/current_labels/selected_query_records.rds')
report<-c('SCRIPT 11d - CURRENT LABEL COLLECTION',
 paste('Seeds queried:',nrow(x)),paste('Candidate seed/label records:',nrow(m)),
 paste('Seeds with selected identity-supported adverse-reaction labels:',n_distinct(selected$drugbank_id)),
 paste('Selected representative labels:',nrow(selected)),paste('Queries with >50 results:',sum(call$search_truncated)),
 'Outputs: current_label_candidates.csv/current_label_selection.csv; outputs/11d_*; raw response caches.',
 'Select up to three label sets per seed, requiring exact primary name and source UNII agreement',
 'with one harmonised ingredient. Other formulations, combinations and ambiguous names remain reviewable.',
 'This is a bounded representative-label supplement, not exhaustive current-label acquisition.',
 'US-label absence does not establish lack of adverse effects or lack of approval elsewhere.',
 'No extracted adverse event is yet an accepted drug/PT association. Next: controlled-vocabulary extraction and review.')
writeLines(report,'outputs/11d_current_label_collection_summary.txt')
cat(paste0('\nSCRIPT 11d RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/11d_collect_current_labels.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
