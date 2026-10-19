# Resolve named biologic context without assigning the conflicting SIDER CID structure.
library(readr);library(dplyr);library(tibble);library(xml2)
if(!file.exists('config/candidate_identity_exceptions.csv'))stop('Run from project root after 14b.')
source('scripts/helpers/stream_drugbank_fields.R');source('scripts/helpers/cached_public_get.R')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
original<-c('data/raw/sider/drug_names.tsv','data/raw/sider/meddra_all_se.tsv.gz',
 'data/curated/candidate_identity_registry.csv','data/curated/candidate_label_selection.csv')
before<-tools::md5sum(original)
manifest<-r('config/database_manifest.csv');xml_path<-manifest$local_file[manifest$source_id=='drugbank'][[1]]
root<-'data/raw/candidate_labels';cache<-file.path(root,'sacrosidase_DrugBank_record.rds')
sig<-paste(normalizePath(xml_path),file.info(xml_path)$size,as.numeric(file.info(xml_path)$mtime),sep='|')
if(file.exists(cache)){
 local<-readRDS(cache);stopifnot(identical(local$signature,sig));rec<-local$record
}else{
 rec<-NULL
 callback<-function(n){
  get<-function(p)xml_text(xml_find_first(n,p))
  if(!identical(get('./drugbank-id[@primary="true"]'),'DB06760'))return(invisible(NULL))
  rec<<-list(drugbank_id='DB06760',name=get('./name'),unii=get('./unii'),cas=get('./cas-number'),
   type=xml_attr(n,'type'),description=get('./description'),aliases=xml_text(xml_find_all(n,'./synonyms/synonym')),
   xrefs=as.character(xml_find_all(n,'./external-identifiers')))
 }
 stream_drugbank_fields(xml_path,c('drugbank-id','name','unii','cas-number','description','synonyms','external-identifiers'),callback)
 stopifnot(!is.null(rec));saveRDS(list(signature=sig,record=rec),cache)
}
stopifnot(tolower(rec$name)=='sacrosidase',rec$unii=='8A7F670F2Y',rec$type=='biotech')
setid<-'d613bb7f-c3f4-462e-81a2-da2347cc4b6b'
url<-paste0('https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/',setid,'.xml')
z<-public_cache_get(url,file.path(root,paste0('spl_',setid,'.rds')),FALSE);stopifnot(z$status==200L)
doc<-read_xml(z$body);xml_ns_strip(doc)
ing<-xml_find_all(doc,'//ingredient[starts-with(@classCode,"ACTI")]/ingredientSubstance')
ids<-unique(na.omit(xml_attr(xml_find_all(ing,'./code'),'code')))
moieties<-unique(na.omit(xml_attr(xml_find_all(ing,'.//activeMoiety/code'),'code')))
names<-unique(xml_text(xml_find_all(ing,'./name')))
stopifnot(length(ids)==1L,rec$unii %in% c(ids,moieties),'sacrosidase' %in% tolower(names),
 xml_attr(xml_find_first(doc,'/document/code'),'code')=='34391-3',
 xml_attr(xml_find_first(doc,'/document/setId'),'root')==setid)
texts<-unique(xml_text(xml_find_all(doc,'//section[code[@code="34084-4"]]//text')))
tables<-unique(vapply(xml_find_all(doc,'//section[code[@code="34084-4"]]//table'),as.character,character(1)))
stopifnot(length(texts)>0L)
labelid<-xml_attr(xml_find_first(doc,'/document/id'),'root');effective<-xml_attr(xml_find_first(doc,'/document/effectiveTime'),'value')
version<-xml_attr(xml_find_first(doc,'/document/versionNumber'),'value');normal_path<-file.path(root,paste0('normalised_',setid,'.rds'))
saveRDS(list(url=url,status=200L,retrieved_utc=z$retrieved_utc,body=list(results=list(list(id=labelid,set_id=setid,
 effective_time=effective,version=version,adverse_reactions=as.list(texts),adverse_reactions_table=as.list(tables))))),normal_path)
row<-tibble(stitch_flat_id='CID100135113',sider_drug_name='sacrosidase',historical_stereo_id='CID000135113',
 rejected_structure_cid='135113',drugbank_id=rec$drugbank_id,unii=rec$unii,cas_number=rec$cas,
 entity_type='biologic_enzyme',mapping_scope='historical_name_context_to_current_UNII_verified_label',
 structural_equivalence=FALSE,historical_CID_conflict_resolved=FALSE,automatic_structure_mapping_blocked=TRUE,
 set_id=setid,label_id=labelid,effective_time=effective,label_version=version,
 ingredient_UNIIs=paste(ids,collapse='|'),ingredient_names=paste(names,collapse='|'),
 adverse_text_chars=nchar(paste(texts,collapse=' ')),adverse_tables=length(tables),cache_file=normal_path,
 selected_for_context_review=TRUE,include_in_small_molecule_structure_screen=FALSE,
 label_url=url,drugbank_url='https://go.drugbank.com/drugs/DB06760',retrieved_utc=z$retrieved_utc)
write_csv(row,'data/curated/candidate_label_named_context_selection.csv',na='')
write_csv(row |> select(-cache_file),'outputs/14c_sacrosidase_resolution_summary.csv',na='')
stopifnot(identical(before,tools::md5sum(original)))
write_csv(tibble(file=original,md5=unname(before)),'outputs/14c_preserved_inputs_summary.csv')
report<-c('SCRIPT 14c - SACROSIDASE NAMED-BIOLOGIC CONTEXT RECOVERY',
 'Verified local DrugBank DB06760: Sacrosidase, biotech, UNII 8A7F670F2Y.',
 'Verified current Sucraid SPL: literal sacrosidase ingredient and the same UNII, single active ingredient, human Rx.',
 paste('Selected named-context label records: 1; adverse-reaction text characters:',row$adverse_text_chars),
 paste('SPL effective time:',effective,'version:',version),
 'Outputs: candidate_label_named_context_selection.csv (local source-context selection); outputs/14c_* summaries.',
 'Original SIDER files, 14a identity registry and 14b strict selection hashes are unchanged.',
 'The invalid CID association is quarantined, not replaced by another guessed CID or fabricated InChIKey.',
 'Historical name context is recovered; historical SIDER chemical/formulation identity remains unresolved.',
 'Keep this enzyme outside the small-molecule structure screen; it may enter a separately labelled side-effect/context comparison.',
 'The recovered label is not a new strict full-key match and must not be pooled silently into the 14b selection.',
 'Current label source text still needs the same PT/context review as other labels. No eEF1A binding claim is made.')
writeLines(report,'outputs/14c_sacrosidase_resolution_summary.txt')
writeLines('PASS: DrugBank name/type/UNII; raw SPL ingredient name/UNII, document type, set ID and AR text; protected input hashes unchanged.',
 'outputs/14c_sacrosidase_validation.txt')
cat(paste0('\nSCRIPT 14c RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/14c_resolve_sacrosidase_context.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
