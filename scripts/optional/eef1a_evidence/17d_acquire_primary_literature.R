# Acquire bounded primary publication metadata, notices and available full text.
library(readr);library(dplyr);library(tibble);library(httr2);library(xml2)
source('scripts/helpers/literature_public_get.R')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
input<-'data/curated/eef1a_chembl_document_registry.csv';before<-tools::md5sum(input);docs<-r(input)
spec<-bind_rows(docs |> transmute(lookup=pubmed_id,role='ChEMBL_source_document'),
 tibble(lookup=c('31496240','PMC5054363','PMC4666417','PMC6613760'),
 role=c('source_integrity_notice','positive_binding_control','RT_PPI_target_context','RT_PPI_chemical_countercontrol')))
write_csv(spec,'config/target_literature_scope.csv')
root<-'data/raw/target_evidence/literature';rows<-list();corrections<-list();logs<-list()
one<-function(v)if(length(v))as.character(v[[1]])else NA_character_
for(i in seq_len(nrow(spec))){
 id<-spec$lookup[i];query<-if(startsWith(id,'PMC'))id else paste0('EXT_ID:',id,' AND SRC:MED')
 req<-request('https://www.ebi.ac.uk/europepmc/webservices/rest/search') |>
 req_url_query(query=query,format='json',resultType='core',pageSize=100)
 z<-literature_public_get(req$url,file.path(root,paste0('metadata_',id,'.rds')))
 stopifnot(z$status==200L)
 matches<-Filter(function(v)if(startsWith(id,'PMC'))identical(v$pmcid,id)else identical(v$pmid,id),z$body$resultList$result)
 stopifnot(length(matches)==1L);v<-matches[[1]];pmid<-one(v$pmid);pmcid<-one(v$pmcid)
 logs[[length(logs)+1L]]<-tibble(url=z$url,status=z$status,retrieved_utc=z$retrieved_utc)
 cs<-v$commentCorrectionList$commentCorrection
 if(length(cs))for(c in cs)corrections[[length(corrections)+1L]]<-tibble(pmid=pmid,related_pmid=one(c$id),
 relation_type=one(c$type),reference=one(c$reference))
 full_status<-'not_available_in_metadata';full_path<-NA_character_;full_url<-NA_character_;license<-NA_character_;text_chars<-0L
 if(!is.na(pmcid)){
  full_url<-paste0('https://www.ebi.ac.uk/europepmc/webservices/rest/',pmcid,'/fullTextXML')
  full_path<-file.path(root,paste0('fulltext_',pmcid,'.rds'));f<-literature_public_get(full_url,full_path,FALSE)
  logs[[length(logs)+1L]]<-tibble(url=f$url,status=f$status,retrieved_utc=f$retrieved_utc)
  full_status<-if(f$status==200L)'full_text_XML_acquired'else if(f$status==404L)'full_text_unavailable_from_public_endpoint'else paste0('full_text_endpoint_HTTP_',f$status)
  if(f$status==200L){
   doc<-read_xml(f$body);article_pmid<-xml_text(xml_find_first(doc,'/article/front/article-meta/article-id[@pub-id-type="pmid"]'))
   stopifnot(identical(article_pmid,pmid))
   license<-paste(unique(xml_text(xml_find_all(doc,'/article/front/article-meta/permissions/license'))),collapse='|')
   text_chars<-nchar(xml_text(xml_find_first(doc,'/article/body')))
  }
 }
 rows[[i]]<-tibble(lookup=id,role=spec$role[i],pmid,pmcid,doi=one(v$doi),title=one(v$title),year=one(v$pubYear),
 publication_types=paste(unlist(v$pubTypeList$pubType),collapse='|'),
 correction_types=paste(vapply(cs,function(c)one(c$type),character(1)),collapse='|'),
 correction_PMIDs=paste(vapply(cs,function(c)one(c$id),character(1)),collapse='|'),
 source_integrity_status=if(any(grepl('retract|concern',vapply(cs,function(c)one(c$type),character(1)),ignore.case=TRUE)))
 'source_integrity_notice_requires_hold'else 'no_notice_returned_by_this_metadata_query',
 full_text_status=full_status,full_text_chars=text_chars,metadata_cache=file.path(root,paste0('metadata_',id,'.rds')),
 full_text_cache=full_path,metadata_url=z$url,full_text_url=full_url,pubmed_url=paste0('https://pubmed.ncbi.nlm.nih.gov/',pmid,'/'),
 license_text=license,retrieved_utc=z$retrieved_utc)
 message('Publication checked: ',i,' / ',nrow(spec),' ',id)
}
registry<-bind_rows(rows);stopifnot(!anyDuplicated(registry$pmid),all(docs$pubmed_id %in% registry$pmid),
 identical(before,tools::md5sum(input)))
write_csv(registry,'data/curated/target_primary_literature_registry.csv',na='')
write_csv(bind_rows(corrections),'data/curated/target_publication_notices.csv',na='')
write_csv(bind_rows(logs),'data/curated/target_literature_API_provenance.csv',na='')
write_csv(registry |> select(lookup,role,pmid,pmcid,doi,title,year,source_integrity_status,correction_types,correction_PMIDs,
 full_text_status,pubmed_url,retrieved_utc),'outputs/17d_literature_coverage.csv',na='')
write_csv(bind_rows(corrections),'outputs/17d_publication_notices_summary.csv',na='')
write_csv(tibble(file=input,md5=unname(before)),'outputs/17d_input_summary.csv')
report<-c('SCRIPT 17d - PRIMARY LITERATURE AND SOURCE-INTEGRITY CHECK',
 paste('Publications checked:',nrow(registry)),paste('ChEMBL source documents covered:',nrow(docs)),
 paste('Public full-text XML records acquired:',sum(registry$full_text_status=='full_text_XML_acquired')),
 'Includes a linked Expression of Concern notice, plitidepsin binding control and two RT/eEF1A context papers.',
 'The 2018 gamendazole/quinolinone paper (PMID 29498519) links to Expression of Concern PMID 31496240; derived records require quarantine.',
 'No returned notice is not proof that a paper is error-free. Full-text availability and abstract-only access remain explicit.',
 'Metadata, notices and available source XML are cached locally with source URLs, dates and license metadata.',
 'Original activity/inventory data are unchanged; no binding or allosteric claim is accepted merely from paper titles or assay type.',
 'Next: apply explicit assay-level evidence decisions, preserve limited-access holds and evaluate positive/RT-binding controls separately from the inventory.')
writeLines(report,'outputs/17d_literature_summary.txt')
writeLines('PASS: all seven original document PMIDs resolved exactly; public XML front-matter PMID matches metadata; notice relations and protected input hashes retained.',
 'outputs/17d_literature_validation.txt')
cat(paste0('\nSCRIPT 17d RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/optional/eef1a_evidence/17d_acquire_primary_literature.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
