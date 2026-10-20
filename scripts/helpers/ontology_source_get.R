# Frozen public ontology/annotation acquisition; never refresh an existing snapshot silently.
ontology_source_get<-function(url,path){
 allowed<-c('https://current.geneontology.org/','https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/',
 'https://bioconductor.org/packages/3.23/data/annotation/src/contrib/',
 'https://fms.alliancegenome.org/download/')
 if(!any(vapply(allowed,function(prefix)startsWith(url,prefix),logical(1))))stop('Unexpected ontology source URL')
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);metadata<-paste0(path,'.metadata.rds')
 if(!file.exists(path)){
  message('Acquiring ',basename(path))
  response<-httr2::request(url)|>httr2::req_timeout(180)|>httr2::req_retry(max_tries=3,retry_on_failure=TRUE)|>
   httr2::req_perform(path=paste0(path,'.partial'))
  stopifnot(httr2::resp_status(response)==200L,file.rename(paste0(path,'.partial'),path))
  saveRDS(list(url=url,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
   last_modified=httr2::resp_header(response,'last-modified'),file_bytes=file.info(path)$size,
   md5=unname(tools::md5sum(path))),metadata)
 }
 stopifnot(file.exists(metadata));m<-readRDS(metadata)
 stopifnot(identical(m$url,url),identical(m$md5,unname(tools::md5sum(path))))
 tibble::tibble(url=url,local_file=path,file_bytes=m$file_bytes,md5=m$md5,
 retrieved_utc=m$retrieved_utc,last_modified=if(length(m$last_modified))m$last_modified else NA_character_)
}
