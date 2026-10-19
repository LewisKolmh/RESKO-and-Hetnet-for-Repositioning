# Cached official HTTPS target-source JSON requests with provenance and pagination support.
target_public_get<-function(url,path){
 if(!grepl('^https://(www[.]ebi[.]ac[.]uk/chembl/api/data/|rest[.]uniprot[.]org/uniprotkb/)',url))stop('Unexpected target API host/path')
 if(file.exists(path)){z<-readRDS(path);stopifnot(identical(z$url,url));return(z)}
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);Sys.sleep(.3)
 resp<-httr2::request(url) |> httr2::req_timeout(60) |> httr2::req_retry(max_tries=3,retry_on_failure=TRUE) |>
 httr2::req_error(is_error=function(r)FALSE) |> httr2::req_perform()
 if(httr2::resp_status(resp)!=200L)stop('Target API HTTP ',httr2::resp_status(resp),': ',url)
 z<-list(url=url,status=200L,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),body=httr2::resp_body_json(resp))
 saveRDS(z,path);z
}
