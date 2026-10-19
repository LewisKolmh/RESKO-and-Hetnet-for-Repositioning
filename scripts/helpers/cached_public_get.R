# Cached GET-only public requests with provenance and explicit HTTP handling.
public_cache_get <- function(url,path,json=TRUE) {
 if(!grepl('^https://(pubchem[.]ncbi[.]nlm[.]nih[.]gov|dailymed[.]nlm[.]nih[.]gov|api[.]fda[.]gov)/',url))stop('Unexpected public API host.')
 if(file.exists(path)){
  x<-readRDS(path);if(!identical(x$url,url))stop('Cache URL changed: ',path)
 }else{
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);Sys.sleep(.35)
  response<-httr2::request(url) |> httr2::req_timeout(60) |>
   httr2::req_retry(max_tries=3,retry_on_failure=TRUE) |>
   httr2::req_error(is_error=function(r)FALSE) |> httr2::req_perform()
  status<-httr2::resp_status(response)
  if(!status %in% c(200L,404L))stop('Public API status ',status,': ',url)
  x<-list(url=url,status=status,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
   body=if(json)httr2::resp_body_json(response)else httr2::resp_body_string(response))
  saveRDS(x,path)
 }
 x
}
