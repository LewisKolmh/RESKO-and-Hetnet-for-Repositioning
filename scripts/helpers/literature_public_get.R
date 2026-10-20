# Official Europe PMC literature requests; raw source caches stay local.
literature_public_get<-function(url,path,json=TRUE){
 if(!startsWith(url,'https://www.ebi.ac.uk/europepmc/webservices/rest/'))stop('Unexpected literature API URL')
 if(file.exists(path)){z<-readRDS(path);stopifnot(identical(z$url,url));return(z)}
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);Sys.sleep(.35)
 response<-httr2::request(url) |> httr2::req_timeout(60) |> httr2::req_retry(max_tries=3,retry_on_failure=TRUE) |>
 httr2::req_error(is_error=function(r)FALSE) |> httr2::req_perform()
 status<-httr2::resp_status(response);if(!status %in% c(200L,404L,500L,502L,503L,504L))stop('Literature API HTTP ',status)
 z<-list(url=url,status=status,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
 body=if(json&&status==200L)httr2::resp_body_json(response)else httr2::resp_body_string(response))
 saveRDS(z,path);z
}
