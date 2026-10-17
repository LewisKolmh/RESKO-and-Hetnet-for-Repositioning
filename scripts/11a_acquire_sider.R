# Acquire the original SIDER 4.1 aggregate tables, retaining raw files and provenance.
library(httr2);library(readr);library(dplyr);library(tibble)
if (!all(dir.exists(c("scripts", "config", "data")))) stop('Run from the project root.')
root<-'data/raw/sider';dir.create(root,recursive=TRUE,showWarnings=FALSE)
files<-c('README','drug_names.tsv','drug_atc.tsv','meddra_all_se.tsv.gz',
 'meddra_freq.tsv.gz','meddra_all_indications.tsv.gz','meddra.tsv.gz')
base<-'https://sideeffects.embl.de/media/download/'
meta_path<-file.path(root,'acquisition_manifest.csv')
old<-if(file.exists(meta_path)) read_csv(meta_path,show_col_types=FALSE,col_types=cols(.default=col_character())) else tibble(file=character())
rows<-list()
for(f in files){
 path<-file.path(root,f);u<-paste0(base,f);downloaded<-FALSE
 if(!file.exists(path)){
   temp<-paste0(path,'.partial')
   request(u)|>req_timeout(120)|>req_retry(max_tries=3,retry_on_failure=TRUE)|>req_perform(path=temp)
   if(endsWith(f,'.gz')){
     con<-file(temp,'rb');magic<-readBin(con,'raw',n=2);close(con)
     if(!identical(magic,as.raw(c(31,139)))) stop('Not gzip: ',f)
   } else if(grepl('<!DOCTYPE|<html',paste(readLines(temp,n=3,warn=FALSE),collapse=''),ignore.case=TRUE))stop('HTML instead of data: ',f)
   if(!file.rename(temp,path)) stop('Could not finalise ',f)
   downloaded<-TRUE
 }
 if(endsWith(f,'.gz')){con<-gzfile(path,'rt');first<-readLines(con,n=2,warn=FALSE);close(con)} else first<-readLines(path,n=2,warn=FALSE)
 if(!length(first)) stop('Empty SIDER input: ',f)
 previous<-old |> filter(file==f)
 md5<-unname(tools::md5sum(path))
 if(nrow(previous)&&!identical(previous$md5[[1]],md5))stop('Cached raw input changed: ',f)
 rows[[f]]<-tibble(file=f,url=u,release_version='4.1',release_date='2015-10-21',
   retrieved_utc=if(nrow(previous))previous$retrieved_utc[[1]] else format(Sys.time(),tz='UTC',usetz=TRUE),
   provenance_note=if(nrow(previous))previous$provenance_note[[1]] else if(downloaded)'downloaded by 11a' else 'downloaded during same-session schema inspection; registered by 11a',
   file_bytes=as.character(file.info(path)$size),md5=md5,local_file=path)
 message('SIDER ready: ',f)
}
meta<-bind_rows(rows);write_csv(meta,meta_path)
manifest<-read_csv('config/database_manifest.csv',col_types=cols(.default=col_character()))
i<-manifest$source_id=='sider';manifest$status[i]<-'aggregate_release_acquired'
manifest$release_version[i]<-'4.1 (2015-10-21)';manifest$retrieved_date[i]<-format(Sys.Date())
manifest$download_url[i]<-base;manifest$local_file[i]<-meta_path
manifest$licence_note[i]<-'Publisher download page states CC BY-SA 4.0 except where otherwise noted; preserve original attribution and per-file notices.'
write_csv(manifest,'config/database_manifest.csv',na='')
report<-c('SCRIPT 11a - SIDER 4.1 ACQUISITION',paste('Raw files acquired/registered:',nrow(meta)),
 paste('Compressed/raw bytes:',sum(as.numeric(meta$file_bytes))),
 'Inputs: original publisher URLs. Outputs: data/raw/sider/ and acquisition_manifest.csv.',
 'Contains names, ATC, aggregate side effects, frequencies, indication mentions and term mapping.',
 'Label-level source files are not acquired in this step; aggregate rows are not independent trials.',
 'SIDER 4.1 was released in 2015. Missing drugs/events do not establish safety or lack of activity.',
 'MedDRA PT and LLT are distinct; UMLS CUIs remain identifiers, not ontology acquisition.',
 'The raw README specifies each column; retain placebo rows and frequency ranges explicitly.',
 'Next: map seed identities with dplyr, preserving flat versus stereo IDs and coverage limitations.')
writeLines(report,'outputs/11a_sider_acquisition_summary.txt')
cat(paste0('\nSCRIPT 11a RUN - ',format(Sys.time(),tz='Europe/Rome',usetz=TRUE),'\n',
 'Script: scripts/11a_acquire_sider.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
