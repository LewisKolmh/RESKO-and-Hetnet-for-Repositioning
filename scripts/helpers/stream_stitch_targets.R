stream_stitch_targets<-function(path,protein_ids){
 con<-gzfile(path,'rt');on.exit(close(con),add=TRUE)
 header<-readLines(con,n=1L,warn=FALSE);stopifnot(identical(header,'chemical\tprotein\texperimental\tprediction\tdatabase\ttextmining\tcombined_score'))
 # Exact tab-delimited protein field, never a substring/gene-name match.
 needles<-paste0('\t',protein_ids,'\t')
 out<-list();n<-0L;blocks<-0L
 repeat{z<-readLines(con,n=100000L,warn=FALSE);if(!length(z))break
  n<-n+length(z);blocks<-blocks+1L;selected<-rep(FALSE,length(z));for(needle in needles)selected<-selected|grepl(needle,z,fixed=TRUE);hit<-z[selected]
  if(length(hit))out[[length(out)+1L]]<-hit
  if(blocks%%25L==0L)message('STITCH rows scanned: ',n)
 }
 list(header=header,lines=unlist(out,use.names=FALSE),rows_scanned=n)
}
