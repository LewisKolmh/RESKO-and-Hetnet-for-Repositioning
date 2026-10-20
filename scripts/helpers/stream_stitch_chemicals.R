# Exact first-field chemical filtering; all protein neighbours and source scores retained.
stream_stitch_chemicals <- function(path, chemical_ids, chunk_lines = 100000L) {
  stopifnot(length(chemical_ids)>0L,!anyNA(chemical_ids))
  con<-gzfile(path,'rt');on.exit(close(con),add=TRUE)
  header<-readLines(con,n=1L,warn=FALSE)
  stopifnot(identical(header,'chemical\tprotein\texperimental\tprediction\tdatabase\ttextmining\tcombined_score'))
  rows<-0L;blocks<-0L;out<-list()
  repeat {
    lines<-readLines(con,n=chunk_lines,warn=FALSE);if(!length(lines))break
    first_tab<-regexpr('\t',lines,fixed=TRUE)
    if(any(first_tab<1L))stop('Malformed STITCH row')
    selected<-substr(lines,1L,first_tab-1L) %in% chemical_ids
    if(any(selected))out[[length(out)+1L]]<-lines[selected]
    rows<-rows+length(lines);blocks<-blocks+1L
    if(blocks%%25L==0L)message('STITCH chemical scan: ',rows,' source rows')
  }
  list(header=header,lines=unlist(out,use.names=FALSE),rows_scanned=rows)
}
