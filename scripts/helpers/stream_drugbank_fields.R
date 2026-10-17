# Stream selected direct-child fields from this indented DrugBank XML export.
# Callback receives one compact XML drug node. Large unselected children are
# discarded before accumulation, avoiding the original NADH-record slowdown.
stream_drugbank_fields <- function(path, fields, callback) {
  con <- file(path,"rt",encoding="UTF-8");on.exit(close(con),add=TRUE)
  pending <- character();inside <- FALSE;keep <- FALSE;n <- 0L
  wanted <- paste0("^  <(",paste(fields,collapse="|"),")([[:space:]>]|/>)")
  retain <- function(segment) {
    starts <- grep("^  <[a-z][a-z0-9-]*([[:space:]>]|/>)",segment,perl=TRUE)
    states <- c(keep,grepl(wanted,segment[starts],perl=TRUE))
    pending <<- c(pending,segment[states[findInterval(seq_along(segment),starts)+1L]])
    keep <<- tail(states,1L)
  }
  repeat {
    chunk <- readLines(con,n=10000,warn=FALSE)
    if (!length(chunk)) break
    starts <- grep("^<drug\\s+type=",chunk,perl=TRUE)
    ends <- grep("^</drug>\\s*$",chunk,perl=TRUE);pos <- 1L
    for (b in sort(c(starts,ends))) {
      if (b %in% starts) {
        if (inside) stop("Nested top-level drug record.")
        pending <- chunk[b];inside <- TRUE;keep <- FALSE;pos <- b+1L
      } else {
        if (!inside) stop("Unmatched drug record end.")
        if (pos<b) retain(chunk[pos:(b-1L)])
        pending <- c(pending,chunk[b])
        doc <- xml2::read_xml(paste('<drugbank xmlns="http://www.drugbank.ca">',
          paste(pending,collapse="\n"),'</drugbank>',sep="\n"))
        xml2::xml_ns_strip(doc)
        callback(xml2::xml_find_first(doc,"./drug"))
        n <- n+1L
        if (n%%1000L==0L) message("Indexed ",n," DrugBank records")
        pending <- character();inside <- FALSE;pos <- b+1L
      }
    }
    if (inside && pos<=length(chunk)) retain(chunk[pos:length(chunk)])
  }
  if (inside) stop("Incomplete XML record.")
  invisible(n)
}
