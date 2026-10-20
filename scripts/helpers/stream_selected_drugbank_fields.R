# Parse requested direct-child fields only for selected primary DrugBank IDs.
# Unselected target sections are discarded before XML construction.
stream_selected_drugbank_fields <- function(path, record_ids, fields, callback,
                                           chunk_lines = 10000L) {
  stopifnot(length(record_ids) > 0L, !anyNA(record_ids), !anyDuplicated(record_ids))
  con <- file(path, 'rt', encoding = 'UTF-8')
  on.exit(close(con), add = TRUE)
  inside <- FALSE; selected <- NA; keep <- FALSE
  pending <- character(); seen <- character(); retained <- character(); n <- 0L
  wanted <- paste0('^  <(', paste(fields, collapse = '|'), ')([[:space:]>]|/>)')
  retain <- function(segment) {
    if (!length(segment)) return(invisible(NULL))
    if (is.na(selected)) {
      ids <- grep('^  <drugbank-id primary="true">[^<]+</drugbank-id>', segment, value = TRUE)
      if (length(ids) != 1L) stop('Expected one primary ID at start of each DrugBank record')
      id <- sub('^  <drugbank-id primary="true">([^<]+)</drugbank-id>.*$', '\\1', ids)
      if (id %in% seen) stop('Repeated primary DrugBank ID')
      seen <<- c(seen, id)
      selected <<- id %in% record_ids
      if (!selected) pending <<- character()
    }
    if (!selected) return(invisible(NULL))
    starts <- grep('^  <[a-z][a-z0-9-]*([[:space:]>]|/>)', segment, perl = TRUE)
    states <- c(keep, grepl(wanted, segment[starts], perl = TRUE))
    pending <<- c(pending, segment[states[findInterval(seq_along(segment), starts) + 1L]])
    keep <<- tail(states, 1L)
  }
  repeat {
    chunk <- readLines(con, n = chunk_lines, warn = FALSE)
    if (!length(chunk)) break
    starts <- grep('^<drug\\s+type=', chunk, perl = TRUE)
    ends <- grep('^</drug>\\s*$', chunk, perl = TRUE)
    pos <- 1L
    for (b in sort(c(starts, ends))) {
      if (b %in% starts) {
        if (inside) stop('Nested top-level DrugBank record')
        inside <- TRUE; selected <- NA; keep <- FALSE
        pending <- chunk[b]; pos <- b + 1L
      } else {
        if (!inside) stop('Unmatched DrugBank record end')
        if (pos < b) retain(chunk[pos:(b - 1L)])
        if (is.na(selected)) stop('Record closed before primary ID')
        if (selected) {
          doc <- xml2::read_xml(paste(c('<drugbank>', pending, chunk[b], '</drugbank>'), collapse = '\n'))
          xml2::xml_ns_strip(doc)
          node <- xml2::xml_find_first(doc, './drug')
          id <- xml2::xml_text(xml2::xml_find_first(node, './drugbank-id[@primary="true"]'))
          if (!id %in% record_ids) stop('Selected record identity mismatch')
          callback(node)
          retained <- c(retained, id)
        }
        n <- n + 1L
        if (n %% 2000L == 0L) message('Scanned ', n, ' records; selected ', length(retained))
        inside <- FALSE; pending <- character(); pos <- b + 1L
      }
    }
    if (inside && pos <= length(chunk)) retain(chunk[pos:length(chunk)])
  }
  if (inside) stop('Incomplete DrugBank XML record')
  if (!setequal(retained, record_ids)) stop('Requested DrugBank records missing from XML')
  list(records_scanned = n, selected_records = length(retained), primary_ids = retained)
}
