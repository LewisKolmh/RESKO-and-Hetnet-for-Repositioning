# Search the six STITCH-linked inventory entities; acquisition is not evidence acceptance.
library(readr)
library(dplyr)
library(tibble)
library(httr2)
source('scripts/helpers/literature_public_get.R')
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
input <- 'outputs/17c_linked_inventory_summary.csv'; before <- tools::md5sum(input)
leads <- r(input)
terms <- tribble(~source_name,~chemical_query,
 'Fluorouracil','(fluorouracil OR "5-fluorouracil" OR "5-FU")',
 'Estradiol','(estradiol OR "17beta-estradiol" OR "17-beta-estradiol")',
 'Doxorubicin','(doxorubicin OR adriamycin)',
 'Dexamethasone','dexamethasone',
 'Testosterone enanthate','("testosterone enanthate" OR "testosterone enantate")',
 'Epirubicin','(epirubicin OR "4-epidoxorubicin")')
stopifnot(setequal(leads$source_name,terms$source_name))
queries <- leads |> select(entity_id,source_name) |> left_join(terms,by='source_name',relationship='one-to-one') |>
 mutate(query=paste0('(EEF1A1 OR EEF1A2 OR eEF1A OR "elongation factor 1 alpha" OR "elongation factor 1-alpha") AND ',chemical_query,
 ' AND FIRST_PDATE:[1900-01-01 TO 2026-10-20]'),
 query_scope='Europe PMC indexed full text and metadata; no human filter, species must be reviewed',
 identity_scope='Exact named compound aliases; no testosterone-parent substitution for enanthate')
write_csv(queries,'config/18a_association_literature_queries.csv')
root <- 'data/raw/target_evidence/association_literature'; hits <- list(); logs <- list()
scalar <- function(v) if(length(v)) as.character(v[[1]]) else NA_character_
for(i in seq_len(nrow(queries))){
 cursor <- '*'; page <- 0L; acquired <- 0L; first_total <- NULL; seen <- character()
 repeat {
  page <- page+1L
  if(page>100L) stop('Pagination safety boundary reached; do not report complete coverage')
  url <- request('https://www.ebi.ac.uk/europepmc/webservices/rest/search') |>
   req_url_query(query=queries$query[i],format='json',resultType='core',pageSize=100,cursorMark=cursor)
  path <- file.path(root,paste0(gsub(':','_',queries$entity_id[i]),'_page',page,'.rds'))
  z <- literature_public_get(url$url,path);stopifnot(z$status==200L)
  total <- as.integer(z$body$hitCount)
  if(is.null(first_total))first_total <- total
  if(total!=first_total)stop('Search hit count changed during pagination; rerun fresh acquisition')
  records <- z$body$resultList$result
  for(v in records){
   key <- paste(scalar(v$source),scalar(v$id),sep=':')
   if(key %in% seen)stop('Duplicate search record across pages')
   seen <- c(seen,key)
   hits[[length(hits)+1L]] <- tibble(entity_id=queries$entity_id[i],source_name=queries$source_name[i],
    literature_key=key,pmid=scalar(v$pmid),pmcid=scalar(v$pmcid),doi=scalar(v$doi),title=scalar(v$title),year=scalar(v$pubYear),
    abstract=scalar(v$abstractText),publication_types=paste(unlist(v$pubTypeList$pubType),collapse='|'),
    notice_types=paste(vapply(v$commentCorrectionList$commentCorrection,function(c)scalar(c$type),character(1)),collapse='|'),
    review_status='unreviewed_search_hit_not_binding_evidence',query=queries$query[i],source_url=z$url,
    retrieved_utc=z$retrieved_utc,cache=path)
  }
  acquired <- acquired+length(records)
  logs[[length(logs)+1L]] <- tibble(entity_id=queries$entity_id[i],source_name=queries$source_name[i],
   page=page,total_hits=total,page_records=length(records),url=z$url,retrieved_utc=z$retrieved_utc,cache=path)
  if(acquired==total)break
  next_cursor <- scalar(z$body$nextCursorMark)
  if(!length(records)||is.na(next_cursor)||identical(next_cursor,cursor)||acquired>total)stop('Incomplete search pagination')
  cursor <- next_cursor
 }
 message(queries$source_name[i],': ',acquired,' search hits, all pages acquired')
}
records <- bind_rows(hits); provenance <- bind_rows(logs)
stopifnot(!anyDuplicated(records[c('entity_id','literature_key')]),identical(before,tools::md5sum(input)))
write_csv(records,'data/curated/association_literature_search_hits.csv',na='')
write_csv(provenance,'data/curated/association_literature_search_provenance.csv')
coverage <- provenance |> group_by(entity_id,source_name) |> summarise(total_hits=first(total_hits),
 acquired_records=sum(page_records),pages=n(),.groups='drop')
stopifnot(nrow(coverage)==6L,all(coverage$total_hits==coverage$acquired_records))
write_csv(coverage,'outputs/18a_search_coverage.csv')
write_csv(tibble(file=input,md5=unname(before)),'outputs/18a_input_summary.csv')
report <- c('SCRIPT 18a - SIX-LEAD PRIMARY LITERATURE SEARCH ACQUISITION',
 paste('Compound/publication search hits:',nrow(records)),paste('Distinct publication keys:',n_distinct(records$literature_key)),
 'All six exact-compound query scopes fully paginated; search cutoff 20 October 2026.',
 'Indexed full-text co-occurrence can include reference lists, housekeeping genes, combination therapies and nonhuman experiments.',
 'Search hits are not evidence of direct binding and do not establish the original source of STITCH 5 edges.',
 'Testosterone enanthate remains distinct from testosterone; epirubicin remains distinct from doxorubicin.',
 'Raw core metadata/abstract caches and source-rich hits remain local; queries and compact coverage are versioned.',
 'No BH adjustment or binding ranking applies to these search counts. Next: primary-source adjudication and provenance tracing.')
writeLines(report,'outputs/18a_search_summary.txt')
writeLines('PASS: six queries, exact compound names, complete pagination, no duplicate entity/publication keys and input hash unchanged.',
 'outputs/18a_search_validation.txt')
cat(paste0('\nSCRIPT 18a RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/optional/eef1a_evidence/18a_search_association_literature.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
