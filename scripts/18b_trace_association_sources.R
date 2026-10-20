# Bounded source tracing: supplementary annotations are not assumed to be STITCH sources.
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(httr2)
library(xml2)
source('scripts/helpers/literature_public_get.R')
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
input <- 'data/curated/association_literature_search_hits.csv'; before <- tools::md5sum(input)
hits <- r(input)
chemical_patterns <- c('Epirubicin'='epirubicin|epidoxorubicin','Fluorouracil'='fluorouracil|5.fu',
 'Estradiol'='estradiol','Doxorubicin'='doxorubicin|adriamycin','Dexamethasone'='dexamethasone',
 'Testosterone enanthate'='testosterone enanthate|testosterone enantate')
gene_pattern <- 'eef1a[12]?|elongation factor.?1.?alpha|elongation factor 1.?α'
triage <- hits |> mutate(title_abstract=paste(title,coalesce(abstract,'')),
 gene_in_title_abstract=str_detect(title_abstract,regex(gene_pattern,ignore_case=TRUE)),
 compound_in_title_abstract=mapply(function(s,n)grepl(chemical_patterns[[n]],s,ignore.case=TRUE),title_abstract,source_name),
 triage_status=if_else(gene_in_title_abstract & compound_in_title_abstract,
 'priority_abstract_review_required','full_text_context_pending_not_rejected'))
priority <- triage |> filter(triage_status=='priority_abstract_review_required')
write_csv(triage |> select(-title_abstract),'data/curated/association_literature_triage.csv',na='')
write_csv(priority |> select(-title_abstract),'data/curated/association_literature_priority.csv',na='')
root <- 'data/raw/target_evidence/association_sources';dir.create(root,recursive=TRUE,showWarnings=FALSE)
# Restricted public-source cache; application errors must never count as source evidence.
public_html <- function(url,path){
 if(!grepl('^(https://drugmap[.]idrblab[.]net/data/dot/details/OT9Z23K5|http://stitch[.]embl[.]de/cgi/)',url))stop('Unexpected source URL')
 if(file.exists(path)){z<-readRDS(path);stopifnot(identical(z$url,url));return(z)}
 response<-request(url)|>req_timeout(45)|>req_retry(max_tries=2,retry_on_failure=TRUE)|>req_error(is_error=function(r)FALSE)|>req_perform()
 z<-list(url=url,status=resp_status(response),retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),body=resp_body_string(response))
 saveRDS(z,path);z
}
drugmap_url <- 'https://drugmap.idrblab.net/data/dot/details/OT9Z23K5'
z <- public_html(drugmap_url,file.path(root,'drugmap_EEF1A2.rds'));stopifnot(z$status==200L)
h <- read_html(z$body)
# Exact reference row numbers on the acquired page; verify compound-row link targets.
spec <- tribble(~entity_id,~source_name,~reference_number,~expected_pmid,
 'seed:DB00997','Doxorubicin','39','29803840',
 'seed:DB00783','Estradiol','41','20106945',
 'seed:DB01234','Dexamethasone','48','16144937')
annotations <- list()
for(i in seq_len(nrow(spec))){
 ref <- xml_find_all(h,paste0('//tr[th[normalize-space(.)="',spec$reference_number[i],'"]]'))
 stopifnot(length(ref)==1L)
 link <- xml_attr(xml_find_first(ref,'.//a[contains(@href,"pubmed")]'),'href')
 pmid <- str_extract(link,'[0-9]{6,}')
 stopifnot(identical(pmid,spec$expected_pmid[i]))
 rows <- xml_find_all(h,'//tr[td]')
 matched <- Filter(function(x){cells<-xml_find_all(x,'./td');length(cells)>1L && any(trimws(xml_text(cells))==spec$source_name[i])},as.list(rows))
 stopifnot(length(matched)==1L)
 row_text <- str_squish(xml_text(matched[[1]]))
 refid <- xml_attr(ref,'id')
 stopifnot(any(xml_attr(xml_find_all(matched[[1]],'.//a'),'href')==paste0('#',refid)),grepl('expression',row_text))
 annotations[[i]] <- tibble(entity_id=spec$entity_id[i],source_name=spec$source_name[i],gene='EEF1A2',
 supplementary_source='DrugMAP',annotation_text=row_text,reference_number=spec$reference_number[i],pmid,
 primary_reference=str_squish(xml_text(ref)),source_url=drugmap_url,retrieved_utc=z$retrieved_utc,
 relation='expression_annotation_not_binding',STITCH_original_source_link_verified=FALSE)
}
annotations <- bind_rows(annotations)
write_csv(annotations,'data/curated/association_supplementary_source_links.csv',na='')
# Reproducible direct legacy-edge lookup. HTTP 200 alone is not success.
edges <- r('data/curated/eef1a_STITCH_inventory_association_evidence.csv')
attempts <- list()
for(i in seq_len(nrow(edges))){
 url <- request('http://stitch.embl.de/cgi/show_edge_details.pl') |>
  req_url_query(node1=paste0('-1.',edges$chemical[i]),node2=edges$protein[i])
 z <- public_html(url$url,file.path(root,paste0('stitch_edge_',i,'.rds')))
 text <- str_squish(xml_text(read_html(z$body)))
 status <- if(z$status!=200L)paste0('HTTP_',z$status)else if(grepl('not provided any input|error',text,ignore.case=TRUE))
 'application_error_not_valid_evidence'else 'unverified_response_requires_manual_review'
 attempts[[i]] <- tibble(entity_id=edges$entity_id[i],chemical=edges$chemical[i],protein=edges$protein[i],
 status=status,http_status=z$status,source_url=z$url,retrieved_utc=z$retrieved_utc,
 original_study_provenance='unresolved_no_verified_edge_study_links')
}
write_csv(bind_rows(attempts),'data/curated/association_STITCH_provenance_attempts.csv')
# Resolve supplementary reference PMIDs exactly; annotations require primary verification.
metadata <- list()
for(id in unique(annotations$pmid)){
 url <- request('https://www.ebi.ac.uk/europepmc/webservices/rest/search') |>
 req_url_query(query=paste0('EXT_ID:',id,' AND SRC:MED'),format='json',resultType='core',pageSize=100)
 z <- literature_public_get(url$url,file.path(root,paste0('metadata_',id,'.rds')));stopifnot(z$status==200L)
 vs <- Filter(function(v)identical(v$pmid,id),z$body$resultList$result);stopifnot(length(vs)==1L);v<-vs[[1]]
 one <- function(v)if(length(v))as.character(v[[1]])else NA_character_
 metadata[[length(metadata)+1L]] <- tibble(pmid=id,pmcid=one(v$pmcid),doi=one(v$doi),title=one(v$title),
 abstract=one(v$abstractText),source_url=z$url,retrieved_utc=z$retrieved_utc)
}
write_csv(bind_rows(metadata),'data/curated/association_supplementary_primary_metadata.csv',na='')
coverage <- triage |> count(source_name,triage_status,name='search_records')
write_csv(coverage,'outputs/18b_triage_coverage.csv')
write_csv(annotations |> select(entity_id,source_name,gene,pmid,relation,STITCH_original_source_link_verified),
 'outputs/18b_supplementary_links_summary.csv')
write_csv(bind_rows(attempts) |> count(status,http_status,name='source_edge_attempts'),'outputs/18b_STITCH_provenance_summary.csv')
stopifnot(nrow(triage)==nrow(hits),nrow(priority)==15L,nrow(annotations)==3L,
 nrow(edges)==13L,identical(before,tools::md5sum(input)))
write_csv(tibble(file=input,md5=unname(before)),'outputs/18b_input_summary.csv')
report <- c('SCRIPT 18b - ASSOCIATION SOURCE TRACING AND LITERATURE TRIAGE',
 paste('Priority title/abstract co-occurrence records:',nrow(priority)),
 paste('Remaining search records explicitly pending full-text context:',nrow(triage)-nrow(priority)),
 'Lexical triage is a reading queue, not a binding classifier; index tokenization can broaden named-alias query matches.',
 'Three DrugMAP EEF1A2 expression annotations traced to exact primary PMIDs: doxorubicin 29803840, estradiol 20106945, dexamethasone 16144937.',
 'These are supplementary expression links, not verified STITCH 5 upstream records or physical-binding measurements.',
 'All 13 direct STITCH edge-details attempts are retained with response validation; no original study provenance established.',
 'Additional interactive probes also returned blank edge popups; failures do not indicate absent evidence.',
 'Outputs: source-rich triage, annotations, primary metadata and provenance attempts locally; compact coverage in outputs/18b_*.',
 'Next: manually adjudicate priority abstracts and distinguish exact compound, human isoform, expression, combination therapy and nonhuman context.')
writeLines(report,'outputs/18b_source_summary.txt')
writeLines('PASS: every search record status-accounted; exact supplementary reference PMIDs/compound-row anchors verified; STITCH HTTP/application status retained; input unchanged.',
 'outputs/18b_source_validation.txt')
cat(paste0('\nSCRIPT 18b RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/18b_trace_association_sources.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
