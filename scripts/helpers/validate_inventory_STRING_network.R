# Independent raw full-file anchor check and inventory/network consistency checks.
library(readr);library(dplyr)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
anchor<-r('outputs/21b_anchor_identity_summary.csv')
path<-'data/raw/string/v12.5/9606.protein.links.detailed.v12.5.txt.gz'
con<-gzfile(path,'rt');header<-readLines(con,n=1L,warn=FALSE)
fields<-strsplit(header,' ',fixed=TRUE)[[1]];fields<-fields[nzchar(fields)]
rows<-0L;selected<-list()
repeat{
 lines<-readLines(con,n=100000L,warn=FALSE);if(!length(lines))break
 keep<-rep(FALSE,length(lines))
 for(id in anchor$current_STRING_id)keep<-keep|startsWith(lines,paste0(id,' '))|grepl(paste0(' ',id,' '),lines,fixed=TRUE)
 if(any(keep))selected[[length(selected)+1L]]<-lines[keep]
 rows<-rows+length(lines)
}
close(con)
# Base R parser is independent of the readr source parser used by acquisition.
fresh<-read.table(text=paste(c(header,unlist(selected,use.names=FALSE)),collapse='\n'),header=TRUE,
 stringsAsFactors=FALSE,check.names=FALSE,colClasses=c('character','character',rep('integer',length(fields)-2L))) |>
 filter(protein1<protein2) |> arrange(protein1,protein2)
full<-readRDS('data/curated/STRING_human_functional_network_v12.5.rds')
expected<-full$edges |> filter(protein1 %in% anchor$current_STRING_id|protein2 %in% anchor$current_STRING_id) |>
 arrange(protein1,protein2)
stopifnot(rows==full$source_directed_rows,
 identical(as.data.frame(fresh),as.data.frame(expected)))
local<-readRDS('data/curated/inventory_STRING_functional_induced_network_v12.5.rds')
links<-r('data/curated/inventory_current_STRING_protein_links.csv')
coverage<-r('data/curated/inventory_current_STRING_entity_coverage.csv')
legacy<-r('data/curated/STITCH_STRING_historical_bridge_coverage.csv')
stitch<-r('data/curated/STITCH_STRING_protein_mapping_evidence.csv')
stopifnot(nrow(coverage)==200L,!anyDuplicated(coverage$entity_id),
 all(links$current_STRING_id %in% local$nodes$protein),
 all(anchor$current_STRING_id %in% local$nodes$protein),
 all(local$edges$protein1 %in% local$nodes$protein),all(local$edges$protein2 %in% local$nodes$protein),
 !anyDuplicated(local$edges[c('protein1','protein2')]),
 all(is.na(stitch$current_STRING_id[stitch$mapping_status!='unique_historical_accession_bridge_projection'])),
 all(legacy$current_node_candidates[legacy$mapping_status=='unique_historical_accession_bridge_projection']=='1'))
# Every induced source edge is retained with original channel values, including lower-confidence rows.
node_ids<-local$nodes$protein
expected_local<-full$edges[full$edges$protein1 %in% node_ids & full$edges$protein2 %in% node_ids,,drop=FALSE]
stopifnot(identical(as.data.frame(expected_local),as.data.frame(local$edges)))
hashes<-r('outputs/21b_input_summary.csv');stopifnot(identical(unname(tools::md5sum(hashes$file)),hashes$md5))
files<-list.files('scripts',pattern='[.]R$',recursive=TRUE,full.names=TRUE)
for(f in files)invisible(parse(f))
validation<-c(paste('PASS: independent complete raw source scan of',rows,'directed rows reproduces every human eEF1A-anchor edge and score channel.'),
 'PASS: induced edges exactly match source endpoints/channel values; nodes retain all mapped proteins/anchors, including isolates.',
 'PASS: all 200 entities retained; ambiguous historical mappings have no assigned current node; unique bridge acceptance and protected input hashes verified.',
 'PASS: all project R scripts parse. Source/current-release mapping and functional-association scope remain explicit.')
cat(paste(validation,collapse='\n'),'\n')
cat(paste0(paste(validation,collapse='\n'),'\n'),file='outputs/21b_harmonisation_validation.txt',append=TRUE)
cat(paste0('\nSCRIPT 21b - INDEPENDENT VALIDATION\n',paste(validation,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
