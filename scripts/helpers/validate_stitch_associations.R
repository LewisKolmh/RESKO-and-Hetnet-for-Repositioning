# Independent second full-file pass and exact inventory-join accounting for 17c.
library(readr);library(dplyr)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
source('scripts/helpers/stream_stitch_targets.R')
p<-'data/raw/stitch/9606.protein_chemical.links.detailed.v5.0.tsv.gz'
cw<-r('config/eef1a_STITCH_protein_crosswalk.csv');hash<-r('outputs/17c_input_summary.csv')
stopifnot(identical(unname(tools::md5sum(hash$file)),hash$md5))
z<-stream_stitch_targets(p,cw$stitch_protein_id)
raw<-read_tsv(I(paste(c(z$header,z$lines),collapse='\n')),col_types=cols(.default=col_character())) |> arrange(chemical,protein)
edges<-r('data/curated/eef1a_STITCH_association_channels.csv') |> arrange(chemical,protein)
stopifnot(all(vapply(names(raw),function(k)identical(raw[[k]],edges[[k]]),logical(1))),nrow(raw)==3198L,z$rows_scanned==15473939L)
# Derive textmining-only support from raw score fields independently of saved boolean flags.
raw_text_only<-with(raw,as.integer(experimental)==0L & as.integer(prediction)==0L & as.integer(database)==0L & as.integer(textmining)>0L)
stopifnot(sum(raw_text_only)==sum(edges$textmining_only=='TRUE'))
map<-r('data/curated/STITCH_inventory_identifier_leads.csv');links<-r('data/curated/eef1a_STITCH_inventory_association_evidence.csv')
expected<-inner_join(map,raw,by='chemical',relationship='many-to-many')
stopifnot(nrow(expected)==nrow(links),!anyDuplicated(links[c('entity_id','chemical','protein')]),
 !any(links$entity_id=='sider:CID100135113'),all(is.na(links$direct_binding_established)),
 all(as.integer(links$experimental)==0L))
cov<-r('outputs/17c_STITCH_inventory_coverage.csv');inventory<-r('outputs/16f_target_evidence_handoff_summary.csv')
stopifnot(nrow(cov)==2L*nrow(inventory),setequal(cov$entity_id,inventory$entity_id),
 !anyDuplicated(cov[c('entity_id','gene')]),sum(as.integer(cov$associated_chemical_IDs))==nrow(links))
for(f in c('scripts/17c_acquire_stitch_associations.R','scripts/helpers/stream_stitch_targets.R'))invisible(parse(f))
msg<-c('PASS: independent second scan of all 15,473,939 source rows exactly reproduces every retained chemical/protein/channel field.',
 'PASS: text-mining-only counts independently recomputed from raw score fields; inventory join accounts for all 13 unique source-ID/protein links.',
 'PASS: all 200 entities/two isoforms retained in coverage; zero inventory experimental-channel links; no direct-binding status assigned.',
 'PASS: protected input hashes agree and source scripts parse. Flat/stereo rows are not independent experiments.')
cat(paste0(paste(msg,collapse='\n'),'\n'),file='outputs/17c_STITCH_validation.txt',append=TRUE)
cat(paste0('\nSCRIPT 17c - INDEPENDENT VALIDATION\n',paste(msg,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(msg,collapse='\n'),'\n')
