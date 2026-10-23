# Verify every structure link and supplied mol block against the raw ChEMBL snapshot.
suppressPackageStartupMessages(library(readr))
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
x <- r('data/curated/RESKO24_chemistry_entity_coverage.csv')
inv <- r('outputs/16f_target_evidence_handoff_summary.csv')
rec <- readRDS('data/curated/RESKO24_ChEMBL_structure_records.rds')
links <- r('data/curated/RESKO24_entity_ChEMBL_structure_links.csv')
stopifnot(nrow(inv)==200L,!anyDuplicated(x$entity_id),!anyDuplicated(rec$chembl_id),
 setequal(x$entity_id[x$scope=='original_inventory'],inv$entity_id),
 all(paste(x$entity_id[x$scope=='original_inventory'],x$expected_inchikey[x$scope=='original_inventory']) %in%
 paste(inv$entity_id,inv$current_inchikey)),
 !any(links$entity_id %in% x$entity_id[x$structure_mapping_blocked=='TRUE']),
 !any(x$historical_equivalence_proven=='TRUE'))
logs <- r('data/curated/RESKO24_ChEMBL_API_provenance.csv')
source <- list()
for(p in logs$cache_path) {
 z <- readRDS(p);stopifnot(z$status==200L,is.null(z$body$page_meta[['next']]))
 for(m in z$body$molecules) source[[m$molecule_chembl_id]] <- m
}
stopifnot(setequal(names(source),rec$chembl_id))
scalar <- function(z) if(length(z)) as.character(z[[1]]) else NA_character_
normal <- function(z) gsub('\r\n','\n',z,fixed=TRUE)
for(i in seq_len(nrow(rec))) {
 z <- rec[i,];m <- source[[z$chembl_id]];s <- m$molecule_structures
 stopifnot(identical(z$standard_inchikey,scalar(s$standard_inchi_key)),
 identical(z$standard_inchi,scalar(s$standard_inchi)),
 identical(z$canonical_smiles,scalar(s$canonical_smiles)),
 identical(z$molfile,scalar(s$molfile)),
 identical(z$parent_chembl_id,scalar(m$molecule_hierarchy$parent_chembl_id)))
 if(!is.na(z$molfile)&&nzchar(z$molfile)) {
  p <- paste0('data/raw/chemistry24/sdf/',z$chembl_id,'.sdf')
  lines <- readLines(p,warn=FALSE)
  stopifnot(tail(lines,1)=='$$$$',
   identical(paste(lines[-length(lines)],collapse='\n'),
    sub('[\r\n]+$','',normal(z$molfile))))
 }
}
for(i in seq_len(nrow(links))) {
 z <- links[i,];match <- rec[rec$chembl_id==z$chembl_id,]
 stopifnot(nrow(match)==1L,z$expected_inchikey==match$standard_inchikey)
}
for(i in seq_len(nrow(x))) {
 z <- x[i,];stopifnot(as.integer(z$ChEMBL_records)==sum(links$entity_id==z$entity_id))
}
a <- readRDS('data/raw/chemistry24/status_before.rds')
b <- readRDS('data/raw/chemistry24/status_after.rds')
stopifnot(identical(a$body$chembl_db_version,b$body$chembl_db_version))
h <- r('outputs/24a_input_hashes_summary.csv')
stopifnot(all(unname(tools::md5sum(h$path))==h$md5))
invisible(lapply(list.files('scripts',pattern='\\.R$',recursive=TRUE,full.names=TRUE),parse))
writeLines(c('PASS: original 200 entities retained; no blocked structure linked; original source keys unchanged.',
 paste('PASS:',nrow(rec),'ChEMBL molecule records match raw source keys, SMILES, InChIs, parent IDs and mol blocks.'),
 'PASS: every entity/record link is an exact full-key match; every SDF reproduces its supplied mol block.',
 'PASS: per-entity coverage counts, release consistency, input/cache hashes and all R script parses.',
 'LIMIT: source consistency is not independent chemical parsing, stereochemical validation or a fingerprint result.'),
 'outputs/24a_structure_validation.txt')
cat(readLines('outputs/24a_structure_validation.txt'),sep='\n')
