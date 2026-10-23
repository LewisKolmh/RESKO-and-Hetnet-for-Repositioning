suppressPackageStartupMessages(library(readr))
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
x <- r('data/curated/RESKO23_candidate_identity_registry.csv')
index <- r('data/curated/sider_compound_index.csv')
q1 <- r('data/curated/RESKO_mechanism_SIDER_primary_group_review_queue.csv')
q2 <- r('data/curated/RESKO_mechanism_SIDER_exploratory_group_review_queue.csv')
expected <- union(q1$flat,q2$flat)
aug <- 'data/curated/RESKO23_augmented_IMiD_group_review_queue.csv'
if(file.exists(aug)) expected <- union(expected,r(aug)$flat)
stopifnot(!anyDuplicated(x$lead_id),
 setequal(x$flat,expected),
 all(paste(x$flat,x$stitch_stereo_id,x$sider_pubchem_cid) %in%
 paste(index$stitch_flat_id,index$stitch_stereo_id,index$sider_pubchem_cid)))
for(i in seq_len(nrow(x))) {
 z <- x[i,]; cid <- z$sider_pubchem_cid
 p <- readRDS(paste0('data/raw/section23_identity/properties_',cid,'.rds'))
 s <- readRDS(paste0('data/raw/section23_identity/synonyms_',cid,'.rds'))
 stopifnot(p$status==200L,s$status %in% c(200L,404L))
 pr <- p$body$PropertyTable$Properties[[1]]
 stopifnot(z$returned_cid==as.character(pr$CID),z$current_inchikey==pr$InChIKey,
 z$current_inchi==pr$InChI,z$molecular_formula==pr$MolecularFormula)
 supported <- if(s$status==404L) FALSE else any(vapply(s$body$InformationList$Information,function(a)
 as.character(a$CID)==z$returned_cid && any(tolower(trimws(as.character(unlist(a$Synonym))))==
 tolower(trimws(z$historical_name))),logical(1)))
 stopifnot((z$historical_name_alias_supported=='TRUE')==supported)
}
links <- r('data/curated/RESKO23_candidate_DrugBank_structure_links.csv')
db <- r('data/curated/drugbank_structure_identity_index.csv')
stopifnot(all(paste(links$record_drugbank_id,links$current_inchikey) %in%
 paste(db$record_drugbank_id,db$record_inchikey)),
 !any(x$historical_formulation_equivalence_proven=='TRUE'))
inv <- r('outputs/16f_target_evidence_handoff_summary.csv')
il <- r('data/curated/RESKO23_existing_inventory_structure_links.csv')
stopifnot(nrow(inv)==200L,all(paste(il$entity_id,il$current_inchikey) %in%
 paste(inv$entity_id,inv$current_inchikey)),!any(il$entity_id=='seed:DB06760'))
h <- r('outputs/23c_input_hashes_summary.csv')
stopifnot(all(unname(tools::md5sum(h$path))==h$md5))
writeLines(c(paste('PASS:',length(unique(x$flat)),'qualifying groups /',nrow(x),'stereo-CID leads retained without pooling.'),
 'PASS: raw cached PubChem properties, names and requested/returned CIDs agree for every lead.',
 'PASS: exact full-key DrugBank and inventory links agree with source; original inventory remains 200.',
 'PASS: historical equivalence never asserted; blocked sacrosidase structure not recovered; input/cache hashes unchanged.'),
 'outputs/23c_identity_validation.txt')
cat(readLines('outputs/23c_identity_validation.txt'),sep='\n')
