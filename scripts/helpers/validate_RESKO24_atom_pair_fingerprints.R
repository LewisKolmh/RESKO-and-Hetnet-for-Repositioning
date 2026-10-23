# Independent checks in R: raw graph counts, bit membership and every Tanimoto cell.
suppressPackageStartupMessages(library(readr))
source('scripts/helpers/prepare_R_chemistry_runtime.R')
source('scripts/helpers/RESKO_atom_pair_chemistry.R')
prepare_R_chemistry_runtime(FALSE)
bits <- readRDS('data/curated/RESKO24_atom_pair_2048_bits.rds')
ap <- readRDS('data/curated/RESKO24_atom_pair_descriptors.rds')
t <- readRDS('data/curated/RESKO24_structure_Tanimoto.rds')
em <- readRDS('data/curated/RESKO24_entity_Tanimoto.rds')
q <- read_csv('data/curated/RESKO24_atom_pair_structure_QC.csv',show_col_types=FALSE)
e <- read_csv('data/curated/RESKO24_atom_pair_entity_coverage.csv',show_col_types=FALSE)
source_coverage <- read_csv('data/curated/RESKO24_chemistry_entity_coverage.csv',show_col_types=FALSE)
v <- read_csv('data/curated/RESKO24_atom_pair_vocabulary.csv',
 col_types=cols(.default=col_character()))$descriptor
stopifnot(ncol(bits)==2048L,!anyDuplicated(rownames(bits)),!anyDuplicated(v),
 identical(colnames(bits),v),all(bits %in% c(0,1)),all(rowSums(bits)>0L))
# Reconstruct features without desc2fp, including its exact numeric-to-string convention.
manual_bits <- t(vapply(rownames(bits),function(id)
 as.integer(v %in% as.character(ap[[id]])),integer(length(v))))
stopifnot(all(bits==manual_bits))
# Full independent set arithmetic, including numerator and denominator.
for(i in seq_len(nrow(bits))) for(j in seq_len(nrow(bits))) {
 a <- which(bits[i,]==1);b <- which(bits[j,]==1)
 numerator <- length(intersect(a,b));denominator <- length(union(a,b))
 stopifnot(t$intersection[i,j]==numerator,t$union[i,j]==denominator,
 abs(t$similarity[i,j]-numerator/denominator)<1e-12)
}
stopifnot(all(is.finite(t$similarity)),all(t$similarity>=0 & t$similarity<=1),
 all(diag(t$similarity)==1),isTRUE(all.equal(t$similarity,t(t$similarity))))
# Read raw V2000 counts/bond endpoints using fixed columns, not the ChemmineR parser.
for(i in seq_len(nrow(q))) {
 z <- q[i,];raw <- readLines(paste0('data/raw/chemistry24/sdf/',z$chembl_id,'.sdf'),warn=FALSE)
 na <- as.integer(substr(raw[4],1,3));nb <- as.integer(substr(raw[4],4,6))
 stopifnot(na==z$atom_count,nb==z$bond_count)
 parent <- seq_len(na)
 root <- function(k) {while(parent[k]!=k) k <- parent[k];k}
 if(nb>0L) for(line in raw[4L+na+seq_len(nb)]) {
  a <- as.integer(substr(line,1,3));b <- as.integer(substr(line,4,6))
  stopifnot(a>=1L,b>=1L,a<=na,b<=na,a!=b)
  ra <- root(a);rb <- root(b);if(ra!=rb) parent[ra] <- rb
 }
 components <- length(unique(vapply(seq_len(na),root,integer(1))))
 if(!is.na(z$components))stopifnot(components==z$components)
 if(z$status=='fingerprint_available')stopifnot(components==1L,na>=2L,z$active_bits>0L)
}
stopifnot(nrow(em)==204L,ncol(em)==204L,nrow(e)==204L,
 identical(rownames(em),source_coverage$entity_id),identical(colnames(em),source_coverage$entity_id))
for(i in seq_len(nrow(e))) {
 z <- e[i,]
 if(z$entity_status!='fingerprint_available') {
  stopifnot(all(is.na(em[z$entity_id,])),all(is.na(em[,z$entity_id])))
 } else {
  ready <- e[e$entity_status=='fingerprint_available',]
  expected <- t$similarity[z$representative_chembl_id,ready$representative_chembl_id]
  stopifnot(all(abs(em[z$entity_id,ready$entity_id]-expected)<1e-12))
 }
}
blocked <- source_coverage$entity_id[source_coverage$structure_mapping_blocked];stopifnot(length(blocked)==1L,all(is.na(em[blocked,])))
# Meaningful arithmetic and disconnected-graph fixtures.
fixture <- rbind(a=c(1,1,0,0),b=c(0,1,1,0),c=c(0,0,1,1),empty=c(0,0,0,0))
f <- RESKO_binary_Tanimoto(fixture)$similarity
stopifnot(abs(f['a','b']-1/3)<1e-12,f['a','c']==0,f['a','a']==1,
 is.na(f['empty','empty']),f['a','empty']==0,
 RESKO_component_count(4L,rbind(c(1,2),c(3,4)))==2L,
 RESKO_component_count(3L,rbind(c(1,2),c(2,3)))==1L,
 inherits(try(RESKO_component_count(2L,matrix(c(1,3),nrow=1)),silent=TRUE),'try-error'))
h <- read_csv('outputs/24b_input_hashes_summary.csv',show_col_types=FALSE)
stopifnot(all(unname(tools::md5sum(h$path))==h$md5))
# Scientific workflow scripts are R; package-native code under ignored data is separate.
files <- list.files('scripts',recursive=TRUE,full.names=TRUE)
stopifnot(!any(grepl('\\.(py|sh|js|ts|ipynb|jl)$',files,ignore.case=TRUE)))
rscripts <- files[grepl('\\.R$',files)];invisible(lapply(rscripts,parse))
writeLines(c(paste('PASS:',nrow(bits),'2048-bit fingerprints independently rebuilt from frozen AP descriptors/vocabulary.'),
 paste('PASS:',nrow(bits)*(nrow(bits)-1L)/2L,'unique structure-pair scores and all diagonal/directed cells verified by set arithmetic.'),
 'PASS: all raw SDF atom/bond counts and graph components agree; no disconnected graph fingerprinted.',
 'PASS: full 204-entity matrix/projection agrees; missing/held rows and diagonals stay NA.',
 'PASS: partial overlap, disjoint, identity, empty-bit, disconnected and invalid-bond fixtures.',
 paste('PASS: protected hashes unchanged;',length(rscripts),'workflow scripts parse as R; no non-R analysis scripts.'),
 'LIMIT: descriptor implementation and finite vocabulary are documented choices; 2D similarity is not identity, affinity or allostery.'),
 'outputs/24b_atom_pair_validation.txt')
cat(readLines('outputs/24b_atom_pair_validation.txt'),sep='\n')
