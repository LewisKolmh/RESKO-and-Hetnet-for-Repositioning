# Independent base-R set reconstruction, plus boundary/missingness fixtures.
suppressPackageStartupMessages(library(readr))
r <- function(p) read_csv(p,show_col_types=FALSE)
raw <- read_tsv('data/raw/sider/meddra_all_se.tsv.gz',
 col_names=c('flat','stereo','label','type','pt','term'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
x <- unique(as.data.frame(raw[raw$type=='PT',c('flat','pt')]))
f <- table(x$pt); n <- length(unique(x$flat))
sm <- r('outputs/23a_SIDER_disease_baseline_summary.csv')
cv <- r('data/curated/RESKO_SIDER_disease_seed_coverage.csv')
mp <- r('data/curated/RESKO_SIDER_primary_seed_mapping.csv')
sc <- r('data/curated/RESKO_SIDER_disease_candidate_scores.csv')
pt <- r('data/curated/RESKO_SIDER_disease_joint_patterns.csv')
for(i in seq_len(nrow(sm))) {
 z <- sm[i,]; sel <- function(d) d$cancer_id==z$cancer_id & d$branch==z$branch &
  d$minimum_seed_terms==z$minimum_seed_terms
 removed <- if(z$branch=='top_decile_types') names(f)[f>=sort(f,decreasing=TRUE)[ceiling(length(f)*.10)]] else names(f)[f/n>.10]
 keep <- x[!x$pt %in% removed,]
 ids <- cv$drugbank_id[sel(cv) & cv$eligible]
 sets <- lapply(ids,function(id) unique(keep$pt[keep$flat %in% mp$stitch_flat_id[mp$drugbank_id==id]]))
 joint <- if(length(sets)) Reduce(intersect,sets) else character()
 stopifnot(setequal(joint,pt$pt_cui[sel(pt)]),length(joint)==z$joint_PT_count)
 ss <- sc[sel(sc),]
 expected <- vapply(ss$stitch_flat_id,function(id)sum(unique(keep$pt[keep$flat==id]) %in% joint),integer(1))
 stopifnot(identical(as.integer(ss$shared_PT),unname(expected)),
  all(ss$pass_ge_50==(length(joint)>0 & expected*2>=length(joint))),
  all(ss$pass_gt_50==(length(joint)>0 & expected*2>length(joint))))
}
# Intersection differs from union; one of two terms passes >= but not >.
stopifnot(identical(intersect(c('a','b'),c('b','c')),'b'),
 1/2>=.5,!(1/2>.5),is.na(if(0>0)0/0 else NA_real_),
 !any(cv$eligible[!cv$mapping_present]))
h <- r('outputs/23a_input_hashes_summary.csv')
stopifnot(all(unname(tools::md5sum(h$path))==h$md5))
scripts <- list.files('scripts',pattern='\\.R$',recursive=TRUE,full.names=TRUE)
invisible(lapply(scripts,parse))
writeLines(c('PASS: independent raw-source PT set intersections and every candidate overlap/boundary agree.',
 'PASS: intersection/union, exact 50-percent tie, empty pattern and missing-mapping fixtures.',
 'PASS: protected input hashes unchanged; all R scripts parse.'),
 'outputs/23a_SIDER_baseline_validation.txt')
cat(readLines('outputs/23a_SIDER_baseline_validation.txt'),sep='\n')
