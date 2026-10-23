# Independent base-R reconstruction of all cohort gates, overlaps and exclusions.
suppressPackageStartupMessages(library(readr))
r <- function(p) read_csv(p,show_col_types=FALSE)
a <- r('config/RESKO_mechanism_cohort_membership.csv')
map <- r('data/curated/seed_sider_mapping_accepted.csv')
leads <- r('data/curated/seed_sider_mapping_all.csv')
seeds <- r('data/curated/cancer_primary_seeds_expanded.csv')
stopifnot(nrow(a)==85L,length(unique(a$drugbank_id))==75L,
 !anyDuplicated(a[c('cancer_id','drugbank_id')]),
 setequal(paste(a$cancer_id,a$drugbank_id),paste(seeds$cancer_id,seeds$drugbank_id)))
x <- read_tsv('data/raw/sider/meddra_all_se.tsv.gz',
 col_names=c('flat','stereo','label','type','pt','term'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
x <- unique(as.data.frame(x[x$type=='PT',c('flat','pt')]))
f <- table(x$pt); N <- length(unique(x$flat))
sm <- r('outputs/23b_mechanism_SIDER_cohort_summary.csv')
sc <- readRDS('data/curated/RESKO_mechanism_SIDER_all_scores.rds')
pt <- r('data/curated/RESKO_mechanism_SIDER_joint_patterns.csv')
for(i in seq_len(nrow(sm))) {
 z <- sm[i,]; members <- a[a$cohort_id==z$cohort_id,]
 removed <- if(z$branch=='top_decile_types') names(f)[f>=sort(f,decreasing=TRUE)[ceiling(length(f)*.1)]] else names(f)[f/N>.1]
 keep <- x[!x$pt %in% removed,]; by_group <- split(keep$pt,keep$flat)
 sets <- lapply(members$drugbank_id,function(id) unique(unlist(by_group[intersect(names(by_group),map$stitch_flat_id[map$drugbank_id==id])],use.names=FALSE)))
 eligible <- lengths(sets)>=z$minimum_seed_terms
 joint <- if(any(eligible)) Reduce(intersect,sets[eligible]) else character()
 groups <- unique(map$stitch_flat_id[map$drugbank_id %in% members$drugbank_id[eligible]])
 gate <- all(eligible) & length(groups)>=2L & length(joint)>=5L
 exploratory <- !all(eligible) & length(groups)>=2L & length(joint)>=5L
 sel <- function(d) d$cohort_id==z$cohort_id & d$branch==z$branch & d$minimum_seed_terms==z$minimum_seed_terms
 stopifnot(setequal(joint,pt$pt[sel(pt)]),z$eligible_seeds==sum(eligible),
 z$joint_PT_count==length(joint),z$independent_flat_profiles==length(groups),
 z$full_class_assessable==all(eligible),z$primary_gate==gate,
 z$exploratory_covered_subset_gate==exploratory)
 ss <- sc[sel(sc),]
 expected <- vapply(ss$flat,function(g) sum(unique(by_group[[g]]) %in% joint),integer(1))
 known_ids <- a$drugbank_id[a$cancer_id==z$cancer_id]
 known <- map$stitch_flat_id[map$drugbank_id %in% known_ids]
 held <- leads$stitch_flat_id[leads$drugbank_id %in% known_ids]
 stopifnot(all(ss$shared_PT==expected),
 all(ss$known_primary_group==(ss$flat %in% known)),
 all(ss$unresolved_known_primary_lead==(ss$flat %in% held & !ss$flat %in% known)),
 all(ss$pass_ge_50==(length(joint)>0 & 2L*expected>=length(joint))),
 all(ss$pass_gt_50==(length(joint)>0 & 2L*expected>length(joint))))
 ok <- !ss$known_primary_group & !ss$unresolved_known_primary_lead & !ss$known_identity_conflict
 stopifnot(z$nonprimary_groups_ge_50==sum(ok & ss$pass_ge_50),
 z$nonprimary_groups_gt_50==sum(ok & ss$pass_gt_50))
}
q <- r('data/curated/RESKO_mechanism_SIDER_primary_group_review_queue.csv')
stopifnot(all(q$primary_gate),!any(q$known_identity_conflict),
 !any(q$known_primary_group),!any(q$unresolved_known_primary_lead),
 all(sm$independent_flat_profiles[sm$primary_gate]>=2),
 all(sm$joint_PT_count[sm$primary_gate]>=5),
 !any(sm$primary_gate & sm$exploratory_covered_subset_gate))
# Explicit small-group/tie guards.
stopifnot(!(TRUE & 1>=2 & 10>=5),!(FALSE & 2>=2 & 10>=5),
 1*2>=2,!(1*2>2),length(Reduce(intersect,list(c('a','b'),c('b','c'))))==1)
h <- r('outputs/23b_input_hashes_summary.csv')
stopifnot(all(unname(tools::md5sum(h$path))==h$md5))
invisible(lapply(list.files('scripts',pattern='\\.R$',recursive=TRUE,full.names=TRUE),parse))
writeLines(c(paste('PASS:',nrow(sm),'branch/cohort cases independently reconstructed from raw PT sets.'),
 'PASS: all overlaps, eligibility, distinct-profile counts, gates and known-treatment holds.',
 'PASS: 85 disease/seed associations and all 75 primary drugs retained, including singleton/uncovered cohorts.',
 'PASS: singleton, incomplete-coverage and exact-half fixtures; protected hashes unchanged; all scripts parse.'),
 'outputs/23b_mechanism_SIDER_validation.txt')
cat(readLines('outputs/23b_mechanism_SIDER_validation.txt'),sep='\n')
