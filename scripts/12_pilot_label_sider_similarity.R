# Descriptive feasibility pilot: recent explicit single-agent profiles versus SIDER flat-group profiles.
# Scores are term-set overlaps, not drug-target evidence, p-values or corrected significance.
library(readr);library(dplyr);library(tidyr);library(stringr)
if(!file.exists('data/curated/current_label_PT_supplement.csv'))stop('Run 11g from the project root first.')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
f<-r('data/curated/current_label_PT_supplement.csv')
seed<-f |> filter(comparison_stratum=='explicit_monotherapy_or_single_agent') |>
 distinct(drugbank_id,drug_name,pt_cui)
seed_sizes<-seed |> count(drugbank_id,drug_name,name='seed_PT_terms') |> filter(seed_PT_terms>=5L)
seed<-seed |> semi_join(seed_sizes,by='drugbank_id')
se<-read_tsv('data/raw/sider/meddra_all_se.tsv.gz',col_names=c('stitch_flat_id','stereo','label_cui','type','pt_cui','term'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
ref<-se |> filter(type=='PT') |> distinct(stitch_flat_id,pt_cui)
ref_sizes<-ref |> count(stitch_flat_id,name='SIDER_PT_terms')
names<-read_tsv('data/raw/sider/drug_names.tsv',col_names=c('stitch_flat_id','SIDER_drug_name'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
# Flat-group profiles are the SIDER comparison unit, not verified exact molecular entities.
shared<-seed |> inner_join(ref,by='pt_cui',relationship='many-to-many') |>
 group_by(drugbank_id,stitch_flat_id) |> summarise(shared_PT_terms=n_distinct(pt_cui),.groups='drop')
scores<-crossing(drugbank_id=seed_sizes$drugbank_id,stitch_flat_id=ref_sizes$stitch_flat_id) |>
 left_join(seed_sizes,by='drugbank_id') |> left_join(ref_sizes,by='stitch_flat_id') |>
 left_join(shared,by=c('drugbank_id','stitch_flat_id')) |>
 mutate(shared_PT_terms=replace_na(shared_PT_terms,0L),
 jaccard=shared_PT_terms/(seed_PT_terms+SIDER_PT_terms-shared_PT_terms),
 seed_term_containment=shared_PT_terms/seed_PT_terms) |> left_join(names,by='stitch_flat_id')
# Exclude known self-label comparisons; unverified alias conflicts do not define self identity.
self<-r('data/curated/seed_sider_mapping_all.csv') |>
 filter(mapping_status=='accepted_full_InChIKey_match'|str_to_lower(drug_name)==str_to_lower(sider_drug_name)) |>
 distinct(drugbank_id,stitch_flat_id)
scores<-scores |> anti_join(self,by=c('drugbank_id','stitch_flat_id')) |>
 group_by(drugbank_id) |> arrange(desc(jaccard),desc(shared_PT_terms),stitch_flat_id,.by_group=TRUE) |>
 mutate(descriptive_rank=row_number()) |> ungroup() |>
 mutate(analysis_scope='descriptive_cross_source_single_agent_pilot',
 interpretation='historical_side_effect_term_overlap_only; no_inhibition_or_significance_claim')
top<-scores |> filter(descriptive_rank<=5L)
write_csv(scores,'data/curated/current_label_SIDER_similarity_pilot.csv',na='')
write_csv(top,'outputs/12_label_sider_top_matches_summary.csv',na='')
write_csv(seed_sizes,'outputs/12_label_sider_profile_coverage.csv')
# Independent set operations check a nonzero overlap calculation and prevent duplicate-row inflation.
if(nrow(top)){
 t<-top[1,];A<-unique(seed$pt_cui[seed$drugbank_id==t$drugbank_id]);B<-unique(ref$pt_cui[ref$stitch_flat_id==t$stitch_flat_id])
 stopifnot(t$shared_PT_terms==length(intersect(A,B)),
 abs(t$jaccard-length(intersect(A,B))/length(union(A,B)))<1e-12)
}
stopifnot(all(scores$jaccard>=0&scores$jaccard<=1),all(scores$seed_term_containment>=0&scores$seed_term_containment<=1),
 !anyDuplicated(scores[c('drugbank_id','stitch_flat_id')]))
report<-c('SCRIPT 12 - DESCRIPTIVE CURRENT-LABEL/SIDER OVERLAP PILOT',
 paste('Explicit monotherapy/single-agent seed profiles with >=5 PTs:',nrow(seed_sizes)),
 paste('SIDER flat-group comparison profiles:',nrow(ref_sizes)),paste('Non-self pairs scored:',nrow(scores)),
 paste('Top-five rows retained:',nrow(top)),
 'Inputs: reviewed supplementary PT features and original SIDER 4.1 PT associations.',
 'Outputs: full scores in data/curated/current_label_SIDER_similarity_pilot.csv;',
 'top-five descriptive summary and input-profile coverage in outputs/12_*.',
 'Jaccard measures shared PTs divided by union PTs. Containment measures shared PTs divided by seed PTs.',
 'This demonstrates vocabulary compatibility for newer seeds; it does not establish candidate efficacy.',
 'Only explicit single-agent/monotherapy lists are included; regimen/mixed-cohort profiles remain separate.',
 'The comparison mixes a bounded common-reaction subset with older fuller labels; profile size,',
 'common toxicities, route/formulation and source-era differences can dominate apparent similarity.',
 'Reference flat groups are not automatically verified chemical identities. Validate any expanded candidate before target joins.',
 'No p-values or BH adjustment are produced: the pilot has no calibrated null model.',
 'Next: review overlap contributions, construct consistently processed current-label comparison profiles,',
 'and integrate target/pathway evidence while retaining the original RESKO baseline.')
writeLines(report,'outputs/12_label_sider_pilot_summary.txt')
cat(paste0('\nSCRIPT 12 RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/12_pilot_label_sider_similarity.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
