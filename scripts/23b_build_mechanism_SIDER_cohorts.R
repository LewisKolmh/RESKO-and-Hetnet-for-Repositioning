# Predeclared mechanism-cohort adaptation; source-supported classes, not patient subtypes.
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr)})
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
inputs <- c('data/curated/cancer_primary_seeds_expanded.csv','data/curated/drugbank_drugs.csv',
 'data/curated/seed_sider_mapping_accepted.csv','data/raw/sider/meddra_all_se.tsv.gz',
 'config/cancers.csv','data/curated/seed_sider_mapping_all.csv',
 'data/raw/sider/drug_names.tsv')
stopifnot(all(file.exists(inputs))); hashes <- tools::md5sum(inputs)
# All 75 primary seeds are assigned, including singleton classes and unavailable drugs.
# Classes are specified independently of SIDER coverage/overlap; never choose best subsets.
classes <- list(
 CDK4_6_inhibitor=c('Abemaciclib','Palbociclib','Ribociclib'),
 PI3K_inhibitor=c('Alpelisib','Inavolisib'), AKT_inhibitor='Capivasertib',
 RAF_MEK_inhibitor='Avutometinib', FAK_inhibitor='Defactinib',
 aromatase_inhibitor=c('Anastrozole','Exemestane','Letrozole'),
 estrogen_receptor_modulator=c('Tamoxifen','Toremifene'),
 estrogen_receptor_degrader=c('Elacestrant','Fulvestrant','Imlunestrant'),
 estrogen_receptor_PROTAC='Vepdegestrant', GnRH_agonist='Goserelin',
 progestin='Megestrol acetate', glucocorticoid_agonist='Dexamethasone',
 glucocorticoid_receptor_antagonist='Relacorilant',
 platinum_DNA_crosslinker=c('Carboplatin','Cisplatin'),
 taxane_microtubule_stabiliser=c('Docetaxel','Paclitaxel'),
 epothilone_microtubule_stabiliser='Ixabepilone',
 vinca_microtubule_inhibitor=c('Vinblastine','Vinorelbine'),
 halichondrin_microtubule_inhibitor='Eribulin',
 anthracycline=c('Doxorubicin','Epirubicin'),
 fluoropyrimidine=c('Capecitabine','Fluorouracil'),
 nucleoside_antimetabolite='Gemcitabine', antifolate='Methotrexate',
 nitrogen_mustard=c('Cyclophosphamide','Melphalan','Melphalan flufenamide'),
 nitrosourea='Carmustine', topoisomerase_I_inhibitor='Topotecan',
 DNA_minor_groove_binder='Trabectedin',
 PARP_inhibitor=c('Niraparib','Olaparib','Rucaparib','Talazoparib'),
 proteasome_inhibitor=c('Bortezomib','Carfilzomib','Ixazomib'),
 cereblon_IMiD=c('Lenalidomide','Pomalidomide'),
 mTOR_inhibitor='Everolimus', XPO1_inhibitor='Selinexor',
 HER2_kinase_inhibitor=c('Lapatinib','Neratinib','Tucatinib'),
 HER2_antibody=c('Trastuzumab','Pertuzumab','Margetuximab'),
 HER2_ADC=c('Trastuzumab deruxtecan','Trastuzumab emtansine'),
 TROP2_ADC=c('Datopotamab deruxtecan','Sacituzumab govitecan'),
 FRalpha_ADC='Mirvetuximab soravtansine', BCMA_ADC='Belantamab mafodotin',
 VEGF_antibody='Bevacizumab', PD1_antibody='Pembrolizumab',
 CD38_antibody=c('Daratumumab','Isatuximab'), SLAMF7_antibody='Elotuzumab',
 BCMA_CAR_T=c('Ciltacabtagene autoleucel','Idecabtagene vicleucel'),
 BCMA_CD3_bispecific=c('Elranatamab','Linvoseltamab','Teclistamab'),
 GPRC5D_CD3_bispecific='Talquetamab')
classification <- bind_rows(lapply(names(classes),function(k)
 tibble(drug_name=classes[[k]],mechanism_class=k)))
seeds <- r(inputs[1]) |> semi_join(r(inputs[5]) |> filter(include=='TRUE'),by='cancer_id') |>
 distinct(cancer_id,cancer_name,drugbank_id,drug_name)
stopifnot(!anyDuplicated(classification$drug_name),
 setequal(classification$drug_name,seeds$drug_name))
db <- r(inputs[2]) |> select(drugbank_id,mechanism_of_action,source_export_date)
assignment <- seeds |> left_join(classification,by='drug_name') |>
 left_join(db,by='drugbank_id') |>
 mutate(cohort_id=paste(cancer_id,mechanism_class,sep='__'),
  classification_basis='analyst mechanism grouping; frozen DrugBank mechanism text retained locally',
  mechanism_source_url=paste0('https://go.drugbank.com/drugs/',drugbank_id))
stopifnot(!anyNA(assignment$mechanism_class),
 all(!is.na(assignment$mechanism_of_action) & nzchar(assignment$mechanism_of_action)))
freeze <- function(x,p) {
 if(file.exists(p)) stopifnot(identical(as.data.frame(r(p)),
 as.data.frame(mutate(x,across(everything(),as.character))))) else write_csv(x,p)
}
freeze(assignment |> select(cancer_id,cohort_id,drugbank_id,drug_name,mechanism_class,
 source_export_date,mechanism_source_url,classification_basis),
 'config/RESKO_mechanism_cohort_membership.csv')
policy <- tibble(parameter=c('adaptation','cohort_formation','minimum_independent_profiles',
 'minimum_joint_PT','missingness','pruning','seed_minimum','candidate_boundaries',
 'known_treatments','identity','mechanism_limits'),
 value=c('mechanism-stratified primary cancer seeds; not original paper fallback',
 'all declared classes crossed with each seed disease; no subset optimisation',
 '2 distinct accepted SIDER flat groups; singleton coverage never a multi-drug pattern',
 '5 retained joint PT for qualifying search',
 'full intended class required for primary gate; incomplete covered subsets exploratory only',
 'top-decile types with all cutoff ties; prevalence >10%; same global SIDER universe',
 '>=3 and >3 post-pruning PT','>=50% and >50%, both reported',
 'exclude accepted primary groups and hold unresolved primary mapping leads; no identity pooling',
 'historical flat groups require fresh identity review; known sacrosidase conflict excluded',
 'broad pharmacological classes do not imply interchangeable indications, patient subtypes or identical mechanisms; ADC payloads may differ'))
freeze(policy,'config/RESKO_mechanism_cohort_policy.csv')
message('Membership and policy frozen before SIDER intersection calculation.')
write_csv(assignment,'data/curated/RESKO_mechanism_cohort_source_evidence.csv')
raw <- read_tsv(inputs[4],col_names=c('flat','stereo','label','type','pt','term'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
profiles <- raw |> filter(type=='PT',!is.na(pt),nzchar(pt)) |> distinct(flat,pt)
freq <- profiles |> count(pt,name='n')
cut <- sort(freq$n,decreasing=TRUE)[ceiling(nrow(freq)*.1)]
map <- r(inputs[3]) |> distinct(drugbank_id,stitch_flat_id)
all_leads <- r(inputs[6]) |> distinct(drugbank_id,stitch_flat_id)
names <- read_tsv(inputs[7],col_names=c('flat','historical_name'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
stopifnot(!anyDuplicated(names$flat))
coverage <- summaries <- patterns <- scores <- list()
for(branch in c('top_decile_types','prevalence_gt_0.10')) {
 removed <- if(branch=='top_decile_types') freq$pt[freq$n>=cut] else
 freq$pt[freq$n/n_distinct(profiles$flat)>.1]
 kept <- profiles |> filter(!pt %in% removed)
 sp <- assignment |> select(cohort_id,drugbank_id) |>
 inner_join(map,by='drugbank_id',relationship='many-to-many') |>
 inner_join(kept,by=c('stitch_flat_id'='flat'),relationship='many-to-many') |>
 distinct(cohort_id,drugbank_id,pt)
 counts <- sp |> count(cohort_id,drugbank_id,name='retained_PT')
 for(minimum in c(3L,4L)) {
  cv <- assignment |> select(cohort_id,cancer_id,mechanism_class,drugbank_id,drug_name) |>
   left_join(counts,by=c('cohort_id','drugbank_id')) |>
   mutate(retained_PT=replace_na(retained_PT,0L),mapping_present=drugbank_id %in% map$drugbank_id,
    eligible=mapping_present & retained_PT>=minimum,branch=branch,minimum_seed_terms=minimum)
  coverage[[length(coverage)+1L]] <- cv
  for(cohort in unique(assignment$cohort_id)) {
   cc <- cv |> filter(cohort_id==cohort); ids <- cc$drugbank_id[cc$eligible]
   groups <- unique(map$stitch_flat_id[map$drugbank_id %in% ids])
   joint <- sp |> filter(cohort_id==cohort,drugbank_id %in% ids) |>
    count(pt) |> filter(n==length(ids),length(ids)>0) |> pull(pt)
   J <- length(joint); cancer <- cc$cancer_id[1]
   full <- all(cc$eligible); multi <- length(groups)>=2L
   gate <- full & multi & J>=5L
   exploratory <- !full & multi & J>=5L
   context <- tibble(cohort_id=cohort,cancer_id=cancer,mechanism_class=cc$mechanism_class[1],
    branch=branch,minimum_seed_terms=minimum,intended_seeds=nrow(cc),eligible_seeds=length(ids),
    independent_flat_profiles=length(groups),missing_mapping=sum(!cc$mapping_present),
    joint_PT_count=J,full_class_assessable=full,primary_gate=gate,
    exploratory_covered_subset_gate=exploratory)
   patterns[[length(patterns)+1L]] <- tibble(cohort_id=cohort,branch=branch,
    minimum_seed_terms=minimum,pt=joint)
   known <- map$stitch_flat_id[map$drugbank_id %in% seeds$drugbank_id[seeds$cancer_id==cancer]]
   held_known <- all_leads$stitch_flat_id[all_leads$drugbank_id %in% seeds$drugbank_id[seeds$cancer_id==cancer]]
   hits <- kept |> filter(pt %in% joint) |> count(flat,name='shared_PT')
   sc <- tibble(flat=sort(unique(profiles$flat))) |> left_join(hits,by='flat') |>
    mutate(shared_PT=replace_na(shared_PT,0L),joint_PT_count=J,
     coverage_fraction=if(J>0)shared_PT/J else NA_real_,
     pass_ge_50=J>0 & 2L*shared_PT>=J,pass_gt_50=J>0 & 2L*shared_PT>J,
     known_primary_group=flat %in% known,
     unresolved_known_primary_lead=flat %in% held_known & !flat %in% known,known_identity_conflict=flat=='CID100135113',
     cohort_id=cohort,branch=branch,minimum_seed_terms=minimum,
     primary_gate=gate,exploratory_covered_subset_gate=exploratory)
   scores[[length(scores)+1L]] <- sc
   summaries[[length(summaries)+1L]] <- context |> mutate(
    nonprimary_groups_ge_50=sum(sc$pass_ge_50 & !sc$known_primary_group & !sc$unresolved_known_primary_lead & !sc$known_identity_conflict),
    nonprimary_groups_gt_50=sum(sc$pass_gt_50 & !sc$known_primary_group & !sc$unresolved_known_primary_lead & !sc$known_identity_conflict))
  }
 }
}
sm <- bind_rows(summaries); sc <- bind_rows(scores) |> left_join(names,by="flat")
write_csv(bind_rows(coverage),'data/curated/RESKO_mechanism_SIDER_seed_coverage.csv')
write_csv(bind_rows(patterns) |> left_join(raw |> filter(type=='PT') |> distinct(pt,term),
 by='pt',relationship='many-to-many'),'data/curated/RESKO_mechanism_SIDER_joint_patterns.csv')
saveRDS(sc,'data/curated/RESKO_mechanism_SIDER_all_scores.rds')
write_csv(sc |> filter(primary_gate,pass_ge_50,!known_primary_group,!unresolved_known_primary_lead,!known_identity_conflict),
 'data/curated/RESKO_mechanism_SIDER_primary_group_review_queue.csv')
write_csv(sc |> filter(exploratory_covered_subset_gate,pass_ge_50,!known_primary_group,!unresolved_known_primary_lead,!known_identity_conflict),
 'data/curated/RESKO_mechanism_SIDER_exploratory_group_review_queue.csv')
write_csv(sm,'outputs/23b_mechanism_SIDER_cohort_summary.csv')
write_csv(sm |> count(branch,minimum_seed_terms,primary_gate,exploratory_covered_subset_gate),
 'outputs/23b_mechanism_SIDER_gate_summary.csv')
stopifnot(identical(hashes,tools::md5sum(inputs)))
protected <- c(inputs,'config/RESKO_mechanism_cohort_membership.csv',
 'config/RESKO_mechanism_cohort_policy.csv')
write_csv(tibble(path=protected,md5=unname(tools::md5sum(protected))),
 'outputs/23b_input_hashes_summary.csv')
print(sm |> filter(primary_gate | exploratory_covered_subset_gate),n=100,width=Inf)
message('23b complete. No p-values, binding claims or original-paper fallback equivalence.')
