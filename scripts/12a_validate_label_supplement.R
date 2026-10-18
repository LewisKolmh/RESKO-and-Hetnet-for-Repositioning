# Final cross-file checks for source-aware current-label supplementation and overlap.
library(readr);library(dplyr);library(stringr);library(tibble)
if(!file.exists('data/curated/current_label_PT_supplement.csv'))stop('Run from the project root.')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
x<-r('data/curated/seed_identity_crosscheck_reviewed.csv')
s<-r('data/curated/current_label_selection_supplemented.csv')
f<-r('data/curated/current_label_PT_supplement.csv')
c<-r('data/curated/current_label_reviewed_positive_clauses.csv')
links<-r('data/curated/current_label_clause_PT_review.csv')
p<-r('data/curated/current_label_SIDER_similarity_pilot.csv')
base<-r('data/curated/seed_sider_PT_binary_long.csv')
stopifnot(nrow(x)==170,!anyDuplicated(x$drugbank_id),all(s$drugbank_id %in% x$drugbank_id),
 all(f$drugbank_id %in% c$drugbank_id),all(f$present=='1'),nrow(base)==5946,
 !anyDuplicated(base[c('drugbank_id','meddra_umls_cui')]),
 nrow(r('data/curated/cancer_drug_review_all.csv'))==238)
stopifnot(!any(f$pt_name %in% c('Breast cancer','Ovarian cancer','Small cell lung cancer','Breast cancer stage IV')),
 all(f$comparison_stratum[str_detect(f$regimen_scope,'combination|mixed')]=='context_only_regimen_or_formulation'))
co<-f |> filter(drugbank_id %in% c('DB12282','DB15254'))
stopifnot(n_distinct(co$evidence_group)==1L,
 setequal(co$pt_cui[co$drugbank_id=='DB12282'],co$pt_cui[co$drugbank_id=='DB15254']))
hot<-links |> filter(alias_norm %in% c('hot flash','hot flashes','hot flush','hot flushes'),
 mapping_review=='accepted_positive_clause_PT_mapping')
stopifnot(nrow(hot)>0L,all(hot$pt_cui=='C0600142'))
# Every pilot seed must have an explicit reviewed single-agent source, not a regimen-only profile.
mono<-f |> filter(comparison_stratum=='explicit_monotherapy_or_single_agent') |> distinct(drugbank_id,pt_cui)
stopifnot(all(p$drugbank_id %in% mono$drugbank_id),all(as.numeric(p$jaccard)>=0 & as.numeric(p$jaccard)<=1))
manifest<-r('data/raw/sider/acquisition_manifest.csv')
stopifnot(all(unname(tools::md5sum(manifest$local_file))==manifest$md5))
for(path in list.files('scripts',pattern='[.]R$',recursive=TRUE,full.names=TRUE))invisible(parse(path))
# Record supplementary sources without replacing the original nine-source RESKO manifest.
sources<-tribble(~source_id,~source_name,~scope,~canonical_url,~version_provenance,
 'pubchem','PubChem','Partial compound identities/synonyms','https://pubchem.ncbi.nlm.nih.gov/','data/raw/identity_api/; outputs/10b_api_provenance.csv',
 'chebi','ChEBI','Cisplatin/carboplatin identity support','https://www.ebi.ac.uk/chebi/','data/curated/platinum_chebi_identity_records.csv',
 'openfda_labels','openFDA drug labels','Bounded seed label queries; up to 50 results/seed','https://open.fda.gov/apis/drug/label/','outputs/11d_current_label_query_provenance.csv; per-label dates/versions',
 'dailymed_labels','DailyMed SPL','Missing-label and explicit active-moiety/route support','https://dailymed.nlm.nih.gov/dailymed/','outputs/11f_DailyMed_query_provenance.csv; cached SPL versions',
 'ema_pepaxti','EMA Pepaxti medicine overview','Named European combination-regimen context','https://www.ema.europa.eu/en/medicines/human/EPAR/pepaxti','data/raw/current_labels/EMA_pepaxti_overview.rds; publication date not inferred from retrieval')
write_csv(sources,'config/supplementary_source_manifest.csv')
report<-c('SCRIPT 12a - LABEL SUPPLEMENT VALIDATION',
 paste('Seeds with selected current US labels:',n_distinct(s$drugbank_id)),
 paste('Selected US label records:',nrow(s)),paste('Reviewed supplementary seed/PT pairs:',nrow(distinct(f,drugbank_id,pt_cui))),
 paste('Reviewed supplementary primary seed profiles:',n_distinct(f$drugbank_id)),
 paste('Pilot seed profiles:',n_distinct(p$drugbank_id)),
 'Passed: label/seed coverage, phenotype-header exclusion, hot-flush CUI resolution,',
 'combination stratum separation, shared co-pack evidence, unique original SIDER features,',
 'SIDER raw-file hashes, original 238 disease-review rows, pilot strata and all R-script parsing.',
 'Supplementary source manifest: config/supplementary_source_manifest.csv.',
 'Original nine-source manifest remains distinct. Source archives and rich curated evidence stay local.',
 'These checks do not establish exhaustive clinical curation, drug-specific causality,',
 'eEF1A binding/inhibition, allostery, HIV-1 RT PPI disruption or a calibrated statistical null.')
writeLines(report,'outputs/12a_label_supplement_validation.txt')
cat(paste0('\nSCRIPT 12a RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/12a_validate_label_supplement.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
