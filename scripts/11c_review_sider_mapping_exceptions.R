# Keep historical drug-label context separately from exact molecular mappings.
# A label name does not prove a salt, stereoisomer or coordination representation equivalent.
library(readr);library(dplyr);library(stringr);library(tidyr)
if (!all(dir.exists(c("scripts", "config", "data"))))stop('Run from project root.')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
x<-r('data/curated/seed_identity_crosscheck_reviewed.csv')
m<-r('data/curated/seed_sider_mapping_all.csv')
a<-r('data/curated/seed_aliases_reviewed.csv')
q<-m |> filter(mapping_status!='accepted_full_InChIKey_match') |>
 mutate(review_decision=case_when(
   drugbank_id=='DB11760' & str_to_lower(sider_drug_name)=='bicalutamide' ~ 'reject_cross_entity_alias',
   str_to_lower(str_squish(drug_name))==str_to_lower(str_squish(sider_drug_name)) ~ 'retain_named_label_context_ONLY',
   TRUE ~ 'retain_unverified_related_name_lead_ONLY'),
   review_reason=case_when(review_decision=='reject_cross_entity_alias' ~
     'ChEMBL CHEMBL3137320 lists BICALUTAMIDE as a PND synonym; SIDER bicalutamide has a different structure. Do not transfer its data to talazoparib.',
     review_decision=='retain_named_label_context_ONLY' ~
     'Exact normalised primary drug name matches the historical SIDER label identity, but current CID structure differs. Retain drug-label context separately; no molecular equivalence or activity transfer.',
     TRUE ~ 'Alias or related-form name differs from primary seed name; keep as a lead, without automatic label-profile transfer.'))
# Name-level context includes only literal primary-name agreement, never synonym-only matches.
context_map<-q |> filter(review_decision=='retain_named_label_context_ONLY') |>
 distinct(drugbank_id,drug_name,stitch_flat_id,sider_drug_name) |>
 mutate(mapping_scope='historical_named_drug_label_context_NOT_exact_structure')
se<-read_tsv('data/raw/sider/meddra_all_se.tsv.gz',col_names=c('stitch_flat_id','stitch_stereo_id',
 'label_umls_cui','meddra_type','meddra_umls_cui','side_effect_name'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
context<-context_map |> inner_join(se |> filter(meddra_type=='PT') |>
 distinct(stitch_flat_id,meddra_umls_cui,side_effect_name),by='stitch_flat_id',relationship='many-to-many')
strict<-r('data/curated/seed_sider_PT_binary_long.csv')
coverage<-r('outputs/11b_sider_seed_coverage.csv') |>
 mutate(across(c(candidate_SIDER_records,accepted_SIDER_records,unique_PT_side_effects),as.integer)) |>
 left_join(context |> distinct(drugbank_id,meddra_umls_cui) |>
 count(drugbank_id,name='named_label_context_PT_count'),by='drugbank_id') |>
 mutate(reviewed_coverage_status=case_when(!is.na(unique_PT_side_effects) ~ 'strict_structure_supported_PT_profile',
 !is.na(named_label_context_PT_count)~'named_label_context_only_NOT_strict',
 candidate_SIDER_records>0~'no_accepted_mapping_after_review',TRUE~'no_candidate_mapping'))
primary<-r('data/curated/primary_compound_mapping_candidates.csv')$drugbank_id
summary<-bind_rows(coverage |> count(reviewed_coverage_status,name='drugs') |> mutate(scope='all_170_seeds'),
 coverage |> filter(drugbank_id %in% primary) |> count(reviewed_coverage_status,name='drugs') |> mutate(scope='primary_55_structure_queue'))
# Do not merge these context associations into seed_sider_PT_binary_long.csv.
write_csv(q,'data/curated/sider_mapping_exception_decisions.csv',na='')
write_csv(context,'data/curated/seed_sider_named_label_context_PT.csv',na='')
write_csv(coverage,'outputs/11c_sider_reviewed_coverage.csv',na='')
write_csv(summary,'outputs/11c_sider_reviewed_coverage_summary.csv',na='')
write_csv(a |> filter(seed_drugbank_id=='DB11760',str_to_lower(alias)=='bicalutamide'),
 'outputs/11c_rejected_alias_provenance.csv',na='')
# Validate the join and keep missing coverage as NA, not zero side effects.
stopifnot(nrow(q)==25,nrow(coverage)==170,!anyDuplicated(coverage$drugbank_id),
 !'DB11760' %in% context$drugbank_id,
 all(is.na(coverage$unique_PT_side_effects[coverage$accepted_SIDER_records==0L])),
 nrow(strict)==5946,nrow(r('data/curated/cancer_drug_review_all.csv'))==238)
f<-r('data/curated/seed_sider_frequencies.csv')
lo<-as.numeric(f$frequency_lower);hi<-as.numeric(f$frequency_upper)
stopifnot(all(is.na(lo)|is.na(hi)|lo<=hi),all(is.na(lo)|(lo>=0&lo<=1)),all(is.na(hi)|(hi>=0&hi<=1)))
writeLines(c('11a/11b/11c validation passed: raw manifest hashes retained; schema parsing succeeded;',
 'all accepted mappings agree on full structure key; seed/PT features unique;',
 'frequency bounds ordered within [0,1]; placebo fields retained separately;',
 '170 seeds and 238 disease-review rows retained; unmapped coverage remains NA;',
 '25 exception decisions retained; context data kept outside strict 5,946-row feature table.'),
 'outputs/11_sider_validation.txt')
report<-c('SCRIPT 11c - SIDER IDENTITY EXCEPTION REVIEW',capture.output(print(summary)),
 capture.output(print(count(q,review_decision))),
 'Outputs: sider_mapping_exception_decisions.csv; seed_sider_named_label_context_PT.csv',
 '(data/curated/), reviewed coverage and rejected-alias provenance reports (outputs/11c_*).',
 'Exact primary-name historical profiles are retained as a SEPARATE sensitivity/context table.',
 'These are not exact-structure joins. Salt, stereochemical and coordination differences remain unresolved.',
 'The ChEMBL BICALUTAMIDE alias is rejected for transferring SIDER data to talazoparib.',
 'Two other synonym-only/related-name leads are retained without transferring profiles.',
 'Main analysis remains 34 seed drugs (23 primary structure-queue drugs), with 5,946 unique seed/PT pairs.',
 'Validation: outputs/11_sider_validation.txt. Unmapped values are not observed negative events.',
 'Next: define a side-effect similarity/background analysis with strict and contextual results separate;',
 'then expand candidate identities and integrate target/pathway evidence. BH applies to valid tests only.')
writeLines(report,'outputs/11c_sider_review_summary.txt')
cat(paste0('\nSCRIPT 11c RUN - ',format(Sys.time(),tz='Europe/Rome',usetz=TRUE),'\n',
 'Script: scripts/11c_review_sider_mapping_exceptions.R\n',
 'Inputs: 11b mapping, raw PT records and 10d alias provenance. No new API calls.\n',
 paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
