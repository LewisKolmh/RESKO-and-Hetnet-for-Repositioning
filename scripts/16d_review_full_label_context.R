# Full-source context review for every retained candidate profile and sacrosidase.
library(readr);library(dplyr);library(stringr);library(tibble);library(xml2)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
paths<-c('data/curated/reviewed_profile_compatibility_registry.csv','data/curated/candidate_label_identity_review.csv',
 'data/curated/candidate_label_named_context_selection.csv')
before<-tools::md5sum(paths);p<-r(paths[1]) |> filter(cohort=='candidate');labels<-r(paths[2]);bio<-r(paths[3])
norm<-function(z)str_to_lower(str_squish(str_replace_all(z,'[\u2010-\u2015]','-')))
rules<-tribble(~entity_id,~regimen_review,~cohort_review,~threshold_operator,~threshold_percent,~comparator_rule,~frequency_limitation,~marker1,~marker2,
 'sider:CID100006116','mixed_formulation_trials_monotherapy_not_established','ESRD_hemodialysis_solid_and_liquid_trials','>',10,
 'open_label_and_placebo_crossover','Summary >10% heading does not describe all listed pooled table rates; nausea 3.6%, vomiting 2.4%, hypercalcemia 12.6%',
 '98 enrolled ESRD hemodialysis patients','Nausea6 (3.6)',
 'sider:CID100040632','single_study_drug_vs_placebo_background_not_fully_specified','IPF_pooled_studies_1_2_3_623_treated_624_placebo','>=',10,
 'more_common_than_placebo','Placebo-filtered list; cannot assume same selection rule as uncontrolled cancer summaries',
 '623 patients received 2403 mg/day of pirfenidone and 624 patients received placebo','more frequent in the pirfenidone than placebo',
 'sider:CID100148121','single_investigational_agent_single_arm_supportive_care_not_excluded','PDX_008_relapsed_refractory_PTCL','>',35,
 'uncontrolled_single_arm','Very high reporting cutoff yields only four reviewed terms; absence of other terms is not absence of effects',
 'Study PDX-008, an open-label, single-arm','Most common adverse reactions (>35%)',
 'sider:CID109912092','concomitant_DMARDs_and_other_medications_allowed','psoriatic_arthritis_PsA_1_2_3','>=',5,
 'placebo_controlled_with_background_treatment','Summary cutoff differs from table cutoff and early/late treatment periods; monotherapy cannot be inferred',
 'Subjects were allowed to receive stable doses of concomitant methotrexate','at least one DMARD (65.0%)',
 'sider:CID110482134','active_control_crossover_with_dietary_background','adult_UCD_short_term_study_1','>=',10,
 'sodium_phenylbutyrate_active_control','Retain adult short-term list; do not pool long-term or pediatric lists or interpret comparator events as test-drug effects',
 'crossover, 4-week study (Study 1)','must be used with dietary protein restriction',
 'sider:CID116126651','explicitly_pooled_monotherapy_and_combination','pooled_B_cell_malignancies','>=',30,
 'mixed_studies_and_regimens','Pooled list is not a single-agent profile; different monotherapy tables are separate evidence',
 'as a single agent at 420 mg orally once daily (475 patients)','in combination with other drugs at 420 mg orally once daily (827 patients)') |>
 mutate(strict_monotherapy_established=FALSE,frequency_cohort_alignment_established=FALSE,
 review_method='assistant_reviewed_full_cached_SPL_AR_and_clinical_context')
write_csv(rules,'config/full_candidate_context_review.csv',na='')
results<-list();evidence<-list()
for(i in seq_len(nrow(rules))){
 a<-rules[i,];b<-p |> filter(entity_id==a$entity_id);stopifnot(nrow(b)==1L)
 l<-labels |> filter(stitch_flat_id==b$stitch_flat_id,label_id==b$label_id);stopifnot(nrow(l)==1L)
 z<-readRDS(l$raw_cache_file);doc<-read_xml(z$body);xml_ns_strip(doc)
 text<-norm(xml_text(doc));markers<-c(a$marker1,a$marker2)
 stopifnot(xml_attr(xml_find_first(doc,'/document/id'),'root')==b$label_id,
 xml_attr(xml_find_first(doc,'/document/setId'),'root')==b$set_id)
 for(m in markers){loc<-str_locate(text,fixed(norm(m)));stopifnot(!is.na(loc[1]))
  evidence[[length(evidence)+1L]]<-tibble(entity_id=a$entity_id,label_id=b$label_id,marker=m,
   source_excerpt=str_sub(text,max(1,loc[1]-120),min(nchar(text),loc[2]+220)),
   raw_cache_file=l$raw_cache_file,raw_cache_md5=unname(tools::md5sum(l$raw_cache_file)))
 }
 results[[i]]<-bind_cols(b |> select(profile_id,entity_id,entity_name,current_record_name,label_id,set_id,label_version,
  effective_time,routes,reviewed_PT_terms),a |> select(-entity_id,-marker1,-marker2)) |>
 mutate(source_url=l$url,retrieved_utc=l$retrieved_utc)
}
# Named biologic remains a separate reference: counts and postmarketing reports, not a common-frequency profile.
stopifnot(nrow(bio)==1L)
z<-readRDS(bio$cache_file);rec<-z$body$results[[1]];text<-norm(paste(unlist(rec$adverse_reactions),collapse=' '))
stopifnot(str_detect(text,fixed('a total of 52 patients')),str_detect(text,fixed('population of uncertain size')),
 str_detect(text,fixed('abdominal pain (4)')),str_detect(text,fixed('dehydration (1)')))
bio_review<-bio |> select(entity_id=stitch_flat_id,label_id,set_id,label_version,effective_time,label_url) |>
 mutate(entity_id=paste0('sider:',entity_id),entity_name='sacrosidase',cohort_review='52_treated_patients_and_postmarketing_reports',
 frequency_interpretation='patient_counts_not_common_frequency_profile',regimen_review='enzyme_replacement_named_context_only',
 comparison_ready=FALSE,small_molecule_screen_eligible=FALSE,
 decision='Retain biologic reference separately; do not map the conflicting CID or pool clinical counts with postmarketing reports')
context<-bind_rows(results)
stopifnot(nrow(context)==6L,!anyDuplicated(context$profile_id),identical(before,tools::md5sum(paths)))
write_csv(context,'data/curated/reviewed_candidate_full_context.csv',na='')
write_csv(bind_rows(evidence),'data/curated/candidate_full_context_source_evidence.csv',na='')
write_csv(bio_review,'data/curated/sacrosidase_context_review.csv',na='')
write_csv(context,'outputs/16d_candidate_context_summary.csv',na='')
write_csv(bio_review,'outputs/16d_biologic_context_summary.csv',na='')
write_csv(tibble(file=paths,md5=unname(before)),'outputs/16d_input_summary.csv')
report<-c('SCRIPT 16d - FULL CANDIDATE SOURCE CONTEXT REVIEW',
 'All six retained candidate profiles and the sacrosidase reference reviewed against cached source labels.',
 'Pirfenidone: placebo-filtered >=10% list from pooled IPF trials; single study drug versus placebo, without proof of absent background treatment.',
 'Pralatrexate: single-arm PDX-008 context; >35% summary cutoff explains the four-term sparse list.',
 'Apremilast: concomitant DMARDs allowed; 65% used at least one DMARD. Its PsA summary is not an isolated monotherapy profile.',
 'Ibrutinib: the selected B-cell list pools monotherapy and combination exposures; do not label it single-agent.',
 'Glycerol phenylbutyrate: adult short-term active-control crossover with dietary background; retain other ages and durations separately.',
 'Calcium acetate: pooled solid/liquid ESRD studies; summary >10% heading is inconsistent with nausea/vomiting pooled table rates below 10%.',
 'Sacrosidase: clinical patient counts plus postmarketing context; not a common-frequency small-molecule profile.',
 'Source markers, document IDs/versions, cache hashes and local evidence excerpts are retained.',
 'No strict candidate monotherapy claim is upgraded; thresholds/comparator rules are documented but cross-drug cohort alignment is not established.',
 'Outputs: reviewed_candidate_full_context.csv, candidate_full_context_source_evidence.csv, sacrosidase_context_review.csv (local); outputs/16d_* summaries.',
 'Next: re-assess all pair barriers with these context findings, then close section 16 with consolidated features, review coverage and an explicit limitations ledger.')
writeLines(report,'outputs/16d_full_context_summary.txt')
writeLines('PASS: six distinct candidate profiles, pinned source IDs/set IDs, all reviewed source markers and sacrosidase count/context markers; input hashes unchanged.',
 'outputs/16d_full_context_validation.txt')
cat(paste0('\nSCRIPT 16d RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/16d_review_full_label_context.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
