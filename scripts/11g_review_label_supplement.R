# Curate bounded, explicitly positive reaction lists for primary compounds lacking strict SIDER profiles.
# Regimen and formulation context are part of the feature identity, not discarded covariates.
library(readr);library(dplyr);library(stringr);library(tibble);library(xml2);library(httr2)
if(!file.exists('data/curated/current_label_PT_dictionary.csv'))stop('Run 11d/11f/11e from the project root first.')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
norm<-function(z)str_to_lower(str_squish(str_replace_all(z,'[\u2010-\u2015]','-')))
s<-r('data/curated/current_label_selection_supplemented.csv');d<-r('data/curated/current_label_PT_dictionary.csv')
x<-r('data/curated/seed_identity_crosscheck_reviewed.csv');coverage<-r('outputs/11c_sider_reviewed_coverage.csv')
primary<-r('data/curated/primary_compound_mapping_candidates.csv')$drugbank_id
# Start/end delimiters are reviewed against cached source paragraphs. No whole-section auto-acceptance.
rules<-tribble(~drugbank_id,~start_regex,~end_regex,~regimen_scope,~required_marker,~forbidden_marker,
 'DB00014','adverse event profile was similar for women','Tumor flare|To report','female_label_context_regimen_unspecified','','',
 'DB00361','Most common adverse reactions','\\(\\s*6\\.1\\s*\\)|To report','label_summary_regimen_unspecified','','',
 'DB00441','most common.{0,90}adverse reactions','\\(\\s*6\\.1\\s*\\)|To report','intravenous_label_context_regimen_unspecified','','INLEXZO',
 'DB00515','Common adverse reactions are','The following adverse reactions have been identified','label_summary_regimen_unspecified','','',
 'DB00539','Most common adverse reactions are','\\(\\s*6\\.1\\s*\\)|To report','label_summary_regimen_unspecified','','',
 'DB00675','most frequent adverse reaction to tamoxifen is','Other adverse reactions','label_summary_regimen_unspecified','','',
 'DB00958','Most common adverse reactions.{0,200}recurrent ovarian cancer.{0,100}as a single agent','\\(\\s*6\\.1\\s*\\)|To report','explicit_single_agent','','',
 'DB01030','Small cell lung cancer','Cervical cancer|To report','small_cell_lung_cancer_label_context','','',
 'DB01229','SUMMARY.{0,25}OF ADVERSE EVENTS IN PATIENTS WITH SOLID TUMORS RECEIVING SINGLE.AGENT PACLITAXEL','Hypersensitivity Reaction','explicit_single_agent_conventional_paclitaxel','','albumin-bound|protein-bound',
 'DB01248','most common adverse reactions across all (?:docetaxel|BEIZRAY) indications','Incidence varies|\\(\\s*6\\s*\\)|To report','mixed_indications_regimen_unspecified','','',
 'DB06374','most common.{0,80}adverse reactions.{0,100}of ORSERDU were','\\(\\s*6\\.1\\s*\\)|To report','single_agent_regulatory_cohort','ORSERDU','',
 'DB08871','most common adverse reactions.{0,50}in metastatic breast cancer were','\\(\\s*6\\.1\\s*\\)|The most common','metastatic_breast_cancer_label_context','','',
 'DB09073','Most common adverse reactions.{0,150}in combination with either letrozole or fulvestrant','\\(\\s*6\\.1\\s*\\)|To report','combination_letrozole_or_fulvestrant','','',
 'DB09074','as a single agent were','\\(\\s*6\\.1\\s*\\)|in combination','explicit_single_agent','','',
 'DB09570','most common adverse reactions','\\(\\s*6\\.1\\s*\\)|To report','combination_lenalidomide_dexamethasone','','',
 'DB11652','most common adverse reactions.{0,60}with TUKYSA in combination with trastuzumab and capecitabine','\\(\\s*6\\.1\\s*\\)|The most common','combination_trastuzumab_capecitabine','','',
 'DB11730','In patients with early breast cancer','\\(\\s*6\\s*\\)|In patients with advanced','combination_endocrine_therapy_early_breast','','',
 'DB11760','Most common adverse reactions.{0,60}as a single agent','\\(\\s*6\\.1\\s*\\)|Most common adverse reactions','explicit_single_agent','','',
 'DB11793','Most common adverse reactions.{0,80}in patients who received ZEJULA were','\\(\\s*6\\.1\\s*\\)|To report','monotherapy_confirmed_in_source_cohort','','',
 'DB11828','NERLYNX as a single agent:','\\(\\s*6\\s*\\)|NERLYNX in combination','explicit_single_agent','','',
 'DB11942','most common adverse reactions.{0,100}multiple myeloma who receive Xd are','\\(\\s*6\\.1\\s*\\)|To report','combination_selinexor_dexamethasone_Xd','','',
 'DB12001','Most common adverse reactions','\\(\\s*6\\s*\\)|To report','mixed_monotherapy_and_combination_cohorts','','',
 'DB12015','Most common adverse reactions','\\(\\s*6\\.1\\s*\\)|To report','combination_fulvestrant','PIQRAY','VIJOICE',
 'DB12218','HR-positive, HER2-negative breast cancer were','\\(\\s*6\\.1\\s*\\)|PTEN-deficient','combination_fulvestrant_breast','','',
 'DB12332','Most common adverse reactions.{0,80}among patients with ovarian cancer were','\\(\\s*6\\.1\\s*\\)|Most common adverse reactions','ovarian_cancer_label_context','','',
 'DB14976','Most common adverse reactions','\\(\\s*6\\s*\\)|To report','combination_nab_paclitaxel','','',
 'DB15275','most common.{0,80}adverse reactions','\\(\\s*6\\.1\\s*\\)|To report','combination_palbociclib_fulvestrant','','',
 'DB19006','Most common.{0,80}adverse reactions with VEPPANU','\\(\\s*6\\.1\\s*\\)|To report','single_agent_regulatory_cohort','','',
 'DB19043','most common.{0,80}adverse reactions.{0,100}INLURIYO monotherapy were','\\(\\s*6\\.1\\s*\\)|The most common','explicit_monotherapy','','')
write_csv(rules,'config/current_label_positive_clause_rules.csv')
clauses<-list();failures<-list()
extract_clause<-function(text,start,end){
 p<-str_locate(text,regex(start,ignore_case=TRUE));if(is.na(p[1]))return(NA_character_)
 tail<-str_sub(text,p[1]);e<-str_locate(str_sub(tail,nchar(str_sub(text,p[1],p[2]))+1L),regex(end,ignore_case=TRUE))
 if(is.na(e[1]))return(NA_character_)
 str_sub(tail,1L,nchar(str_sub(text,p[1],p[2]))+e[1]-1L)
}
for(i in seq_len(nrow(rules))){
 rule<-rules[i,];options<-s |> filter(drugbank_id==rule$drugbank_id);found<-FALSE
 for(j in seq_len(nrow(options))){
  m<-options[j,];z<-readRDS(m$cache_file);rec<-Filter(function(a)identical(a$id,m$label_id),z$body$results)[[1]]
  text<-str_squish(paste(unlist(rec$adverse_reactions,use.names=FALSE),collapse=' '))
  if(nzchar(rule$required_marker)&&!str_detect(text,regex(rule$required_marker,ignore_case=TRUE)))next
  if(nzchar(rule$forbidden_marker)&&str_detect(text,regex(rule$forbidden_marker,ignore_case=TRUE)))next
  clause<-extract_clause(text,rule$start_regex,rule$end_regex)
  if(is.na(clause)||nchar(clause)>6000L)next
  clauses[[length(clauses)+1L]]<-tibble(drugbank_id=m$drugbank_id,drug_name=m$drug_name,set_id=m$set_id,
   label_id=m$label_id,effective_time=m$effective_time,source='FDA_SPL_via_openFDA_or_DailyMed',source_url=m$label_url,
   label_identity_status=m$label_identity_status,regimen_scope=rule$regimen_scope,positive_clause=clause,
   evidence_group=paste0('SPL:',m$set_id),review_method='assistant_reviewed_positive_clause_with_explicit_boundaries')
  found<-TRUE;break
 }
 if(!found)failures[[length(failures)+1L]]<-tibble(drugbank_id=rule$drugbank_id,reason='No selected source matched both reviewed clause boundaries; not auto-promoted')
}
# A co-pack supports regimen context for both ingredients, never independent single-drug attribution.
setid<-'65889527-1406-45df-8a85-e2690dcf427d'
z<-readRDS(paste0('data/raw/current_labels/dailymed/spl_',setid,'.rds'));doc<-read_xml(z$body);xml_ns_strip(doc)
text<-str_squish(paste(xml_text(xml_find_all(doc,'//section[code[@code="34084-4"]]//text')),collapse=' '))
clause<-extract_clause(text,'most common.{0,80}adverse reactions','\\(\\s*6\\.1\\s*\\)|To report')
moieties<-xml_attr(xml_find_all(doc,'//ingredient[starts-with(@classCode,"ACTI")]//activeMoiety/code'),'code')
for(id in c('DB12282','DB15254')){
 seed<-x |> filter(drugbank_id==id);stopifnot(seed$unii %in% moieties,!is.na(clause))
 clauses[[length(clauses)+1L]]<-tibble(drugbank_id=id,drug_name=seed$drug_name,set_id=setid,
 label_id=xml_attr(xml_find_first(doc,'/document/id'),'root'),effective_time=xml_attr(xml_find_first(doc,'/document/effectiveTime'),'value'),
 source='DailyMed_combination_product_SPL',source_url=paste0('https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=',setid),
 label_identity_status='explicit_seed_moiety_in_two_ingredient_combination',
 regimen_scope='combination_avutometinib_defactinib_CO_PACK_SHARED_EVIDENCE',positive_clause=clause,
 evidence_group=paste0('SPL:',setid),review_method='assistant_reviewed_regimen_context_not_individual_drug_attribution')
}
# Regulatory overview supplement for the remaining primary drug without a current US-label match.
url<-'https://www.ema.europa.eu/en/medicines/human/EPAR/pepaxti'
path<-'data/raw/current_labels/EMA_pepaxti_overview.rds'
if(file.exists(path))z<-readRDS(path)else{
 resp<-request(url)|>req_timeout(60)|>req_retry(max_tries=3,retry_on_failure=TRUE)|>req_perform()
 z<-list(url=url,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),body=resp_body_string(resp));saveRDS(z,path)}
text<-str_squish(xml_text(xml_find_first(read_html(z$body),'//body')))
stopifnot(str_detect(text,regex('active substance melphalan flufenamide',ignore_case=TRUE)))
clause<-extract_clause(text,'The most common side effects with Pepaxti','The most frequent serious side effects|Pepaxti must')
stopifnot(!is.na(clause))
clauses[[length(clauses)+1L]]<-tibble(drugbank_id='DB16627',drug_name='Melphalan flufenamide',set_id=NA_character_,label_id='EMA_PEPAXTI_OVERVIEW',
 effective_time=NA_character_,source='EMA_regulatory_medicine_overview',source_url=url,
 label_identity_status='regulatory_named_active_substance_context',regimen_scope='combination_dexamethasone_EU_overview',
 positive_clause=clause,evidence_group='EMA:PEPAXTI:OVERVIEW',
 review_method='assistant_reviewed_regulatory_common_reaction_paragraph; overview may predate current detailed SmPC')
cl<-bind_rows(clauses) |> mutate(clause_id=sprintf('POSCLAUSE%03d',row_number()))
# Match only the listed reactions after the cohort/formulation prefix.
reaction_body<-function(text,id){
 if(id=='DB01229'){
  p<-str_locate(text,regex('Neutropenia',ignore_case=TRUE));return(if(is.na(p[1]))NA_character_ else str_sub(text,p[1]))
 }
 if(id=='DB11828'){
  p<-str_locate(text,regex('NERLYNX as a single agent:',ignore_case=TRUE));return(if(is.na(p[1]))NA_character_ else str_sub(text,p[2]+1L))
 }
 p<-str_locate(text,regex('\\b(?:are|were|is|included)\\b\\s*:?',ignore_case=TRUE))
 if(is.na(p[1]))NA_character_ else str_sub(text,p[2]+1L)
}
stopifnot(!str_detect(reaction_body('Most common adverse reactions in breast cancer were nausea and vomiting.','test'),
 regex('breast cancer',ignore_case=TRUE)))
cl$reaction_list_text<-mapply(reaction_body,cl$positive_clause,cl$drugbank_id,USE.NAMES=FALSE)
stopifnot(!anyNA(cl$reaction_list_text))
terms<-unique(d$alias_norm);terms<-terms[order(nchar(terms),decreasing=TRUE)]
pattern<-paste0('(?<![\\p{L}\\p{N}])(?:',paste(str_escape(terms),collapse='|'),')(?![\\p{L}\\p{N}])')
hits<-list()
for(i in seq_len(nrow(cl))){
 text<-norm(cl$reaction_list_text[i]);matches<-unique(str_extract_all(text,regex(pattern))[[1]])
 if(length(matches))hits[[i]]<-tibble(clause_id=cl$clause_id[i],alias_norm=matches)
}
# Resolve the literal hot-flush variants using the supplied MedDRA term-CUI rows,
# rather than also assigning the broader SIDER label-level 'Menopausal symptoms' mapping.
meddra<-read_tsv('data/raw/sider/meddra.tsv.gz',col_names=c('cui','kind','meddra_id','term'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
hot_aliases<-c('hot flash','hot flashes','hot flush','hot flushes')
proof<-meddra |> filter(norm(term) %in% hot_aliases)
stopifnot(nrow(proof)>0L,all(proof$cui=='C0600142'))
linked<-bind_rows(hits) |> inner_join(d,by='alias_norm',relationship='many-to-many') |>
 left_join(cl,by='clause_id') |>
 mutate(manual_mapping_evidence=if_else(alias_norm %in% hot_aliases & pt_cui=='C0600142',
   'MedDRA supplied literal term-CUI mapping resolves Hot flush; broader Menopausal symptoms alternative held',NA_character_),
 mapping_review=case_when(!is.na(manual_mapping_evidence)~'accepted_positive_clause_PT_mapping',
 ambiguous_dictionary_mapping=='TRUE'~'hold_ambiguous_dictionary_mapping',
 broad_term=='TRUE'~'hold_broad_term',TRUE~'accepted_positive_clause_PT_mapping'))
accepted<-linked |> filter(mapping_review=='accepted_positive_clause_PT_mapping')
features<-accepted |> distinct(drugbank_id,drug_name,pt_cui,pt_name,source,source_url,effective_time,
 regimen_scope,evidence_group,label_identity_status) |> mutate(present=1L,
 feature_layer='current_regulatory_label_supplement_not_original_SIDER')
# Only explicit single-agent/monotherapy clauses enter the narrow comparison subset.
features<-features |> mutate(comparison_stratum=if_else(str_detect(regimen_scope,'single_agent|monotherapy_confirmed|explicit_monotherapy'),
 'explicit_monotherapy_or_single_agent','context_only_regimen_or_formulation'))
write_csv(cl,'data/curated/current_label_reviewed_positive_clauses.csv',na='')
write_csv(linked,'data/curated/current_label_clause_PT_review.csv',na='')
write_csv(features,'data/curated/current_label_PT_supplement.csv',na='')
write_csv(features,'outputs/11g_current_label_PT_features.csv',na='')
write_csv(if(length(failures))bind_rows(failures)else tibble(drugbank_id=character(),reason=character()),'outputs/11g_clause_followup_queue.csv')
summary<-coverage |> filter(drugbank_id %in% primary) |>
 left_join(features |> distinct(drugbank_id,pt_cui) |> count(drugbank_id,name='supplement_PT_terms'),by='drugbank_id')
write_csv(summary,'outputs/11g_supplement_coverage.csv',na='')
stopifnot(nrow(summary)==55,all(features$present==1L),all(features$pt_cui %in% d$pt_cui),
 !any(features$pt_name %in% c('Breast cancer','Ovarian cancer','Small cell lung cancer','Breast cancer stage IV')),
 !anyDuplicated(features[c('drugbank_id','pt_cui','regimen_scope','evidence_group')]),
 length(unique(features$evidence_group[features$drugbank_id %in% c('DB12282','DB15254')]))==1L)
report<-c('SCRIPT 11g - CONTEXT-REVIEWED REGULATORY LABEL SUPPLEMENT',
 paste('Primary seed drugs with reviewed positive clauses:',n_distinct(cl$drugbank_id)),
 paste('Primary seed drugs with supplementary PT features:',n_distinct(features$drugbank_id)),
 paste('Distinct supplementary seed/PT pairs:',nrow(distinct(features,drugbank_id,pt_cui))),
 paste('Narrow single-agent/monotherapy profiles:',n_distinct(features$drugbank_id[features$comparison_stratum=='explicit_monotherapy_or_single_agent'])),
 paste('Clause-boundary followups:',length(failures)),
 'Inputs: identity-supported FDA/DailyMed labels; explicit co-pack moiety evidence; named EMA Pepaxti overview.',
 'Outputs: reviewed positive clauses, per-alias mapping decisions and separate current-label PT supplement (data/curated/).',
 'Versioned compact feature export: outputs/11g_current_label_PT_features.csv; coverage: outputs/11g_supplement_coverage.csv.',
 'Context decisions and delimiters: config/current_label_positive_clause_rules.csv.',
 'Only explicit positive reaction lists are promoted, with unambiguous existing dictionary mappings.',
 'Cohort, indication and formulation prefixes are removed before matching; full source clauses remain available.',
 'Literal hot-flush variants are manually resolved to C0600142 using the supplied MedDRA term-CUI rows; broader alternatives remain held.',
 'Other broad/ambiguous matches remain held. This is an assistant-reviewed common-reaction subset, not exhaustive adverse-event curation.',
 'Combination and mixed-cohort features are labelled context only; individual causality is not inferred.',
 'The avutometinib/defactinib co-pack uses ONE shared evidence group, not two independent confirmations.',
 'Pepaxti is an EU overview supplement, not a current US label; its publication date is not inferred from retrieval.',
 'No incidence values are pooled or extracted. Different formulations/routes remain source-labelled.',
 'Original SIDER tables remain unchanged. Missing features are not negative observations.',
 'Next: run a descriptive term-overlap pilot with explicit monotherapy profiles, keeping context-only results separate.')
writeLines(report,'outputs/11g_regulatory_supplement_summary.txt')
cat(paste0('\nSCRIPT 11g RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/11g_review_label_supplement.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
