# Paper-aligned completion audit; source acquisition alone is not method replication.
library(readr);library(dplyr);library(tibble);library(httr2)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
inputs<-c('config/database_manifest.csv','config/RESKO_evidence_scope.csv',
 'outputs/16f_target_evidence_handoff_summary.csv','outputs/18c_six_lead_summary.csv')
before<-tools::md5sum(inputs)
manifest<-r(inputs[1]);scope<-r(inputs[2]);inventory<-r(inputs[3]);leads<-r(inputs[4])
stopifnot(scope$value[scope$setting=='require_published_eef1a_affinity_for_RESKO_inclusion']=='FALSE',
 nrow(inventory)==200L,nrow(leads)==6L)
url<-'https://sure.sunderland.ac.uk/9701/1/DrugRepositionMcGarry.pdf'
path<-'data/raw/method_reference/DrugRepositionMcGarry_accepted_2018.pdf'
dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
if(!file.exists(path)){
 temporary<-paste0(path,'.partial')
 z<-request(url)|>req_timeout(60)|>req_retry(max_tries=3,retry_on_failure=TRUE)|>req_perform(path=temporary)
 stopifnot(resp_status(z)==200L,file.rename(temporary,path))
}
con<-file(path,'rb');magic<-readBin(con,'raw',n=5);close(con);stopifnot(rawToChar(magic)=='%PDF-')
source<-tibble(title='RESKO: Repositioning drugs by using side effects and knowledge from ontologies',
 doi='10.1016/j.knosys.2018.06.017',source_url=url,version='Accepted manuscript deposited by authors',
 local_file=path,md5=unname(tools::md5sum(path)),checked_date='2026-10-20',
 usage='Methods reference only; source PDF remains local, no republication')
write_csv(source,'config/RESKO_method_reference.csv')
steps<-tribble(~step,~paper_section,~required_component,~project_status,~next_work,
 '01','3.1','Disease UMLS identifiers and conventional disease-treatment cohort','partial',
 'Validate cancer/subtype CUI mappings and provenance; current treatment review used indication text, not complete UMLS-driven extraction',
 '02','3.1 / 3.2 / Algorithm 1','SIDER drug/PT profiles and joint side-effect pattern per disease','partial',
 'Reproduce original disease-cohort SIDER baseline separately from label supplements and exploratory pairwise Jaccard',
 '03','3.1 / Figure 3','Remove common side effects using declared interpretation','method_ambiguity',
 'Text describes top 10 percent most common types, Figure 3 says prevalence above 10 percent of drugs; run labelled sensitivity branches',
 '04','3.2 / Algorithm 1','Joint pattern coverage and candidate retrieval at 50 percent threshold','partial',
 'Define intersection denominator and threshold boundary; reproduce primary conventional-seed baseline per cancer, keep investigational tiers separate',
 '05','3.1 / 3.5','ChEMBL structures, 2048-bit atom-pair fingerprints and Tanimoto','partial_identity_only',
 'Complete seed/candidate structure records and implement paper-matched fingerprints; identity checks do not provide similarity matrices',
 '06','3.1 / 3.3','DrugBank/STITCH candidate targets and human STRING protein networks','partial_target_subset_only',
 'Extract all selected seed/candidate protein associations and human PPI context, not only direct eEF1A neighbours',
 '07','3.1','Reactome human protein-pathway annotations','not_acquired',
 'Acquire release/source provenance, map all network proteins and construct pathway matrices',
 '08','3.4','GO ontology plus human gene annotations and Wang similarity','not_acquired',
 'Acquire ontology and annotation releases separately, propagate accepted relationships and implement semantic similarity',
 '09','3.6','DO ontology plus gene-disease annotations and hypergeometric enrichment','not_acquired',
 'Acquire ontology and compatible human disease-gene annotations; specify gene universe and valid test family',
 '10','3.3','Network closeness and clustering statistics','not_implemented',
 'Construct reproducible networks and define disconnected-node/edge-confidence handling',
 '11','3.7','Block-matrix integration / constrained matrix factorisation','not_implemented',
 'Verify dimensions, objective, constraints, hyperparameters and optimisation details before coding; do not replace with an arbitrary weighted sum',
 '12','3.8','Latent-network evaluation, candidate ranking and literature validation','not_implemented',
 'Define validation labels/splits before assessment, prevent source leakage, distinguish unknown from true negative labels',
 '13','User-requested extension','Benjamini-Hochberg adjusted enrichment p-values','pending_valid_test_family',
 'Apply p.adjust(method="BH") to a prespecified family of valid tests; retain raw p-values, effects and universe',
 '14','User-requested target follow-up','Human EEF1A1/2 binding, selectivity, RT-PPI and allostery annotation','bounded_followup_complete_not_exhaustive',
 'Record available evidence only; missing affinity never excludes a RESKO candidate or prevents integration completion')
steps<-steps |> mutate(method_doi=source$doi,paper_source=url,reviewed_on='2026-10-20')
write_csv(steps,'config/RESKO_method_completion_checklist.csv')
write_csv(steps,'outputs/19_method_completion_summary.csv')
write_csv(manifest |> select(source_id,source_name,status,release_version),
 'outputs/19_database_coverage.csv',na='')
ambiguities<-tribble(~issue,~paper_locations,~handling,
 'Common-side-effect pruning','3.1 versus Figure 3','Keep top-decile term removal and above-10-percent drug-prevalence removal separate sensitivity settings',
 'Minimum seed side effects','3.2 prose versus Algorithm 1 lines 4/13','Record >3 versus >=3 explicitly and report boundary sensitivity',
 'Candidate similarity boundary','Methods summary versus Algorithm 1','Report >50 percent and >=50 percent separately at exact ties; retain numerator/denominator',
 'Fallback seed reduction','Figure 3 versus Algorithm 1','Do not choose a favourable half of the cohort after viewing candidates; predeclare any fallback and preserve full-cohort result',
 'Matrix factorisation reproducibility','3.7 / equation 7','Obtain/verify implementation details; disclose adaptations and do not claim bit-for-bit reproduction',
 'Cancer application versus original Alzheimer example','Disease-forward workflow','Describe cancer runs as a method adaptation; reproduce an Alzheimer reference run if original data and settings are available')
write_csv(ambiguities,'config/RESKO_method_ambiguities.csv')
write_csv(ambiguities,'outputs/19_method_ambiguities_summary.csv')
stopifnot(identical(before,tools::md5sum(inputs)))
write_csv(tibble(file=inputs,md5=unname(before)),'outputs/19_input_summary.csv')
report<-c('SCRIPT 19 - PAPER-ALIGNED RESKO COMPLETION AUDIT',
 'Method reference: McGarry et al. (2018), DOI 10.1016/j.knosys.2018.06.017; accepted manuscript acquired and hashed locally.',
 'Nine original source registrations audited; UMLS, STRING, Reactome, GO and DO are still not acquired in the manifest.',
 'DrugBank and SIDER data are available; ChEMBL and STITCH integration is partial, not full candidate-network completion.',
 '14 method/extension checklist components and six interpretation/reproducibility issues documented.',
 'Published eEF1A affinity is an optional annotation, never a candidate inclusion gate. All six association leads and original 200 entities retained.',
 'Core unfinished work includes disease-wise SIDER pattern baseline, atom-pair chemical similarity, broad human protein/pathway/ontology integration, matrix factorisation and validated ranking.',
 'Current-label supplements, cancer seed expansion and target follow-up remain labelled adaptations; STITCH is already an original paper source.',
 'No claim that RESKO is complete or that exploratory Jaccard/association evidence reproduces the published integrated ranking.',
 'Outputs: method checklist, original database coverage, ambiguity register, input hashes and validation in outputs/19_* and config/RESKO_*.',
 'Next core acquisition: Reactome human protein-pathway mappings, then broad seed/candidate targets and remaining human network/ontology sources.')
writeLines(report,'outputs/19_method_audit_summary.txt')
writeLines(c('PASS: source PDF magic/header and checksum, original source manifest retained, scope forbids affinity exclusion.',
 'PASS: all 200 entities/six leads preserved; protected input hashes unchanged; method gaps and paper ambiguities explicit.'),
 'outputs/19_method_audit_validation.txt')
notes<-readLines('Notes.txt',warn=FALSE);start<-grep('^CURRENT CHECKPOINT',notes)[1];end<-grep('^PROJECT AIM',notes)[1]
stopifnot(start<end);writeLines(c(notes[seq_len(start-1)],'CURRENT CHECKPOINT - RETURN TO CORE RESKO; SCRIPT 19 METHOD AUDIT COMPLETE',report,'',notes[end:length(notes)]),'Notes.txt')
cat(paste0('\nSCRIPT 19 RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/maintenance/19_audit_RESKO_method_completion.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
