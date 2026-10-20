# Acquire complete official mapping files and retain the human pathway layer.
library(readr);library(dplyr);library(tibble);library(httr2)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
root<-'data/raw/reactome';dir.create(root,recursive=TRUE,showWarnings=FALSE)
version_url<-'https://reactome.org/ContentService/data/database/version'
get_version<-function(){
 z<-request(version_url)|>req_timeout(30)|>req_retry(max_tries=3,retry_on_failure=TRUE)|>req_perform()
 v<-trimws(resp_body_string(z));stopifnot(grepl('^[0-9]+$',v));v
}
files<-c('UniProt2Reactome.txt','UniProt2Reactome_All_Levels.txt','ReactomePathways.txt','ReactomePathwaysRelation.txt')
# Reuse the frozen snapshot on rerun; do not silently mix it with a newer release.
manifest_path<-file.path(root,'acquisition_manifest.csv')
if(file.exists(manifest_path)){
 acquired<-r(manifest_path)
 stopifnot(setequal(acquired$file,files),all(file.exists(acquired$local_file)),
 all(unname(tools::md5sum(acquired$local_file))==acquired$md5))
}else{
 before_version<-get_version();rows<-list()
 for(f in files){
  url<-paste0('https://reactome.org/download/current/',f);path<-file.path(root,f);temporary<-paste0(path,'.partial')
  message('Downloading ',f)
  response<-request(url)|>req_timeout(300)|>req_retry(max_tries=3,retry_on_failure=TRUE)|>req_perform(path=temporary)
  stopifnot(resp_status(response)==200L,file.info(temporary)$size>1000,file.rename(temporary,path))
  rows[[length(rows)+1L]]<-tibble(file=f,local_file=path,url,
   reported_database_version=before_version,version_url,
   retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
   last_modified=resp_header(response,'last-modified'),
   file_bytes=as.character(file.info(path)$size),md5=unname(tools::md5sum(path)),
   release_scope='Current mapping snapshot; version service checked before and after acquisition; file hashes freeze actual contents')
 }
 after_version<-get_version();stopifnot(identical(before_version,after_version))
 acquired<-bind_rows(rows);write_csv(acquired,manifest_path,na='')
}
input<-c('outputs/16f_target_evidence_handoff_summary.csv','config/eef1a_target_registry.csv')
before<-tools::md5sum(input)
read_mapping<-function(file){
 z<-read_tsv(file.path(root,file),col_names=c('uniprot','pathway_id','pathway_url','pathway_name','evidence_code','species'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
 stopifnot(nrow(problems(z))==0L,ncol(z)==6L,!anyNA(z$uniprot),!anyNA(z$species))
 list(source_rows=nrow(z),human=z |> filter(species=='Homo sapiens') |> distinct())
}
lowest<-read_mapping(files[1]);all_levels<-read_mapping(files[2])
pathways<-read_tsv(file.path(root,files[3]),col_names=c('pathway_id','pathway_name','species'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
relations<-read_tsv(file.path(root,files[4]),col_names=c('parent_pathway_id','child_pathway_id'),
 col_types=cols(.default=col_character()),quote='',progress=FALSE)
stopifnot(nrow(problems(pathways))==0L,nrow(problems(relations))==0L)
human_pathways<-pathways |> filter(species=='Homo sapiens') |> distinct()
human_relations<-relations |> filter(parent_pathway_id %in% human_pathways$pathway_id,child_pathway_id %in% human_pathways$pathway_id) |> distinct()
stopifnot(!anyDuplicated(human_pathways$pathway_id),
 all(all_levels$human$pathway_id %in% human_pathways$pathway_id),
 all(lowest$human$pathway_id %in% human_pathways$pathway_id),
 all(grepl('^R-HSA-',all_levels$human$pathway_id)),
 nrow(anti_join(lowest$human |> distinct(uniprot,pathway_id),
 all_levels$human |> distinct(uniprot,pathway_id),by=c('uniprot','pathway_id')))==0L)
write_csv(lowest$human,'data/curated/Reactome_human_lowest_level_annotations.csv')
write_csv(all_levels$human,'data/curated/Reactome_human_all_level_annotations.csv')
write_csv(human_pathways,'data/curated/Reactome_human_pathways.csv')
write_csv(human_relations,'data/curated/Reactome_human_pathway_hierarchy.csv')
targets<-r(input[2]);stopifnot(all(c('P68104','Q05639') %in% targets$uniprot))
target_context<-all_levels$human |> inner_join(targets |> select(gene,uniprot),by='uniprot',relationship='many-to-one') |>
 mutate(annotation_scope='Pathway participation, not compound binding or direction of drug effect')
write_csv(target_context,'data/curated/Reactome_eef1a_human_pathway_context.csv')
write_csv(tibble(layer=c('lowest_level','all_levels'),
 source_rows=c(lowest$source_rows,all_levels$source_rows),
 human_annotation_rows=c(nrow(lowest$human),nrow(all_levels$human)),
 distinct_human_proteins=c(n_distinct(lowest$human$uniprot),n_distinct(all_levels$human$uniprot)),
 distinct_human_pathways=c(n_distinct(lowest$human$pathway_id),n_distinct(all_levels$human$pathway_id))),
 'outputs/20a_mapping_coverage.csv')
write_csv(targets |> select(gene,uniprot) |> left_join(target_context |> count(gene,uniprot,name='all_level_annotation_rows'),
 by=c('gene','uniprot'),relationship='one-to-one') |> mutate(all_level_annotation_rows=coalesce(all_level_annotation_rows,0L)),
 'outputs/20a_target_context_summary.csv')
write_csv(acquired |> select(-local_file),'outputs/20a_source_summary.csv',na='')
write_csv(tibble(human_pathway_nodes=nrow(human_pathways),human_hierarchy_edges=nrow(human_relations)),
 'outputs/20a_hierarchy_summary.csv')
manifest<-r('config/database_manifest.csv');k<-manifest$source_id=='reactome';stopifnot(sum(k)==1L)
manifest$status[k]<-'human_protein_pathway_mappings_and_hierarchy_acquired_integration_pending'
manifest$release_version[k]<-paste0('Version service ',unique(acquired$reported_database_version),'; frozen mapping snapshot hashes')
manifest$retrieved_date[k]<-'2026-10-20';manifest$download_url[k]<-'https://reactome.org/download/'
manifest$local_file[k]<-manifest_path
write_csv(manifest,'config/database_manifest.csv',na='')
# Preserve the historical script-19 audit; update the live checklist only.
checklist<-r('config/RESKO_method_completion_checklist.csv');k<-checklist$step=='07'
checklist$project_status[k]<-'human_mappings_and_hierarchy_acquired_integration_pending'
checklist$next_work[k]<-'Join broad selected drug-protein associations to acquired human pathway mappings; construct pathway matrices and valid enrichment tests'
write_csv(checklist,'config/RESKO_method_completion_checklist.csv')
stopifnot(identical(before,tools::md5sum(input)))
write_csv(tibble(file=input,md5=unname(before)),'outputs/20a_input_summary.csv')
report<-c('SCRIPT 20a - CORE RESKO REACTOME ACQUISITION',
 paste('Reported Reactome database version:',unique(acquired$reported_database_version)),
 paste('Official source files acquired:',nrow(acquired)),
 paste('Human lowest-level protein/pathway annotation rows:',nrow(lowest$human)),
 paste('Human all-level protein/pathway annotation rows:',nrow(all_levels$human)),
 paste('Human pathway nodes:',nrow(human_pathways)),paste('Human parent/child pathway edges:',nrow(human_relations)),
 'Lowest-level and all-level annotations are separate; ancestor-expanded rows are not independent experiments.',
 'Human pathway membership does not establish drug binding, inhibition, allostery or treatment direction.',
 'All-human mapping layer is acquired for later broad candidate-network integration, not restricted to eEF1A.',
 'Current data are a versioned/hash-frozen adaptation, not the original 2018 release. Source file modification/retrieval dates are distinct.',
 'Raw source files and full human mapping tables remain local; compact coverage/provenance/hierarchy reports in outputs/20a_*.',
 'No enrichment or BH test performed yet; selected drug-protein sets and an explicit background are required.',
 'Next: broad drug-protein extraction and human network/GO/DO integration; the complete RESKO analysis remains unfinished.')
writeLines(report,'outputs/20a_acquisition_summary.txt')
writeLines(c('PASS: official downloads/checksums and version checks; all files parse without problems.',
 'PASS: human pathway referential integrity; all lowest-level protein/pathway pairs occur in all-level mappings.',
 'PASS: protected inventory/target input hashes unchanged; no binding-affinity gate or statistical ranking introduced.'),
 'outputs/20a_acquisition_validation.txt')
notes<-readLines('Notes.txt',warn=FALSE);start<-grep('^CURRENT CHECKPOINT',notes)[1];end<-grep('^PROJECT AIM',notes)[1]
stopifnot(start<end);writeLines(c(notes[seq_len(start-1)],'CURRENT CHECKPOINT - CORE RESKO CONTINUES; REACTOME SCRIPT 20a COMPLETE',report,'',notes[end:length(notes)]),'Notes.txt')
cat(paste0('\nSCRIPT 20a RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/20a_acquire_reactome.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
