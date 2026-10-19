# Resolve exact human eEF1A isoforms and acquire source activity records for review.
library(readr);library(dplyr);library(tibble);library(httr2);library(stringr)
if(!file.exists('outputs/16f_target_evidence_handoff_summary.csv'))stop('Run from project root after section 16.')
source('scripts/helpers/target_public_get.R')
root<-'data/raw/target_evidence';dir.create(root,recursive=TRUE,showWarnings=FALSE)
expected<-tibble(gene=c('EEF1A1','EEF1A2'),uniprot=c('P68104','Q05639'),taxon=9606L)
write_csv(expected,'config/eef1a_target_accessions.csv')
one<-function(z)if(length(z))as.character(z[[1]])else NA_character_
logs<-list();targets<-list();activities<-list();coverage<-list();stitch<-list()
get<-function(url,label){z<-target_public_get(url,file.path(root,paste0(label,'.rds')))
 logs[[length(logs)+1L]]<<-tibble(url=z$url,status=z$status,retrieved_utc=z$retrieved_utc);z}
status<-get('https://www.ebi.ac.uk/chembl/api/data/status.json','chembl_status')
# STITCH's HTTPS certificate is checked normally. HTTP public lookups are recorded explicitly, never treated as protected transport.
probe<-function(url,label,gene,expected_string_ids){
 path<-file.path(root,paste0('stitch_probe_',label,'.rds'))
 if(file.exists(path)){z<-readRDS(path);stopifnot(identical(z$url,url))}else{
  result<-tryCatch(request(url) |> req_timeout(25) |> req_error(is_error=function(r)FALSE) |> req_perform(),error=function(e)e)
  z<-if(inherits(result,'error'))list(url=url,status=NA_integer_,error=conditionMessage(result),body=NA_character_)else
   list(url=url,status=resp_status(result),error=NA_character_,body=resp_body_string(result))
  z$retrieved_utc<-format(Sys.time(),tz='UTC',usetz=TRUE);saveRDS(z,path)
 }
 rows<-tibble();if(!is.na(z$status)&&z$status==200L&&!is.na(z$body)&&startsWith(trimws(z$body),'stringId\t'))
  rows<-read_tsv(I(z$body),col_types=cols(.default=col_character()),show_col_types=FALSE)
 preferred<-if('preferredName' %in% names(rows))paste(unique(rows$preferredName),collapse='|')else NA_character_
 id<-if('stringId' %in% names(rows))paste(unique(rows$stringId),collapse='|')else NA_character_
 name_ok<-nrow(rows)==1L&&!is.na(preferred)&&identical(preferred,gene)&&rows$ncbiTaxonId[[1]]=='9606'
 accession_ok<-name_ok&&id %in% expected_string_ids
 tibble(gene,url,status=z$status,error=z$error,retrieved_utc=z$retrieved_utc,transport=ifelse(startsWith(url,'https'),'HTTPS','HTTP_public_query'),
 returned_preferred_name=preferred,returned_string_id=id,preferred_gene_and_taxon_match=name_ok,
 uniprot_STRING_crossreference_match=accession_ok,
 disposition=if(accession_ok)'verified_target_identifier_only'else if(name_ok)'review_historical_identifier_crosswalk'else 'rejected_or_unavailable_target_lookup')
}
for(i in seq_len(nrow(expected))){
 e<-expected[i,];u<-get(paste0('https://rest.uniprot.org/uniprotkb/',e$uniprot,'.json'),paste0('uniprot_',e$uniprot))
 gene_names<-vapply(u$body$genes,function(g)one(g$geneName$value),character(1))
 stopifnot(u$body$primaryAccession==e$uniprot,u$body$organism$taxonId==e$taxon,e$gene %in% gene_names)
 refs<-u$body$uniProtKBCrossReferences
 string_ids<-vapply(Filter(function(a)identical(a$database,'STRING'),refs),function(a)one(a$id),character(1))
 z<-get(paste0('https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession=',e$uniprot,'&limit=100'),paste0('chembl_target_',e$uniprot))
 stopifnot(is.null(z$body$page_meta$`next`),length(z$body$targets)==z$body$page_meta$total_count)
 chosen<-Filter(function(a)identical(as.integer(a$tax_id),9606L)&&identical(a$target_type,'SINGLE PROTEIN')&&
 length(a$target_components)==1L&&identical(a$target_components[[1]]$accession,e$uniprot),z$body$targets)
 stopifnot(length(chosen)==1L);t<-chosen[[1]];tid<-t$target_chembl_id
 targets[[i]]<-tibble(gene=e$gene,uniprot=e$uniprot,taxon=9606L,chembl_target_id=tid,preferred_name=t$pref_name,
 target_type=t$target_type,uniprot_entry_version=one(u$body$entryAudit$entryVersion),
 sequence_version=one(u$body$entryAudit$sequenceVersion),sequence_length=u$body$sequence$length,
 current_STRING_ids=paste(string_ids,collapse='|'),target_identity_status='human_single_component_exact_accession',
 uniprot_url=u$url,chembl_target_url=z$url,retrieved_utc=z$retrieved_utc)
 url<-paste0('https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=',tid,'&limit=1000');page<-0L;total<-NA_integer_;seen<-character();a_list<-list()
 repeat{
  if(url %in% seen)stop('Repeated activity page');seen<-c(seen,url);page<-page+1L
  a<-get(url,paste0('chembl_activity_',tid,'_page',page));if(page==1L)total<-as.integer(a$body$page_meta$total_count)
  if(total>10000L)stop('Activity count exceeds declared acquisition bound; increase explicitly rather than truncate')
  for(v in a$body$activities){
   fields<-c('activity_id','molecule_chembl_id','molecule_pref_name','parent_molecule_chembl_id','assay_chembl_id',
    'assay_description','assay_type','assay_variant_accession','assay_variant_mutation','target_chembl_id',
    'target_organism','target_tax_id','standard_type','standard_relation','standard_value','standard_units','standard_flag',
    'pchembl_value','activity_comment','data_validity_comment','potential_duplicate','document_chembl_id','document_journal','document_year',
    'bao_format','bao_label','src_id')
   row<-as_tibble(setNames(lapply(fields,function(f)one(v[[f]])),fields))
   a_list[[length(a_list)+1L]]<-row |> mutate(queried_gene=e$gene,queried_uniprot=e$uniprot,source_page_url=a$url,retrieved_utc=a$retrieved_utc)
  }
  next_page<-a$body$page_meta$`next`;if(is.null(next_page))break
  url<-if(startsWith(next_page,'https://'))next_page else paste0('https://www.ebi.ac.uk',next_page)
 }
 act<-bind_rows(a_list);stopifnot(nrow(act)==total,all(act$target_chembl_id==tid),!anyDuplicated(act$activity_id))
 activities[[i]]<-act
 coverage[[i]]<-tibble(gene=e$gene,uniprot=e$uniprot,chembl_target_id=tid,activity_rows=nrow(act),
 unique_molecules=n_distinct(act$molecule_chembl_id),unique_assays=n_distinct(act$assay_chembl_id),pages=page,complete_target_query=TRUE)
 for(identifier in c(e$gene,e$uniprot))stitch[[length(stitch)+1L]]<-probe(
  paste0('http://stitch.embl.de/api/tsv/resolve?identifier=',identifier,'&species=9606'),identifier,e$gene,string_ids)
 if(i==1L)stitch[[length(stitch)+1L]]<-probe(
  paste0('https://stitch.embl.de/api/tsv/resolve?identifier=',e$uniprot,'&species=9606'),'HTTPS_check',e$gene,string_ids)
 message('Target acquired: ',e$gene,'; ',nrow(act),' ChEMBL activity rows')
}
targets<-bind_rows(targets);a<-bind_rows(activities);sp<-bind_rows(stitch)
# An assay B label can describe effects on tRNA binding, not direct compound/protein binding.
a<-a |> mutate(source_evidence_status='unreviewed_activity_record',
 review_route=case_when(str_detect(str_to_lower(assay_description),'surface plasmon|isothermal titration|dissociation|binding affinity')~'candidate_physical_binding_method_requires_assay_review',
 str_detect(str_to_lower(assay_description),'trna|translation|protein synthesis|nucleotide exchange')~'functional_or_partner_binding_effect_requires_review',
 TRUE~'assay_context_review_required'),direct_eef1a_binding_established=NA,
 isoform_selectivity_established=NA,allostery_established=NA,RT_PPI_disruption_established=NA)
write_csv(targets,'config/eef1a_target_registry.csv',na='')
write_csv(a,'data/curated/eef1a_chembl_activity_records.csv',na='')
write_csv(bind_rows(logs),'data/curated/target_source_API_provenance.csv',na='')
write_csv(targets,'outputs/17a_target_identity_summary.csv',na='')
write_csv(bind_rows(coverage),'outputs/17a_chembl_activity_coverage.csv',na='')
write_csv(a |> count(queried_gene,assay_type,review_route,name='activity_rows'),'outputs/17a_evidence_review_summary.csv',na='')
write_csv(sp,'outputs/17a_STITCH_resolution_summary.csv',na='')
# Source status is partial; target-resolution responses do not count as STITCH interaction acquisition.
manifest<-read_csv('config/database_manifest.csv',col_types=cols(.default=col_character()))
manifest$status[manifest$source_id=='chembl']<-'partial_identity_and_eef1a_target_activities_acquired'
manifest$release_version[manifest$source_id=='chembl']<-one(status$body$chembl_db_version)
manifest$status[manifest$source_id=='stitch']<-'identifier_resolution_probed_associations_not_acquired'
manifest$retrieved_date[manifest$source_id %in% c('chembl','stitch')]<-format(Sys.Date())
write_csv(manifest,'config/database_manifest.csv',na='')
report<-c('SCRIPT 17a - TARGET IDENTITY AND SOURCE ACQUISITION',
 paste('Resolved human isoforms:',nrow(targets)),paste('ChEMBL activity rows acquired:',nrow(a)),
 paste('Unique target-assay pairs:',nrow(a |> distinct(target_chembl_id,assay_chembl_id))),
 paste('Verified STITCH identifier responses:',sum(sp$uniprot_STRING_crossreference_match)),
 'Exact human UniProt accessions are required in single-component ChEMBL targets; isoforms remain separate.',
 'Activity pagination is exhaustive for these two target IDs; target families, other species and indirectly annotated records are outside this query scope.',
 'Raw ChEMBL records and UniProt source responses remain cached locally with URLs, dates and source versions.',
 'STITCH resolution is validated against gene/taxon and UniProt STRING cross-references. Wrong gene resolutions are rejected.',
 'STITCH HTTPS certificate validation is not bypassed; public HTTP lookups have explicit transport provenance.',
 'No STITCH chemical-protein associations are acquired by identifier resolution alone.',
 'Assay type B is not proof of direct compound binding; partner-binding and translation effects remain separate review routes.',
 'No inventory compound has been promoted to eEF1A binding, isoform selectivity, allostery or RT-PPI disruption.',
 'Outputs: target registry (config/); eef1a_chembl_activity_records.csv and API provenance (local); outputs/17a_* coverage/identity/resolution reports.',
 'Next: retrieve complete assay and ligand records, cross-check inventory identities, and review physical-binding versus functional evidence; acquire STITCH associations only after valid target resolution.')
writeLines(report,'outputs/17a_target_acquisition_summary.txt')
writeLines('PASS: exact human accessions and one-component target identity, exhaustive activity page totals, unique activity IDs, and checked STITCH returned identities.',
 'outputs/17a_target_acquisition_validation.txt')
cat(paste0('\nSCRIPT 17a RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/17a_acquire_target_records.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
