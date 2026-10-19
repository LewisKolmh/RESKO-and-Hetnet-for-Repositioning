# Expand assay/ligand provenance and perform exact identity joins; source claims remain review candidates.
library(readr);library(dplyr);library(tibble);library(stringr)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
source('scripts/helpers/target_public_get.R')
paths<-c('data/curated/eef1a_chembl_activity_records.csv','outputs/16f_target_evidence_handoff_summary.csv',
 'data/curated/seed_identity_crosscheck_reviewed.csv','config/eef1a_target_registry.csv',
 'outputs/17a_STITCH_resolution_summary.csv','data/raw/target_evidence/uniprot_P68104.rds','data/raw/target_evidence/uniprot_Q05639.rds')
before<-tools::md5sum(paths);a<-r(paths[1]);inventory<-r(paths[2]);seeds<-r(paths[3]);targets<-r(paths[4])
one<-function(z)if(length(z))as.character(z[[1]])else NA_character_
root<-'data/raw/target_evidence';logs<-list()
get<-function(resource,id){url<-paste0('https://www.ebi.ac.uk/chembl/api/data/',resource,'/',id,'.json')
 z<-target_public_get(url,file.path(root,paste0('chembl_',resource,'_',id,'.rds')))
 logs[[length(logs)+1L]]<<-tibble(url=z$url,retrieved_utc=z$retrieved_utc,status=z$status);z}
molecules<-list();assays<-list();docs<-list()
ids<-sort(unique(a$molecule_chembl_id))
for(i in seq_along(ids)){
 z<-get('molecule',ids[i]);v<-z$body;stopifnot(v$molecule_chembl_id==ids[i])
 molecules[[i]]<-tibble(molecule_chembl_id=ids[i],molecule_name=one(v$pref_name),molecule_type=one(v$molecule_type),
 standard_inchikey=one(v$molecule_structures$standard_inchi_key),standard_inchi=one(v$molecule_structures$standard_inchi),
 canonical_smiles=one(v$molecule_structures$canonical_smiles),parent_chembl_id=one(v$molecule_hierarchy$parent_chembl_id),
 source_url=z$url,retrieved_utc=z$retrieved_utc)
 if(i%%5L==0L)message('Ligand records: ',i,' / ',length(ids))
}
ids<-sort(unique(a$assay_chembl_id))
for(i in seq_along(ids)){
 z<-get('assay',ids[i]);v<-z$body;stopifnot(v$assay_chembl_id==ids[i])
 fields<-c('assay_chembl_id','assay_type','description','target_chembl_id','confidence_score','relationship_type',
 'assay_organism','assay_tax_id','assay_cell_type','assay_tissue','assay_subcellular_fraction','assay_test_type',
 'bao_format','document_chembl_id','src_id')
 assays[[i]]<-as_tibble(setNames(lapply(fields,function(f)one(v[[f]])),fields)) |>
 mutate(source_url=z$url,retrieved_utc=z$retrieved_utc)
 if(i%%5L==0L)message('Assay records: ',i,' / ',length(ids))
}
ids<-sort(unique(na.omit(a$document_chembl_id)))
for(i in seq_along(ids)){
 z<-get('document',ids[i]);v<-z$body;stopifnot(v$document_chembl_id==ids[i])
 fields<-c('document_chembl_id','pubmed_id','doi','title','journal','year','doc_type')
 docs[[i]]<-as_tibble(setNames(lapply(fields,function(f)one(v[[f]])),fields)) |>
 mutate(source_url=z$url,retrieved_utc=z$retrieved_utc)
}
m<-bind_rows(molecules);s<-bind_rows(assays);d<-bind_rows(docs)
# Match current full representations only. No parent/salt/connectivity/name mergers.
links<-inventory |> filter(!is.na(current_inchikey),structure_mapping_blocked=='FALSE') |>
 inner_join(m |> filter(!is.na(standard_inchikey)),by=c('current_inchikey'='standard_inchikey'),relationship='many-to-many') |>
 select(entity_id,origin,source_name,current_record_name,drugbank_ids,current_inchikey,identity_status,
 molecule_chembl_id,molecule_name,parent_chembl_id) |> distinct()
# Source identifier agreements are a separate review channel, not a substitute for structure validation.
id_links<-seeds |> select(drugbank_id,drug_name,chembl_id,inchikey) |>
 inner_join(m,by=c('chembl_id'='molecule_chembl_id'),relationship='many-to-many') |>
 transmute(entity_id=paste0('seed:',drugbank_id),drug_name,molecule_chembl_id=chembl_id,
 source_inchikey=inchikey,chembl_inchikey=standard_inchikey,
 identity_status=case_when(is.na(source_inchikey)|is.na(chembl_inchikey)~'identifier_agreement_structure_missing',
 source_inchikey==chembl_inchikey~'identifier_and_full_key_agreement',TRUE~'identifier_agreement_full_key_disagreement_review'))
# Pin full assay metadata alongside each activity; conflicting isoform descriptions cannot silently pass.
evidence<-a |> left_join(s |> rename(assay_record_target=target_chembl_id,assay_record_description=description,
 assay_record_type=assay_type,assay_record_document=document_chembl_id,assay_record_source_url=source_url,
 assay_record_retrieved_utc=retrieved_utc,assay_src_id=src_id,assay_bao_format=bao_format),by='assay_chembl_id') |>
 left_join(m |> select(molecule_chembl_id,molecule_name,standard_inchikey,parent_chembl_id,molecule_type),by='molecule_chembl_id') |>
 mutate(isoform_description_conflict=case_when(queried_gene=='EEF1A2'~str_detect(str_to_lower(assay_record_description),'eef1a1|ef1.alpha.?1|elongation factor 1.alpha 1'),
 queried_gene=='EEF1A1'~str_detect(str_to_lower(assay_record_description),'eef1a2|ef1.alpha.?2|elongation factor 1.alpha 2'),TRUE~FALSE),
 source_assay_target_agreement=target_chembl_id==assay_record_target,
 evidence_review_status=if_else(isoform_description_conflict,'hold_isoform_description_target_conflict','pending_primary_assay_review'))
stopifnot(nrow(evidence)==nrow(a),!anyDuplicated(evidence$activity_id),all(evidence$source_assay_target_agreement),
 all(links$current_inchikey %in% na.omit(m$standard_inchikey)),!any(links$entity_id=='sider:CID100135113'),
 identical(before,tools::md5sum(paths)))
inventory_evidence<-links |> inner_join(evidence,by=c('molecule_chembl_id','molecule_name','parent_chembl_id'),relationship='many-to-many')
# Resolve the STITCH alternative protein ID via independently retrieved current UniProt Ensembl cross-references.
sp<-r('outputs/17a_STITCH_resolution_summary.csv') |> filter(preferred_gene_and_taxon_match=='TRUE',transport=='HTTP_public_query')
crosswalk<-list()
for(i in seq_len(nrow(sp))){
 t<-targets |> filter(gene==sp$gene[i]);u<-readRDS(file.path(root,paste0('uniprot_',t$uniprot,'.rds')))$body
 ensembl<-Filter(function(v)identical(v$database,'Ensembl'),u$uniProtKBCrossReferences)
 protein_ids<-unlist(lapply(ensembl,function(v)vapply(Filter(function(z)identical(z$key,'ProteinId'),v$properties),function(z)one(z$value),character(1))),use.names=FALSE)
 mapped<-sp$returned_string_id[i] %in% paste0('9606.',sub('[.][0-9]+$','',protein_ids))
 canonical<-sp$returned_string_id[i] %in% strsplit(t$current_STRING_ids,'|',fixed=TRUE)[[1]]
 crosswalk[[i]]<-tibble(gene=sp$gene[i],uniprot=t$uniprot,stitch_protein_id=sp$returned_string_id[i],
 canonical_STRING_crossreference_match=canonical,current_UniProt_Ensembl_protein_match=mapped,
 crosswalk_supported=canonical||mapped,sequence_variant_equivalence_established=FALSE,
 source_query_url=sp$url[i],uniprot_source_url=paste0('https://rest.uniprot.org/uniprotkb/',t$uniprot,'.json'))
}
crosswalk<-bind_rows(crosswalk) |> distinct(gene,stitch_protein_id,.keep_all=TRUE)
stopifnot(all(crosswalk$crosswalk_supported),setequal(crosswalk$gene,c('EEF1A1','EEF1A2')))
write_csv(crosswalk,'config/eef1a_STITCH_protein_crosswalk.csv')
write_csv(crosswalk,'outputs/17b_STITCH_crosswalk_summary.csv')
coverage<-inventory |> select(entity_id,origin,source_name,current_inchikey,structure_mapping_blocked) |>
 left_join(links |> count(entity_id,name='exact_target_ligand_links'),by='entity_id') |>
 mutate(exact_target_ligand_links=tidyr::replace_na(exact_target_ligand_links,0L),
 source_coverage_status=case_when(structure_mapping_blocked=='TRUE'~'structure_quarantined_named_context_only',
 is.na(current_inchikey)~'structure_key_missing',exact_target_ligand_links>0~'exact_structure_source_activity_for_review',
 TRUE~'not_found_in_exact_two_target_ChEMBL_query_NOT_negative_evidence'))
write_csv(coverage,'outputs/17b_inventory_target_coverage.csv',na='')
write_csv(m,'data/curated/eef1a_chembl_ligand_registry.csv',na='')
write_csv(s,'data/curated/eef1a_chembl_assay_registry.csv',na='')
write_csv(d,'data/curated/eef1a_chembl_document_registry.csv',na='')
write_csv(evidence,'data/curated/eef1a_chembl_enriched_activity_evidence.csv',na='')
write_csv(links,'data/curated/eef1a_inventory_exact_ligand_links.csv',na='')
write_csv(id_links,'data/curated/eef1a_inventory_identifier_review.csv',na='')
write_csv(inventory_evidence,'data/curated/eef1a_inventory_source_activity_evidence.csv',na='')
write_csv(bind_rows(logs),'data/curated/target_assay_ligand_API_provenance.csv',na='')
write_csv(m |> select(-standard_inchi,-canonical_smiles),'outputs/17b_ligand_identity_summary.csv',na='')
write_csv(links,'outputs/17b_inventory_ligand_match_summary.csv',na='')
write_csv(id_links,'outputs/17b_identifier_crosscheck_summary.csv',na='')
write_csv(evidence |> count(queried_gene,evidence_review_status,confidence_score,name='activity_rows'),
 'outputs/17b_assay_review_coverage.csv',na='')
write_csv(d |> select(document_chembl_id,pubmed_id,doi,year,doc_type,source_url),
 'outputs/17b_source_documents_summary.csv',na='')
write_csv(evidence |> filter(isoform_description_conflict) |> select(activity_id,queried_gene,queried_uniprot,target_chembl_id,
 assay_chembl_id,molecule_chembl_id,evidence_review_status),'outputs/17b_isoform_conflicts_summary.csv',na='')
write_csv(tibble(file=paths,md5=unname(before)),'outputs/17b_input_summary.csv')
report<-c('SCRIPT 17b - FULL ASSAY/LIGAND EXPANSION AND INVENTORY IDENTITY LINKS',
 paste('Ligands:',nrow(m)),paste('Assay records:',nrow(s)),paste('Source documents:',nrow(d)),
 paste('Enriched activity rows:',nrow(evidence)),paste('Inventory entities with exact full-key ligand matches:',n_distinct(links$entity_id)),
 paste('Inventory source-activity rows linked:',nrow(inventory_evidence)),
 paste('Activity rows held for isoform-description conflicts:',sum(evidence$isoform_description_conflict)),
 paste('STITCH target gene/protein crosswalks supported by current UniProt references:',nrow(crosswalk)),
 'EEF1A1 alternative STITCH protein ID is present in current UniProt Ensembl cross-references; canonical sequence-variant equivalence is not asserted.',
 'All available assay fields and complete cached source JSON are retained; confidence score and organism are not dropped.',
 'Exact full-key links attach source activities, not automatically confirmed binding or inhibitory effects. ChEMBL-ID-only links stay separate.',
 'Parent/child compounds, undefined stereochemistry, salts and source-name aliases are not silently merged.',
 'Assay/target identity contradictions are held, not reassigned between EEF1A1 and EEF1A2 by text alone.',
 'Primary documents, assay design, activity relations/units and controls must be reviewed before physical binding, affinity or selectivity is asserted.',
 'No RT-PPI disruption, allostery or translation-sparing claim is created from these activity rows.',
 'Missing exact-target ChEMBL records is not negative evidence; target families, other annotations and literature remain outside this bounded query.',
 'Outputs: eef1a_chembl_* registries/evidence and eef1a_inventory_* links locally; compact identity/conflict/document/coverage reports in outputs/17b_*.',
 'Next: validate STITCH protein crosswalks and association channels; review source documents and add positive-control literature before prioritising binding candidates.')
writeLines(report,'outputs/17b_assay_ligand_summary.txt')
writeLines('PASS: unique source assay/ligand/activity IDs; exact source-target agreement; full-key inventory links; protected inputs unchanged and sacrosidase wrong structure excluded.',
 'outputs/17b_assay_ligand_validation.txt')
cat(paste0('\nSCRIPT 17b RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',
 'Script: scripts/17b_expand_assays_and_match_ligands.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
