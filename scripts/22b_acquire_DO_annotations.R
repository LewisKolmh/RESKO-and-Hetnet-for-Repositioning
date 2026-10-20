# DO hierarchy plus an explicitly versioned compiled human gene-association context layer.
library(readr);library(dplyr);library(tibble);library(DBI);library(RSQLite);library(igraph)
source('scripts/helpers/ontology_source_get.R');source('scripts/helpers/read_ontology_obo.R')
sources<-bind_rows(
 ontology_source_get('https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo','data/raw/do/doid.obo'),
 ontology_source_get('https://bioconductor.org/packages/3.23/data/annotation/src/contrib/HDO.db_1.0.0.tar.gz','data/raw/do/HDO.db_1.0.0.tar.gz'))
do<-read_ontology_obo(sources$local_file[1],'DOID');active<-do$terms$term_id[!do$terms$is_obsolete]
safe_edges<-do$edges |> filter(relation=='is_a',child %in% active,parent %in% active)
g<-graph_from_data_frame(safe_edges |> select(child,parent),directed=TRUE,vertices=tibble(name=active));stopifnot(is_dag(g))
# Inspect package data/build provenance without installing it or executing its build script.
archive<-sources$local_file[2];wanted<-c('HDO.db/inst/extdata/HDO.sqlite','HDO.db/inst/extdata/get_sqlite.r',
 'HDO.db/DESCRIPTION','HDO.db/man/HDOGENE.Rd')
stopifnot(all(wanted %in% untar(archive,list=TRUE)))
untar(archive,files=wanted,exdir='data/raw/do/HDO_source')
path<-'data/raw/do/HDO_source/HDO.db/inst/extdata/HDO.sqlite'
connection<-dbConnect(SQLite(),path,flags=SQLITE_RO)
annotations<-as_tibble(dbReadTable(connection,'do_gene'))
metadata<-as_tibble(dbReadTable(connection,'metadata'));map_metadata<-as_tibble(dbReadTable(connection,'map_metadata'))
dbDisconnect(connection)
build_code<-readLines('data/raw/do/HDO_source/HDO.db/inst/extdata/get_sqlite.r',warn=FALSE)
stopifnot(any(grepl('olddata',build_code,fixed=TRUE)),any(grepl('DO2EG',build_code,fixed=TRUE)),
 any(grepl('rbind(anno, anno_old)',build_code,fixed=TRUE)))
source_date<-map_metadata$source_date[map_metadata$map_name=='GENE'];stopifnot(length(source_date)==1L)
human_gene_ids<-AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,keytype='ENTREZID')
a<-annotations |> transmute(source_row_id=row_number(),source_term_id=as.character(doid),entrez_id=as.character(gene)) |>
 left_join(resolve_ontology_ids(unique(annotations$doid),do),by='source_term_id',relationship='many-to-one') |>
 mutate(annotation_status=case_when(!entrez_id %in% human_gene_ids~'no_current_human_Entrez_key_hold',
 !term_status %in% c('active_primary_ID','active_documented_alternative_ID')~term_status,
 TRUE~'accepted_compiled_human_gene_disease_context'),
 evidence_layer='HDO_db_compiled_legacy_context',source_package='HDO.db 1.0.0',
 annotation_metadata_date=source_date,ontology_release=do$data_version,
 original_row_provenance='Alliance-derived and older DO2EG inputs merged by package build; per-row origin not retained',
 original_qualifiers='not_retained_in_compiled_map',original_evidence_codes='not_retained_in_compiled_map',
 original_references='not_retained_in_compiled_map',causal_or_direct_disease_evidence_established=FALSE)
accepted<-a |> filter(annotation_status=='accepted_compiled_human_gene_disease_context') |>
 left_join(do$terms |> select(term_id,name),by='term_id',relationship='many-to-one')
stopifnot(nrow(a)==nrow(annotations),!anyNA(a$annotation_status),all(accepted$entrez_id %in% human_gene_ids))
saveRDS(do,'data/curated/DO_ontology_snapshot.rds')
write_csv(do$terms,'data/curated/DO_term_registry.csv',na='')
write_csv(do$xrefs,'data/curated/DO_term_crossreferences.csv')
write_csv(do$edges,'data/curated/DO_all_source_relations.csv')
write_csv(safe_edges,'data/curated/DO_safe_annotation_relations.csv')
write_csv(a,'data/curated/DO_compiled_annotation_status_evidence.csv',na='')
write_csv(accepted,'data/curated/DO_compiled_human_context_annotations.csv',na='')
write_csv(metadata,'data/raw/do/HDO_database_metadata.csv')
write_csv(map_metadata,'data/raw/do/HDO_map_source_metadata.csv')
sources<-sources |> mutate(source_role=c('current_DO_hierarchy','compiled_gene_annotation_snapshot'),
 source_version=c(do$data_version,paste0('HDO.db 1.0.0; annotation metadata date ',source_date)))
write_csv(sources,'data/raw/do/acquisition_manifest.csv',na='')
write_csv(sources |> select(-local_file),'outputs/22b_source_summary.csv',na='')
write_csv(a |> count(annotation_status,name='source_annotation_rows'),'outputs/22b_annotation_coverage.csv')
write_csv(tibble(ontology_release=do$data_version,terms=nrow(do$terms),active_terms=length(active),safe_is_a_edges=nrow(safe_edges),
 compiled_source_rows=nrow(a),accepted_context_rows=nrow(accepted),human_gene_ids=n_distinct(accepted$entrez_id),
 annotated_disease_terms=n_distinct(accepted$term_id),compiled_annotation_metadata_date=source_date),
 'outputs/22b_ontology_summary.csv')
policy<-tibble(component=c('current_Alliance_human_gene_annotations','HDO_compiled_gene_term_map','positive_causal_training_labels'),
 status=c('Earlier official endpoint GET timed out with zero bytes; current raw annotations not acquired',
 'Frozen HDO.db 1.0.0 snapshot used as legacy context; metadata date 20240723; mixed historical row origins',
 'Not established from compiled map: qualifiers/evidence/references unavailable'),
 workflow_use=c('Coverage limitation remains open; do not imply current raw source acquisition',
 'Source-aware disease/function context only; exact human Entrez and current DO identifier checks',
 'Do not use these rows as independently verified causal positives or evidence of compound binding'))
write_csv(policy,'config/DO_annotation_source_policy.csv');write_csv(policy,'outputs/22b_source_limitations_summary.csv')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
m<-r('config/database_manifest.csv');k<-m$source_id=='do'
m$status[k]<-'DO_hierarchy_and_compiled_legacy_human_gene_context_acquired_current_Alliance_raw_unavailable'
m$release_version[k]<-paste(do$data_version,'; HDO.db1.0.0 annotation metadata 20240723')
m$retrieved_date[k]<-'2026-10-20';m$download_url[k]<-'https://bioconductor.org/packages/3.23/data/annotation/html/HDO.db.html'
m$local_file[k]<-'data/raw/do/acquisition_manifest.csv';write_csv(m,'config/database_manifest.csv',na='')
report<-c('SCRIPT 22b - DISEASE ONTOLOGY AND VERSIONED HUMAN GENE CONTEXT',
 paste('Ontology release:',do$data_version),paste('Ontology terms:',nrow(do$terms)),
 paste('Compiled gene/disease source rows:',nrow(a)),paste('Accepted human context rows:',nrow(accepted)),
 paste('Annotation metadata date:',source_date),
 'Current raw Alliance annotation download timed out; current raw source acquisition remains incomplete.',
 'HDO.db 1.0.0 build merges Alliance-derived and older DO2EG mappings without per-row source, qualifier, evidence or reference fields.',
 'Compiled rows are legacy association context, not newly adjudicated causal/positive disease evidence. Source limitations are explicit.',
 'Human Entrez keys checked against versioned org.Hs.eg.db; current DO primary/alternative identifiers checked, obsolete/unknown terms held.',
 'All source ontology relations retained; active is_a hierarchy verified acyclic for later propagation.',
 'No enrichment/BH tests or drug-binding/allostery claims. Outputs: source-rich data locally and compact coverage/source limitations in outputs/22b_*.')
writeLines(report,'outputs/22b_acquisition_summary.txt')
writeLines('PASS: source hashes, package SQLite read-only extraction, mixed-provenance build audit, ontology DAG, all-row status and human/current-term checks.',
 'outputs/22b_acquisition_validation.txt')
write_csv(tibble(package='org.Hs.eg.db',version=as.character(packageVersion('org.Hs.eg.db'))),
 'outputs/22b_mapping_package_summary.csv')
cat(paste0('\nSCRIPT 22b RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/22b_acquire_DO_annotations.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
