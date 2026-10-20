# GO hierarchy and human annotation evidence; qualifiers and proteoforms remain explicit.
library(readr);library(dplyr);library(tibble);library(igraph)
source('scripts/helpers/ontology_source_get.R');source('scripts/helpers/read_ontology_obo.R')
sources<-bind_rows(
 ontology_source_get('https://current.geneontology.org/ontology/go-basic.obo','data/raw/go/go-basic.obo'),
 ontology_source_get('https://current.geneontology.org/annotations/gaf/HUMAN-uniprot.gaf.gz','data/raw/go/HUMAN-uniprot.gaf.gz'))
go<-read_ontology_obo(sources$local_file[1],'GO')
stopifnot(nrow(go$terms)>30000L)
active<-go$terms$term_id[!go$terms$is_obsolete]
# Only is_a/part_of are the declared annotation-propagation relationships; retain others as source edges.
safe_edges<-go$edges |> filter(relation %in% c('is_a','part_of'),child %in% active,parent %in% active)
g<-graph_from_data_frame(safe_edges |> select(child,parent),directed=TRUE,vertices=tibble(name=active))
stopifnot(is_dag(g))
gaf_path<-sources$local_file[2];con<-gzfile(gaf_path,'rt');header<-character()
repeat{line<-readLines(con,n=1L,warn=FALSE);if(!length(line)||!startsWith(line,'!'))break;header<-c(header,line)};close(con)
gaf_version<-sub('^!go-version: ','',grep('^!go-version: ',header,value=TRUE))
stopifnot(length(gaf_version)==1L,grepl(sub('^releases/','',go$data_version),gaf_version,fixed=TRUE))
columns<-c('db','object_id','symbol','qualifier','source_term_id','references','evidence_code','with_from',
 'aspect','object_name','synonyms','object_type','taxon','annotation_date','assigned_by','annotation_extension','gene_product_form_id')
a<-read_tsv(gaf_path,comment='!',col_names=columns,col_types=cols(.default=col_character()),quote='',progress=FALSE)
stopifnot(nrow(problems(a))==0L,ncol(a)==17L)
a<-a |> mutate(source_row_id=row_number(),primary_taxon=sub('[|].*$','',taxon),
 negative_NOT=grepl('(^|[|])NOT([|]|$)',coalesce(qualifier,'')),
 no_data_evidence=evidence_code=='ND',
 annotation_uniprot=case_when(is.na(gene_product_form_id)|gene_product_form_id==''~object_id,
 startsWith(gene_product_form_id,'UniProtKB:')~sub('^UniProtKB:','',gene_product_form_id),TRUE~NA_character_)) |>
 left_join(resolve_ontology_ids(unique(a$source_term_id),go),by='source_term_id',relationship='many-to-one') |>
 mutate(annotation_status=case_when(primary_taxon!='taxon:9606'~'nonhuman_primary_taxon_hold',
 db!='UniProtKB'~'non_UniProt_source_hold',negative_NOT~'negative_NOT_not_positive_feature',
 no_data_evidence~'ND_no_positive_annotation',!term_status %in% c('active_primary_ID','active_documented_alternative_ID')~term_status,
 is.na(annotation_uniprot)~'unresolved_gene_product_form_hold',TRUE~'accepted_positive_human_annotation'),
 evidence_layer=case_when(evidence_code=='IEA'~'electronic',
 evidence_code %in% c('EXP','IDA','IPI','IMP','IGI','IEP','HTP','HDA','HMP','HGI','HEP')~'experimental_annotation',
 TRUE~'other_evidence_code'),ontology_release=go$data_version)
stopifnot(nrow(a)>100000L,!anyNA(a$annotation_status),!anyDuplicated(a$source_row_id))
accepted<-a |> filter(annotation_status=='accepted_positive_human_annotation') |>
 left_join(go$terms |> select(term_id,name,namespace),by='term_id',relationship='many-to-one')
aspects<-c(F='molecular_function',P='biological_process',C='cellular_component')
stopifnot(all(accepted$namespace==unname(aspects[accepted$aspect])),!any(accepted$negative_NOT),
 all(accepted$primary_taxon=='taxon:9606'))
saveRDS(go,'data/curated/GO_ontology_snapshot.rds')
write_csv(go$terms,'data/curated/GO_term_registry.csv',na='')
write_csv(go$edges,'data/curated/GO_all_source_relations.csv')
write_csv(safe_edges,'data/curated/GO_safe_annotation_relations.csv')
write_csv(a,'data/curated/GO_human_annotation_status_evidence.csv',na='')
write_csv(accepted,'data/curated/GO_human_positive_annotation_evidence.csv',na='')
writeLines(header,'data/raw/go/HUMAN_annotation_header.txt')
sources<-sources |> mutate(ontology_release=go$data_version)
write_csv(sources,'data/raw/go/acquisition_manifest.csv',na='')
write_csv(sources |> select(-local_file),'outputs/22a_source_summary.csv',na='')
write_csv(a |> count(annotation_status,name='source_annotation_rows'),'outputs/22a_annotation_coverage.csv')
write_csv(accepted |> group_by(namespace,evidence_layer) |> summarise(annotation_rows=n(),
 distinct_protein_forms=n_distinct(annotation_uniprot),distinct_terms=n_distinct(term_id),.groups='drop'),
 'outputs/22a_evidence_summary.csv')
write_csv(tibble(ontology_release=go$data_version,terms=nrow(go$terms),active_terms=length(active),
 all_source_relations=nrow(go$edges),safe_is_a_part_of_relations=nrow(safe_edges),
 source_annotation_rows=nrow(a),accepted_positive_rows=nrow(accepted)), 'outputs/22a_ontology_summary.csv')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
m<-r('config/database_manifest.csv');k<-m$source_id=='go'
m$status[k]<-'ontology_and_human_GAF_acquired_qualifiers_evidence_forms_preserved';m$release_version[k]<-go$data_version
m$retrieved_date[k]<-'2026-10-20';m$download_url[k]<-'https://geneontology.org/docs/download-go-annotations/downloads/'
m$local_file[k]<-'data/raw/go/acquisition_manifest.csv';write_csv(m,'config/database_manifest.csv',na='')
report<-c('SCRIPT 22a - GO ONTOLOGY AND HUMAN GAF ANNOTATIONS',paste('Ontology release:',go$data_version),
 paste('Ontology terms:',nrow(go$terms)),paste('Human GAF source rows:',nrow(a)),paste('Accepted positive human annotation rows:',nrow(accepted)),
 'GO ontology and GAF headers agree on ontology release. All source rows status-accounted.',
 'NOT qualifiers and ND evidence do not become positive features; obsolete/unknown terms and unresolved protein forms are held.',
 'Exact protein-form identifiers, evidence codes, references, assignment dates and providers retained. Electronic and other evidence layers remain separate.',
 'All source relations retained; only active is_a/part_of edges are declared safe for later annotation propagation and verified acyclic.',
 'No enrichment/BH tests, semantic-similarity ranking or binding inference performed at source acquisition.',
 'Outputs: full ontology/annotation evidence locally; source, status, aspect/evidence coverage and validation in outputs/22a_*.')
writeLines(report,'outputs/22a_acquisition_summary.txt')
writeLines('PASS: source checksums, ontology DAG, matching release headers, 17-column GAF parsing, all-row statuses, human/qualifier/aspect checks.',
 'outputs/22a_acquisition_validation.txt')
cat(paste0('\nSCRIPT 22a RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/22a_acquire_GO_annotations.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
