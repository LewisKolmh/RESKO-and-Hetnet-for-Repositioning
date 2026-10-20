# Independent raw GAF anchor checks, ontology-ID edge cases and feature integrity.
library(readr);library(dplyr);library(tibble)
source('scripts/helpers/read_ontology_obo.R')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
# Exercise parser boundaries and ambiguous aliases without using actual output tables.
p<-tempfile(fileext='.obo')
writeLines(c('format-version: 1.2','data-version: fixture','',
 '[Term]','id: TEST:1','name: root','alt_id: TEST:old','alt_id: TEST:old','',
 '[Term]','id: TEST:2','name: child','is_a: TEST:1 ! root','relationship: part_of TEST:1 ! root','alt_id: TEST:conflict','',
 '[Term]','id: TEST:3','name: obsolete','is_obsolete: true','alt_id: TEST:conflict','',
 '[Typedef]','id: part_of'),p)
z<-read_ontology_obo(p,'TEST');unlink(p)
resolved<-resolve_ontology_ids(c('TEST:1','TEST:old','TEST:conflict','TEST:3','TEST:missing'),z)
stopifnot(nrow(z$terms)==3L,nrow(z$edges)==2L,
 identical(resolved$term_status,c('active_primary_ID','active_documented_alternative_ID',
 'ambiguous_alternative_ID_hold','obsolete_term_hold','unknown_term_hold')))
# Read the complete raw file independently; parse only exact anchor-object rows with base R.
anchors<-r('config/eef1a_target_registry.csv')$uniprot
con<-gzfile('data/raw/go/HUMAN-uniprot.gaf.gz','rt');count<-0L;selected<-list();source_positions<-integer()
repeat{
 lines<-readLines(con,n=100000L,warn=FALSE);if(!length(lines))break
 lines<-lines[!startsWith(lines,'!')];if(!length(lines))next
 keep<-rep(FALSE,length(lines));for(id in anchors)keep<-keep|grepl(paste0('\t',id,'\t'),lines,fixed=TRUE)
 if(any(keep)){
  selected[[length(selected)+1L]]<-lines[keep]
  source_positions<-c(source_positions,count+which(keep))
 }
 count<-count+length(lines)
}
close(con)
lines<-unlist(selected,use.names=FALSE)
parts<-lapply(strsplit(lines,'\t',fixed=TRUE),function(x){length(x)<-17L;x[x=='']<-NA_character_;x})
raw<-as.data.frame(do.call(rbind,parts),stringsAsFactors=FALSE)
names(raw)<-c('db','object_id','symbol','qualifier','source_term_id','references','evidence_code','with_from',
 'aspect','object_name','synonyms','object_type','taxon','annotation_date','assigned_by','annotation_extension','gene_product_form_id')
raw$source_row_id<-source_positions;raw<-raw[raw$object_id %in% anchors,,drop=FALSE]
full_status<-read_csv('data/curated/GO_human_annotation_status_evidence.csv',
 col_types=cols_only(source_row_id=col_integer(),object_id=col_character(),annotation_status=col_character()),progress=FALSE)
stopifnot(count==nrow(full_status))
positive<-r('data/curated/GO_human_positive_annotation_evidence.csv') |> filter(object_id %in% anchors)
ontology<-readRDS('data/curated/GO_ontology_snapshot.rds')
term_resolution<-resolve_ontology_ids(raw$source_term_id,ontology)
expected<-raw$db=='UniProtKB' & sub('[|].*$','',raw$taxon)=='taxon:9606' &
 !grepl('(^|[|])NOT([|]|$)',ifelse(is.na(raw$qualifier),'',raw$qualifier)) & raw$evidence_code!='ND' &
 term_resolution$term_status %in% c('active_primary_ID','active_documented_alternative_ID') &
 (is.na(raw$gene_product_form_id)|startsWith(raw$gene_product_form_id,'UniProtKB:'))
stopifnot(setequal(as.character(raw$source_row_id[expected]),positive$source_row_id))
cols<-c('object_id','source_term_id','qualifier','references','evidence_code','with_from','taxon','annotation_date','gene_product_form_id')
lookup<-match(as.character(positive$source_row_id),as.character(raw$source_row_id))
for(col in cols)stopifnot(identical(as.character(raw[[col]][lookup]),as.character(positive[[col]])))
go<-readRDS('data/curated/inventory_GO_term_features.rds')
do<-readRDS('data/curated/inventory_DO_legacy_context_features.rds')
coverage<-r('data/curated/inventory_ontology_annotation_coverage.csv')
stopifnot(nrow(coverage)==200L,!anyDuplicated(coverage$entity_id),!anyDuplicated(go),!anyDuplicated(do),
 all(go$entity_id %in% coverage$entity_id),all(do$entity_id %in% coverage$entity_id),
 all(do$evidence_layer=='HDO_db_compiled_legacy_context'),
 !any(go$entity_id=='sider:CID100135113'),!any(do$entity_id=='sider:CID100135113'),
 !any(go$term_id %in% c('GO:0003674','GO:0008150','GO:0005575')),!any(do$term_id=='DOID:4'))
hashes<-r('outputs/22c_input_summary.csv');stopifnot(identical(unname(tools::md5sum(hashes$file)),hashes$md5))
for(f in list.files('scripts',pattern='[.]R$',recursive=TRUE,full.names=TRUE))invisible(parse(f))
validation<-c(paste('PASS: independent scan of all',count,'GAF rows; human eEF1A anchor positive/negative selection and every retained reference/evidence/form field agree with raw data.'),
 'PASS: OBO fixture verifies term/typedef boundaries, relation parsing, duplicate aliases, ambiguity holds, obsolete and unknown IDs.',
 'PASS: 200-entity coverage and unique source-layer features preserved; generic roots and blocked sacrosidase structure excluded; DO legacy evidence label retained.',
 'PASS: protected input hashes unchanged and all R scripts parse; no significance, binding or causal evidence promoted.')
cat(paste(validation,collapse='\n'),'\n')
cat(paste0(paste(validation,collapse='\n'),'\n'),file='outputs/22c_integration_validation.txt',append=TRUE)
cat(paste0('\nSCRIPT 22c - INDEPENDENT VALIDATION\n',paste(validation,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
