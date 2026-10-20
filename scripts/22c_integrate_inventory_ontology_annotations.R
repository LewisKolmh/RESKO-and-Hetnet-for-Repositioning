# Join GO protein annotations and separately labelled DO gene-context projections to inventory nodes.
library(readr);library(dplyr);library(tidyr);library(tibble)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
inputs<-c('outputs/16f_target_evidence_handoff_summary.csv','data/curated/inventory_current_STRING_protein_links.csv',
 'data/curated/STRING_current_direct_UniProt_aliases.csv','data/raw/string/v12.5/9606.protein.aliases.v12.5.txt.gz',
 'data/curated/GO_human_positive_annotation_evidence.csv','data/curated/DO_compiled_human_context_annotations.csv',
 'data/curated/GO_term_registry.csv','data/curated/DO_term_registry.csv','config/cancers.csv',
 'config/eef1a_target_registry.csv','data/curated/DO_term_crossreferences.csv')
before<-tools::md5sum(inputs)
inventory<-r(inputs[1]);links<-r(inputs[2]) |>
 distinct(entity_id,current_STRING_id,source_layer,evidence_stratum)
current_alias<-r(inputs[3]);go<-r(inputs[5]);do<-r(inputs[6]);go_terms<-r(inputs[7]);do_terms<-r(inputs[8])
# GO requires an exact protein/form accession to one current node, never a base-ID or gene-symbol fallback.
counts<-current_alias |> count(uniprot,name='current_node_candidates')
go_map<-tibble(annotation_uniprot=sort(unique(go$annotation_uniprot))) |>
 left_join(counts,by=c('annotation_uniprot'='uniprot'),relationship='one-to-one') |>
 mutate(current_node_candidates=coalesce(current_node_candidates,0L),mapping_status=case_when(
 current_node_candidates==0L~'no_exact_current_protein_form_match',current_node_candidates>1L~'ambiguous_current_protein_nodes_hold',
 TRUE~'unique_exact_current_protein_form_match'))
go_accept<-current_alias |> inner_join(go_map |> filter(current_node_candidates==1L),
 by=c('uniprot'='annotation_uniprot'),relationship='many-to-one') |>
 select(annotation_uniprot=uniprot,current_STRING_id=protein,current_alias_sources=alias_sources)
stopifnot(!anyDuplicated(go_accept$annotation_uniprot))
go_node_evidence<-go |> inner_join(go_accept,by='annotation_uniprot',relationship='many-to-one')
# Retain qualifiers, refs, forms and evidence in provenance; features are presence-only, not independent experiments.
go_node_features<-go_node_evidence |>
 filter(!term_id %in% c('GO:0003674','GO:0008150','GO:0005575')) |>
 distinct(current_STRING_id,term_id,namespace,evidence_layer,annotation_relation=qualifier)
# DO compilation uses Entrez genes. Current direct HGNC-supplied Entrez aliases avoid symbol matching.
a<-read_tsv(inputs[4],col_types=cols(.default=col_character()),quote='',progress=FALSE)
stopifnot(nrow(problems(a))==0L);names(a)<-c('current_STRING_id','entrez_id','alias_sources')
gene_candidates<-a |> filter(alias_sources=='Ensembl_HGNC_entrez_id',grepl('^[0-9]+$',entrez_id)) |>
 distinct(current_STRING_id,entrez_id,alias_sources)
gene_counts<-gene_candidates |> count(current_STRING_id,name='entrez_candidates')
node_universe<-r('data/curated/STRING_human_protein_registry_v12.5.csv') |> select(current_STRING_id=protein)
gene_coverage<-node_universe |> left_join(gene_counts,by='current_STRING_id',relationship='one-to-one') |>
 mutate(entrez_candidates=coalesce(entrez_candidates,0L),mapping_status=case_when(entrez_candidates==0L~'no_direct_Entrez_alias',
 entrez_candidates>1L~'ambiguous_Entrez_gene_hold',TRUE~'unique_direct_Entrez_gene_alias'))
accepted_genes<-gene_candidates |> inner_join(gene_coverage |> filter(entrez_candidates==1L),
 by='current_STRING_id',relationship='many-to-one')
stopifnot(!anyDuplicated(accepted_genes$current_STRING_id))
do_node_evidence<-do |> inner_join(accepted_genes,by='entrez_id',relationship='many-to-many') |>
 mutate(projection_scope='Gene-level compiled association projected to protein node; not isoform-specific causal evidence')
do_node_features<-do_node_evidence |> filter(term_id!='DOID:4') |>
 distinct(current_STRING_id,term_id,evidence_layer)
go_entity_features<-links |> inner_join(go_node_features,by='current_STRING_id',relationship='many-to-many') |>
 distinct(entity_id,source_layer,evidence_stratum,term_id,namespace,evidence_layer,annotation_relation)
do_entity_features<-links |> inner_join(do_node_features,by='current_STRING_id',relationship='many-to-many') |>
 distinct(entity_id,source_layer,evidence_stratum,term_id,evidence_layer)
# Direct/current GO and compiled legacy DO remain separate matrices/evidence layers.
saveRDS(go_entity_features,'data/curated/inventory_GO_term_features.rds',compress='gzip')
saveRDS(do_entity_features,'data/curated/inventory_DO_legacy_context_features.rds',compress='gzip')
write_csv(go_map,'data/curated/GO_current_protein_form_mapping_coverage.csv')
write_csv(gene_coverage,'data/curated/STRING_Entrez_gene_mapping_coverage.csv')
write_csv(accepted_genes,'data/curated/STRING_direct_Entrez_gene_links.csv')
write_csv(go_node_features,'data/curated/STRING_GO_protein_term_features.csv')
write_csv(do_node_features,'data/curated/STRING_DO_legacy_gene_context_features.csv')
saveRDS(go_node_evidence,'data/curated/STRING_GO_source_annotation_evidence.rds',compress='gzip')
saveRDS(do_node_evidence,'data/curated/STRING_DO_legacy_source_context_evidence.rds',compress='gzip')
# Disease names select exact active DO concepts; crossrefs do not constitute full UMLS acquisition.
cancers<-r(inputs[9]) |> filter(include=='TRUE')
cancer_map<-cancers |> mutate(lookup=tolower(cancer_name)) |>
 left_join(do_terms |> filter(is_obsolete=='FALSE') |> transmute(lookup=tolower(name),doid=term_id,DO_name=name),
 by='lookup',relationship='one-to-one') |> select(-lookup)
stopifnot(nrow(cancer_map)==3L,!anyNA(cancer_map$doid))
xrefs<-r(inputs[11])
cancer_xrefs<-cancer_map |> left_join(xrefs,by=c('doid'='term_id'),relationship='one-to-many') |>
 mutate(crossref_scope='Ontology cross-reference only; not proof of licensed UMLS acquisition or disease-treatment equivalence')
write_csv(cancer_map,'config/cancer_DO_concepts.csv')
write_csv(cancer_xrefs,'data/curated/cancer_DO_source_crossreferences.csv',na='')
write_csv(cancer_map,'outputs/22c_cancer_concept_summary.csv')
# Term count coverage keeps all original entities and source/role strata.
covered<-function(features,label)features |> distinct(entity_id,term_id) |> count(entity_id,name=label)
coverage<-inventory |> select(entity_id,origin,source_name,seed_tiers,cancers) |>
 left_join(covered(go_entity_features,'GO_terms'),by='entity_id',relationship='one-to-one') |>
 left_join(covered(do_entity_features,'DO_legacy_context_terms'),by='entity_id',relationship='one-to-one') |>
 mutate(across(c(GO_terms,DO_legacy_context_terms),~coalesce(.x,0L)),
 GO_coverage=if_else(GO_terms>0L,'annotations_available','no_accepted_annotation_links'),
 DO_coverage=if_else(DO_legacy_context_terms>0L,'compiled_legacy_context_available','no_accepted_context_links'),
 missing_coverage_is_not_negative=TRUE,affinity_required=FALSE)
anchors<-r(inputs[10]) |> select(gene,uniprot) |>
 inner_join(current_alias,by='uniprot',relationship='one-to-one') |> select(gene,uniprot,current_STRING_id=protein)
anchor_summary<-anchors |>
 left_join(go_node_features |> group_by(current_STRING_id) |> summarise(GO_terms=n_distinct(term_id),.groups='drop'),
 by='current_STRING_id',relationship='one-to-one') |>
 left_join(do_node_features |> group_by(current_STRING_id) |> summarise(DO_legacy_context_terms=n_distinct(term_id),.groups='drop'),
 by='current_STRING_id',relationship='one-to-one') |>
 mutate(across(c(GO_terms,DO_legacy_context_terms),~coalesce(.x,0L)))
stopifnot(nrow(coverage)==200L,!anyDuplicated(coverage$entity_id),nrow(anchors)==2L,
 all(go_node_evidence$annotation_status=='accepted_positive_human_annotation'),
 !any(go_node_evidence$negative_NOT=='TRUE'),all(go_node_evidence$primary_taxon=='taxon:9606'),
 all(go_entity_features$entity_id %in% inventory$entity_id),all(do_entity_features$entity_id %in% inventory$entity_id),
 !any(go_entity_features$entity_id=='sider:CID100135113'),
 all(go_entity_features$term_id %in% go_terms$term_id[go_terms$is_obsolete=='FALSE']),
 all(do_entity_features$term_id %in% do_terms$term_id[do_terms$is_obsolete=='FALSE']),
 identical(before,tools::md5sum(inputs)))
write_csv(coverage,'data/curated/inventory_ontology_annotation_coverage.csv')
write_csv(coverage |> count(origin,GO_coverage,DO_coverage,name='inventory_entities'),'outputs/22c_entity_coverage.csv')
write_csv(go_map |> count(mapping_status,name='GO_protein_forms'),'outputs/22c_GO_mapping_coverage.csv')
write_csv(gene_coverage |> count(mapping_status,name='STRING_registry_nodes'),'outputs/22c_gene_mapping_coverage.csv')
write_csv(anchor_summary,'outputs/22c_anchor_annotation_summary.csv')
layer_summary<-bind_rows(
 go_entity_features |> group_by(source_layer,evidence_stratum) |> summarise(ontology_layer='GO_positive_annotations',
 inventory_entities=n_distinct(entity_id),distinct_terms=n_distinct(term_id),feature_rows=n(),.groups='drop'),
 do_entity_features |> group_by(source_layer,evidence_stratum) |> summarise(ontology_layer='DO_compiled_legacy_context',
 inventory_entities=n_distinct(entity_id),distinct_terms=n_distinct(term_id),feature_rows=n(),.groups='drop'))
write_csv(layer_summary,'outputs/22c_source_layer_summary.csv')
write_csv(tibble(file=inputs,md5=unname(before)),'outputs/22c_input_summary.csv')
checklist<-r('config/RESKO_method_completion_checklist.csv')
checklist$project_status[checklist$step=='08']<-'GO_ontology_human_evidence_and_inventory_term_features_available_Wang_pending'
checklist$next_work[checklist$step=='08']<-'Declare ontology propagation/semantic similarity and matrix design; retain GO qualifier/evidence/form layers'
checklist$project_status[checklist$step=='09']<-'current_DO_hierarchy_and_legacy_gene_context_integrated_current_raw_source_unavailable'
checklist$next_work[checklist$step=='09']<-'Retain source-date/qualifier limits; define human annotation universe and valid enrichment family before hypergeometric/BH analysis'
write_csv(checklist,'config/RESKO_method_completion_checklist.csv')
report<-c('SCRIPT 22c - HUMAN INVENTORY GO/DO ANNOTATION INTEGRATION',
 paste('Inventory entries retained:',nrow(coverage)),paste('Entries with GO term features:',sum(coverage$GO_terms>0L)),
 paste('Entries with DO compiled legacy context:',sum(coverage$DO_legacy_context_terms>0L)),
 paste('Current protein nodes with GO features:',n_distinct(go_node_features$current_STRING_id)),
 paste('Current protein nodes with DO legacy context:',n_distinct(do_node_features$current_STRING_id)),
 'GO uses exact current protein-form accession aliases with unique-node mapping; NOT/ND and unknown/obsolete terms stay excluded from positive features.',
 'Positive GO annotation relation qualifiers remain explicit in node/entity feature tables; differently qualified relationships are not silently equated.',
 'DO uses unique direct Entrez gene aliases; gene-level compiled contexts are not isoform-specific causal evidence and remain separate from GO evidence.',
 'Source protein roles, STITCH chemical-identity strata, evidence types, seed tiers and cancers retained. All 200 original entries remain in coverage.',
 'Three cancer concepts mapped by exact active DO names; ontology crossrefs do not constitute full UMLS acquisition.',
 'Feature tables are deduplicated presence-only associations, not independent experiment counts. Generic ontology roots excluded from features.',
 'Hierarchy DAGs are available for later propagation; no new ancestor-expanded observations or semantic/enrichment ranking claimed here.',
 'Current Alliance raw annotations remain unavailable; DO legacy snapshot provenance limitations explicitly retained.',
 'No drug-binding/allostery inference, no affinity filter, no p-values/BH adjustment at this annotation stage.',
 'Outputs: source-rich node evidence and entity feature RDS locally; coverage, anchor, cancer concept and source-layer reports in outputs/22c_*.',
 'Next: paper-aligned feature/ontology-similarity matrices and the remaining chemistry/SIDER-baseline prerequisites; define enrichment background before testing.')
writeLines(report,'outputs/22c_integration_summary.txt')
writeLines('PASS: exact/unique protein-form and Entrez mapping, current active ontology IDs, source/identity strata and all 200 entries preserved; protected input hashes unchanged.',
 'outputs/22c_integration_validation.txt')
notes<-readLines('Notes.txt',warn=FALSE);start<-grep('^CURRENT CHECKPOINT',notes)[1];end<-grep('^PROJECT AIM',notes)[1]
stopifnot(start<end);writeLines(c(notes[seq_len(start-1)],'CURRENT CHECKPOINT - GO AND VERSIONED DO CONTEXT ANNOTATION STAGES 22a/22b/22c COMPLETE',report,'',notes[end:length(notes)]),'Notes.txt')
cat(paste0('\nSCRIPT 22c RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/22c_integrate_inventory_ontology_annotations.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
