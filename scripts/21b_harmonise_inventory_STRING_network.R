# Core RESKO: explicit protein bridges and a source-aware human network layer.
library(readr);library(dplyr);library(tidyr);library(tibble)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
inputs<-c('outputs/16f_target_evidence_handoff_summary.csv',
 'data/curated/inventory_DrugBank_human_protein_links.csv',
 'data/curated/inventory_broad_STITCH_association_evidence.csv',
 'data/curated/STITCH_historical_UniProt_alias_review.csv',
 'data/curated/STRING_human_protein_registry_v12.5.csv',
 'data/raw/string/v12.5/9606.protein.aliases.v12.5.txt.gz',
 'data/curated/STRING_human_functional_network_v12.5.rds','config/eef1a_target_registry.csv')
before<-tools::md5sum(inputs)
inventory<-r(inputs[1]);db<-r(inputs[2]);stitch<-r(inputs[3]);legacy_alias<-r(inputs[4]);info<-r(inputs[5])
targets<-r(inputs[8]);stopifnot(nrow(inventory)==200L,!anyDuplicated(inventory$entity_id))
settings_path<-'config/STRING_network_projection_settings.csv'
if(!file.exists(settings_path))write_csv(tibble(setting=c('primary_combined_score','sensitivity_combined_scores','expansion_scope',
 'current_alias_sources','legacy_bridge_policy','affinity_required'),
 value=c('700','400|900','Induced network over accepted inventory protein nodes plus two human eEF1A anchors; no neighbour expansion',
 'UniProt_AC|Ensembl_UniProt|Ensembl_HGNC_uniprot_ids',
 'All direct historical UniProt aliases bridged to current aliases; accept only one distinct current node per legacy node; no gene-name or same-ENSP fallback',
 'FALSE'),basis='Declared project adaptation; retain full published network and explicit sensitivity layers, not a reproduction claim'),settings_path)
settings<-r(settings_path);get_setting<-function(k)settings$value[match(k,settings$setting)]
primary<-as.integer(get_setting('primary_combined_score'));sensitivity<-as.integer(strsplit(get_setting('sensitivity_combined_scores'),'|',fixed=TRUE)[[1]])
stopifnot(length(primary)==1L,!is.na(primary),primary>0L,primary<=1000L,
 all(sensitivity>0L&sensitivity<=1000L),get_setting('affinity_required')=='FALSE')
aliases<-read_tsv(inputs[6],col_types=cols(.default=col_character()),quote='',progress=FALSE)
stopifnot(nrow(problems(aliases))==0L,ncol(aliases)==3L);names(aliases)<-c('protein','uniprot','alias_sources')
uniprot_pattern<-'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$'
allowed_sources<-strsplit(get_setting('current_alias_sources'),'|',fixed=TRUE)[[1]]
stopifnot(all(allowed_sources %in% c('UniProt_AC','Ensembl_UniProt','Ensembl_HGNC_uniprot_ids')))
source_pattern<-paste0('(^| )(',paste(allowed_sources,collapse='|'),')( |$)')
current_alias<-aliases |> filter(grepl(uniprot_pattern,uniprot),grepl(source_pattern,alias_sources)) |>
 group_by(protein,uniprot) |> summarise(alias_sources=paste(sort(unique(alias_sources)),collapse='|'),.groups='drop')
stopifnot(all(current_alias$protein %in% info$protein))
write_csv(current_alias,'data/curated/STRING_current_direct_UniProt_aliases.csv')
accession_counts<-current_alias |> count(uniprot,name='current_node_candidates')
accession_coverage<-tibble(uniprot=sort(unique(c(db$uniprot,targets$uniprot)))) |>
 left_join(accession_counts,by='uniprot',relationship='one-to-one') |>
 mutate(current_node_candidates=coalesce(current_node_candidates,0L),mapping_status=case_when(
 current_node_candidates==0L~'no_direct_current_accession_match',current_node_candidates>1L~'ambiguous_current_node_hold',
 TRUE~'unique_exact_current_accession_alias'))
accepted_accessions<-current_alias |> inner_join(accession_coverage |> filter(current_node_candidates==1L),
 by='uniprot',relationship='many-to-one')
stopifnot(!anyDuplicated(accepted_accessions$uniprot))
db_review<-db |> left_join(accession_coverage,by='uniprot',relationship='many-to-one') |>
 left_join(accepted_accessions |> select(uniprot,current_STRING_id=protein,current_alias_sources=alias_sources),
 by='uniprot',relationship='many-to-one')
db_current<-db_review |> filter(mapping_status=='unique_exact_current_accession_alias')
# Do not equate identical Ensembl text across releases with sequence/isoform identity.
legacy_candidates<-legacy_alias |> select(legacy_protein=protein,uniprot,historical_alias_sources=alias_sources) |>
 inner_join(current_alias |> rename(current_STRING_id=protein,current_alias_sources=alias_sources),
 by='uniprot',relationship='many-to-many') |> distinct()
legacy_counts<-legacy_candidates |> group_by(legacy_protein) |> summarise(
 bridge_accessions=n_distinct(uniprot),current_node_candidates=n_distinct(current_STRING_id),.groups='drop')
legacy_coverage<-tibble(legacy_protein=sort(unique(stitch$protein))) |>
 left_join(legacy_alias |> count(protein,name='historical_direct_alias_rows'),by=c('legacy_protein'='protein'),relationship='one-to-one') |>
 left_join(legacy_counts,by='legacy_protein',relationship='one-to-one') |>
 mutate(across(c(historical_direct_alias_rows,bridge_accessions,current_node_candidates),~coalesce(.x,0L)),
 mapping_status=case_when(historical_direct_alias_rows==0L~'no_direct_historical_UniProt_alias',
 current_node_candidates==0L~'historical_alias_no_current_node_match',
 current_node_candidates>1L~'ambiguous_legacy_to_current_nodes_hold',
 TRUE~'unique_historical_accession_bridge_projection'),sequence_equivalence_established=FALSE)
accepted_legacy<-legacy_candidates |> inner_join(legacy_coverage |> filter(current_node_candidates==1L),
 by='legacy_protein',relationship='many-to-one') |> distinct(legacy_protein,current_STRING_id)
stopifnot(!anyDuplicated(accepted_legacy$legacy_protein))
stitch_review<-stitch |> rename(legacy_protein=protein) |>
 left_join(legacy_coverage,by='legacy_protein',relationship='many-to-one') |>
 left_join(accepted_legacy,by='legacy_protein',relationship='many-to-one')
stitch_current<-stitch_review |> filter(mapping_status=='unique_historical_accession_bridge_projection')
anchor<-targets |> select(gene,uniprot) |>
 left_join(accepted_accessions |> select(uniprot,current_STRING_id=protein),by='uniprot',relationship='one-to-one') |>
 left_join(info |> select(current_STRING_id=protein,preferred_name),by='current_STRING_id',relationship='one-to-one')
stopifnot(nrow(anchor)==2L,!anyNA(anchor$current_STRING_id),all(anchor$preferred_name==anchor$gene),
 all(grepl('^9606[.]',anchor$current_STRING_id)))
# Separate source/role/compound identity layers; no cross-source confidence pooling.
node_links<-bind_rows(
 db_current |> transmute(entity_id,current_STRING_id,source_layer='DrugBank',
 evidence_stratum=paste0('DrugBank_role:',protein_role),source_protein_id=uniprot,
 mapping_basis='unique_exact_current_accession_alias'),
 stitch_current |> transmute(entity_id,current_STRING_id,source_layer='STITCH_5',
 evidence_stratum=compound_identity_layer,source_protein_id=legacy_protein,
 mapping_basis='unique_historical_accession_bridge_not_sequence_equivalence')) |> distinct()
stopifnot(!anyNA(node_links$current_STRING_id),all(node_links$current_STRING_id %in% info$protein),
 all(node_links$entity_id %in% inventory$entity_id),!any(node_links$entity_id=='sider:CID100135113'))
network<-readRDS(inputs[7]);stopifnot(network$source_release=='12.5',network$taxon==9606L)
nodes<-sort(unique(c(node_links$current_STRING_id,anchor$current_STRING_id)))
induced<-network$edges |> filter(protein1 %in% nodes,protein2 %in% nodes)
node_table<-info |> filter(protein %in% nodes) |>
 left_join(anchor |> select(protein=current_STRING_id,eef1a_anchor=gene),by='protein',relationship='one-to-one') |>
 mutate(is_eef1a_anchor=!is.na(eef1a_anchor))
stopifnot(nrow(node_table)==length(nodes))
# Keep the complete published-confidence induced network, including explicit isolated nodes.
saveRDS(list(source_release=network$source_release,taxon=9606L,
 network_type=network$network_type,scope=get_setting('expansion_scope'),nodes=node_table,
 edges=induced,primary_combined_score=primary,sensitivity_combined_scores=sensitivity,
 source_input_md5=before,identity_note='STITCH historical bridges remain labelled projections; sources are not independent experiments'),
 'data/curated/inventory_STRING_functional_induced_network_v12.5.rds',compress='gzip')
source_floor<-min(network$edges$combined_score)
cutoffs<-sort(unique(c(source_floor,primary,sensitivity)))
network_summary<-bind_rows(lapply(cutoffs,function(cutoff){
 e<-induced |> filter(combined_score>=cutoff);touched<-unique(c(e$protein1,e$protein2))
 tibble(combined_score_cutoff=cutoff,layer=if(cutoff==primary)'declared_primary'else if(cutoff==source_floor)'all_published_source_edges'else 'declared_sensitivity',
 inventory_and_anchor_nodes=length(nodes),induced_undirected_edges=nrow(e),nodes_with_edges=length(touched),
 nodes_retained_without_edges=sum(!nodes %in% touched))
}))
anchor_summary<-bind_rows(lapply(seq_len(nrow(anchor)),function(i){
 id<-anchor$current_STRING_id[i]
 full<-network$edges |> filter((protein1==id|protein2==id)&combined_score>=primary)
 local<-induced |> filter((protein1==id|protein2==id)&combined_score>=primary)
 tibble(gene=anchor$gene[i],uniprot=anchor$uniprot[i],current_STRING_id=id,combined_score_cutoff=primary,
 full_human_network_associations=nrow(full),induced_inventory_network_associations=nrow(local),
 interpretation='Functional protein association degree; not number of compounds binding eEF1A')
}))
entity_counts<-node_links |> group_by(entity_id) |> summarise(
 accepted_current_protein_nodes=n_distinct(current_STRING_id),
 DrugBank_current_nodes=n_distinct(current_STRING_id[source_layer=='DrugBank']),
 STITCH_projected_current_nodes=n_distinct(current_STRING_id[source_layer=='STITCH_5']),.groups='drop')
coverage<-inventory |> select(entity_id,origin,source_name,seed_tiers,cancers) |>
 left_join(entity_counts,by='entity_id',relationship='one-to-one') |>
 mutate(across(c(accepted_current_protein_nodes,DrugBank_current_nodes,STITCH_projected_current_nodes),~coalesce(.x,0L)),
 network_mapping_status=if_else(accepted_current_protein_nodes>0L,'accepted_current_protein_links','no_accepted_current_protein_links'),
 missing_coverage_is_not_negative=TRUE,affinity_required=FALSE)
stopifnot(nrow(coverage)==200L,!anyDuplicated(coverage$entity_id),
 identical(before,tools::md5sum(inputs)),all(induced$protein1<induced$protein2),
 !anyDuplicated(induced[c('protein1','protein2')]))
write_csv(accession_coverage,'data/curated/DrugBank_STRING_accession_mapping_coverage.csv')
write_csv(db_review,'data/curated/DrugBank_STRING_protein_mapping_evidence.csv',na='')
write_csv(legacy_candidates,'data/curated/STITCH_STRING_historical_bridge_candidates.csv',na='')
write_csv(legacy_coverage,'data/curated/STITCH_STRING_historical_bridge_coverage.csv')
write_csv(stitch_review,'data/curated/STITCH_STRING_protein_mapping_evidence.csv',na='')
write_csv(node_links,'data/curated/inventory_current_STRING_protein_links.csv')
write_csv(coverage,'data/curated/inventory_current_STRING_entity_coverage.csv')
write_csv(anchor,'outputs/21b_anchor_identity_summary.csv')
write_csv(anchor_summary,'outputs/21b_anchor_network_summary.csv')
write_csv(network_summary,'outputs/21b_network_coverage.csv')
write_csv(coverage |> count(origin,network_mapping_status,name='inventory_entities'),'outputs/21b_entity_coverage.csv')
write_csv(accession_coverage |> count(mapping_status,name='requested_accessions'),'outputs/21b_DrugBank_mapping_coverage.csv')
write_csv(legacy_coverage |> count(mapping_status,name='legacy_protein_ids'),'outputs/21b_STITCH_mapping_coverage.csv')
write_csv(node_links |> group_by(source_layer,evidence_stratum) |> summarise(inventory_entities=n_distinct(entity_id),
 distinct_current_proteins=n_distinct(current_STRING_id),distinct_entity_protein_links=n(),.groups='drop'),
 'outputs/21b_source_layer_summary.csv')
write_csv(tibble(file=inputs,md5=unname(before)),'outputs/21b_input_summary.csv')
manifest<-r('config/database_manifest.csv');manifest$status[manifest$source_id=='string']<-
 'human_v12.5_network_and_inventory_harmonisation_acquired_source_layers_separate'
write_csv(manifest,'config/database_manifest.csv',na='')
checklist<-r('config/RESKO_method_completion_checklist.csv')
checklist$project_status[checklist$step=='06']<-'DrugBank_STITCH_STRING_network_layers_available_missing_ambiguous_bridges_explicit'
checklist$next_work[checklist$step=='06']<-'Use source-aware human protein/network layers in matrix integration; retain identity/coverage holds'
checklist$project_status[checklist$step=='10']<-'human_network_ready_topology_statistics_pending'
checklist$next_work[checklist$step=='10']<-'Compute declared network statistics with disconnected-node handling after feature/background design; no binding interpretation'
write_csv(checklist,'config/RESKO_method_completion_checklist.csv')
report<-c('SCRIPT 21b - INVENTORY PROTEIN HARMONISATION AND HUMAN NETWORK PREPARATION',
 paste('Inventory entities retained:',nrow(coverage)),paste('Entities with accepted current STRING protein links:',sum(coverage$accepted_current_protein_nodes>0L)),
 paste('DrugBank entities with current protein links:',sum(coverage$DrugBank_current_nodes>0L)),
 paste('STITCH entities with projected current protein links:',sum(coverage$STITCH_projected_current_nodes>0L)),
 paste('Distinct current inventory-associated protein nodes:',n_distinct(node_links$current_STRING_id)),
 paste('Legacy STITCH protein IDs with unique accession bridge:',nrow(accepted_legacy)),
 paste('All-published-confidence induced network edges:',nrow(induced)),paste('Declared primary combined-score threshold:',primary),
 'Both human EEF1A anchors map by exact UniProt accession and agree with current STRING gene labels.',
 'Current exact-accession links and historical accession-bridge projections stay separate; no gene-name, BLAST-only or identical-ENSP fallback.',
 'Ambiguous/missing proteins are retained in coverage. Isoform suffixes are not stripped; historical bridges do not establish sequence equivalence.',
 'Complete source and all-published-confidence induced network retained; threshold 700 and sensitivity 400/900 are explicit project adaptations, not significance tests.',
 'Induced nodes are accepted inventory proteins plus anchors, with no external-neighbour expansion. Nodes without edges remain present.',
 'DrugBank roles and STITCH chemical-identity strata remain separate. Shared source reports do not become independent experiments.',
 'Network degree/proximity cannot establish compound binding, inhibition, allostery or RT-PPI disruption. No affinity exclusion or BH correction.',
 'Outputs: source-rich bridges/coverage/network locally; compact source-layer, anchor and confidence-sensitivity summaries in outputs/21b_*.',
 'Next: GO and DO human annotations, followed by declared feature matrices, network statistics, integration and validation.')
writeLines(report,'outputs/21b_harmonisation_summary.txt')
writeLines(c('PASS: explicit current/historical alias provenance, unique-node acceptance, ambiguity holds, exact human anchor labels.',
 'PASS: all 200 entities retained; network nodes include isolates; canonical unique edge keys and protected input hashes verified.'),
 'outputs/21b_harmonisation_validation.txt')
notes<-readLines('Notes.txt',warn=FALSE);start<-grep('^CURRENT CHECKPOINT',notes)[1];end<-grep('^PROJECT AIM',notes)[1]
stopifnot(start<end);writeLines(c(notes[seq_len(start-1)],'CURRENT CHECKPOINT - HUMAN STRING ACQUISITION/HARMONISATION 21a/21b COMPLETE',report,'',notes[end:length(notes)]),'Notes.txt')
cat(paste0('\nSCRIPT 21b RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/21b_harmonise_inventory_STRING_network.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
