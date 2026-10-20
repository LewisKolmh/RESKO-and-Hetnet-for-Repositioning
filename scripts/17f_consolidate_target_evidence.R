# Consolidate bounded section 17 without turning missing coverage into inactivity.
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
paths <- c('outputs/16f_target_evidence_handoff_summary.csv',
 'data/curated/eef1a_reviewed_activity_evidence.csv',
 'data/curated/eef1a_STITCH_reviewed_inventory_evidence.csv',
 'data/curated/eef1a_inventory_exact_ligand_links.csv',
 'data/curated/eef1a_positive_binding_controls.csv',
 'data/curated/target_primary_literature_registry.csv')
hashes <- tools::md5sum(paths)
inventory <- r(paths[1]); activities <- r(paths[2]); stitch <- r(paths[3]); links <- r(paths[4]); controls <- r(paths[5]); lit <- r(paths[6])
stopifnot(nrow(inventory)==200L,!anyDuplicated(inventory$entity_id),nrow(links)==0L,
 nrow(activities)==39L,n_distinct(activities$assay_chembl_id)==23L,
 n_distinct(activities$molecule_chembl_id)==20L,
 !any(inventory$entity_id=='seed:DB04977'),
 all(activities$target_tax_id=='9606'),all(activities$queried_gene %in% c('EEF1A1','EEF1A2')))
associations <- stitch |> group_by(entity_id,gene) |> summarise(
 STITCH_source_links=n(),experimental_channel_links=sum(experimental_channel_present=='TRUE'),
 STITCH_mapping_relations=paste(sort(unique(mapping_relation)),collapse='|'),.groups='drop')
coverage <- crossing(inventory,gene=c('EEF1A1','EEF1A2')) |>
 left_join(associations,by=c('entity_id','gene'),relationship='one-to-one') |>
 mutate(target_species='Homo sapiens',target_taxon=9606L,
 target_uniprot=if_else(gene=='EEF1A1','P68104','Q05639'),
 STITCH_source_links=coalesce(STITCH_source_links,0L),
 experimental_channel_links=coalesce(experimental_channel_links,0L),
 exact_ChEMBL_ligand_matches=0L,
 human_direct_binding_conclusion='not_established_in_bounded_sources',
 human_isoform_selectivity_conclusion='not_established_in_bounded_sources',
 human_RT_PPI_disruption_conclusion='not_established_in_bounded_sources',
 human_allostery_conclusion='not_established_in_bounded_sources',
 translation_sparing_conclusion='not_established_in_bounded_sources',
 interpretation='Missing source/identity coverage is not evidence of inactivity; association scores are not binding evidence')
stopifnot(nrow(coverage)==400L, sum(coverage$STITCH_source_links)==13L,
 sum(coverage$experimental_channel_links)==0L,
 identical(sort(unique(coverage$entity_id)),sort(inventory$entity_id)))
write_csv(coverage,'data/curated/section17_human_target_evidence_coverage.csv',na='')
write_csv(coverage |> count(gene,human_direct_binding_conclusion,STITCH_source_links,name='inventory_entities'),
 'outputs/17f_human_target_coverage.csv')
# Human scope is mandatory. Rabbit binding is retained only as contextual reference.
controls <- controls |> mutate(human_target_qualification='does_not_qualify_as_purified_human_binding',
 target_scope='Human EEF1A1 and EEF1A2; nonhuman evidence is contextual only')
write_csv(controls,'data/curated/eef1a_positive_binding_controls.csv')
write_csv(controls,'outputs/17f_control_scope_summary.csv')
references <- activities |> group_by(molecule_chembl_id,standard_inchikey,molecule_name) |>
 summarise(source_records=n(),review_layers=paste(sort(unique(review_layer)),collapse='|'),
 inventory_exact_match=FALSE,role='noninventory_reference_requires_review_not_ranked_candidate',.groups='drop')
write_csv(references,'data/curated/section17_noninventory_reference_ligands.csv',na='')
# Source access and evidence gaps are explicit completion boundaries.
limitations <- tribble(~issue,~status,~next_action,
 'Human isoform-specific direct binding','Not established for the 200 inventory entities','Acquire verified human EEF1A1/EEF1A2 binding evidence before qualification',
 'SPR source PMID 21599001','Six reported records; full methodology/probe identity pending','Retrieve and review full methods and supporting information',
 'PMID 29498519','Eleven source activities quarantined; Expression of Concern','Retain hold until source integrity can be resolved',
 'STITCH upstream study identifiers','Unresolved in detailed export; six inventory association leads only','Trace original supporting studies; do not count overlapping sources independently',
 'Primary full-text access','Only two of eleven scoped papers acquired as full XML','Recover remaining full texts; endpoint errors are access failures',
 'RT/eEF1A context literature','PPI mechanism context acquired; chemical countercontrol full-text unavailable','Verify compound-binding partner from primary methods before using as control',
 'Human selectivity/allostery/translation sparing','Not established by this bounded section','Keep separate evidence requirements; modelling alone cannot qualify allostery',
 'Database completeness','ChEMBL exact-target and STITCH legacy queries are bounded','Continue remaining original source coverage and targeted human literature searches',
 'Statistical ranking','No valid enrichment test family defined in section 17','Specify null/background and test family before BH adjustment')
write_csv(limitations,'outputs/17f_limitations_summary.csv')
# All protected inputs except the intentionally enriched control registry must be unchanged.
stopifnot(identical(hashes[-5],tools::md5sum(paths[-5])))
write_csv(tibble(file=paths,input_md5=unname(hashes),output_md5=unname(tools::md5sum(paths))),
 'outputs/17f_input_summary.csv')
report <- c('SCRIPT 17f - HUMAN TARGET EVIDENCE CONSOLIDATION; BOUNDED SECTION 17 COMPLETE',
 'Scope: human EEF1A1 (P68104) and EEF1A2 (Q05639); nonhuman evidence cannot qualify a human binder.',
 '200 original inventory entities retained as 400 human entity/isoform coverage rows.',
 '39 ChEMBL records / 23 assays / 20 noninventory ligands / 7 source publications; zero exact inventory matches.',
 '11 records quarantined for Expression of Concern; six reported SPR records await full methods; one isoform contradiction held.',
 'STITCH: 3198 target association rows; 13 inventory links across six drugs, zero inventory experimental-channel links.',
 'Plitidepsin remains a separate contextual binding reference, not a seed: rabbit eEF1A2 Kd 79 +/- 30 nM does not qualify purified-human binding.',
 'No inventory compound qualifies for human direct binding plus selective/allosteric RT-PPI inhibition on evidence reviewed here. This is not evidence that no such compound exists.',
 'Full-text, upstream STITCH provenance and human selectivity gaps remain explicit; completion is bounded source accounting, not exhaustive literature adjudication.',
 'Outputs: full human coverage and noninventory references locally; compact coverage, limitations, control scope, hashes and validation in outputs/17f_*.',
 'Next: section 18 should resolve human binding/SPR identity gaps and expand original-source coverage before any candidate ranking or structural allostery claims.')
writeLines(report,'outputs/17f_section_summary.txt')
writeLines(c('PASS: 200 original identities and tiers retained, 400 unique human entity/isoform rows.',
 'PASS: 39 source activities / 23 assays / 20 ligands accounted for; zero exact inventory links.',
 'PASS: 13 STITCH inventory links accounted for; zero inventory experimental-channel links.',
 'PASS: nonhuman positive control remains outside inventory and cannot qualify human binding; protected inputs unchanged.'),
 'outputs/17f_section_validation.txt')
notes <- readLines('Notes.txt',warn=FALSE)
checkpoint_start <- grep('^CURRENT CHECKPOINT',notes)[1]
checkpoint_end <- grep('^PROJECT AIM',notes)[1]
stopifnot(!is.na(checkpoint_start),!is.na(checkpoint_end),checkpoint_end>checkpoint_start)
notes <- c(notes[seq_len(checkpoint_start-1)],report,'',notes[checkpoint_end:length(notes)])
writeLines(notes,'Notes.txt')
cat(paste0('\nSCRIPT 17f RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
