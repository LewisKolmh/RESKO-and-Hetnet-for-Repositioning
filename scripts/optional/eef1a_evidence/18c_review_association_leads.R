# Explicit abstract-level review. Pending full-text records remain pending.
library(readr)
library(dplyr)
library(tibble)
r <- function(p) read_csv(p,col_types=cols(.default=col_character()))
paths <- c('data/curated/association_literature_triage.csv',
 'data/curated/association_literature_priority.csv',
 'data/curated/association_supplementary_primary_metadata.csv',
 'outputs/17c_linked_inventory_summary.csv','outputs/16f_target_evidence_handoff_summary.csv')
before <- tools::md5sum(paths)
triage<-r(paths[1]);priority<-r(paths[2]);supp<-r(paths[3]);leads<-r(paths[4]);inventory<-r(paths[5])
# Reviewed against the cached primary abstracts, not inferred from search ranks.
decisions <- tribble(~literature_key,~review_layer,~species_context,~gene_context,~interpretation,
 'AGR:IND609269770','duplicate_index_context','Human CRC cell models and in vivo study; exact in vivo species not resolved from abstract','EEF1A2',
 'Same normalized title and abstract as MED:40311992; do not count as independent support',
 'MED:40311992','resistance_pathway_context','Human CRC cell models and in vivo study; exact in vivo species not resolved from abstract','EEF1A2',
 'UGDH/EEF1A2/PRDX1 pathway implicated in 5-FU resistance; no compound-eEF1A affinity experiment reported in abstract; differs from original STITCH EEF1A1 link',
 'MED:41711048','secondary_review_context','Mixed literature review','EEF1A2',
 'Review of lactylation and chemotherapy resistance; not a new direct-binding experiment',
 'MED:36172672','nonhuman_comparator','Hamster','Eef1a1',
 'Atractylodes extract study includes 5-FU comparator; Eef1a1 expression context does not establish human binding',
 'MED:19561488','indirect_expression_context','Human HEp-2 cell line as reported','EEF1A1',
 '5-FU decreases NATH; EEF1A1 protein increases after NATH knockdown, not a measured direct 5-FU-eEF1A1 interaction',
 'MED:39919587','different_compound_nonhuman','Rat Leydig cells','EEF1A1|EEF1A2',
 'Bisphenol B exposure alters steroid output and protein markers; estradiol is an outcome, not evidence of estradiol-eEF1A binding',
 'PMC:PMC6284925','different_form_nonhuman_reference_gene','Beef cattle','EEF1A2',
 'Estradiol benzoate/trenbolone implants; EEF1A2 used as reference gene; ester cannot be merged with estradiol binding evidence',
 'MED:25458689','nonhuman_reference_gene','Rana nigromaculata frog','elongation factor 1 alpha 1',
 'Hormone exposure/reference-gene validation; not a drug-binding assay',
 'MED:26525398','nonhuman_reproductive_context','Beef cows','EEF1A1',
 'Estrus and hormone-protocol gene-expression context; not human binding or isolated estradiol mechanism',
 'MED:17288598','different_compound_nonhuman_reference_gene','Fathead minnow fish','ef1a unspecified human isoform',
 '17alpha-ethinylestradiol changes reference-gene expression; distinct compound and species from estradiol/human eEF1A',
 'MED:12101787','different_compound_nonhuman_reference_gene','Japanese medaka fish','elongation factor 1 alpha unspecified human isoform',
 '17alpha-ethinyloestradiol/nonylphenol exposure; elongation factor is PCR internal standard, not binding target',
 'MED:19487283','human_expression_context','Human cancer cell lines; xenograft context','EEF1A1',
 'Doxorubicin-associated senescence accompanies eEF1A1 reduction; expression effect, not affinity or inhibition measurement; does not verify original EEF1A2 STITCH edge',
 'MED:34357135','combination_therapy_review','Clinical/preclinical review','eEF1A isoforms',
 'eEF1A-targeting claims concern plitidepsin; dexamethasone co-treatment does not inherit those claims',
 'MED:32406881','combination_therapy_review','Clinical/preclinical review','EEF1A2',
 'Plitidepsin target discussion and dexamethasone combination; no dexamethasone-eEF1A2 binding measurement',
 'MED:30111169','combination_therapy_review','Clinical/preclinical review','EEF1A2',
 'Plitidepsin binding discussed alongside dexamethasone combination trial; no binding attribution to dexamethasone')
stopifnot(setequal(priority$literature_key,decisions$literature_key),!anyDuplicated(decisions$literature_key))
write_csv(decisions,'config/18c_abstract_review_decisions.csv')
review <- priority |> left_join(decisions,by='literature_key',relationship='many-to-one') |>
 mutate(review_basis='Cached primary abstract or identified secondary-review abstract; full methods not reviewed',
 human_direct_binding='not_established_by_abstract_review',human_allostery='not_established',
 human_RT_PPI_disruption='not_established',review_date='2026-10-20',
 STITCH_original_source_verified=FALSE)
# Assert duplicated AGR/MED abstract/title before grouping their evidence.
normalise <- function(x)gsub('[^[:alnum:]]','',tolower(x))
dup <- review |> filter(literature_key %in% c('AGR:IND609269770','MED:40311992'))
stopifnot(nrow(dup)==2L,n_distinct(normalise(dup$title))==1L,n_distinct(normalise(dup$abstract))==1L)
review <- review |> mutate(publication_evidence_group=if_else(literature_key=='AGR:IND609269770','MED:40311992',literature_key))
write_csv(review,'data/curated/association_priority_abstract_review.csv',na='')
# Abstracts establish study modality but do not expose the gene-level supplementary rows.
supp_review <- supp |> mutate(
 review_basis='Primary abstract plus exact reference link in supplementary DrugMAP annotation',
 primary_modality=case_when(pmid=='29803840'~'Human cardiac-microtissue doxorubicin transcriptomics',
 pmid=='20106945'~'Human HepG2/HepaRG chemical-exposure transcriptomics',
 pmid=='16144937'~'Human ovarian-cancer glucocorticoid apoptosis and expression study'),
 EEF1A2_specific_annotation_verification='pending_primary_gene_level_table_or_data; EEF1A2 is not named in acquired abstract',
 binding_conclusion='not_established',STITCH_original_source_verified=FALSE)
stopifnot(all(!grepl('EEF1A2|eEF1A2',supp_review$abstract)))
write_csv(supp_review,'data/curated/association_supplementary_primary_review.csv',na='')
# All source hits remain status-accounted. No silent full-text rejection.
statuses <- triage |> left_join(review |> select(entity_id,literature_key,review_layer,publication_evidence_group),
 by=c('entity_id','literature_key'),relationship='one-to-one') |>
 mutate(review_status=if_else(is.na(review_layer),'full_text_context_pending_not_rejected','abstract_reviewed_full_methods_pending'))
write_csv(statuses,'data/curated/association_literature_review_status.csv',na='')
summary <- tribble(~source_name,~finding,~qualification,~next_action,
 'Fluorouracil','Human-cell indirect EEF1A1 expression context and EEF1A2 resistance-pathway context',
 'No direct human binding established; resistance does not imply drug inhibition of the resistance protein',
 'Review 19561488 and 40311992 full methods; retain A1/A2 distinction',
 'Doxorubicin','Human-cell EEF1A1 senescence/expression context; supplementary EEF1A2 transcriptomics annotation',
 'No direct binding established; original STITCH isoform link still unresolved',
 'Review 19487283 full methods and gene-level data for 29803840',
 'Estradiol','Supplementary EEF1A2 expression annotation linked to human transcriptomics paper',
 'EEF1A2-specific row not verified from primary abstract; nonhuman and other-estrogen hits excluded as binding support',
 'Recover gene-level data for 20106945; preserve estradiol/ethinylestradiol/ester identities',
 'Dexamethasone','Supplementary expression annotation; priority abstracts describe plitidepsin co-treatment',
 'Plitidepsin binding cannot be assigned to dexamethasone; primary EEF1A2-specific annotation pending',
 'Recover gene-level data for 16144937 and original STITCH provenance',
 'Testosterone enanthate','Five indexed co-occurrence hits; none passes title/abstract co-mention triage',
 'Full-text context pending; no conclusion of inactivity or equivalence with testosterone',
 'Review five full-text contexts using exact ester identity',
 'Epirubicin','Merged-group STITCH context; no title/abstract co-mention passes current lexical triage',
 'Doxorubicin findings cannot be transferred to epirubicin; indexed full-text hits remain pending',
 'Resolve historical merged-group source and full-text contexts before chemical attribution')
summary <- leads |> select(entity_id,source_name,target_genes) |> left_join(summary,by='source_name',relationship='one-to-one') |>
 mutate(human_direct_binding_status='not_established_in_bounded_review',
 human_allosteric_RT_PPI_status='not_established_in_bounded_review',original_STITCH_study_status='unresolved')
stopifnot(nrow(summary)==6L,!anyNA(summary$finding),nrow(statuses)==2769L,
 sum(statuses$review_status=='abstract_reviewed_full_methods_pending')==15L,
 nrow(inventory)==200L,identical(before,tools::md5sum(paths)))
write_csv(summary,'outputs/18c_six_lead_summary.csv')
write_csv(review |> select(source_name,literature_key,pmid,doi,review_layer,species_context,gene_context,interpretation,publication_evidence_group),
 'outputs/18c_abstract_review_summary.csv',na='')
write_csv(statuses |> count(source_name,review_status,name='source_records'),'outputs/18c_review_coverage.csv')
write_csv(tibble(file=paths,md5=unname(before)),'outputs/18c_input_summary.csv')
report <- c('SCRIPT 18c - BOUNDED SIX-LEAD ABSTRACT REVIEW',
 '15 priority index records reviewed as 14 evidence groups; duplicate AGR/MED title and abstract verified, not independent experiments.',
 '2754 indexed full-text-context records remain explicitly pending. Lexical triage is not a full-text exclusion decision.',
 'Fluorouracil: indirect human EEF1A1 expression context and a human-cell EEF1A2 resistance pathway; neither establishes drug binding.',
 'Doxorubicin: human EEF1A1 senescence-associated reduction; differs from original EEF1A2 association and does not establish binding.',
 'Estradiol/dexamethasone/doxorubicin: supplementary EEF1A2 expression annotations traced to primary references, but gene-level rows remain unverified.',
 'Dexamethasone combination-treatment abstracts do not transfer plitidepsin target binding to dexamethasone.',
 'Testosterone enanthate and epirubicin remain unresolved; parent/ester and anthracycline identities are not pooled.',
 'Zero six-lead compounds qualify as direct human eEF1A binders or selective allosteric RT-PPI inhibitors on this bounded abstract review.',
 'Original STITCH upstream study provenance remains unresolved; no claim of exhaustive absence or completion of full-text review.',
 'No inventory additions, seed-tier changes, binding score or BH adjustment. Original 200-entity inventory protected.',
 'Outputs: detailed review/status tables locally; six-drug findings, review coverage and validation in outputs/18c_*.',
 'Next: retrieve primary gene-level tables and full methods, review the five testosterone-enanthate contexts and recover high-value full-text sources before expanding prioritisation.')
writeLines(report,'outputs/18c_review_summary.txt')
writeLines(c('PASS: all 2769 entity/publication search records accounted for; 15 priority rows manually adjudicated.',
 'PASS: duplicate AGR/MED abstracts verified and grouped; original six leads and 200-entity inventory preserved.',
 'PASS: human A1/A2 distinction, compound identity differences, pending gene-level expression annotations and unresolved provenance retained.',
 'PASS: protected input hashes unchanged; no direct-binding or allostery claim promoted.'),'outputs/18c_review_validation.txt')
notes <- readLines('Notes.txt',warn=FALSE);start<-grep('^SCRIPT 17f - HUMAN TARGET',notes)[1];end<-grep('^PROJECT AIM',notes)[1]
if(!is.na(start)&&!is.na(end)&&start<end){notes<-c(notes[seq_len(start-1)],
 'CURRENT CHECKPOINT - SECTION 18a/18b/18c COMPLETE; FULL-TEXT AND ORIGINAL STITCH PROVENANCE REVIEW OPEN',
 report,'',notes[end:length(notes)]);writeLines(notes,'Notes.txt')}
cat(paste0('\nSCRIPT 18c RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/optional/eef1a_evidence/18c_review_association_leads.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
