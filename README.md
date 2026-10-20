# RESKO and HetNet for Repositioning

Repository: https://github.com/LewisKolmh/RESKO-and-Hetnet-for-Repositioning

Local folder: `RESKO HetNet eEF1A_Cancer` (formerly `RESKONew`).

R workflow for disease-forward drug repurposing around eEF1A, using ovarian cancer,
breast cancer and multiple myeloma as starting contexts. Direct eEF1A binding,
inhibition, HIV-1 RT interaction disruption, allostery and translation-sparing
selectivity are separate evidence requirements; none is established by side-effect
similarity alone.

## Running the project

Open this folder in VS Code and use it as the R working directory. Scripts check
for project contents rather than requiring a particular folder name. Downloaded
source files and API responses stay in `data/raw/`; curated tables are under
`data/curated/`, with coverage and audit reports under `outputs/`.

Run scripts in order and read `Notes.txt` for the actual run history and results.
Early scripts 00–02 are currently empty files: their console-run setup is described
in the notes and must be restored before claiming complete end-to-end reproduction.
DrugBank requires the user's own permitted XML export, registered in the local
manifest; the original export is not part of this repository.

```r
# From the project root, after the previously completed seed and identity steps:
source("scripts/11d_collect_current_labels.R")
source("scripts/11f_resolve_label_moieties.R")
source("scripts/11e_extract_current_label_terms.R")
source("scripts/11g_review_label_supplement.R")
source("scripts/12_pilot_label_sider_similarity.R")
source("scripts/12a_validate_label_supplement.R")
source("scripts/13_audit_side_effect_similarity.R")
source("scripts/14a_resolve_candidate_identities.R")
source("scripts/14b_collect_candidate_labels.R")
source("scripts/helpers/validate_candidate_labels.R")
source("scripts/14c_resolve_sacrosidase_context.R")
source("scripts/15_extract_comparable_label_terms.R")
source("scripts/helpers/validate_comparable_labels.R")
source("scripts/16a_review_candidate_summaries.R")
source("scripts/16b_review_seed_summaries.R")
source("scripts/helpers/validate_seed_summary_review.R")
source("scripts/16c_audit_profile_compatibility.R")
source("scripts/helpers/validate_profile_compatibility.R")
source("scripts/16d_review_full_label_context.R")
source("scripts/16e_reassess_context_and_overlap.R")
source("scripts/16f_finalize_section16.R")
source("scripts/17a_acquire_target_records.R")
source("scripts/17b_expand_assays_and_match_ligands.R")
source("scripts/helpers/validate_target_source_records.R")
source("scripts/17c_acquire_stitch_associations.R")
source("scripts/helpers/validate_stitch_associations.R")
```

## Data layers

- Cancer seed tiers and indication decisions remain separate.
- Original SIDER 4.1 features are the historical baseline.
- Historical named-label context is separate from exact-structure SIDER mappings.
- Current openFDA labels supplement coverage; dictionary hits require contextual
  review before they become supplementary drug/PT features.
- Unknown coverage is missing, not a negative adverse-event observation.
- Current-label and SIDER-derived profiles must not be silently pooled. A later
  comparison needs compatible processing for candidate drugs as well as seeds.

The original RESKO source framework is retained: UMLS, DrugBank, SIDER, ChEMBL,
STITCH, STRING, Reactome, GO and Disease Ontology. PubChem, ChEBI and current FDA
labels are supplementary sources. Acquisition of one table does not imply full
acquisition of its source database.

## Section completion and GitHub

After each completed section: run appropriate checks, update `Notes.txt` with
inputs, outputs, counts, meaning and limitations, then commit and push the changes
to the user-designated repository. Preserve remote history and report the commit.
Downloaded raw data, cached API responses and temporary R session files are excluded
from Git; reproducible scripts, configuration, compact result tables and notes are versioned. Source-rich curated tables remain local;
curated code and approved compact feature exports are versioned.
The designated remote is `LewisKolmh/RESKO-and-Hetnet-for-Repositioning`, branch `main`.

## Current-label supplement results

The representative US-label collection covers 79 seeds (155 label records).
A separate reviewed common-reaction subset supplies 303 drug/PT pairs for the
32 primary compounds without strict SIDER mappings. Combination, mixed-cohort
and formulation contexts remain explicit. Eight profiles with explicit single-agent
context and at least five terms enter the descriptive overlap pilot; its ranks
are not validated target-binding predictions or statistical significance.

Compact features and coverage/validation reports are versioned in `outputs/`.
Full source evidence and API caches remain local. The original SIDER layer is
unchanged. See `Notes.txt` for counts, decisions, limitations and run history.

Sources: [SIDER](https://sideeffects.embl.de/download/),
[openFDA labels](https://open.fda.gov/apis/drug/label/),
[DailyMed](https://dailymed.nlm.nih.gov/dailymed/), and the specific
[EMA Pepaxti overview](https://www.ema.europa.eu/en/medicines/human/EPAR/pepaxti).
Preserve source attribution and publisher terms when reusing data.

## Shared-term and rank audit (script 13)

Script 13 retains the raw overlap scores, adds smoothed inverse-frequency weighted
Jaccard, and runs a declared sensitivity analysis removing terms seen in at least
20% of SIDER profiles. Too-small reduced profiles remain unranked. Settings are
versioned in `config/side_effect_audit_settings.csv`; score ties remain visible.

Of 88 leading baseline pairs including ties, 73 share only common terms. Thirty
SIDER groups enter the descriptive identity-review queue; eleven have some
lower-frequency overlap. These are unverified reference groups, not eEF1A-binding
candidates. Common-only matches are explicitly flagged. No p-values or BH
adjustment are produced by this descriptive analysis.

Term contributions, rank comparisons, stability, reference term prevalence and
the identity-review queue are exported as `outputs/13_*summary.csv`. Input hashes
and validation/session reports are retained. Full tables stay local. Next, resolve
chemical identities and construct comparable current-label profiles before
interpreting the ranks or joining target evidence.

## Candidate identity cross-check (script 14a)

All 30 queued SIDER groups have current PubChem structure keys; 19 groups match
DrugBank by full InChIKey, recovering source aliases and external identifiers.
These are current-record links, not confirmed historical SIDER identities.
Sixteen historical names have exact alias support; the remaining names require
review. The sacrosidase lead resolves to an unrelated small-molecule record and
is explicitly blocked from automatic label mapping in
`config/candidate_identity_exceptions.csv`. One group already matches a seed key.

API responses and DrugBank enrichment are cached locally; reruns reuse caches.
Compact results are in `outputs/14a_*`; source-rich tables remain local. Script 14b follows with ingredient/active-moiety validation, keeping unresolved
identities and combination/formulation contexts separate.

## Current candidate labels (script 14b)

The bounded DailyMed collection retains 51 distinct SPL sets across 22 SIDER
groups (54 group/label records). Eleven groups have 23 single-ingredient records
with adverse-reaction text supported by source-record UNII and ingredient or
explicit active-moiety name agreement. These remain inputs for context review;
they are not accepted PT profiles or evidence of eEF1A activity.

Searches use UNIIs first, names as fallback, one API page per query and at most
three recently published sets per group. Thirteen searches hit a panel/page limit.
The sacrosidase conflict is blocked; unmatched, combination and nonhuman records
stay outside selection. Coverage gaps are missing data. Current clinical label
attribution does not establish historical SIDER identity or salt equivalence.

Coverage, source crosswalk, input hashes and validation are in `outputs/14b_*`.
Raw XML, normalised AR text/tables and source-rich review tables stay local. The
validation helper independently checks selected XML records and text. Next is
script 15: process candidate and seed labels comparably and review PT context.
API reference: [DailyMed SPL search](https://dailymed.nlm.nih.gov/dailymed/webservices-help/v2/spls_api.cfm).

## Sacrosidase exception (script 14c)

Script 14c verifies DrugBank DB06760 (a biotech enzyme) and the Sucraid SPL
ingredient name/UNII, recovering one current label in a separate named-context
selection: `data/curated/candidate_label_named_context_selection.csv`. This is
clinical label attribution, not a repair of the historical CID or a strict
structure match. The incorrect small-molecule lead remains blocked, the original
SIDER files and strict selections remain unchanged, and sacrosidase stays outside
the small-molecule structure screen. Script 15 must retain this context distinction.
The source and decision summary are in `outputs/14c_*`.

## Comparable label extraction (script 15)

Script 15 uses the same dictionary, normalisation, boundaries and sentence/table
processing for all 210 collected entity/label associations (206 distinct label
documents). It processes 155 seed associations, all 54 candidate review records
and the separately verified sacrosidase label. Unresolved candidate identities
and named-biologic context remain explicitly separate from supported identities.

The extraction yields 59,500 candidate evidence rows. Every term retains a source
segment, label version, cohort, identity layer and context flags. All matches
require review; none becomes a confirmed adverse effect solely through a literal
dictionary match. Repeated labels and table/narrative mentions are not independent
evidence. Full source text and review queues remain local.

`outputs/15_*` contains coverage, review workload, input hashes and validation.
The validation helper checks every match against its source span/dictionary and
exports full-inventory coverage, including drugs without collected labels.
Historical SIDER data and previously reviewed features remain separate. Next is
bounded source-context review before constructing comparable PT features or ranks.

## Candidate summary context review (script 16a)

All 39 candidate common-summary passages (234 queued evidence rows) have explicit
batch decisions. Six reviewed exact reaction-list spans yield 38 entity/PT/context
features. Rules and four exact-PT disambiguation decisions are versioned in
`config/candidate_positive_list_review_rules.csv` and
`config/candidate_PT_mapping_review.csv`. All other surrounding content stays
outside those accepted spans. Full source text/review evidence remains local.

Technetium-pentetate radiopharmaceutical summaries and carfilzomib combination
summaries are held. Duplicate source lists, distinct populations and discontinuation
endpoints are not silently pooled. Regimen/population qualifiers remain attached;
all profiles have `comparison_ready = FALSE`. This is a bounded candidate batch,
not completion of seed, table, narrative or named-biologic review. The compact
feature export is `outputs/16a_candidate_PT_features.csv`; coverage and validation
are in `outputs/16a_*`. Next is seed summary review using the same approach.

## Seed representative summary review (script 16b)

Seventy representative seed passages were inspected using pinned source keys and
exact list boundaries. Sixty-four source-context profiles yield 634 reviewed
entity/PT features. Six representatives are held for source formatting, comparator
repair or endpoint/grade-specific review. The other 2,848 queued matches from
nonrepresentative passages remain pending; they are not labelled reviewed.

Features retain source versions, route, population/regimen qualifiers, seed tiers
and cancer membership. Intravesical, ocular, albumin-bound and non-cancer transplant
contexts remain explicit. All profiles have `comparison_ready = FALSE`; missing
regimen details are not inferred single-agent. Four exact PT decisions and the
previously reviewed hot-flush mapping are explicit; other ambiguities stay held.

Rules and holds are in `config/seed_*review*.csv` and
`config/seed_summary_context_holds.csv`. Compact features, coverage, workflow counts
and validation are in `outputs/16b_*`; full source/review evidence stays local. The
validation helper independently checks accepted source spans, PT provenance and
original routes. Next is a compatibility audit across seed and candidate identity,
route, regimen and profile coverage before ranking.

## Profile compatibility audit (script 16c)

All 384 seed/candidate combinations from 64 seed profiles and six candidate
profiles are retained with explicit context constraints. The declared strict
sensitivity stratum requires explicit single-agent context, a shared single label
route, at least five reviewed PTs, no same-key self-comparison and reviewed
frequency/cohort alignment. These constraints define this analysis stratum;
different routes are not universally incompatible for exploratory comparison.

Four seed profiles and no candidate summaries explicitly establish single-agent
context. Ten profiles fall below the declared term minimum, 159 pairs share a
single route, and frequency/cohort alignment remains pending. Consequently no
strict pair is ready. A single active ingredient does not establish trial
monotherapy. No missing context is turned into zero similarity or negative
eEF1A evidence. Prior readiness flags are not used as the audit decision rule.

`outputs/16c_*` contains profile/pair audits, actionable profile follow-up, all 200
inventory entries, checksums and validation. Exact source-list prefixes stay local
in `data/curated/profile_comparison_source_context.csv`. Settings are versioned
in `config/profile_compatibility_settings.csv`. Next is 16d: inspect full candidate
source cohorts/regimens and threshold definitions before reassessing compatibility.

## Section 16 complete: source context and handoff (16d–16f)

16d verifies all six candidate full-label contexts plus sacrosidase. The selected
apremilast trials allow concomitant DMARDs; the ibrutinib summary pools mono and
combination treatment. Pirfenidone is a placebo-filtered >=10% IPF list; the
pralatrexate >35% summary is sparse. Calcium acetate has conflicting summary/table
frequency presentation. Glycerol phenylbutyrate retains adult short-term active
control context; sacrosidase retains clinical-count/postmarketing context.

16e keeps all 384 pairs with explicit limitations. No original strict comparison
is ready. A declared, relaxed descriptive stratum permits one olaparib/pirfenidone
comparison with the same nominal cutoff, route and adequate list size: 8 shared
terms across a 23-term union (Jaccard 0.3478). All eight terms are common under the
original SIDER prevalence rule. This does not provide an eEF1A-specific signal,
statistical significance or a validated target ranking. Supportive/background
therapy uncertainty is not confused with absence of study-arm information.

16f consolidates 672 reviewed entity/PT/context features (64 seed and six candidate
profiles). A local ledger accounts for every one of the 59,500 extracted evidence
rows: 675 accepted mentions support the feature export; the other rows retain
explicit held, rejected, pending or other statuses. Completion refers to the
bounded review, audit and handoff; it does not claim exhaustive manual review of
all narratives and tables. Compact results are in `outputs/16f_*`; source-rich
evidence remains local. Final validation checks all six stage input manifests,
source provenance, feature keys, independent aligned-pair arithmetic and parsing.

All 200 original entities remain available for independent target-evidence search.
The handoff keeps EEF1A1/EEF1A2 binding and downstream activity fields missing until
assessed; missing side-effect coverage is not a negative filter. Sacrosidase's
wrong CID remains quarantined. Next is section 17 target-evidence integration,
including STITCH/ChEMBL, while preserving the original nine-source framework.

## Target evidence acquisition (17a–17b; section 17 ongoing)

Human EEF1A1 (P68104/CHEMBL1795120) and EEF1A2 (Q05639/CHEMBL1795122) are
resolved against UniProt and exact single-component ChEMBL targets. The complete
query returns 39 activity records: 28 for EEF1A1 and 11 for EEF1A2, representing
23 assays, 20 ligands and seven source documents. Full source JSON, assay metadata,
confidence/organism fields and versions stay local. UniProt is supplementary
identity evidence; the original nine-source framework is retained.

No source ligand exactly matches an inventory full InChIKey. All 200 inventory
entries retain explicit source-coverage status; no match is not a negative binding
result. One assay description names EEF1A1 but is assigned to EEF1A2 and is held
for review. Assay type B can describe an effect on tRNA binding rather than direct
compound binding. Neither database annotation nor an affinity-named assay is
automatically accepted as direct binding, allostery or RT-PPI disruption.

The EEF1A1 STITCH gene-name query is ambiguous; accession lookup plus current
UniProt Ensembl cross-references supports its alternate protein ID. EEF1A2 matches
current STRING/Ensembl cross-references. Both checked crosswalks are in
`config/eef1a_STITCH_protein_crosswalk.csv`; sequence-variant equivalence remains
unasserted. HTTPS certificate validation failed at STITCH, and public HTTP
transport is explicit in acquisition provenance; TLS validation is not disabled.
No chemical-protein associations are acquired merely by resolving a target ID.

Compact reports are in `outputs/17a_*` and `outputs/17b_*`. Exact source records
remain under `data/raw/target_evidence/` and source-rich curated evidence stays
local. The independent helper checks every ligand key/assay target and the full
inventory join. Stage 17c below adds the detailed STITCH associations; primary-source and
positive-control review follows.

References: [ChEMBL web services](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services),
[EEF1A1 UniProt](https://www.uniprot.org/uniprotkb/P68104/entry),
[EEF1A2 UniProt](https://www.uniprot.org/uniprotkb/Q05639/entry), and
[STITCH 5 source paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4702904/).

## STITCH detailed association acquisition (17c complete)

The official complete human detailed STITCH 5.0 file is downloaded locally
(84,697,010 compressed bytes; source modified 2016-08-16). This is legacy,
unsupported data. Streaming all 15,473,939 source rows retains 3,198 exact
eEF1A protein-ID association rows: 1,651 EEF1A1 and 1,547 EEF1A2, involving
1,769 distinct chemical IDs. No API neighbourhood cap or additional score
threshold is applied. Raw source/caches remain outside Git.

Experimental, prediction, database, textmining and combined scores remain separate
on the original 0–1000 scale. There are 3,069 text-mining-only rows and 36 rows
with experimental-channel support across the whole target subset. Those are
association rows, not independent experiments: merged/stereo representations may
repeat support, and this file does not separate direct from transferred evidence.
Scores are neither affinities nor p-values and do not establish binding or
isoform selectivity.

Thirteen source-ID/protein links touch six inventory drugs: fluorouracil, estradiol,
doxorubicin, dexamethasone, testosterone enanthate and epirubicin. None has
experimental-channel support. Five have stereo-ID links with existing current
identity support; epirubicin is merged-group context only. Original STITCH 5
chemical InChIKeys are not separately revalidated. Namespace translation follows
the cached publisher README; flat identifiers never become molecular structures.

Coverage retains all 200 entities against both isoforms (400 rows). Missing
identifiers/records remain missing coverage, not negative activity. Sacrosidase's
incorrect CID stays excluded. Public HTTP download transport is recorded; HTTPS
validation is not disabled. Publisher README/download documentation, detailed
source rows and the acquisition manifest stay local under `data/raw/stitch/` and
`data/curated/`. Compact source, channel and inventory reports are in
`outputs/17c_*`. The validation helper independently rescans the entire source
and checks every retained field, join count and input hash.

Next: 17d primary-source/provenance review of experimental/database-supported
leads, historical chemical forms and positive controls before binding
prioritisation. STITCH imports can overlap ChEMBL/DrugBank evidence and must not
be counted as independent confirmation.

### Section 17 completion: human target evidence

The bounded acquisition/review/consolidation stages 17a–17f are complete. The
required targets are **human EEF1A1 (P68104)** and **human EEF1A2 (Q05639)**.
Nonhuman experiments are contextual references and cannot qualify human binding.

The 200 original inventory entities remain unchanged, with 400 human
entity/isoform coverage rows. Exact-target ChEMBL acquisition yielded 39 activity
records, 23 assays and 20 ligands; none exactly matches the inventory. Eleven
records are quarantined because their source paper has an Expression of Concern;
six reported SPR records await full-method/probe-identity review. An isoform
annotation conflict remains held. STITCH contributes association-only leads for
six inventory drugs, with no experimental-channel inventory links.

Plitidepsin is a separate reference, not an added cancer seed. The reviewed
radioligand affinity experiment used rabbit muscle eEF1A2; its human-cell evidence
and modelled binding site do not establish purified-human affinity, isoform
selectivity, experimental allostery or translation sparing.

No inventory compound meets the combined human binding/selectivity/allosteric
RT-PPI requirements in this bounded review. Missing coverage is not inactivity.
Full-text access gaps and unresolved original STITCH study provenance remain open;
see `outputs/17f_limitations_summary.csv`. No statistical test family was defined,
so section 17 does not apply BH to association scores or affinity values.

Run `scripts/17e_review_target_evidence.R` after 17d, then
`scripts/17f_consolidate_target_evidence.R`. Full source-rich tables remain local;
compact results and validations are in `outputs/17e_*` and `outputs/17f_*`.
Section 18 should resolve human binding and identity gaps before ranking candidates.

### Section 18: association provenance follow-up

Scripts 18a–18c search and review the six STITCH-linked inventory drugs without
promoting associations into binding hits. Six fully paginated Europe PMC searches
returned 2,769 entity/publication index records. Title/abstract co-mention triage
selected 15 records (14 evidence groups after checking an AGR/MED duplicate).
The other 2,754 records remain pending full-text context, not rejected.

The bounded abstract review identifies human expression/resistance context for
fluorouracil and doxorubicin, not direct compound–eEF1A binding. Isoforms remain
separate: doxorubicin's eEF1A1 expression context does not verify its original
EEF1A2 STITCH edge. Three supplementary DrugMAP expression annotations were traced
to primary PMIDs, but the EEF1A2-specific gene-level data still need verification.
These supplementary references are not established as STITCH 5 upstream sources.
Direct legacy edge-page attempts returned application errors; HTTP 200 responses
were not accepted as evidence. Combination-treatment reviews cannot transfer
plitidepsin binding claims to dexamethasone. Epirubicin and doxorubicin remain
distinct, as do testosterone enanthate and testosterone.

Findings are in `outputs/18c_six_lead_summary.csv`; exact queries, source response
statuses, pending review coverage and input hashes are retained. Original human
target requirements and all 200 inventory identities/seed tiers remain unchanged.

### Core RESKO completion scope (20 October 2026)

Published human eEF1A binding affinity is an optional evidence annotation, not a
requirement for retaining a RESKO candidate or completing the computational
workflow. Missing affinity remains unknown. The core workflow is tracked against
McGarry et al. (2018), DOI `10.1016/j.knosys.2018.06.017`, in
`config/RESKO_method_completion_checklist.csv` and `outputs/19_*`.

RESKO is not yet fully reproduced: broad candidate targets/PPI networks,
Reactome, GO/DO annotations, chemical fingerprints, matrix integration and final
validation/ranking remain unfinished. The original paper already includes
STITCH. Current-label supplements and BH correction are explicitly labelled
adaptations/extensions. Conflicting paper threshold descriptions and missing
implementation detail are recorded in `config/RESKO_method_ambiguities.csv`;
these must not be hidden behind a claim of exact replication.

### Script 20a: Reactome human pathway layer acquired

Four official Reactome mapping/hierarchy files are frozen locally with source
URLs, retrieval/file-modification dates and checksums. The database-version
service reported version 97 before and after acquisition. There are 54,699
human lowest-level annotation rows, 161,425 all-level rows, 2,883 human pathway
nodes and 2,899 hierarchy edges. Lowest/all-level mappings remain separate to
avoid counting inherited memberships as independent observations.

The human mapping layer covers the database broadly for later seed/candidate
protein joins. Acquisition is complete; candidate-network integration and
pathway enrichment remain pending. This is a current-data adaptation of the
paper, not recovery of its original 2018 release. Full source-rich mappings stay
local; compact source, coverage, hierarchy and validation reports are in
`outputs/20a_*`. No affinity requirement or enrichment p-value is introduced.

### Script 20b: broad DrugBank protein roles joined to human pathways

The 200-entity inventory is retained. Reviewed identifiers connect 190 entities
to 189 DrugBank records; ten entities have no reviewed DrugBank identifier.
Selective streaming scanned the full XML export and extracted 1,096 role entries,
with 2,249 source article references. The new reader skips unselected target
sections before XML parsing; a boundary fixture and an independent full-record
comparison validate the extraction.

There are 164 entities with accepted human protein links, spanning 393 exact
UniProt accessions across targets, enzymes, carriers and transporters. Target-role
links alone cover 156 entities and 251 distinct human proteins. Roles remain
separate: metabolism/transport proteins are not automatically therapeutic targets.
Human joins require explicit taxonomy 9606, consistent organism labels and exact
accessions; isoform suffixes are not stripped. Named-biologic identity context
and sacrosidase's historical wrong structure remain separated/excluded from
strict structure-linked candidate evidence.

Reactome all-level features cover 161 inventory entities and 1,461 pathways
across all roles. Target-role pathway features cover 151 entities and 1,337
pathways. Lowest-level and ancestor-expanded evidence remain separate. Pathway
membership and source action annotations do not establish drug inhibition,
binding affinity, selectivity or allostery. Missing coverage remains unknown.

Full joins, references, record metadata and role-specific features are local under
`data/curated/`; compact source, role, entity and pathway counts plus input hashes
and validation are in `outputs/20b_*`. Published affinity is never an inclusion
gate. Broad STITCH chemical/protein associations, historical protein crosswalks,
STRING and GO/DO integration remain open; the core RESKO analysis is unfinished.
