# RESKO and HetNet for Repositioning

An R workflow applying the disease-forward RESKO method to ovarian cancer,
breast cancer and multiple myeloma, with human **EEF1A1 (P68104)** and
**EEF1A2 (Q05639)** as the subsequent target-interest layer.

Method reference: [McGarry et al., 2018 — RESKO](https://doi.org/10.1016/j.knosys.2018.06.017).
Repository: https://github.com/LewisKolmh/RESKO-and-Hetnet-for-Repositioning

## Scientific workflow

1. Identify cancer-treatment seeds and retain primary, investigational and
   documented-use tiers separately.
2. Resolve drug identities and acquire SIDER side-effect profiles. Reviewed
   current-label evidence supplements coverage as a separate data layer.
3. Identify side-effect-derived candidates and assess profile comparability.
4. Collect chemical structures, DrugBank/STITCH protein associations, human
   STRING networks, Reactome pathways and GO/DO annotations.
5. Construct the paper-aligned similarity/network/ontology matrices, integrate
   them and validate a candidate ranking.
6. Annotate available human eEF1A evidence and prioritise hypotheses for testing.

**Published eEF1A binding affinity is optional evidence, not a candidate inclusion
filter.** Side-effect similarity, protein association, pathway membership, direct
binding, inhibition, isoform selectivity, RT-PPI disruption and allostery are
separate claims. Missing evidence does not demonstrate inactivity.

The full RESKO integration and final ranking are not yet complete. Cancer cohorts,
current-label supplements, modern database releases and BH correction are explicit
adaptations; STITCH is already a source in the original paper. Remaining work and
paper ambiguities are tracked in:

- `config/RESKO_method_completion_checklist.csv`
- `config/RESKO_method_ambiguities.csv`

## Running in VS Code

Open the project folder and use it as the R working directory. Current local
folder: `RESKO HetNet eEF1A_Cancer` (formerly `RESKONew`). Scripts identify the
project by its contents rather than folder name.

Use [scripts/README.md](scripts/README.md) for the stage order and dependencies.
Do not blindly rerun every script: some review stages preserve decisions and
some source acquisitions require existing local caches or licensed exports.

The three zero-byte setup placeholders were removed. A fresh environment needs
R packages, cancer configuration and the local database manifest established
before running the collection stages. The original console setup is recorded in
`Notes.txt`; full fresh-checkout reproduction remains unfinished.

Core protein/pathway stages, after their documented prerequisites:

```r
source("scripts/20a_acquire_reactome.R")
source("scripts/20b_integrate_drugbank_proteins_pathways.R")
source("scripts/helpers/validate_drugbank_pathway_integration.R")
source("scripts/20c_integrate_broad_STITCH_associations.R")
source("scripts/21a_acquire_human_STRING.R")
source("scripts/21b_harmonise_inventory_STRING_network.R")
source("scripts/helpers/validate_inventory_STRING_network.R")
source("scripts/22a_acquire_GO_annotations.R")
source("scripts/22b_acquire_DO_annotations.R")
source("scripts/22c_integrate_inventory_ontology_annotations.R")
source("scripts/helpers/validate_inventory_ontology_annotations.R")
```

Optional target-specific evidence collection lives in
`scripts/optional/eef1a_evidence/`. The paper-method completion audit is in
`scripts/maintenance/`. These scripts remain available because they document
relevant evidence and reproducibility checks; they are not affinity gates in the
core RESKO workflow. Historical Notes entries retain the paths used at the time.

## Data and evidence boundaries

- The integration inventory currently contains 170 cancer seed drugs and 30
  side-effect-derived review entries. These are not a completed ranked shortlist.
- Preserve exact structural matches, named-biologic context, historical merged
  chemical groups and unresolved identities separately. In particular,
  sacrosidase's incorrect historical small-molecule structure remains excluded.
- SIDER 4.1 is the historical baseline. Current-label profiles use separate,
  reviewed contexts; do not silently pool different evidence sources or eras.
- DrugBank targets, enzymes, carriers and transporters remain separate.
  Pharmacokinetic proteins are not automatically therapeutic targets.
- STITCH scores retain their evidence channels. They are association confidence,
  not affinities or p-values; overlapping database reports are not independent
  experiments. Historical protein aliases are not sequence-equivalence proof.
- Reactome lowest-level and ancestor-expanded mappings remain separate.
- Apply Benjamini–Hochberg only to a defined family of valid statistical tests
  with a stated null and background, retaining raw p-values and effect sizes.

## Results and provenance

`Notes.txt` records results by script with counts, interpretation and limitations.
Compact summaries, coverage, source hashes and validation reports are in
`outputs/`. Source-rich curated tables stay local in `data/curated/`; raw exports
and request caches stay local in `data/raw/`.

Script 20a acquired four official Reactome files: 54,699 human lowest-level
protein/pathway annotation rows, 161,425 all-level rows and 2,883 pathways.

Script 20b retained all 200 inventory entries and extracted 1,096 DrugBank
protein-role entries with 2,249 source article references. Accepted human links
cover 164 entries and 393 distinct proteins across all roles; target-role links
cover 156 entries and 251 proteins. Reactome features cover 161 entries and
1,461 pathways across all roles. See `outputs/20b_*` for role-specific counts.

DrugBank requires the user's own permitted XML export, registered in the local
manifest. Original framework sources are UMLS, DrugBank, SIDER, ChEMBL, STITCH,
STRING, Reactome, GO and Disease Ontology; PubChem, ChEBI and current FDA labels
are supplements. A downloaded table does not imply an entire source is integrated.

After each completed stage, validate the results, update `Notes.txt`, review the
Git diff, and commit/push code, configuration and approved compact reports to
`main`. Raw/cached data, source-rich tables and session files remain outside Git.

Script 20c scanned all 15,473,939 human STITCH source rows. It retained 37,701
selected chemical/protein rows, yielding 38,940 inventory/source-ID/protein links
for 116 inventory entries. There are 9,438 distinct legacy human protein IDs.
Historical direct Ensembl/UniProt aliases support 6,369 unique current-Reactome
projections; eight proteins have multiple current accession matches and remain
held. Other missing mappings remain explicit. Compound identity strata and source
score channels are separate; these are association/context results, not a new
list of proven binding hits. See `outputs/20c_*` for coverage and validation.

## Human STRING network layer (scripts 21a–21b)

The complete human STRING 12.5 functional-association release is acquired locally:
19,699 registry proteins, 10,826,054 directed source rows and 5,413,027 unique
undirected edges. Reciprocal rows were checked for exact score-channel equality
before deduplication. Source URLs, checksums, release and retrieval/file dates are
retained; full source data are not committed.

Current direct UniProt aliases link 164 DrugBank inventory entries into the
network. Explicit historical-accession bridges project STITCH associations for
116 entries, with 181 entries covered across the two sources. All 200 entries
remain in coverage. There are 8,946 unique legacy protein bridges, 184 ambiguous
bridges held and 308 missing/unsupported legacy mappings. Gene names and equal
Ensembl identifier text across releases are not used as mapping fallbacks.

The induced inventory network contains 8,969 current protein nodes and 2,011,287
edges at the published source floor of 150. All nodes and full source scores are
retained. The declared primary cutoff of 700 yields 113,712 edges; cutoffs 400
and 900 are separate sensitivity layers. These are documented project adaptations,
not significance tests or claims of exact original-paper replication.

Human EEF1A1 and EEF1A2 have separately verified accession/gene-label anchors.
Historical projections, DrugBank protein roles and STITCH chemical identity
strata remain distinct. Functional associations and network proximity do not
establish direct binding, inhibition, allostery or RT-PPI disruption. GO/DO
annotations, feature matrices, topology statistics and final RESKO ranking remain
unfinished. Compact reports are in `outputs/21a_*` and `outputs/21b_*`.

## GO and Disease Ontology annotation layers (scripts 22a–22c)

GO ontology release 2026-07-26 and its matching human GAF source are frozen
locally with URLs/checksums and release headers. Of 910,476 source annotation
rows, 860,943 are accepted positive human annotations. Negative NOT qualifiers,
ND evidence and non-UniProt objects remain explicitly status-accounted rather
than becoming positive protein features. Exact protein-form accessions, relation
qualifiers, references, dates, providers and evidence-code layers are preserved.
Only active is_a/part_of hierarchy relations are declared for later GO propagation.

DO hierarchy release 2026-09-30 is acquired. The current raw Alliance human
annotation endpoint timed out; **current raw DO gene annotation acquisition
remains incomplete**. The available annotation layer uses a frozen HDO.db 1.0.0
compiled snapshot with metadata date 2024-07-23: 76,544 of 77,019 gene/term rows
pass current human Entrez and DO-term checks. Its build merges Alliance-derived
and older DO2EG mappings without per-row qualifiers, evidence codes or references.
It is retained as **compiled legacy context**, not adjudicated causal evidence.
See `config/DO_annotation_source_policy.csv` for the explicit use boundary.

GO and DO legacy-context features each cover 181 of the 200 inventory entries;
all 200 remain tracked. Current protein-form mappings and gene-level disease
projections stay separate. Cancer concepts map by exact active DO names to
DOID:1612 (breast cancer), DOID:2394 (ovarian cancer) and DOID:9538 (multiple
myeloma); DO cross-references do not imply full UMLS acquisition.

Source-rich annotations and feature RDS files stay local. Compact source/status,
protein/gene mapping, evidence-layer, cancer-concept and validation reports are
in `outputs/22a_*`, `outputs/22b_*` and `outputs/22c_*`. Generic ontology roots
are excluded from features. Wang similarity, ancestor-expanded matrices,
valid enrichment/BH testing and final RESKO ranking remain unfinished. Annotation
membership does not establish compound binding, inhibition or allostery.

## Completed SIDER candidate stage (scripts 23a–23f)

Section 23 preserves the broad cancer baseline, predeclared mechanism cohorts,
cutoff sensitivities, current structure reconciliation and known-treatment context.
Final identity accounting covers 40 historical groups / 47 stereo-CID leads.
Historical forms are never silently pooled. A label-verified omitted thalidomide
myeloma indication prompted a separate corrected-cohort analysis; the original
IMiD results remain sensitivity context. Full run order is in scripts/README.md.

Five current-representation leads pass all strict setting checks without an
outstanding form/name-support hold: memantine and aripiprazole (breast aromatase
cohort); busulfan, methotrexate and bortezomib (ovarian nitrogen-mustard cohort).
These are chemistry-annotation proposals, not eEF1A inhibitors, a final RESKO
ranking or proven novel indications. Threshold-sensitive, known-treatment and
identity-held contexts remain in the final local handoff. The original 200-entry
inventory is unchanged. See outputs/23f_section_completion_summary.txt.

Next is the paper-matched atom-pair fingerprint/Tanimoto layer. Seed completeness,
current raw DO coverage, remaining matrices/statistics, valid enrichment/BH and
model validation remain explicit unfinished parts of the overall workflow.
