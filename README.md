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
source("scripts/11e_extract_current_label_terms.R")
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
