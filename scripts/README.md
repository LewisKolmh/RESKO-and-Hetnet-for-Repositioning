# R workflow guide

Run scripts from the project root, regardless of the script's subdirectory.
This guide describes the roles of the code; `Notes.txt` records completed runs,
counts and limitations. The complete RESKO model and final ranking are unfinished.

| Main stages | Purpose |
|---|---|
| 03–04 | Register the permitted local DrugBank XML and extract drug records |
| 05–08 | Identify, review and tier cancer-treatment seeds |
| 09–10d | Extract identifiers, verify chemical identities and preserve distinct forms |
| 11a–11c | Acquire historical SIDER data and resolve seed mappings |
| 11d–11g | Supplement missing coverage with reviewed current-label evidence |
| 12–13 | Inspect side-effect overlap and sensitivity; exploratory comparisons, not final RESKO ranking |
| 14a–14c | Resolve candidate identities, labels and historical mapping conflicts |
| 15–16f | Build and review compatible side-effect profiles and the integration inventory |
| 20a | Acquire human Reactome mapping and pathway hierarchy |
| 20b | Extract broad DrugBank protein roles and join accepted human proteins to Reactome |
| 20c | Extract all STITCH neighbours for accepted inventory chemical IDs and review historical protein/pathway projections |
| 21a | Acquire the complete human STRING functional network, protein registry and aliases |
| 21b | Harmonise current/legacy protein links and prepare the source-aware human inventory network |
| 22a | Acquire GO hierarchy and human GAF evidence, retaining qualifiers, references and protein forms |
| 22b | Acquire current DO hierarchy and audit/version the available compiled human gene-context snapshot |
| 22c | Join ontology annotations to inventory proteins and preserve evidence/source strata |
| 23a | Disease-wise SIDER-only joint-pattern baseline, threshold sensitivities and explicit coverage limits |

The existing inventory is an integration/review cohort, not a completed
paper-matched candidate ranking. The disease-wise SIDER joint-pattern baseline,
chemistry, broader network and ontology layers and final model still need work.

## Optional human eEF1A evidence annotation

`scripts/optional/eef1a_evidence/` contains 17a–17f and 18a–18d.
These acquire target records, review source integrity and annotate published
binding/association context. Run in their original numbered order if this layer
needs refreshing. Missing measured affinity does not exclude a RESKO candidate.
17a/17b generated target identity/crosswalk inputs used by the optional 17c
eEF1A slice. Core 20c independently derives chemical leads from the earlier
identity tables and does not require running the optional target studies.

## Maintenance and validation

`scripts/maintenance/19_audit_RESKO_method_completion.R` checks progress against
the published method. Helpers contain request/cache utilities, streaming readers
and evidence checks; they are necessary support code, not competing analyses.
The live method checklist and ambiguity register are in `config/RESKO_*`.

The three empty setup placeholders were removed. Install the required R packages
and supply the cancer configuration and local manifests before a fresh run.
Historical notes and output hashes retain their original paths so previous results
remain auditable; active code references use the current script locations.

After 20c, run `scripts/helpers/validate_broad_STITCH_integration.R` to independently
rescan the full source and verify all retained chemical/protein/score fields.

After 21b, run `scripts/helpers/validate_inventory_STRING_network.R` to independently
check raw eEF1A-anchor source edges, all induced source values and mapping holds.
The primary confidence cutoff and sensitivity cutoffs are declared in
`config/STRING_network_projection_settings.csv`; functional associations do not
imply physical binding. Source/role/identity layers remain separate.

After 22c, run `scripts/helpers/validate_inventory_ontology_annotations.R` for raw
GAF anchor verification, ontology-ID ambiguity fixtures and feature checks. DO
compiled-context provenance limits are in `config/DO_annotation_source_policy.csv`;
current raw Alliance gene annotations remain an acquisition gap. Neither ontology
annotation nor missing annotation is a binding-activity result.

After 23a, run `scripts/helpers/validate_disease_SIDER_baseline.R`. No current full-cohort cancer pattern is assessable; covered-subset results do not meet the five-term flowchart minimum. Do not promote the two-term ovarian review queue to a final shortlist.
