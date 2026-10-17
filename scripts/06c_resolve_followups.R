# Batch 005: resolve the six deferred associations using targeted source review.
# Findings checked on 2026-10-16. Evidence summaries below are paraphrases.
library(readr)
library(dplyr)
library(tibble)
if (!all(dir.exists(c("scripts", "config", "data")))) stop("Run from the project root.")

followups <- tribble(
 ~cancer_id, ~drugbank_id, ~use_role, ~decision, ~regulatory_status, ~evidence_url, ~evidence_basis, ~research_stage, ~evidence_excerpt, ~review_notes,
 "breast_cancer", "DB05867", "diagnostic", "exclude", "not_applicable",
 "https://pubmed.ncbi.nlm.nih.gov/16322892/", "primary_imaging_study", "diagnostic_research",
 "Paraphrase: Oliva et al. studied technetium-labelled 14F7 to demonstrate tumour-associated ganglioside expression in human breast cancer.",
 "The breast association supports radioimmunodetection/imaging, not an established treatment indication. Do not transfer antitumour claims about the antibody or myeloma to this breast imaging construct.",
 "breast_cancer", "DB05531", "unrelated_mention", "exclude", "not_applicable",
 "https://clinicaltrials.gov/study/NCT00237354?rank=1&term=PV-10&viewType=Table", "primary_trial_registry_identity_check", "identity_mismatch",
 "Paraphrase: NCT00237354 identifies the intralesional PV-10 intervention as a 10% rose bengal disodium preparation.",
 "Do not accept the apparent PV10/PV-10 synonym collision as evidence for alpha-pyrrolidinononanophenone. Trial intervention is rose bengal (different compound/DrugBank identity). Exclude this association for unreliable mapping; do not silently replace the compound or infer lack of biological activity.",
 "breast_cancer", "DB06713", "investigational_treatment", "include_investigational", "investigational",
 "https://pubmed.ncbi.nlm.nih.gov/12711003/", "primary_study_abstract", "preclinical_mechanism",
 "Paraphrase: Pasqualini et al. reported inhibition of estrone sulfatase by norelgestromin in T-47D and MCF-7 breast cancer cells.",
 "Retain as preclinical mechanistic research only. Enzyme effects in cancer cells do not establish tumour regression, approved breast treatment, or eEF1A engagement. Do not combine this evidence tier with clinical-treatment seeds.",
 "breast_cancer", "DB06548", "investigational_treatment", "include_investigational", "investigational",
 "https://pmc.ncbi.nlm.nih.gov/articles/PMC4273397/", "primary_preclinical_study", "preclinical",
 "Paraphrase: Yuen et al. examined nitrogen-containing bisphosphonates, including minodronic acid, in HER-expressing cancer-cell models that included MB231 breast cells.",
 "Direct cell-based anticancer investigation is distinguishable from bone-protective use. Retain in preclinical research tier only; no established clinical breast treatment or eEF1A binding is inferred.",
 "multiple_myeloma", "DB06548", "investigational_treatment", "include_investigational", "investigational",
 "https://pubmed.ncbi.nlm.nih.gov/11137564/", "primary_study_abstract", "preclinical",
 "Paraphrase: Takahashi et al. tested YM529/minodronic acid in human myeloma cell lines and observed growth inhibition and apoptosis.",
 "The source establishes an in-vitro antitumour investigation rather than only skeletal management. It is not clinical efficacy, a current trial, or eEF1A target evidence.",
 "ovarian_cancer", "DB05292", "investigational_treatment", "include_investigational", "investigational",
 "https://www.sec.gov/Archives/edgar/vprr/0305/03057960.pdf", "historical_sponsor_filing_and_DrugBank_description", "historical_clinical",
 "Paraphrase: A historical Medarex filing identifies Osidem/IDM-1 as MDX-210 combined with activated killer cells and reports ovarian-cancer clinical development.",
 "DrugBank XML description independently matches the MAK-cell/bispecific-antibody identity. Retain historical cell-therapy research; no current active programme or ovarian marketing approval claimed. A cell product is unlikely to be suitable for later small-molecule direct-binding analysis."
)
path <- "data/curated/cancer_drug_review_decisions.csv"
review <- read_csv(path, col_types=cols(.default=col_character()))
backup <- "outputs/cancer_drug_review_before_batch_005.csv"
fields <- setdiff(names(followups), c("cancer_id", "drugbank_id"))
for (i in seq_len(nrow(followups))) {
  idx <- which(review$cancer_id==followups$cancer_id[i] &
               review$drugbank_id==followups$drugbank_id[i])
  if (length(idx)!=1L) stop("Follow-up key missing or duplicated.")
  if (review$decision[idx]!="defer" &&
      !all(vapply(fields, function(f) identical(review[[f]][idx], followups[[f]][i]), logical(1)))) {
    stop("Refusing to overwrite a conflicting completed review.")
  }
}
if (!file.exists(backup) && !file.copy(path, backup)) stop("Backup failed.")
for (i in seq_len(nrow(followups))) {
  idx <- which(review$cancer_id==followups$cancer_id[i] &
               review$drugbank_id==followups$drugbank_id[i])
  for (field in fields) review[[field]][idx] <- followups[[field]][i]
  review$jurisdiction[idx] <- NA_character_
  review$reviewed_by[idx] <- "Codex assistant (targeted source review)"
  review$reviewed_on[idx] <- "2026-10-16"
  review$review_batch[idx] <- "005"
}
write_csv(followups, "outputs/cancer_drug_review_batch_005_evidence.csv", na="")
write_csv(review, path, na="")
source("scripts/06_review_cancer_drugs.R")

counts <- reviewed |> count(cancer_name, decision, name="associations")
write_csv(reviewed |> count(review_batch, cancer_name, decision, name="associations"),
          "outputs/cancer_drug_all_batches_summary.csv")
unique_drugs <- reviewed |> group_by(decision) |>
  summarise(associations=n(), unique_drugs=n_distinct(drugbank_id), .groups="drop")
write_csv(unique_drugs, "outputs/cancer_drug_review_unique_counts.csv")
stopifnot(nrow(reviewed)==217L, all(reviewed$review_status=="reviewed"),
          !any(reviewed$decision %in% c("pending", "defer")))

report <- c(
 "SECTION 06 - ASSOCIATION SCREENING COMPLETE (2026-10-16)",
 "All 217 original associations / 189 drugs have an explicit review disposition.",
 "There are no pending or deferred association-screening decisions.",
 "", capture.output(print(unique_drugs)), "", capture.output(print(counts)),
 "",
 "Interpretation:",
 "Primary: official cancer-use corroboration, with jurisdiction retained.",
 "Documented-use: treatment mentioned but current approval not independently verified.",
 "Investigational: historical and preclinical research, not necessarily active trials.",
 "Excluded: non-treatment uses, false mentions or an unreliable compound mapping.",
 "This completes screening of the original 217 matches, not a comprehensive audit",
 "of all cancer drugs or every current label/trial. Evidence tiers are explicit.",
 "NCI list corroboration is weaker than a complete formulation-specific label review.",
 "Do not use absence from a list as proof of no approval or biological activity.",
 "Candidate discovery used DrugBank; NCI/FDA/EMA/UK sources corroborate use only.",
 "The original text search can miss treatments (e.g. shared-list phrasing for paclitaxel",
 "and ovarian cancer); audit seed coverage before describing the set as exhaustive.",
 "Biologics, cell products and conjugates remain identified by drug_type; assess",
 "modality/identity suitability before direct small-molecule target analyses.",
 "No eEF1A binding, inhibition, allostery or PPI selectivity has been established.",
 "", "Next: audit seed coverage and modality, then obtain SIDER and map identifiers.",
 "Preserve primary, documented and investigational sets separately for sensitivity analyses."
)
writeLines(report, "outputs/06_review_completion_summary.txt")
cat(paste0(
 "\nSCRIPT 06c - Batch 005 targeted follow-ups and completion summary\n",
 "Script: scripts/06c_resolve_followups.R\n",
 "Six follow-ups resolved: four research associations accepted; imaging and\n",
 "a likely PV10 identifier collision excluded. Prior deferred decisions remain\n",
 "auditable in the pre-batch backup and batch 002-004 evidence exports.\n",
 "Outputs: outputs/cancer_drug_review_batch_005_evidence.csv;\n",
 "outputs/cancer_drug_review_before_batch_005.csv; updated script 06 exports;\n",
 "outputs/cancer_drug_review_unique_counts.csv;\n",
 "outputs/06_review_completion_summary.txt.\n",
 "Result explanation: all original matches have screening dispositions, with\n",
 "source certainty and research stage preserved. Not all are clinical treatments.\n",
 paste(capture.output(print(unique_drugs)), collapse="\n"), "\n",
 "No unreviewed/deferred rows remain. No eEF1A activity conclusions were drawn.\n",
 "Next: seed-coverage/modality audit, then SIDER acquisition and identifier mapping.\n"
), file="Notes.txt", append=TRUE)
message("Section 06 association screening complete; summary and Notes.txt updated.")
