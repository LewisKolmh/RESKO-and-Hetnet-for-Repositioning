# Apply the first evidence-reviewed batch and regenerate script 06 outputs.
# Evidence was inspected on 2026-10-16, not fetched automatically by this script.
# DailyMed labels supplement the original DrugBank discovery records.
library(readr)
library(dplyr)
library(tibble)

if (!all(dir.exists(c("scripts", "config", "data")))) stop("Run from the project root.")

batch <- tribble(
  ~cancer_id, ~drugbank_id, ~label_id, ~evidence_excerpt, ~review_notes,
  "breast_cancer", "DB12001", "be4bc0de-0fdc-4d46-8d25-be43c79e6a06",
  "Paraphrase: US label includes treatment of HR-positive/HER2-negative breast cancer in specified early and advanced settings.",
  "Include cancer-directed use. Indication is subtype/setting-specific; endocrine combinations and prior-treatment conditions apply. Label sections 1.1-1.2; not a claim of efficacy in all breast cancers.",
  "breast_cancer", "DB09073", "e0e6412f-50b4-4fd4-9364-62818d121a07",
  "Paraphrase: US label includes HR-positive/HER2-negative advanced or metastatic breast cancer with specified endocrine partners.",
  "Include on the basis of labelled HR-positive/HER2-negative combination use, with treatment-history restrictions. Section 1.1; this selected indication is not an exhaustive label summary.",
  "breast_cancer", "DB09074", "741ff3e3-dc1a-45a6-84e5-2481b27131aa",
  "Paraphrase: US label includes selected germline BRCA-mutated, HER2-negative early high-risk and metastatic breast cancers after chemotherapy.",
  "Include cancer-directed treatment with germline BRCA, HER2 and treatment-history restrictions. Sections 1.4-1.5. Maintain a separate ovarian association for this same drug.",
  "ovarian_cancer", "DB09074", "741ff3e3-dc1a-45a6-84e5-2481b27131aa",
  "Paraphrase: US label includes maintenance of BRCA-mutated advanced epithelial ovarian cancer responding to first-line platinum chemotherapy.",
  "Include cancer-directed maintenance therapy; biomarker, histology, stage and chemotherapy-response restrictions apply. Section 1.1 supports this association. Other labelled maintenance settings are not exhaustively represented here.",
  "ovarian_cancer", "DB11793", "c15c7b7e-4b7f-4489-bbbc-884caeee0669",
  "Paraphrase: US capsule label describes maintenance in advanced epithelial ovarian cancer following response to first-line platinum chemotherapy.",
  "Include niraparib active moiety mapped from the niraparib tosylate product. Section 1.1, capsule label revised May 2024. Establishes a labelled ovarian maintenance use; not an exhaustive or formulation-independent eligibility summary.",
  "multiple_myeloma", "DB00188", "b605b541-4bfd-db3e-9710-11c27f6dee6c",
  "Paraphrase: US prescribing information lists treatment of adults with multiple myeloma.",
  "Include cancer-directed treatment. Generic bortezomib injection label, section 1.1; match active moiety to DrugBank record.",
  "multiple_myeloma", "DB00480", "5fa97bf5-28a2-48f1-8955-f56012d296be",
  "Paraphrase: US label includes adult multiple myeloma with dexamethasone and maintenance after autologous stem-cell transplantation.",
  "Include cancer-directed treatment and post-transplant maintenance. Section 1.1. The supportive-care search flag may reflect another indication and does not invalidate this myeloma use."
) |>
  mutate(
    evidence_url = paste0("https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=", label_id),
    review_status = "reviewed",
    use_role = "cancer_directed_treatment",
    regulatory_status = "approved_for_this_cancer",
    decision = "include_primary",
    jurisdiction = "US",
    reviewed_by = "Codex assistant (source review)",
    reviewed_on = "2026-10-16"
  )

path <- "data/curated/cancer_drug_review_decisions.csv"
worksheet <- read_csv(path, col_types = cols(.default = col_character()))
keys <- c("cancer_id", "drugbank_id")
stopifnot(!anyDuplicated(batch[keys]), !anyDuplicated(worksheet[keys]))
if (nrow(anti_join(batch, worksheet, by = keys)) > 0L) {
  stop("A batch association is absent from the original worksheet.")
}
fields <- c("evidence_url", "evidence_excerpt", "review_status", "use_role",
            "regulatory_status", "decision", "jurisdiction", "reviewed_by",
            "reviewed_on", "review_notes")

# Preserve the previous worksheet and refuse to replace conflicting reviews.
backup <- "outputs/cancer_drug_review_before_batch_001.csv"
if (!file.exists(backup)) {
  if (!file.copy(path, backup)) stop("Could not back up the review worksheet.")
}
for (i in seq_len(nrow(batch))) {
  index <- which(worksheet$cancer_id == batch$cancer_id[i] &
                   worksheet$drugbank_id == batch$drugbank_id[i])
  if (worksheet$review_status[index] != "pending") {
    same <- all(vapply(fields, function(field) {
      identical(worksheet[[field]][index], batch[[field]][i])
    }, logical(1)))
    if (!same) stop("Existing review conflicts with batch: ", batch$drugbank_id[i])
  }
}
for (i in seq_len(nrow(batch))) {
  index <- which(worksheet$cancer_id == batch$cancer_id[i] &
                   worksheet$drugbank_id == batch$drugbank_id[i])
  for (field in fields) worksheet[[field]][index] <- batch[[field]][i]
}

write_csv(batch, "outputs/cancer_drug_review_batch_001_evidence.csv", na = "")
write_csv(worksheet, path, na = "")
source("scripts/06_review_cancer_drugs.R")
cat(
  paste0(
    "\nSCRIPT 06a - Evidence review batch 001 applied successfully\n",
    "Script: scripts/06a_apply_review_batch_001.R\n",
    "Evidence checked: 2026-10-16; US labels on DailyMed, reached via NCI.\n",
    "Seven associations / six drugs accepted: abemaciclib, palbociclib,\n",
    "olaparib (breast and ovarian separately), niraparib, bortezomib, lenalidomide.\n",
    "Source summaries are paraphrases. Label restrictions are retained in notes.\n",
    "This is an initial review batch, not a representative or final seed set.\n",
    paste0("Pending associations after validation: ", sum(reviewed$review_status == "pending"), ".\n"),
    "Evidence: outputs/cancer_drug_review_batch_001_evidence.csv.\n",
    "Backup: outputs/cancer_drug_review_before_batch_001.csv.\n",
    "Output explanation: batch evidence traces the accepted associations to labels;\n",
    "the backup preserves pre-batch decisions. Six drugs yield seven associations\n",
    "because olaparib is included for two cancers. These are treatment seeds,\n",
    "not eEF1A activity predictions or a complete cancer-treatment inventory.\n",
    "Next: continue evidence review, then assess SIDER identifier/coverage.\n",
    "eEF1A binding and inhibition remain unassessed.\n"
  ), file = "Notes.txt", append = TRUE
)
