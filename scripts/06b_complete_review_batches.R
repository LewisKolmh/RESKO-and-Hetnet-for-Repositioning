# Complete the remaining association-screening batches, not a claim that
# every historical trial or current formulation has been exhaustively audited.
# Curated decisions use the supplied DrugBank export plus official sources
# inspected on 2026-10-16. No web search is silently repeated on reruns.
library(readr)
library(dplyr)
library(stringr)
library(tibble)

if (!all(dir.exists(c("scripts", "config", "data")))) stop("Run from the project root.")
path <- "data/curated/cancer_drug_review_decisions.csv"
review <- read_csv(path, col_types = cols(.default = col_character()))
keys <- c("cancer_id", "drugbank_id")
stopifnot(!anyDuplicated(review[keys]))
for (field in c("evidence_basis", "research_stage", "review_batch")) {
  if (!field %in% names(review)) review[[field]] <- NA_character_
}
review <- review |>
  mutate(
    evidence_basis = if_else(review_status == "reviewed" & is.na(evidence_basis),
                             "individual_US_label", evidence_basis),
    research_stage = if_else(review_status == "reviewed" & is.na(research_stage),
                             "labelled_treatment", research_stage),
    review_batch = if_else(review_status == "reviewed" & is.na(review_batch),
                           "001", review_batch)
  )

ids <- function(x) strsplit(x, " +")[[1]]
nci <- c(
  breast_cancer = "https://www.cancer.gov/about-cancer/treatment/drugs/breast",
  ovarian_cancer = "https://www.cancer.gov/about-cancer/treatment/drugs/ovarian",
  multiple_myeloma = "https://www.cancer.gov/about-cancer/treatment/drugs/multiple-myeloma"
)
# These are manually checked active-moiety mappings, not fuzzy name matches.
# Supportive-care entries on the NCI lists are deliberately not in these sets.
primary_nci <- list(
  breast_cancer = ids("DB12015 DB01217 DB01101 DB12218 DB00531 DB16410 DB01248 DB00997 DB06374 DB00445 DB08871 DB01590 DB00990 DB00947 DB00441 DB00014 DB19043 DB15275 DB04845 DB01259 DB01006 DB14967 DB00351 DB00563 DB11828 DB01229 DB09037 DB06366 DB11730 DB12893 DB11760 DB00675 DB00539 DB00072 DB14962 DB05773 DB11652 DB19006 DB00570"),
  ovarian_cancer = ids("DB00112 DB00958 DB00515 DB00531 DB00997 DB00441 DB01042 DB12489 DB14976 DB12332 DB01030"),
  multiple_myeloma = ids("DB08889 DB00262 DB16738 DB00531 DB09331 DB00997 DB06317 DB15395 DB16665 DB14811 DB09570 DB01042 DB08910 DB11942 DB16678 DB16655")
)
documented <- list(
  breast_cancer = ids("DB00357 DB06719 DB00286 DB14655 DB02187 DB09381 DB00783 DB01185 DB08905 DB06710 DB00984 DB04575 DB00894 DB13944"),
  ovarian_cancer = ids("DB00488 DB01181 DB00361"),
  multiple_myeloma = character()
)

# Hand-curated special cases. Blank jurisdiction means no jurisdiction claimed.
special <- tribble(
 ~cancer_id, ~drugbank_id, ~use_role, ~decision, ~regulatory_status, ~jurisdiction, ~evidence_url, ~evidence_basis, ~research_stage, ~review_notes,
 "breast_cancer", "DB05867", "uncertain", "defer", "unverified", "",
 "https://go.drugbank.com/drugs/DB05867", "DrugBank_indication_and_mechanism", "role_verification_needed",
 "The indication template says treatment but the named technetium-labelled antibody may be an imaging construct; target binding alone does not establish therapy. Verify the original imaging/treatment study and compound identity before seed inclusion.",
 "breast_cancer", "DB00112", "cancer_directed_treatment", "include_primary", "approved_for_this_cancer", "EU",
 "https://www.ema.europa.eu/en/documents/product-information/avastin-epar-product-information_en.pdf", "EU_SmPC", "labelled_treatment",
 "EU SmPC section 4.1 includes metastatic breast cancer with paclitaxel or selected capecitabine use. Record EU scope; do not infer US breast approval.",
 "breast_cancer", "DB00361", "cancer_directed_treatment", "include_primary", "approved_for_this_cancer", "UK",
 "https://www.medicines.org.uk/emc/product/15312/smpc", "UK_SmPC", "labelled_treatment",
 "UK oral vinorelbine SmPC section 4.1 supports advanced breast cancer after anthracycline-containing regimen failure or unsuitability. Do not transfer that licence to ovarian cancer.",
 "breast_cancer", "DB01616", "investigational_treatment", "include_investigational", "investigational", "",
 "https://go.drugbank.com/drugs/DB01616", "DrugBank_indication", "preclinical",
 "The breast-cancer evidence is a cell-based MG132 sensitisation investigation, not an established treatment or verified clinical trial.",
 "breast_cancer", "DB06202", "prevention", "exclude", "unverified", "",
 "https://go.drugbank.com/drugs/DB06202", "DrugBank_indication", "historical_unspecified",
 "This export describes breast-risk reduction, not treatment of established breast cancer. This classification is specific to the supplied evidence, not a denial of later treatment research.",
 "breast_cancer", "DB01050", "prevention", "exclude", "unverified", "",
 "https://go.drugbank.com/drugs/DB01050", "DrugBank_indication", "historical_unspecified",
 "The breast-cancer mention concerns investigational prophylaxis, not a treatment indication.",
 "breast_cancer", "DB00481", "prevention", "exclude", "not_applicable", "",
 "https://go.drugbank.com/drugs/DB00481", "DrugBank_indication", "prevention",
 "Risk reduction of invasive breast cancer is prevention; do not use as an established breast-treatment seed.",
 "breast_cancer", "DB06713", "uncertain", "defer", "unverified", "",
 "https://go.drugbank.com/drugs/DB06713", "DrugBank_indication", "hypothesis_only",
 "Only a potential mechanism for breast use is stated. No treatment or study evidence sufficient for a seed; retain as unresolved rather than a negative.",
 "breast_cancer", "DB05531", "uncertain", "defer", "unverified", "",
 "https://go.drugbank.com/drugs/DB05531", "DrugBank_indication", "identity_verification_needed",
 "The broad indication template alone does not resolve the unusual compound name or establish a credible breast-treatment investigation. Verify compound identity and primary study before inclusion.",
 "breast_cancer", "DB06548", "uncertain", "defer", "unverified", "",
 "https://go.drugbank.com/drugs/DB06548", "DrugBank_indication", "role_verification_needed",
 "The indication template does not distinguish anticancer investigation from bone-directed supportive use. Do not infer direct antitumour activity from the disease name.",
 "multiple_myeloma", "DB06548", "uncertain", "defer", "unverified", "",
 "https://go.drugbank.com/drugs/DB06548", "DrugBank_indication", "role_verification_needed",
 "The indication template does not distinguish anticancer investigation from bone-directed supportive use. Keep outside accepted seeds until a primary study establishes the intended role.",
 "ovarian_cancer", "DB05292", "uncertain", "defer", "unverified", "",
 "https://go.drugbank.com/drugs/DB05292", "DrugBank_indication", "identity_verification_needed",
 "The entire indication is only 'Ovarian cancer'. Without a supporting use description or primary study this is not enough to accept an investigational treatment seed.",
 "ovarian_cancer", "DB02187", "unrelated_mention", "exclude", "not_applicable", "",
 "https://go.drugbank.com/drugs/DB02187", "DrugBank_indication", "not_applicable",
 "False match spans primary ovarian failure and a separate breast-cancer indication. No ovarian cancer treatment is stated.",
 "ovarian_cancer", "DB13956", "unrelated_mention", "exclude", "not_applicable", "",
 "https://go.drugbank.com/drugs/DB13956", "DrugBank_indication", "not_applicable",
 "False match spans primary ovarian failure and a prostate-cancer indication. No ovarian malignancy treatment is stated.",
 "ovarian_cancer", "DB11735", "unrelated_mention", "exclude", "not_applicable", "",
 "https://go.drugbank.com/drugs/DB11735", "DrugBank_indication", "not_applicable",
 "The ovarian-cancer mention is a proposed dietary association, not cancer treatment.",
 "ovarian_cancer", "DB15254", "cancer_directed_treatment", "include_primary", "approved_for_this_cancer", "US",
 "https://www.fda.gov/drugs/resources-information-approved-drugs/fda-grants-accelerated-approval-combination-avutometinib-and-defactinib-kras-mutated-recurrent-low", "FDA_approval_notice", "labelled_treatment",
 "FDA accelerated approval of the avutometinib/defactinib combination for KRAS-mutated recurrent low-grade serous ovarian cancer after prior systemic therapy; not general ovarian use.",
 "ovarian_cancer", "DB12282", "cancer_directed_treatment", "include_primary", "approved_for_this_cancer", "US",
 "https://www.fda.gov/drugs/resources-information-approved-drugs/fda-grants-accelerated-approval-combination-avutometinib-and-defactinib-kras-mutated-recurrent-low", "FDA_approval_notice", "labelled_treatment",
 "FDA combination approval with avutometinib, not defactinib monotherapy. KRAS-mutated recurrent low-grade serous ovarian cancer after prior systemic therapy.",
 "ovarian_cancer", "DB09037", "cancer_directed_treatment", "include_primary", "approved_for_this_cancer", "US",
 "https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-pembrolizumab-paclitaxel-platinum-resistant-epithelial-ovarian-fallopian-tube-or", "FDA_approval_notice", "labelled_treatment",
 "FDA 2026 approval with paclitaxel, with or without bevacizumab: PD-L1 CPS >=1 platinum-resistant epithelial ovarian/fallopian/peritoneal cancer after one or two prior regimens. Not broad monotherapy approval.",
 "ovarian_cancer", "DB05109", "cancer_directed_treatment", "include_primary", "approved_for_this_cancer", "EU",
 "https://www.ema.europa.eu/en/medicines/human/EPAR/yondelis", "EMA_EPAR", "labelled_treatment",
 "EMA authorisation covers relapsed platinum-sensitive ovarian cancer with pegylated liposomal doxorubicin. Orphan designation alone in the DrugBank text would not establish approval.",
 "multiple_myeloma", "DB17447", "cancer_directed_treatment", "include_primary", "approved_for_this_cancer", "US",
 "https://www.fda.gov/drugs/resources-information-approved-drugs/fda-grants-accelerated-approval-linvoseltamab-gcpt-relapsed-or-refractory-multiple-myeloma", "FDA_approval_notice", "labelled_treatment",
 "FDA accelerated approval for relapsed/refractory adult myeloma after at least four prior lines including specified drug classes.",
 "multiple_myeloma", "DB15719", "cancer_directed_treatment", "include_primary", "approved_for_this_cancer", "US",
 "https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-belantamab-mafodotin-blmf-relapsed-or-refractory-multiple-myeloma", "FDA_approval_notice", "labelled_treatment",
 "FDA 2025 approval with bortezomib/dexamethasone after at least two prior lines including specified classes. Historical withdrawn group must not erase this subsequent combination approval.",
 "multiple_myeloma", "DB01234", "cancer_directed_treatment", "include_primary", "approved_for_this_cancer", "US",
 "https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=6b9cf32c-6d60-4853-bbfd-2304ffc2c129&type=display", "individual_US_label", "labelled_treatment",
 "Hemady oral dexamethasone is labelled with other anti-myeloma products. Do not classify every corticosteroid use as supportive care or transfer this to otic/ophthalmic products.",
 "multiple_myeloma", "DB06603", "cancer_directed_treatment", "include_documented", "unverified", "EU",
 "https://www.ema.europa.eu/en/medicines/human/EPAR/farydak", "EMA_search_and_DrugBank", "documented_treatment",
 "EMA search confirms combination myeloma indication but current full product page retrieval failed. Retain documented-use tier; do not claim independently verified current US approval.",
 "multiple_myeloma", "DB16627", "cancer_directed_treatment", "include_primary", "approved_for_this_cancer", "EU",
 "https://www.ema.europa.eu/en/medicines/human/EPAR/pepaxti", "EMA_EPAR", "labelled_treatment",
 "EU Pepaxti authorisation is separate from US Pepaxto withdrawal. Dexamethasone combination and prior-treatment restrictions apply; do not treat a CHMP variation recommendation as final EC approval."
)

supportive <- tribble(
  ~cancer_id, ~drugbank_id, ~reason,
  "breast_cancer", "DB06643", "Prevention of bone loss during aromatase-inhibitor therapy, not treatment of the breast tumour.",
  "breast_cancer", "DB00380", "Reduces doxorubicin-associated cardiomyopathy or treats extravasation; not cancer-directed treatment.",
  "breast_cancer", "DB00282", "Treatment of osteolytic bone complications/hypercalcaemia; separate from cancer-directed treatment.",
  "multiple_myeloma", "DB06643", "Prevents skeletal-related events; do not equate skeletal protection with direct anti-myeloma activity.",
  "multiple_myeloma", "DB14939", "Stem-cell mobilisation for collection and transplantation; not direct myeloma treatment.",
  "multiple_myeloma", "DB06809", "Stem-cell mobilisation for collection and transplantation; not direct myeloma treatment.",
  "multiple_myeloma", "DB00282", "Bone-directed management of osteolytic myeloma lesions, not established direct tumour treatment.",
  "multiple_myeloma", "DB00399", "Bone-directed/hypercalcaemia management; keep outside the cancer-directed seed set despite NCI listing.",
  "multiple_myeloma", "DB05273", "DrugBank mechanism describes irradiation of metastatic bone sites for pain relief. The indication template alone does not establish an anti-myeloma treatment regimen.",
  "ovarian_cancer", "DB01143", "Reduces cumulative cisplatin renal toxicity; ovarian mention is supportive care, not tumour treatment."
)
diagnostic <- tribble(
  ~cancer_id, ~drugbank_id, ~reason,
  "breast_cancer", "DB15690", "PET imaging of ER-positive lesions as adjunct to biopsy.",
  "breast_cancer", "DB13967", "Sentinel lymph-node mapping dye before biopsy.",
  "breast_cancer", "DB18731", "Intraoperative fluorescence detection of residual cancerous tissue, not treatment.",
  "breast_cancer", "DB09397", "Radiotracer for lymph-node localisation/staging, not cancer treatment.",
  "breast_cancer", "DB09266", "Lymphatic mapping and sentinel-node localisation, not cancer treatment.",
  "ovarian_cancer", "DB15413", "Intraoperative optical identification of malignant lesions, not cancer treatment."
)

pending <- review |> filter(review_status == "pending")
if (nrow(pending) == 0L) {
  message("No unreviewed rows; preserving existing decisions.")
  source("scripts/06_review_cancer_drugs.R")
} else {
  plan <- pending
  for (field in c("evidence_basis", "research_stage", "review_batch")) plan[[field]] <- NA_character_
  for (i in seq_len(nrow(plan))) {
    r <- plan[i, ]; cancer <- r$cancer_id; id <- r$drugbank_id
    s <- special |> filter(cancer_id == cancer, drugbank_id == id)
    b <- supportive |> filter(cancer_id == cancer, drugbank_id == id)
    d <- diagnostic |> filter(cancer_id == cancer, drugbank_id == id)
    # Defaults are metadata only. Unclassified records stop, never auto-accept.
    plan$evidence_url[i] <- r$source_url
    plan$evidence_excerpt[i] <- paste0("Supplied DrugBank export: ", r$indication)
    plan$jurisdiction[i] <- NA_character_
    plan$evidence_basis[i] <- "DrugBank_indication"
    if (nrow(s)) {
      for (field in names(s)[!names(s) %in% keys]) plan[[field]][i] <- s[[field]][1]
    } else if (nrow(b) || nrow(d)) {
      plan$use_role[i] <- if (nrow(b)) "supportive_care" else "diagnostic"
      plan$decision[i] <- "exclude"
      plan$regulatory_status[i] <- "not_applicable"
      plan$research_stage[i] <- "non_treatment_use"
      plan$review_notes[i] <- if (nrow(b)) b$reason else d$reason
    } else if (id %in% primary_nci[[cancer]]) {
      plan$use_role[i] <- "cancer_directed_treatment"
      plan$decision[i] <- "include_primary"
      plan$regulatory_status[i] <- "approved_for_this_cancer"
      plan$jurisdiction[i] <- "US"
      plan$evidence_url[i] <- paste(nci[[cancer]], r$source_url, sep = " | ")
      plan$evidence_basis[i] <- "NCI_cancer_list_plus_DrugBank"
      plan$research_stage[i] <- "labelled_treatment"
      plan$review_notes[i] <- paste(
        "Active-moiety mapping checked against NCI cancer-specific treatment list.",
        "Full DrugBank indication preserves subtype, formulation, combination and",
        "treatment-history context. NCI list corroboration is not an exhaustive",
        "current product-label audit; no general use across every tumour subtype is claimed."
      )
    } else if (id %in% documented[[cancer]]) {
      plan$use_role[i] <- "cancer_directed_treatment"
      plan$decision[i] <- "include_documented"
      plan$regulatory_status[i] <- "unverified"
      plan$research_stage[i] <- "documented_historical_or_off_label_use"
      plan$review_notes[i] <- paste(
        "Export explicitly describes cancer-directed use, including palliative",
        "hormonal treatment where applicable. Current cancer-specific approval,",
        "market availability and jurisdiction not independently verified.",
        "Retain in a separate sensitivity-analysis set, not verified primary seeds."
      )
    } else if (str_detect(r$indication, regex("investigat|early investigation", TRUE)) ||
               (!str_detect(r$groups, "(^|\\|)approved(\\||$)") &&
                str_detect(r$groups, "investigational|experimental") &&
                str_detect(r$indication, regex("treatment", TRUE)))) {
      plan$use_role[i] <- "investigational_treatment"
      plan$decision[i] <- "include_investigational"
      plan$regulatory_status[i] <- "investigational"
      plan$research_stage[i] <- "historical_unspecified"
      plan$review_notes[i] <- paste(
        "Export explicitly records investigation in this cancer, not an approved",
        "cancer-specific treatment. Retain as historical research evidence only:",
        "trial phase, efficacy, clinical versus preclinical status and current",
        "development/enrolment are unverified. This is a separate exploratory set."
      )
    } else {
      stop("No curated disposition for ", cancer, " / ", id, " / ", r$drug_name)
    }
    plan$review_status[i] <- "reviewed"
    plan$reviewed_by[i] <- "Codex assistant (source-screening review)"
    plan$reviewed_on[i] <- "2026-10-16"
    plan$review_batch[i] <- c(breast_cancer="002", multiple_myeloma="003",
                              ovarian_cancer="004")[[cancer]]
  }

  # Preserve the source worksheet and the exact batch decisions for audit.
  backup <- "outputs/cancer_drug_review_before_batches_002_004.csv"
  if (file.exists(backup)) stop("Batch backup already exists; inspect before applying again.")
  if (!file.copy(path, backup)) stop("Could not back up decisions.")
  for (batch in unique(plan$review_batch)) {
    write_csv(plan |> filter(review_batch == batch),
              paste0("outputs/cancer_drug_review_batch_", batch, "_evidence.csv"), na="")
  }
  fields <- c("review_status", "use_role", "regulatory_status", "decision",
              "evidence_url", "evidence_excerpt", "jurisdiction", "reviewed_by",
              "reviewed_on", "review_notes", "evidence_basis", "research_stage", "review_batch")
  # Only previously pending keys are replaced; batch 001 is preserved.
  updated <- review |> anti_join(plan, by=keys) |> bind_rows(plan)
  stopifnot(nrow(updated)==nrow(review), !anyDuplicated(updated[keys]))
  write_csv(updated |> arrange(cancer_name, drug_name), path, na="")
  source("scripts/06_review_cancer_drugs.R")
  batch_summary <- updated |> count(review_batch, cancer_name, decision, name="associations")
  write_csv(batch_summary, "outputs/cancer_drug_all_batches_summary.csv")
  cat(paste0(
    "\nSCRIPT 06b - Remaining review batches applied (2026-10-16)\n",
    "Script: scripts/06b_complete_review_batches.R\n",
    "Batches 002/003/004 cover breast/myeloma/ovarian pending associations.\n",
    "Outputs: batch evidence CSVs, pre-batch backup, all-batches summary, and\n",
    "the validated script 06 seed/exclusion/deferred exports.\n",
    "Output explanation: every original association has a screening disposition.\n",
    "Primary seeds have official cancer-use corroboration; documented-use seeds\n",
    "have unverified current approval; investigations may be historical/preclinical.\n",
    "Deferred records remain unresolved and are not accepted seeds or negatives.\n",
    "Source evidence tiers remain explicit; this is not a complete label/trial audit.\n",
    "No binding, inhibition, allostery or translation-selectivity claims were assessed.\n"
  ), file="Notes.txt", append=TRUE)
}
