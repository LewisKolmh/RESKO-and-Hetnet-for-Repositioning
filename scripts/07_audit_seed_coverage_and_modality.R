# 07: Audit missed indication mentions and prepare modality-aware seed tables.
# This is a local audit: no new associations are automatically accepted.
# Broad mention rules deliberately favour sensitivity and require review.
library(readr)
library(dplyr)
library(stringr)
library(tibble)

if (!all(dir.exists(c("scripts", "config", "data")))) stop("Run from the project root.")
drugs <- read_csv("data/curated/drugbank_drugs.csv",
                  col_types=cols(.default=col_character()))
review <- read_csv("data/curated/cancer_drug_review_decisions.csv",
                   col_types=cols(.default=col_character()))
cancers <- read_csv("config/cancers.csv", show_col_types=FALSE) |> filter(include)
rules <- read_csv("config/coverage_audit_terms.csv",
                  col_types=cols(.default=col_character())) |>
  semi_join(cancers, by="cancer_id")
keys <- c("cancer_id", "drugbank_id")
stopifnot(!anyDuplicated(drugs$drugbank_id), !anyDuplicated(review[keys]),
          all(cancers$cancer_id %in% rules$cancer_id),
          all(review$review_status=="reviewed"))

# 1. Independent, deliberately broad indication-mention screen.
# A mention of normal ovarian physiology can be a false positive. The review
# queue keeps its complete indication so that these can be distinguished.
broad <- lapply(seq_len(nrow(rules)), function(i) {
  pattern <- regex(rules$pattern[i], ignore_case=TRUE)
  drugs |>
    filter(str_detect(coalesce(indication, ""), pattern)) |>
    mutate(
      cancer_id=rules$cancer_id[i],
      audit_pattern=rules$pattern[i],
      matched_text=str_extract(indication, pattern),
      context_excerpt=str_extract(
        str_replace_all(indication, "[\r\n]+", " "),
        regex(paste0(".{0,100}(?:", rules$pattern[i], ").{0,150}"), TRUE)
      )
    )
}) |> bind_rows() |>
  left_join(cancers |> select(cancer_id, cancer_name), by="cancer_id")

# Compare against ALL previous decisions, including exclusions. Do not recycle
# already-reviewed false positives into a new candidate queue.
gaps <- broad |>
  anti_join(review, by=keys) |>
  mutate(
    source_url=paste0("https://go.drugbank.com/drugs/", drugbank_id),
    audit_status="unreviewed_mention",
    interpretation="Mention only; cancer treatment and target activity unverified"
  ) |>
  select(cancer_id, cancer_name, drugbank_id, drug_name, drug_type, groups,
         matched_text, context_excerpt, indication, audit_pattern, source_url,
         source_export_date, audit_status, interpretation) |>
  arrange(cancer_name, drug_name)
stopifnot(!anyDuplicated(gaps[keys]), nrow(inner_join(gaps, review, by=keys))==0L)

coverage <- cancers |> select(cancer_id, cancer_name) |>
  left_join(broad |> count(cancer_id, name="broad_mention_associations"), by="cancer_id") |>
  left_join(review |> count(cancer_id, name="previously_reviewed_associations"), by="cancer_id") |>
  left_join(gaps |> count(cancer_id, name="additional_mentions_to_review"), by="cancer_id") |>
  mutate(across(ends_with("associations"), ~coalesce(.x, 0L)),
         additional_mentions_to_review=coalesce(additional_mentions_to_review, 0L))

# 2. Combine accepted tiers for ANNOTATION only; preserve their distinctions.
# Use the reviewed worksheet, not a drug-level join that could transfer a
# breast decision to an ovarian association for the same drug.
accepted <- review |>
  filter(decision %in% c("include_primary", "include_investigational", "include_documented"))

annotate_modality <- function(table) {
  table |>
    mutate(
      modality_hint=case_when(
        str_detect(drug_name, regex("deruxtecan|emtansine|govitecan|mafodotin|soravtansine|SGN-15", TRUE)) ~
          "conjugate_name_hint",
        str_detect(drug_name, regex("autoleucel|vicleucel|dendritic|IDM-1", TRUE)) ~
          "cell_product_name_hint",
        str_detect(drug_name, regex("vaccine|ancovacivec|ladenovec|MVA-BN|INGN", TRUE)) ~
          "vaccine_or_gene_product_name_hint",
        drug_type=="biotech" ~ "DrugBank_biotech_requires_subtyping",
        drug_type=="small molecule" ~ "DrugBank_small_molecule_requires_structure_check",
        TRUE ~ "unclassified_requires_review"
      ),
      structure_review_route=if_else(
        modality_hint=="DrugBank_small_molecule_requires_structure_check",
        "prioritise_structure_and_identifier_mapping",
        "manual_modality_and_identity_review"
      ),
      direct_binding_eligibility="not_established",
      modality_evidence="DrugBank type plus explicit name-based heuristic; not structure-validated"
    )
}

annotated <- annotate_modality(accepted)
unique_seeds <- annotated |>
  group_by(drugbank_id, drug_name, drug_type, modality_hint,
           structure_review_route, direct_binding_eligibility) |>
  summarise(
    seed_tiers=paste(sort(unique(decision)), collapse=" | "),
    cancers=paste(sort(unique(cancer_name)), collapse=" | "),
    cancer_count=n_distinct(cancer_id),
    .groups="drop"
  ) |>
  arrange(structure_review_route, drug_name)
primary <- annotated |> filter(decision=="include_primary")
primary_mapping_queue <- primary |>
  filter(structure_review_route=="prioritise_structure_and_identifier_mapping")
primary_modality_queue <- primary |>
  filter(structure_review_route=="manual_modality_and_identity_review")
modality_summary <- annotated |>
  group_by(decision, modality_hint, structure_review_route) |>
  summarise(associations=n(), unique_drugs=n_distinct(drugbank_id), .groups="drop")

stopifnot(nrow(annotated)==nrow(accepted),
          nrow(primary_mapping_queue)+nrow(primary_modality_queue)==nrow(primary),
          !anyDuplicated(unique_seeds$drugbank_id))
dir.create("outputs", showWarnings=FALSE)
write_csv(gaps, "outputs/07_additional_indication_mentions.csv", na="")
write_csv(coverage, "outputs/07_seed_coverage_summary.csv")
write_csv(annotated, "data/curated/seed_associations_modality_audit.csv", na="")
write_csv(unique_seeds, "data/curated/unique_seed_drugs_modality_audit.csv", na="")
write_csv(primary_mapping_queue, "outputs/07_primary_structure_mapping_queue.csv", na="")
write_csv(primary_modality_queue, "outputs/07_primary_modality_review_queue.csv", na="")
write_csv(modality_summary, "outputs/07_seed_modality_summary.csv")

report <- c(
  "SCRIPT 07 - COVERAGE AND MODALITY AUDIT",
  paste("DrugBank records:", nrow(drugs)),
  paste("Records with indication text:", sum(!is.na(drugs$indication) & nzchar(drugs$indication))),
  paste("Additional unreviewed mention associations:", nrow(gaps)),
  paste("Unique drugs in accepted tiers:", nrow(unique_seeds)),
  paste("Primary structure-mapping queue:", nrow(primary_mapping_queue), "associations /",
        n_distinct(primary_mapping_queue$drugbank_id), "unique drugs"),
  paste("Primary manual-modality queue:", nrow(primary_modality_queue), "associations /",
        n_distinct(primary_modality_queue$drugbank_id), "unique drugs"),
  "", capture.output(print(coverage)), "",
  "Output interpretation:",
  "07_additional_indication_mentions.csv is a new REVIEW QUEUE, not new accepted seeds.",
  "Broad organ mentions include benign physiology, prevention and adverse effects.",
  "The audit measures gaps relative to the original screen, not recall against a gold standard.",
  "Missing indication text cannot be searched; no claim of exhaustive drug coverage is made.",
  "Modality hints are provisional. DrugBank small-molecule/biotech labels can be coarse",
  "or wrong and do not establish size, permeability, binding or suitability for this project.",
  "Peptides and macrocycles are not automatically discarded. Conjugate payloads must not",
  "inherit parent drug evidence without separately resolving their identities.",
  "Both queues retain all primary associations. No modality-based biological negative is assigned.",
  "Source review decisions and seed tiers are unchanged. No eEF1A activity is predicted.",
  "", "NEXT: review additional mentions and verify modality/structures before final seed selection;",
  "then map the retained drug identifiers to SIDER and report coverage by evidence tier."
)
writeLines(report, "outputs/07_audit_summary.txt")
print(coverage)
message("Additional mention associations to review: ", nrow(gaps))
message("Primary structure-mapping queue: ", n_distinct(primary_mapping_queue$drugbank_id), " drugs")
message("Primary manual-modality queue: ", n_distinct(primary_modality_queue$drugbank_id), " drugs")

# Every successful run records outputs AND what the results mean.
cat(paste0(
  "\nSCRIPT 07 RUN - ", format(Sys.time(), tz="Europe/Rome", usetz=TRUE), "\n",
  "Script: scripts/07_audit_seed_coverage_and_modality.R\n",
  "Inputs: DrugBank drug CSV, section 06 worksheet, cancer config,\n",
  "config/coverage_audit_terms.csv.\n",
  paste(report, collapse="\n"), "\n",
  "Other outputs: outputs/07_seed_coverage_summary.csv;\n",
  "outputs/07_seed_modality_summary.csv;\n",
  "data/curated/seed_associations_modality_audit.csv;\n",
  "data/curated/unique_seed_drugs_modality_audit.csv;\n",
  "outputs/07_primary_structure_mapping_queue.csv;\n",
  "outputs/07_primary_modality_review_queue.csv.\n"
), file="Notes.txt", append=TRUE)
message("Audit complete; Notes.txt updated. Section 06 decisions were not modified.")
