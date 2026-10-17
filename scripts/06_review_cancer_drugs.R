# 06: Prepare a review worksheet, validate decisions, and export reviewed seeds.
# Run from the project root. Existing review decisions are never overwritten.
# After editing the worksheet, rerun this script to validate and export them.

library(readr)
library(dplyr)
library(stringr)

if (!all(dir.exists(c("scripts", "config", "data")))) {
  stop("Run from the project root containing scripts/, config/ and data/.")
}

input_path <- "data/curated/cancer_drug_candidates_for_review.csv"
review_path <- "data/curated/cancer_drug_review_decisions.csv"
notes_path <- "Notes.txt"
candidates <- read_csv(input_path, col_types = cols(.default = col_character()))
keys <- c("cancer_id", "drugbank_id")

if (anyDuplicated(candidates[keys]) || anyNA(candidates[keys])) {
  stop("Candidate keys must be complete and unique.")
}

roles <- c("cancer_directed_treatment", "investigational_treatment", "prevention",
           "supportive_care", "diagnostic", "unrelated_mention", "uncertain")
decisions <- c("pending", "include_primary", "include_investigational",
               "include_documented", "defer", "exclude")
statuses <- c("pending", "reviewed")
regulatory <- c("unverified", "approved_for_this_cancer", "investigational",
                "not_approved_for_this_cancer", "not_applicable")

if (!file.exists(review_path)) {
  # These broad flags only prioritise review; they do not decide eligibility.
  worksheet <- candidates |>
    select(-review_status, -review_notes) |>
    mutate(
      review_hint = case_when(
        str_detect(indication, regex("diagnos|imaging|scintigraph|detect", TRUE)) ~
          "Check for diagnostic use",
        supportive_language_present == "TRUE" ~ "Check for supportive care",
        prevention_language_present == "TRUE" ~ "Check for prevention",
        str_detect(groups, "investigational") ~ "Check investigational status",
        TRUE ~ "Check cancer-directed indication and approval"
      ),
      review_status = "pending",
      use_role = "uncertain",
      regulatory_status = "unverified",
      decision = "pending",
      evidence_url = NA_character_,
      evidence_excerpt = NA_character_,
      jurisdiction = NA_character_,
      reviewed_by = NA_character_,
      reviewed_on = NA_character_,
      review_notes = NA_character_
    )
  write_csv(worksheet, review_path, na = "")
  message("Created worksheet. All decisions are pending: ", review_path)
} else {
  message("Reading existing decisions without overwriting them: ", review_path)
}

review <- read_csv(review_path, col_types = cols(.default = col_character()))
required <- c(keys, "review_status", "use_role", "regulatory_status", "decision",
              "evidence_url", "evidence_excerpt", "jurisdiction", "reviewed_by",
              "reviewed_on", "review_notes")
if (!all(required %in% names(review))) stop("Required review columns are missing.")
if (anyDuplicated(review[keys]) || anyNA(review[keys]) ||
    nrow(anti_join(candidates, review, by = keys)) > 0L ||
    nrow(anti_join(review, candidates, by = keys)) > 0L) {
  stop("Review worksheet must contain exactly the original candidate keys.")
}

allowed <- function(x, choices) all(!is.na(x) & x %in% choices)
if (!allowed(review$review_status, statuses) || !allowed(review$use_role, roles) ||
    !allowed(review$decision, decisions) ||
    !allowed(review$regulatory_status, regulatory)) {
  stop("Invalid review values. See outputs/cancer_drug_review_instructions.txt.")
}

# Rejoin decisions to the original evidence. Edits to copied names/indications
# in the worksheet cannot silently replace the source records in seed outputs.
extra_fields <- intersect(c("evidence_basis", "research_stage", "review_batch"), names(review))
reviewed <- candidates |>
  select(-review_status, -review_notes) |>
  left_join(review |> select(all_of(c(required, extra_fields))), by = keys)

blank <- function(x) is.na(x) | !nzchar(trimws(x))
completed <- reviewed$review_status == "reviewed"
if (any(completed & (blank(reviewed$evidence_url) |
                      blank(reviewed$evidence_excerpt) |
                      blank(reviewed$reviewed_by) |
                      blank(reviewed$reviewed_on) |
                      blank(reviewed$review_notes)))) {
  stop("Reviewed rows need evidence URL, excerpt, reviewer, date and rationale.")
}
if (any(!blank(reviewed$reviewed_on) &
        !str_detect(coalesce(reviewed$reviewed_on, ""), "^\\d{4}-\\d{2}-\\d{2}$"))) {
  stop("Use YYYY-MM-DD for reviewed_on.")
}
if (any((!completed & reviewed$decision != "pending") |
        (completed & reviewed$decision == "pending"))) {
  stop("Pending rows need pending decisions; reviewed rows need a decision.")
}
primary <- reviewed$decision == "include_primary"
investigational <- reviewed$decision == "include_investigational"
if (any(primary & (reviewed$use_role != "cancer_directed_treatment" |
                   reviewed$regulatory_status != "approved_for_this_cancer" |
                   blank(reviewed$jurisdiction)))) {
  stop("Primary seeds require cancer-directed use, cancer-specific approval and jurisdiction.")
}
if (any(investigational & (reviewed$use_role != "investigational_treatment" |
                          reviewed$regulatory_status != "investigational"))) {
  stop("Investigational seeds require investigational use and regulatory status.")
}
if (any(reviewed$decision == "include_documented" &
        reviewed$use_role != "cancer_directed_treatment")) {
  stop("Documented-use seeds require a cancer-directed treatment association.")
}

dir.create("outputs", showWarnings = FALSE)
summary_table <- reviewed |>
  count(cancer_id, cancer_name, review_status, decision, name = "associations")
write_csv(summary_table, "outputs/cancer_drug_review_summary.csv")
write_csv(reviewed |> filter(decision == "include_primary"),
          "data/curated/cancer_treatment_seeds.csv", na = "")
write_csv(reviewed |> filter(decision == "include_investigational"),
          "data/curated/cancer_investigational_seeds.csv", na = "")
write_csv(reviewed |> filter(decision == "include_documented"),
          "data/curated/cancer_documented_use_seeds.csv", na = "")
write_csv(reviewed |> filter(decision == "defer"),
          "data/curated/cancer_drug_deferred.csv", na = "")
write_csv(reviewed |> filter(decision == "exclude"),
          "data/curated/cancer_drug_excluded.csv", na = "")

instructions <- c(
  "SCRIPT 06 REVIEW WORKSHEET",
  paste("Edit:", review_path),
  "Preserve every drugbank_id/cancer_id pair and all original evidence.",
  "Review the FULL indication; review hints are not final classifications.",
  "A DrugBank approved group does not establish approval for this cancer.",
  "Do not infer eEF1A binding or inhibition from a cancer indication.",
  "One row may have several uses: record the cancer-directed use if supported,",
  "and explain additional uses or uncertainty in review_notes.",
  paste("use_role:", paste(roles, collapse = ", ")),
  paste("regulatory_status:", paste(regulatory, collapse = ", ")),
  paste("decision:", paste(decisions, collapse = ", ")),
  "For reviewed rows, fill evidence_url, evidence_excerpt (brief supporting text),",
  "reviewed_by, reviewed_on (YYYY-MM-DD), review_notes (decision rationale).",
  "For primary seeds, verify cancer-specific approval from a regulatory source",
  "and record its jurisdiction (e.g. US or UK). Keep investigations separate.",
  "include_documented retains cancer treatment use without verified current",
  "cancer-specific approval; it is separate from the primary seed set.",
  "Use defer for a completed review that cannot establish a reliable disposition.",
  "Deferred records remain unresolved, not biological negatives or accepted seeds.",
  "Investigational entries may be historical/preclinical, not active clinical trials.",
  "Mark review_status as reviewed only after checking the evidence.",
  "Rerun source('scripts/06_review_cancer_drugs.R') after saving your edits.",
  "Reruns preserve the worksheet and regenerate validated seed exports.",
  "Zero-row seed exports mean no candidates are accepted yet, not no treatments exist."
)
writeLines(instructions, "outputs/cancer_drug_review_instructions.txt")

print(summary_table)
message("Pending associations: ", sum(!completed))
message("Primary seed associations: ", sum(primary))
message("Investigational seed associations: ", sum(investigational))
message("Documented-use associations: ", sum(reviewed$decision == "include_documented"))
message("Deferred associations: ", sum(reviewed$decision == "defer"))
message("Excluded associations: ", sum(reviewed$decision == "exclude"))

# Log only after the script's validation and exports succeed.
log_entry <- c(
  "", paste("SCRIPT 06 RUN -", format(Sys.time(), tz = "Europe/Rome", usetz = TRUE)),
  "Script: scripts/06_review_cancer_drugs.R",
  paste("Input:", input_path),
  paste("Worksheet (preserved on reruns):", review_path),
  paste("Associations:", nrow(reviewed), "| Pending:", sum(!completed),
        "| Primary:", sum(primary), "| Investigational:", sum(investigational)),
  "Outputs: outputs/cancer_drug_review_summary.csv;",
  "outputs/cancer_drug_review_instructions.txt;",
  "data/curated/cancer_treatment_seeds.csv;",
  "data/curated/cancer_investigational_seeds.csv.",
  "Additional outputs: cancer_documented_use_seeds.csv, cancer_drug_deferred.csv,",
  "cancer_drug_excluded.csv (all under data/curated/).",
  paste("Documented-use:", sum(reviewed$decision == "include_documented"),
        "| Deferred:", sum(reviewed$decision == "defer"),
        "| Excluded:", sum(reviewed$decision == "exclude")),
  paste("Result explanation: primary seeds represent", sum(primary),
        "accepted drug-cancer associations across",
        n_distinct(reviewed$drugbank_id[primary]), "unique drugs."),
  "The summary measures review progress; pending associations are unresolved,",
  "not negative results. Empty seed files mean none accepted in that category.",
  "Review hints are automated; treatment and approval decisions require evidence.",
  "No eEF1A binding assessment has been performed."
)
cat(paste0(paste(log_entry, collapse = "\n"), "\n"),
    file = notes_path, append = TRUE)
message("Run recorded in ", notes_path)
