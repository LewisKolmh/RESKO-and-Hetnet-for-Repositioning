# 05: Find provisional cancer-drug associations in DrugBank indications.
# Run from the project root with source("scripts/05_identify_cancer_drugs.R").
# No association is automatically accepted as a treatment or approval claim.

library(readr)
library(dplyr)
library(stringr)
library(tibble)

if (!all(dir.exists(c("scripts", "config", "data")))) {
  stop("Run this script from the project root containing scripts/, config/ and data/.")
}

drugs <- read_csv(
  "data/curated/drugbank_drugs.csv",
  col_types = cols(.default = col_character())
)
cancers <- read_csv("config/cancers.csv", show_col_types = FALSE) |>
  filter(include)
terms <- read_csv(
  "config/cancer_search_terms.csv",
  col_types = cols(.default = col_character())
) |>
  semi_join(cancers, by = "cancer_id")

stopifnot(
  nrow(drugs) > 0L,
  nrow(cancers) > 0L,
  !anyNA(drugs$drugbank_id),
  !anyDuplicated(drugs$drugbank_id),
  !anyDuplicated(cancers$cancer_id),
  !anyDuplicated(terms$term_id),
  all(cancers$cancer_id %in% terms$cancer_id)
)

output_path <- "data/curated/cancer_drug_candidates_for_review.csv"
if (file.exists(output_path)) {
  stop("Review file already exists. Preserve any manual decisions before rerunning.")
}

# Search indication text only. Mentions in descriptions or mechanisms are
# not evidence of a treatment indication and do not enter this search.
matches <- lapply(seq_len(nrow(terms)), function(i) {
  pattern <- regex(terms$pattern[[i]], ignore_case = TRUE)
  hits <- drugs |>
    filter(str_detect(coalesce(indication, ""), pattern)) |>
    transmute(
      cancer_id = terms$cancer_id[[i]],
      drugbank_id,
      term_id = terms$term_id[[i]],
      matched_text = str_extract(indication, pattern)
    )
  hits
}) |>
  bind_rows()

# Several search patterns may find the same association. Collapse these
# before joining, so there is exactly one row per drug and cancer.
matches <- matches |>
  group_by(cancer_id, drugbank_id) |>
  summarise(
    matched_terms = paste(sort(unique(term_id)), collapse = " | "),
    matched_text = paste(unique(matched_text), collapse = " | "),
    .groups = "drop"
  )

candidates <- matches |>
  left_join(cancers |> select(cancer_id, cancer_name), by = "cancer_id") |>
  left_join(drugs, by = "drugbank_id") |>
  mutate(
    source_url = paste0("https://go.drugbank.com/drugs/", drugbank_id),
    review_status = "pending",
    # These flags help review full indications; they are not classifications
    # and may refer to another disease mentioned in the same indication.
    prevention_language_present = str_detect(
      indication, regex("\\b(prevent\\w*|prophyl\\w*|risk reduction)\\b", TRUE)
    ),
    supportive_language_present = str_detect(
      indication,
      regex("\\b(nausea|vomiting|neutropenia|anemia|anaemia|hypercalcemia|hypercalcaemia|skeletal|bone pain)\\b", TRUE)
    ),
    review_notes = NA_character_
  ) |>
  select(
    cancer_id, cancer_name, drugbank_id, drug_name, drug_type, groups,
    matched_terms, matched_text, indication, source_url, source_export_date,
    prevention_language_present, supportive_language_present,
    review_status, review_notes
  ) |>
  arrange(cancer_name, drug_name)

stopifnot(
  nrow(candidates) == nrow(distinct(candidates, cancer_id, drugbank_id)),
  !anyNA(candidates$drug_name)
)

summary_table <- cancers |>
  select(cancer_id, cancer_name) |>
  left_join(
    candidates |> count(cancer_id, name = "candidate_drugs"),
    by = "cancer_id"
  ) |>
  mutate(candidate_drugs = coalesce(candidate_drugs, 0L))

shared_drugs <- candidates |>
  group_by(drugbank_id, drug_name) |>
  summarise(
    cancer_count = n_distinct(cancer_id),
    cancers = paste(sort(unique(cancer_name)), collapse = " | "),
    .groups = "drop"
  ) |>
  filter(cancer_count > 1L) |>
  arrange(desc(cancer_count), drug_name)

dir.create("outputs", showWarnings = FALSE)
write_csv(candidates, output_path, na = "")
write_csv(summary_table, "outputs/cancer_drug_match_counts.csv")
write_csv(shared_drugs, "outputs/drugs_matching_multiple_cancers.csv")

print(summary_table)
message("Unique candidate drugs: ", n_distinct(candidates$drugbank_id))
message("Drug-cancer associations to review: ", nrow(candidates))
message("Drugs matching multiple cancers: ", nrow(shared_drugs))
message("Review table saved: ", output_path)
