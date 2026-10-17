# 08: Review the 21 broad-screen mentions and subtype the primary biotech queue.
# Original 05/06 decisions and outputs are preserved; expanded files are separate.
# External indication checks were performed on 2026-10-16, not on every rerun.
library(readr)
library(dplyr)
library(stringr)
library(tibble)
if (!all(dir.exists(c("scripts", "config", "data")))) stop("Run from the project root.")

read_table <- function(path) read_csv(path, col_types=cols(.default=col_character()))
old <- read_table("data/curated/cancer_drug_review_decisions.csv")
gaps <- read_table("outputs/07_additional_indication_mentions.csv")
keys <- c("cancer_id", "drugbank_id")
stopifnot(nrow(gaps)==21L, !anyDuplicated(gaps[keys]),
          nrow(inner_join(old, gaps, by=keys))==0L)

review_path <- "data/curated/audit_mention_review_decisions.csv"
if (!file.exists(review_path)) {
  benign_ids <- c("DB04468", "DB00269", "DB01406", "DB05382", "DB00050",
                  "DB09126", "DB00882", "DB00286", "DB09066", "DB09381",
                  "DB00783", "DB04574", "DB00066", "DB06785", "DB13178")
  allowed_ids <- c(benign_ids, "DB00544", "DB01229", "DB00305", "DB06703", "DB09156")
  stopifnot(all(gaps$drugbank_id %in% allowed_ids))
  decisions <- gaps |>
    mutate(
      review_status="reviewed", decision="exclude", use_role="unrelated_mention",
      regulatory_status="not_applicable", jurisdiction=NA_character_,
      evidence_url=source_url,
      evidence_excerpt=paste0("Supplied DrugBank indication: ", indication),
      evidence_basis="DrugBank_indication", research_stage="non_cancer_mention",
      review_batch="006_audit_extension", reviewed_on="2026-10-16",
      reviewed_by="Codex assistant (audit source review)",
      review_notes=paste(
        "This indication passage describes benign breast disease, normal ovarian",
        "function/failure, fertility or another non-target use. It does not support",
        "treatment of this cancer. This is a passage-specific exclusion, not proof",
        "that the compound has no cancer research or eEF1A activity."
      )
    )
  for (i in seq_len(nrow(decisions))) {
    id <- decisions$drugbank_id[i]
    if (id %in% c("DB06703", "DB09156")) {
      decisions$use_role[i] <- "diagnostic"
      decisions$research_stage[i] <- "diagnostic_use"
      decisions$review_notes[i] <- "The breast indication is MRI/contrast mammography imaging, not cancer-directed treatment."
    }
    if (id %in% c("DB00544", "DB01229", "DB00305")) {
      decisions$use_role[i] <- "cancer_directed_treatment"
      decisions$research_stage[i] <- "documented_treatment"
      decisions$decision[i] <- "include_documented"
      decisions$regulatory_status[i] <- "unverified"
      decisions$review_notes[i] <- paste(
        "Mitomycin's supplied indication includes female breast malignancy.",
        "The NCI summary does not establish this cancer-specific approval; absence",
        "is not proof of nonapproval. Keep documented use separate from primary seeds."
      )
      decisions$evidence_url[i] <- paste(decisions$source_url[i],
        "https://www.cancer.gov/about-cancer/treatment/drugs/mitomycin", sep=" | ")
    }
    if (id=="DB00544") {
      decisions$decision[i] <- "include_primary"
      decisions$regulatory_status[i] <- "approved_for_this_cancer"
      decisions$jurisdiction[i] <- "US"
      decisions$research_stage[i] <- "labelled_treatment"
      decisions$evidence_basis[i] <- "NCI_individual_summary_plus_DrugBank"
      decisions$evidence_url[i] <- "https://www.cancer.gov/about-cancer/treatment/drugs/fluorouracil"
      decisions$review_notes[i] <- "NCI corroborates breast-cancer use of injectable fluorouracil. Do not transfer the indication to topical preparations. The original regex missed the shared cancer list."
    }
    if (id=="DB01229") {
      decisions$decision[i] <- "include_primary"
      decisions$regulatory_status[i] <- "approved_for_this_cancer"
      decisions$jurisdiction[i] <- "US"
      decisions$research_stage[i] <- "labelled_treatment"
      decisions$evidence_basis[i] <- "NCI_individual_summary_plus_DrugBank"
      decisions$evidence_url[i] <- "https://www.cancer.gov/about-cancer/treatment/drugs/paclitaxel"
      decisions$review_notes[i] <- "NCI corroborates advanced ovarian use of conventional paclitaxel with cisplatin first-line or alone after prior treatment. Albumin-bound formulation indications are distinct."
    }
  }
  write_csv(decisions, review_path, na="")
} else {
  message("Preserving existing audit review decisions.")
}

audit <- read_table(review_path)
stopifnot(nrow(audit)==nrow(gaps), !anyDuplicated(audit[keys]),
          nrow(anti_join(gaps,audit,by=keys))==0L,
          nrow(anti_join(audit,gaps,by=keys))==0L,
          all(audit$review_status=="reviewed"),
          all(audit$decision %in% c("include_primary", "include_documented", "exclude")))
# Only validated decision columns are carried back onto immutable source evidence.
decision_fields <- c("review_status", "decision", "use_role", "regulatory_status",
  "jurisdiction", "evidence_url", "evidence_excerpt", "evidence_basis", "research_stage",
  "review_batch", "reviewed_on", "reviewed_by", "review_notes")
for (field in decision_fields) {
  if (any(is.na(audit[[field]]) | !nzchar(trimws(audit[[field]]))) && field!="jurisdiction") {
    stop("Incomplete audit decision field: ", field)
  }
}
if (any(audit$decision=="include_primary" &
        (audit$regulatory_status!="approved_for_this_cancer" | is.na(audit$jurisdiction)))) {
  stop("Primary additions need cancer-specific approval evidence and jurisdiction.")
}
additions <- gaps |>
  select(cancer_id, cancer_name, drugbank_id, drug_name, drug_type, groups,
         matched_text, indication, source_url, source_export_date) |>
  left_join(audit |> select(all_of(c(keys,decision_fields))), by=keys)
expanded <- bind_rows(old, additions) |> arrange(cancer_name,drug_name)
stopifnot(nrow(expanded)==238L, !anyDuplicated(expanded[keys]))

# Extract descriptions for the 20 flagged primary drug IDs only.
# Retain only requested descriptions, including physical line breaks, while
# discarding the large interaction sections. No whole drug record is buffered.
queue <- read_table("outputs/07_primary_modality_review_queue.csv") |>
  distinct(drugbank_id,drug_name)
manifest <- read_table("config/database_manifest.csv")
xml_path <- manifest$local_file[manifest$source_id=="drugbank"]
extract_descriptions <- function() {
  con <- file(xml_path,open="rt",encoding="UTF-8")
  on.exit(close(con),add=TRUE)
  current_id <- NA_character_; result <- list()
  capturing <- FALSE; pending_text <- character()
  finish_description <- function() {
    result[[length(result)+1L]] <<- tibble(
      drugbank_id=current_id,
      description=xml2::xml_text(xml2::read_xml(paste(pending_text,collapse="\n")))
    )
    capturing <<- FALSE
    pending_text <<- character()
  }
  repeat {
    lines <- readLines(con,n=10000,warn=FALSE)
    if (!length(lines)) break
    events <- grep('^  <drugbank-id primary="true">|^  <description[ >]|</description>',lines)
    position <- 1L
    for (j in events) {
      if (capturing) {
        if (position<=j) pending_text <- c(pending_text,lines[position:j])
        if (grepl("</description>",lines[j],fixed=TRUE)) {
          finish_description()
        } else {
          stop("Unexpected element before description closed for ",current_id)
        }
        position <- j+1L
      } else if (startsWith(lines[j], '  <drugbank-id primary="true">')) {
        current_id <- sub('^  <drugbank-id primary="true">([^<]+)</drugbank-id>.*$', '\\1',lines[j])
      } else if (grepl('^  <description[ >]',lines[j]) &&
                 !is.na(current_id) && current_id %in% queue$drugbank_id) {
        pending_text <- lines[j]
        capturing <- TRUE
        if (grepl("</description>",lines[j],fixed=TRUE)) finish_description()
        position <- j+1L
      }
    }
    if (capturing && position<=length(lines)) {
      pending_text <- c(pending_text,lines[position:length(lines)])
    }
    if (length(result)==nrow(queue)) break
  }
  if (capturing) stop("Incomplete source description at end of XML.")
  bind_rows(result)
}
message("Reading source descriptions for the primary modality queue...")
descriptions <- extract_descriptions()
stopifnot(nrow(descriptions)==nrow(queue), !anyDuplicated(descriptions$drugbank_id))
classes <- tribble(
 ~drugbank_id, ~reviewed_modality,
 "DB00112","monoclonal_antibody", "DB14967","monoclonal_antibody",
 "DB09037","monoclonal_antibody", "DB06366","monoclonal_antibody",
 "DB00072","monoclonal_antibody", "DB09331","monoclonal_antibody",
 "DB06317","monoclonal_antibody", "DB14811","monoclonal_antibody",
 "DB15395","bispecific_antibody", "DB17447","bispecific_antibody",
 "DB16678","bispecific_antibody", "DB16655","bispecific_antibody",
 "DB16410","antibody_drug_conjugate", "DB12893","antibody_drug_conjugate",
 "DB14962","antibody_drug_conjugate", "DB05773","antibody_drug_conjugate",
 "DB15719","antibody_drug_conjugate", "DB12489","antibody_drug_conjugate",
 "DB16738","CAR_T_cell_product", "DB16665","CAR_T_cell_product"
)
modality <- queue |> left_join(classes,by="drugbank_id") |>
  left_join(descriptions,by="drugbank_id") |>
  mutate(
    modality_evidence="Description in supplied DrugBank XML export",
    source_url=paste0("https://go.drugbank.com/drugs/",drugbank_id),
    direct_compound_route="retain_as_disease_seed; separate_from_direct_small_molecule_candidates",
    payload_note=if_else(reviewed_modality=="antibody_drug_conjugate",
      "Free payload is a different entity; do not inherit parent indication or target evidence.",
      "No isolated small-molecule entity inferred from this product."),
    reviewed_on="2026-10-16"
  )
stopifnot(!anyNA(modality$reviewed_modality),!anyNA(modality$description))

primary <- expanded |> filter(decision=="include_primary") |>
  left_join(modality |> select(drugbank_id,reviewed_modality,direct_compound_route),by="drugbank_id") |>
  mutate(
    reviewed_modality=coalesce(reviewed_modality,"structure_check_pending"),
    direct_compound_route=coalesce(direct_compound_route,"prioritise_structure_and_identifier_mapping")
  )
stopifnot(nrow(primary)==sum(expanded$decision=="include_primary"))
write_csv(expanded,"data/curated/cancer_drug_review_all.csv",na="")
for (tier in c("primary","investigational","documented")) {
  write_csv(expanded |> filter(decision==paste0("include_",tier)),
    paste0("data/curated/cancer_",tier,"_seeds_expanded.csv"),na="")
}
write_csv(modality,"data/curated/primary_modality_review.csv",na="")
write_csv(primary,"data/curated/primary_seed_associations_routed.csv",na="")
write_csv(primary |> filter(reviewed_modality=="structure_check_pending"),
          "outputs/08_primary_structure_mapping_queue.csv",na="")
summary <- expanded |> group_by(decision) |>
  summarise(associations=n(),unique_drugs=n_distinct(drugbank_id),.groups="drop")
write_csv(summary,"outputs/08_expanded_seed_summary.csv")
write_csv(audit |> count(cancer_name,decision,name="associations"),
          "outputs/08_audit_review_summary.csv")
write_csv(modality |> count(reviewed_modality,name="unique_drugs"),
          "outputs/08_modality_review_summary.csv")
report <- c(
 "SCRIPT 08 - AUDIT EXTENSION AND PRIMARY MODALITY REVIEW",
 "21 audit mentions reviewed: 2 primary additions, 1 documented-use addition, 18 exclusions.",
 "Original 217 decisions remain in their original file; consolidated review has 238 rows.",
 "Expanded seed files, not original 06 exports, are the inputs for subsequent mapping.",
 capture.output(print(summary)), "",
 "20 primary product identities: 8 monoclonal antibodies, 4 bispecific antibodies,",
 "6 antibody-drug conjugates and 2 CAR-T cell products, based on source descriptions.",
 "They remain disease-treatment seeds, but are separated from direct small-molecule candidates.",
 "55 distinct primary drugs remain in the structure/identifier queue; this is not a",
 "validated list of cell-permeable compounds or evidence of eEF1A binding.",
 "The wider investigational/documented tiers still need their own modality/identity checks.",
 "All reviews are evidence-specific; exclusions are not target-activity negatives.",
 "Next: extract stable identifiers/structures from DrugBank, obtain SIDER and audit mapping coverage."
)
writeLines(report,"outputs/08_review_summary.txt")
print(summary)
cat(paste0("\nSCRIPT 08 RUN - ",format(Sys.time(),tz="Europe/Rome",usetz=TRUE),"\n",
 "Script: scripts/08_review_audit_mentions_and_modality.R\n",
 paste(report,collapse="\n"),"\n",
 "Outputs: data/curated/audit_mention_review_decisions.csv; cancer_drug_review_all.csv;\n",
 "cancer_primary_seeds_expanded.csv; cancer_investigational_seeds_expanded.csv;\n",
 "cancer_documented_seeds_expanded.csv; primary_modality_review.csv;\n",
 "primary_seed_associations_routed.csv. Summaries and mapping queue: outputs/08_*.\n"
),file="Notes.txt",append=TRUE)
message("Script 08 complete; expanded outputs saved and Notes.txt updated.")
