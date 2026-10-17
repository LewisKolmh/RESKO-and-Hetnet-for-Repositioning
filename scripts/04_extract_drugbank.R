# scripts/04_extract_drugbank.R
# Extract drug identity, groups and indications.
# Run from the project root.

library(readr)
library(dplyr)
library(tibble)

if (!all(dir.exists(c("scripts", "config", "data")))) {
  stop("Run this script from the project root containing scripts/, config/ and data/.")
}

# Install the XML parser if needed.
if (!requireNamespace("xml2", quietly = TRUE)) {
  install.packages(
    "xml2",
    repos = "https://cloud.r-project.org"
  )
}

manifest <- read_csv(
  "config/database_manifest.csv",
  col_types = cols(.default = col_character())
)

drugbank_record <- manifest |>
  filter(source_id == "drugbank")

if (nrow(drugbank_record) != 1L) {
  stop("Expected exactly one DrugBank manifest record.")
}

input_path <- drugbank_record$local_file[[1]]
output_path <- "data/curated/drugbank_drugs.csv"
temporary_path <- paste0(output_path, ".partial")

if (is.na(input_path) || !file.exists(input_path)) {
  stop("Registered DrugBank XML file was not found.")
}

if (file.exists(output_path) || file.exists(temporary_path)) {
  stop(
    "An extraction output already exists. ",
    "Review it before rerunning: ",
    output_path, " or ", temporary_path
  )
}

# Parse one complete drug record.
parse_drug <- function(lines) {

  # Restore the namespace inherited from the original XML root.
  text <- paste(
    '<drugbank xmlns="http://www.drugbank.ca">',
    paste(lines, collapse = "\n"),
    "</drugbank>",
    sep = "\n"
  )

  document <- xml2::read_xml(text)
  xml2::xml_ns_strip(document)

  drug <- xml2::xml_find_first(document, "./drug")

  value <- function(xpath) {
    xml2::xml_text(xml2::xml_find_first(drug, xpath))
  }

  primary_id <- value("./drugbank-id[@primary='true']")

  if (is.na(primary_id) || !nzchar(primary_id)) {
    stop("A drug record has no primary DrugBank ID.")
  }

  groups <- xml2::xml_text(
    xml2::xml_find_all(drug, "./groups/group")
  )

  tibble(
    drugbank_id = primary_id,
    drug_name = value("./name"),
    drug_type = xml2::xml_attr(drug, "type"),
    groups = paste(groups, collapse = "|"),
    indication = value("./indication"),
    mechanism_of_action = value("./mechanism-of-action"),
    cas_number = value("./cas-number"),
    source_export_date = drugbank_record$export_date[[1]]
  )
}

extract_drugs <- function() {

  connection <- file(input_path, open = "rt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)

  pending <- character()
  inside_drug <- FALSE
  batch <- list()
  record_count <- 0L
  wrote_header <- FALSE

  # This export indents direct children of <drug> by two spaces.
  # Discard unneeded children BEFORE accumulating a drug record: some
  # interaction sections occupy hundreds of MB and millions of lines.
  keep_field <- FALSE
  wanted <- paste0(
    "^  <(drugbank-id|name|groups|indication|mechanism-of-action|cas-number)",
    "([[:space:]>]|/>)"
  )

  retain_fields <- function(segment) {
    child_starts <- grep("^  <[a-z][a-z0-9-]*([[:space:]>]|/>)",
                         segment, perl = TRUE)
    states <- c(keep_field, grepl(wanted, segment[child_starts], perl = TRUE))
    selected <- states[findInterval(seq_along(segment), child_starts) + 1L]
    keep_field <<- tail(states, 1L)
    pending <<- c(pending, segment[selected])
  }

  flush_batch <- function() {
    if (length(batch) == 0L) return(invisible(NULL))

    write_csv(
      bind_rows(batch),
      temporary_path,
      append = wrote_header,
      col_names = !wrote_header,
      na = ""
    )

    wrote_header <<- TRUE
    batch <<- list()
    message("Extracted ", record_count, " drugs")
  }

  repeat {
    chunk <- readLines(connection, n = 10000, warn = FALSE)

    if (length(chunk) == 0L) break

    # Your export puts top-level drug boundaries on their own lines.
    starts <- grep("^<drug\\s+type=", chunk, perl = TRUE)
    ends <- grep("^</drug>\\s*$", chunk, perl = TRUE)
    boundaries <- sort(c(starts, ends))
    position <- 1L

    for (boundary in boundaries) {

      if (boundary %in% starts) {

        if (inside_drug) {
          stop("Unexpected nested top-level drug boundary.")
        }

        inside_drug <- TRUE
        pending <- chunk[boundary]
        keep_field <- FALSE
        position <- boundary + 1L

      } else {

        if (!inside_drug) {
          stop("Found a drug end without a matching start.")
        }

        if (position < boundary) {
          retain_fields(chunk[position:(boundary - 1L)])
        }
        pending <- c(pending, chunk[boundary])
        row <- parse_drug(pending)

        record_count <- record_count + 1L
        batch[[length(batch) + 1L]] <- row

        if (record_count <= 5L || record_count %% 25L == 0L) {
          message("Parsed drug ", record_count, ": ", row$drug_name[[1]])
        }

        pending <- character()
        inside_drug <- FALSE
        position <- boundary + 1L

        if (length(batch) >= 250L) {
          flush_batch()
        }
      }
    }

    if (inside_drug && position <= length(chunk)) {
      retain_fields(chunk[position:length(chunk)])
    }
  }

  if (inside_drug) {
    stop("The XML ended with an incomplete drug record.")
  }

  flush_batch()

  if (record_count == 0L) {
    stop("No drug records were extracted.")
  }

  record_count
}

message("Starting extraction. This may take several minutes.")
record_count <- extract_drugs()

# Check the compact output before giving it its final filename.
drugs <- read_csv(
  temporary_path,
  col_types = cols(.default = col_character())
)

stopifnot(
  nrow(drugs) == record_count,
  !anyNA(drugs$drugbank_id),
  !anyDuplicated(drugs$drugbank_id),
  !anyNA(drugs$drug_name)
)

if (!file.rename(temporary_path, output_path)) {
  stop("Could not rename the validated extraction output.")
}

summary_table <- drugs |>
  count(drug_type, name = "drug_count")

print(summary_table)

message("Total drugs extracted: ", nrow(drugs))
message(
  "Drugs with nonempty indication text: ",
  sum(!is.na(drugs$indication) & nzchar(drugs$indication))
)
message("Saved: ", output_path)
