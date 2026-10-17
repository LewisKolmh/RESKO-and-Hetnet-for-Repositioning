# scripts/03_register_drugbank.R
# Check the DrugBank export and register its location.
# This reads only the first five lines, not the full database.

library(readr)
library(dplyr)
library(stringr)

if (!all(dir.exists(c("scripts", "config", "data")))) {
  stop("Run this script from the project root containing scripts/, config/ and data/.")
}

drugbank_path <- "/Users/lewis/Downloads/DrugBank full database.xml"
manifest_path <- "config/database_manifest.csv"

if (!file.exists(drugbank_path)) {
  stop("DrugBank XML file not found: ", drugbank_path)
}

if (!file.exists(manifest_path)) {
  stop("Database manifest not found. Run script 02 first.")
}

# Inspect the header without loading the full XML.
header <- readLines(
  drugbank_path,
  n = 5,
  warn = FALSE,
  encoding = "UTF-8"
)

root_line <- header[
  str_detect(header, "^\\s*<drugbank\\b")
]

if (length(root_line) != 1L) {
  stop("Expected DrugBank root element was not found.")
}

xml_version <- str_match(
  root_line,
  '\\bversion="([^"]+)"'
)[1, 2]

drugbank_export_date <- str_match(
  root_line,
  '\\bexported-on="([^"]+)"'
)[1, 2]

if (is.na(xml_version) || is.na(drugbank_export_date)) {
  stop("Could not extract the XML version or export date.")
}

file_details <- file.info(drugbank_path)

manifest <- read_csv(
  manifest_path,
  col_types = cols(.default = col_character())
)

if (sum(manifest$source_id == "drugbank", na.rm = TRUE) != 1L) {
  stop("Expected exactly one DrugBank row in the manifest.")
}

# Add provenance fields if they do not already exist.
for (column in c("export_date", "file_bytes")) {
  if (!column %in% names(manifest)) {
    manifest[[column]] <- NA_character_
  }
}

manifest <- manifest |>
  mutate(
    status = if_else(
      source_id == "drugbank",
      "local_file_registered",
      status
    ),
    release_version = if_else(
      source_id == "drugbank",
      paste0("XML header version ", xml_version),
      release_version
    ),
    export_date = if_else(
      source_id == "drugbank",
     .env$drugbank_export_date,
     .data$export_date
    ),
    file_bytes = if_else(
      source_id == "drugbank",
      as.character(file_details$size),
      file_bytes
    ),
    local_file = if_else(
      source_id == "drugbank",
      normalizePath(drugbank_path, mustWork = TRUE),
      local_file
    )
  )

# retrieved_date stays blank: we do not know when you downloaded it.
# Licence details also remain unverified.
write_csv(manifest, manifest_path, na = "")

print(
  manifest |>
    filter(source_id == "drugbank") |>
    select(
      source_id,
      status,
      release_version,
      export_date,
      file_bytes,
      local_file
    ),
  width = Inf
)

message("DrugBank export registered successfully.")