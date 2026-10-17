# 09: Extract source-reported chemical identities for the EXPANDED seed sets.
# Cross-references are leads for mapping, not independently confirmed identities.
library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(xml2)
if (!all(dir.exists(c("scripts", "config", "data")))) stop("Run from the project root.")
read_table <- function(p) read_csv(p,col_types=cols(.default=col_character()))
review <- read_table("data/curated/cancer_drug_review_all.csv")
accepted <- review |> filter(decision %in% c("include_primary","include_investigational","include_documented"))
seeds <- accepted |> group_by(drugbank_id,drug_name,drug_type) |>
  summarise(seed_tiers=paste(sort(unique(decision)),collapse=" | "),
            cancers=paste(sort(unique(cancer_name)),collapse=" | "),.groups="drop")
stopifnot(!anyDuplicated(seeds$drugbank_id))
manifest <- read_table("config/database_manifest.csv") |> filter(source_id=="drugbank")
stopifnot(nrow(manifest)==1L,file.exists(manifest$local_file))
source_path <- manifest$local_file[[1]]
prefix <- "data/raw/drugbank/09_seed_"
cache_paths <- paste0(prefix,c("identity.csv","properties.csv","external_ids.csv","cache_metadata.csv"))
source_signature <- paste(normalizePath(source_path),file.info(source_path)$size,
                         as.numeric(file.info(source_path)$mtime),sep="|")
target_signature <- paste(sort(seeds$drugbank_id),collapse="|")

extract_selected <- function() {
  con <- file(source_path,open="rt",encoding="UTF-8")
  on.exit(close(con),add=TRUE)
  pending <- character();inside <- FALSE;keep_field <- FALSE
  identities <- list();properties <- list();external <- list();n <- 0L
  wanted <- "^  <(drugbank-id|name|description|cas-number|unii|calculated-properties|external-identifiers)([[:space:]>]|/>)"
  retain <- function(segment) {
    starts <- grep("^  <[a-z][a-z0-9-]*([[:space:]>]|/>)",segment,perl=TRUE)
    states <- c(keep_field,grepl(wanted,segment[starts],perl=TRUE))
    selected <- states[findInterval(seq_along(segment),starts)+1L]
    keep_field <<- tail(states,1L)
    pending <<- c(pending,segment[selected])
  }
  parse_selected <- function() {
    id_line <- grep('^  <drugbank-id primary="true">',pending,value=TRUE)
    if (length(id_line)!=1L) stop("Invalid top-level primary DrugBank ID.")
    id <- sub('^  <drugbank-id primary="true">([^<]+)</drugbank-id>.*$','\\1',id_line)
    if (!id %in% seeds$drugbank_id) return(invisible(NULL))
    doc <- read_xml(paste('<drugbank xmlns="http://www.drugbank.ca">',
                          paste(pending,collapse="\n"),'</drugbank>',sep="\n"))
    xml_ns_strip(doc)
    drug <- xml_find_first(doc,"./drug")
    val <- function(path) xml_text(xml_find_first(drug,path))
    n <<- n+1L
    identities[[n]] <<- tibble(drugbank_id=id,source_drug_name=val("./name"),
      source_drug_type=xml_attr(drug,"type"),cas_number=val("./cas-number"),
      unii=val("./unii"),description=val("./description"))
    p <- xml_find_all(drug,"./calculated-properties/property")
    properties[[n]] <<- tibble(drugbank_id=rep(id,length(p)),
      kind=xml_text(xml_find_first(p,"./kind")),
      value=xml_text(xml_find_first(p,"./value")),
      method=xml_text(xml_find_first(p,"./source")))
    e <- xml_find_all(drug,"./external-identifiers/external-identifier")
    external[[n]] <<- tibble(drugbank_id=rep(id,length(e)),
      resource=xml_text(xml_find_first(e,"./resource")),
      identifier=xml_text(xml_find_first(e,"./identifier")))
    if (n<=3L || n%%25L==0L) message("Extracted chemical identity ",n," / ",nrow(seeds),": ",id)
  }
  repeat {
    chunk <- readLines(con,n=10000,warn=FALSE)
    if (!length(chunk)) break
    starts <- grep("^<drug\\s+type=",chunk,perl=TRUE)
    ends <- grep("^</drug>\\s*$",chunk,perl=TRUE)
    position <- 1L
    for (boundary in sort(c(starts,ends))) {
      if (boundary %in% starts) {
        if (inside) stop("Unexpected drug boundary.")
        inside <- TRUE;pending <- chunk[boundary];keep_field <- FALSE
        position <- boundary+1L
      } else {
        if (!inside) stop("Unmatched drug end.")
        if (position<boundary) retain(chunk[position:(boundary-1L)])
        pending <- c(pending,chunk[boundary]);parse_selected()
        pending <- character();inside <- FALSE;position <- boundary+1L
      }
    }
    if (inside && position<=length(chunk)) retain(chunk[position:length(chunk)])
  }
  if (inside) stop("Incomplete XML drug record.")
  identity <- bind_rows(identities)
  stopifnot(nrow(identity)==nrow(seeds),!anyDuplicated(identity$drugbank_id),
            setequal(identity$drugbank_id,seeds$drugbank_id))
  list(identity=identity,properties=bind_rows(properties),external=bind_rows(external))
}

if (all(file.exists(cache_paths))) {
  metadata <- read_table(cache_paths[4])
  if (nrow(metadata)!=1L || metadata$source_signature!=source_signature ||
      metadata$target_signature!=target_signature) {
    stop("Chemical cache differs from source/seed set. Preserve it and create a new version before extraction.")
  }
  message("Using validated local extraction cache.")
  identity <- read_table(cache_paths[1]);properties <- read_table(cache_paths[2]);external <- read_table(cache_paths[3])
} else {
  if (any(file.exists(cache_paths))) stop("Incomplete chemical extraction cache; inspect before rerunning.")
  message("Streaming XML for ",nrow(seeds)," accepted seed identities...")
  result <- extract_selected()
  identity <- result$identity;properties <- result$properties;external <- result$external
  write_csv(identity,cache_paths[1],na="")
  write_csv(properties,cache_paths[2],na="")
  write_csv(external,cache_paths[3],na="")
  write_csv(tibble(source_signature=source_signature,target_signature=target_signature),cache_paths[4])
}
stopifnot(nrow(identity)==nrow(seeds),setequal(identity$drugbank_id,seeds$drugbank_id))

# A one-to-many property/resource is flagged rather than arbitrarily selected.
single <- function(values) {
  values <- unique(values[!is.na(values) & nzchar(trimws(values))])
  if (length(values)==1L) values else NA_character_
}
property_names <- c(smiles="SMILES",inchi="InChI",inchikey="InChIKey",
                   molecular_formula="Molecular Formula",molecular_weight="Molecular Weight",
                   traditional_iupac_name="Traditional IUPAC Name")
resources <- c(pubchem_cid="PubChem Compound",pubchem_sid="PubChem Substance",chembl_id="ChEMBL")
registry <- seeds |> left_join(identity,by="drugbank_id")
for (field in names(property_names)) {
  values <- properties |> filter(kind==property_names[[field]]) |>
    group_by(drugbank_id) |> summarise(content=single(value),.groups="drop")
  names(values)[2] <- field
  registry <- registry |> left_join(values,by="drugbank_id")
}
for (field in names(resources)) {
  values <- external |> filter(resource==resources[[field]]) |>
    group_by(drugbank_id) |> summarise(content=single(identifier),.groups="drop")
  names(values)[2] <- field
  registry <- registry |> left_join(values,by="drugbank_id")
}
all_drug_names <- read_table("data/curated/drugbank_drugs.csv") |>
  transmute(other_drugbank_id=drugbank_id,other_drug_name=drug_name,
            name_key=str_to_lower(str_squish(drug_name)))
name_collisions <- registry |>
  filter(!is.na(traditional_iupac_name)) |>
  mutate(name_key=str_to_lower(str_squish(traditional_iupac_name))) |>
  inner_join(all_drug_names,by="name_key") |>
  filter(drugbank_id!=other_drugbank_id) |>
  select(drugbank_id,drug_name,traditional_iupac_name,other_drugbank_id,other_drug_name)
collision_notes <- name_collisions |> group_by(drugbank_id) |>
  summarise(source_name_anomaly=paste0(
    "Traditional-name field matches another DrugBank record: ",
    paste(paste0(other_drug_name," [",other_drugbank_id,"]"),collapse=" | "),
    ". Requires independent identity review; do not silently repair."),.groups="drop")
registry <- registry |> left_join(collision_notes,by="drugbank_id") |> mutate(
  molecular_weight_numeric=suppressWarnings(as.numeric(molecular_weight)),
  source_export_date=manifest$export_date[[1]],
  pubchem_cid_format_ok=!is.na(pubchem_cid) & str_detect(pubchem_cid,"^[1-9][0-9]*$"),
  chembl_id_format_ok=!is.na(chembl_id) & str_detect(chembl_id,"^CHEMBL[0-9]+$"),
  inchikey_format_ok=!is.na(inchikey) & str_detect(inchikey,"^[A-Z]{14}-[A-Z]{10}-[A-Z]$"),
  smiles_present=!is.na(smiles) & nzchar(smiles),
  salt_or_mixture_flag=str_detect(coalesce(smiles,""),fixed(".")),
  peptide_or_macrocycle_text_flag=str_detect(coalesce(description,""),regex("peptide|depsipeptide|macrocycl",TRUE)),
  high_molecular_weight_flag=!is.na(molecular_weight_numeric) & molecular_weight_numeric>1000,
  identity_validation_status="source_reported_not_independently_confirmed"
)

resource_conflicts <- external |> filter(!is.na(identifier),nzchar(identifier)) |>
  distinct(drugbank_id,resource,identifier) |> group_by(drugbank_id,resource) |>
  filter(n()>1L) |> ungroup()
shared_identifiers <- external |>
  filter(resource %in% unname(resources),!is.na(identifier),nzchar(identifier)) |>
  distinct(drugbank_id,resource,identifier) |> group_by(resource,identifier) |>
  filter(n_distinct(drugbank_id)>1L) |> ungroup()
property_conflicts <- properties |>
  filter(kind %in% unname(property_names),!is.na(value),nzchar(value)) |>
  distinct(drugbank_id,kind,value) |> group_by(drugbank_id,kind) |>
  filter(n()>1L) |> ungroup()
structure_queue <- read_table("outputs/08_primary_structure_mapping_queue.csv") |>
  distinct(drugbank_id)
primary_chemistry <- registry |> semi_join(structure_queue,by="drugbank_id")
stopifnot(nrow(primary_chemistry)==nrow(structure_queue))

write_csv(registry,"data/curated/seed_chemical_registry.csv",na="")
write_csv(external,"data/curated/seed_external_identifiers_long.csv",na="")
write_csv(primary_chemistry,"data/curated/primary_compound_mapping_candidates.csv",na="")
write_csv(resource_conflicts,"outputs/09_resource_conflicts.csv",na="")
write_csv(shared_identifiers,"outputs/09_shared_external_identifiers.csv",na="")
write_csv(property_conflicts,"outputs/09_structure_property_conflicts.csv",na="")
write_csv(registry |> filter(!is.na(source_name_anomaly)),"outputs/09_source_name_anomalies.csv",na="")
write_csv(name_collisions,"outputs/09_cross_record_name_collisions.csv",na="")
coverage <- bind_rows(
  registry |> mutate(scope="all_accepted_seed_drugs"),
  primary_chemistry |> mutate(scope="primary_structure_queue")
) |> group_by(scope) |> summarise(
  drugs=n(),with_smiles=sum(smiles_present),with_valid_cid_format=sum(pubchem_cid_format_ok),
  with_valid_chembl_format=sum(chembl_id_format_ok),with_valid_inchikey_format=sum(inchikey_format_ok),
  source_name_anomalies=sum(!is.na(source_name_anomaly)),.groups="drop")
write_csv(coverage,"outputs/09_identifier_coverage.csv")
report <- c(
 "SCRIPT 09 - SOURCE CHEMICAL IDENTITY EXTRACTION",capture.output(print(coverage,width=Inf)),"",
 "Result meaning: cross-references and structures are extracted, not independently verified.",
 "CID and SID are separate namespaces; never use a Substance ID as a Compound ID.",
 "Identifier format checks establish syntax only, not molecular identity.",
 "Shared identifiers may reflect active-moiety/formulation relationships or source errors.",
 "Conflicting values are left unselected and reported; no salts are stripped or records merged.",
 "Source traditional-name anomalies are retained and flagged, not silently corrected.",
 "SMILES have not been sanitised or compared to InChI using a chemistry toolkit.",
 "Missing structures in biologics or investigational products are not negative target results.",
 "Peptide/macrocycle and mass flags guide review; they are not exclusion criteria.",
 "ChEMBL/PubChem cross-references in DrugBank do not mean those databases were acquired.",
 "No direct binding, inhibition, allostery, permeability or therapeutic selectivity is established.",
 "Next: independently check identity leads, obtain SIDER and join its compound IDs with explicit coverage reporting."
)
writeLines(report,"outputs/09_identifier_summary.txt")
print(coverage,width=Inf)
cat(paste0("\nSCRIPT 09 RUN - ",format(Sys.time(),tz="Europe/Rome",usetz=TRUE),"\n",
 "Script: scripts/09_extract_seed_identifiers.R\n",
 "Inputs: expanded review, primary structure queue, registered DrugBank XML.\n",
 "Raw extraction cache: data/raw/drugbank/09_seed_*.csv.\n",
 "Outputs: data/curated/seed_chemical_registry.csv; seed_external_identifiers_long.csv;\n",
 "primary_compound_mapping_candidates.csv. Conflict, anomaly and coverage reports: outputs/09_*.\n",
 paste(report,collapse="\n"),"\n"),file="Notes.txt",append=TRUE)
message("Script 09 complete; source identity registry and Notes.txt saved.")
