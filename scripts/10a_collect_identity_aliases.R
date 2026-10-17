# Recover alternative names and additional DrugBank records by FULL InChIKey.
# A name match or shared connectivity block never automatically merges records.
library(readr);library(dplyr);library(stringr);library(tibble);library(xml2)
if (!all(dir.exists(c("scripts", "config", "data")))) stop("Run from the project root.")
source("scripts/helpers/stream_drugbank_fields.R")
read_table <- function(p) read_csv(p,col_types=cols(.default=col_character()))
seeds <- read_table("data/curated/seed_chemical_registry.csv")
drugs <- read_table("data/curated/drugbank_drugs.csv")
manifest <- read_table("config/database_manifest.csv") |> filter(source_id=="drugbank")
path <- manifest$local_file[[1]]
cache <- "data/raw/drugbank/10_identity_alias_cache.rds"
signature <- paste(normalizePath(path),file.info(path)$size,
                   as.numeric(file.info(path)$mtime),paste(sort(seeds$drugbank_id),collapse="|"),sep="|")
if (file.exists(cache)) {
  x <- readRDS(cache)
  if (!identical(signature,x$signature)) stop("Alias cache inputs changed; preserve and version cache first.")
  message("Using local identity/alias cache.")
} else {
  index_rows <- list();alias_rows <- list();xref_rows <- list();n <- 0L;selected <- 0L
  keys <- unique(na.omit(seeds$inchikey))
  one <- function(v) { v <- unique(v[!is.na(v)&nzchar(v)]); if(length(v)==1L) v else NA_character_ }
  process <- function(d) {
    val <- function(p) xml_text(xml_find_first(d,p))
    id <- val("./drugbank-id[@primary='true']");name <- val("./name")
    key <- one(xml_text(xml_find_all(d,"./calculated-properties/property[kind='InChIKey']/value")))
    n <<- n+1L
    index_rows[[n]] <<- tibble(record_drugbank_id=id,record_name=name,record_inchikey=key)
    if (!(id %in% seeds$drugbank_id || (!is.na(key) && key %in% keys))) return(invisible(NULL))
    selected <<- selected+1L
    syn <- xml_text(xml_find_all(d,"./synonyms/synonym"))
    brands <- xml_text(xml_find_all(d,"./international-brands/international-brand/name"))
    alias_rows[[selected]] <<- bind_rows(
      tibble(record_drugbank_id=id,alias=name,alias_type="primary_name"),
      tibble(record_drugbank_id=rep(id,length(syn)),alias=syn,alias_type="DrugBank_synonym"),
      tibble(record_drugbank_id=rep(id,length(brands)),alias=brands,alias_type="brand_label_not_entity_equivalence")
    )
    e <- xml_find_all(d,"./external-identifiers/external-identifier")
    xref_rows[[selected]] <<- tibble(record_drugbank_id=rep(id,length(e)),
      resource=xml_text(xml_find_first(e,"./resource")),identifier=xml_text(xml_find_first(e,"./identifier")))
  }
  message("Indexing full-structure keys and collecting selected names/cross-references...")
  stream_drugbank_fields(path,c("drugbank-id","name","calculated-properties",
    "synonyms","international-brands","external-identifiers"),process)
  x <- list(signature=signature,index=bind_rows(index_rows),
            aliases=bind_rows(alias_rows),external_ids=bind_rows(xref_rows))
  stopifnot(nrow(x$index)==nrow(drugs),!anyDuplicated(x$index$record_drugbank_id))
  saveRDS(x,cache)
}
index <- x$index
exact <- seeds |> filter(!is.na(inchikey),str_detect(inchikey,"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")) |>
  select(seed_drugbank_id=drugbank_id,inchikey) |>
  inner_join(index,by=c("inchikey"="record_inchikey"),relationship="many-to-many") |>
  transmute(seed_drugbank_id,record_drugbank_id,record_name,
    relationship=if_else(seed_drugbank_id==record_drugbank_id,"original_seed_record","source_full_inchikey_match"))
self <- seeds |> transmute(seed_drugbank_id=drugbank_id,record_drugbank_id=drugbank_id,
                           record_name=drug_name,relationship="original_seed_record")
links <- bind_rows(exact,self) |> distinct(seed_drugbank_id,record_drugbank_id,.keep_all=TRUE)
aliases <- links |> inner_join(x$aliases,by="record_drugbank_id",relationship="many-to-many") |>
  filter(!is.na(alias),nzchar(trimws(alias))) |>
  mutate(source="DrugBank supplied XML",alias_validation="source_record_label_not_individually_structure_verified") |>
  distinct(seed_drugbank_id,record_drugbank_id,alias,alias_type,.keep_all=TRUE)
extra_data <- links |> left_join(drugs,by=c("record_drugbank_id"="drugbank_id"))
extra_ids <- links |> inner_join(x$external_ids,by="record_drugbank_id",relationship="many-to-many")
stopifnot(all(seeds$drugbank_id %in% links$seed_drugbank_id))
write_csv(index,"data/curated/drugbank_structure_identity_index.csv",na="")
write_csv(links,"data/curated/seed_equivalent_drugbank_records.csv",na="")
write_csv(aliases,"data/curated/seed_drugbank_aliases.csv",na="")
write_csv(extra_data,"data/curated/seed_equivalent_record_data.csv",na="")
write_csv(extra_ids,"data/curated/seed_equivalent_record_identifiers.csv",na="")
report <- c(
  paste("Indexed DrugBank records:",nrow(index)),
  paste("Seed-to-record links:",nrow(links)),
  paste("Additional exact-key record links:",sum(links$relationship!="original_seed_record")),
  paste("Alias rows with provenance:",nrow(aliases)),
  "Result meaning: full source InChIKeys link exact-structure records, retaining",
  "their names, indications, mechanisms and cross-references as separate evidence rows.",
  "Original seed decisions are not transferred to other indications or formulations.",
  "Brand labels and synonyms are search leads, not proof of equivalent formulations.",
  "No charge/stereo/salt stripping or name-only merges are performed.",
  "Full InChIKey agreement still requires independent-source checking; standard",
  "InChI representations do not preserve every tautomeric representation.",
  "Next: ChEMBL/PubChem structure checks, synonym retrieval and conflict reporting."
)
writeLines(report,"outputs/10a_alias_collection_summary.txt")
cat(paste0("\nSCRIPT 10a RUN - ",format(Sys.time(),tz="Europe/Rome",usetz=TRUE),"\n",
 "Script: scripts/10a_collect_identity_aliases.R\n",
 "Outputs: data/curated/drugbank_structure_identity_index.csv; seed_equivalent_drugbank_records.csv;\n",
 "seed_drugbank_aliases.csv; seed_equivalent_record_data.csv; seed_equivalent_record_identifiers.csv.\n",
 paste(report,collapse="\n"),"\n"),file="Notes.txt",append=TRUE)
cat(paste(report,collapse="\n"),"\n")
