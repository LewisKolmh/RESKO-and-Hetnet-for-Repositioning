# Cross-check DrugBank structure identities and recover record-linked synonyms.
# Full InChIKey comparisons only: connectivity/name matches do NOT merge entities.
library(readr);library(dplyr);library(stringr);library(tibble);library(httr2);library(jsonlite)
if (!all(dir.exists(c("scripts", "config", "data")))) stop("Run from the project root.")
read_table <- function(p) read_csv(p,col_types=cols(.default=col_character()))
registry <- read_table("data/curated/seed_chemical_registry.csv")
local_aliases <- read_table("data/curated/seed_drugbank_aliases.csv")
cache_dir <- "data/raw/identity_api";dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
request_log <- list()
scalar <- function(x) if(is.null(x)||!length(x)) NA_character_ else as.character(x[[1]])
fetch <- function(url,label) {
  path <- file.path(cache_dir,paste0(label,".rds"))
  if (file.exists(path)) {
    x <- readRDS(path)
    if (!identical(x$url,url)) stop("Cached request changed: ",label)
  } else {
    Sys.sleep(0.25) # Public service requests are deliberately rate limited.
    response <- request(url) |> req_user_agent("RESKO academic identity audit") |>
      req_timeout(45) |> req_retry(max_tries=3,max_seconds=120,retry_on_failure=TRUE) |>
      req_error(is_error=function(resp) FALSE) |> req_perform()
    status <- resp_status(response)
    if (!status %in% c(200L,404L)) stop("API request failed (",status,"): ",url)
    body <- resp_body_json(response,simplifyVector=FALSE)
    x <- list(url=url,status=status,retrieved_utc=format(Sys.time(),tz="UTC",usetz=TRUE),body=body)
    saveRDS(x,path)
  }
  request_log[[length(request_log)+1L]] <<- tibble(url=x$url,status=x$status,
    retrieved_utc=x$retrieved_utc,cache_file=path)
  x
}
chunks <- function(x,n=25L) split(x,ceiling(seq_along(x)/n))
chembl_base <- "https://www.ebi.ac.uk/chembl/api/data/"
pubchem_base <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"

# 1. ChEMBL molecule records include synonyms, structural identity and hierarchy.
chembl_ids <- sort(unique(na.omit(registry$chembl_id)))
molecules <- list()
group <- chunks(chembl_ids,10L)
for (i in seq_along(group)) {
  url <- paste0(chembl_base,"molecule.json?molecule_chembl_id__in=",
                paste(group[[i]],collapse=","),"&limit=",length(group[[i]]))
  r <- fetch(url,paste0("chembl_small_batch_",i))
  if (r$status!=200L) stop("ChEMBL batch unavailable.")
  if (!is.null(r$body$page_meta[["next"]])) stop("Unexpected ChEMBL pagination.")
  molecules <- c(molecules,r$body$molecules)
  message("ChEMBL batch ",i," / ",length(group)," fetched")
}
chembl <- bind_rows(lapply(molecules,function(m) tibble(
  chembl_id=scalar(m$molecule_chembl_id),chembl_pref_name=scalar(m$pref_name),
  chembl_inchikey=scalar(m$molecule_structures$standard_inchi_key),
  chembl_inchi=scalar(m$molecule_structures$standard_inchi),
  chembl_smiles=scalar(m$molecule_structures$canonical_smiles),
  chembl_parent_id=scalar(m$molecule_hierarchy$parent_chembl_id),
  chembl_molecule_type=scalar(m$molecule_type)
)))
stopifnot(!anyDuplicated(chembl$chembl_id))
chembl_synonyms <- bind_rows(lapply(molecules,function(m) {
  syn <- m$molecule_synonyms
  bind_rows(tibble(chembl_id=scalar(m$molecule_chembl_id),alias=scalar(m$pref_name),alias_type="preferred_name"),
    bind_rows(lapply(syn,function(s) tibble(chembl_id=scalar(m$molecule_chembl_id),
      alias=scalar(s$molecule_synonym),alias_type=scalar(s$syn_type)))))
}))
joined <- registry |> left_join(chembl,by="chembl_id")

# 2. PubChem source CID leads, followed by exact-key recovery where needed.
props <- "InChIKey,InChI,MolecularFormula,MolecularWeight,IUPACName,Title"
parse_props <- function(r) {
  if (r$status!=200L) return(tibble(cid=character(),pc_inchikey=character(),
    pc_inchi=character(),pc_formula=character(),pc_weight=character(),
    pc_iupac_name=character(),pc_title=character()))
  bind_rows(lapply(r$body$PropertyTable$Properties,function(p) tibble(
    cid=scalar(p$CID),pc_inchikey=scalar(p$InChIKey),pc_inchi=scalar(p$InChI),
    pc_formula=scalar(p$MolecularFormula),pc_weight=scalar(p$MolecularWeight),
    pc_iupac_name=scalar(p$IUPACName),pc_title=scalar(p$Title))))
}
cid_leads <- sort(unique(registry$pubchem_cid[registry$pubchem_cid_format_ok=="TRUE"]))
groups <- chunks(cid_leads,20L);pc_tables <- list()
for (i in seq_along(groups)) {
  r <- fetch(paste0(pubchem_base,"cid/",paste(groups[[i]],collapse=","),
                     "/property/",props,"/JSON"),paste0("pc_source_cids_",i))
  if (r$status==404L) {
    for (cid in groups[[i]]) pc_tables[[length(pc_tables)+1L]] <- parse_props(
      fetch(paste0(pubchem_base,"cid/",cid,"/property/",props,"/JSON"),paste0("pc_cid_",cid)))
  } else pc_tables[[length(pc_tables)+1L]] <- parse_props(r)
}
pc <- bind_rows(pc_tables) |> distinct(cid,.keep_all=TRUE)
for (i in seq_len(nrow(joined))) {
  key <- joined$inchikey[i];cid <- joined$pubchem_cid[i]
  if (is.na(key) || !grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$",key)) next
  source_matches <- !is.na(cid) && any(pc$cid==cid & pc$pc_inchikey==key,na.rm=TRUE)
  if (!source_matches) {
    r <- fetch(paste0(pubchem_base,"inchikey/",key,"/property/",props,"/JSON"),
               paste0("pc_key_",joined$drugbank_id[i]))
    pc <- bind_rows(pc,parse_props(r)) |> distinct(cid,.keep_all=TRUE)
  }
  if (i%%25L==0L) message("PubChem identity leads processed: ",i," / ",nrow(joined))
}
resolve <- function(key,cid) {
  if (is.na(key)) return(NA_character_)
  matches <- unique(pc$cid[!is.na(pc$pc_inchikey) & pc$pc_inchikey==key])
  if (!is.na(cid) && cid %in% matches) cid else if(length(matches)==1L) matches else NA_character_
}
joined$resolved_pubchem_cid <- vapply(seq_len(nrow(joined)),function(i)
  resolve(joined$inchikey[i],joined$pubchem_cid[i]),character(1))
joined <- joined |> left_join(pc,by=c("resolved_pubchem_cid"="cid")) |>
  mutate(
    chembl_key_status=case_when(
      is.na(inchikey) ~ "no_source_structure_key",
      is.na(chembl_inchikey) ~ "no_ChEMBL_structure_returned",
      inchikey==chembl_inchikey ~ "full_key_match",
      substr(inchikey,1,14)==substr(chembl_inchikey,1,14) ~ "connectivity_only_NOT_merged",
      TRUE ~ "structure_conflict_NOT_merged"
    ),
    pubchem_key_status=case_when(
      is.na(inchikey) ~ "no_source_structure_key",
      !is.na(resolved_pubchem_cid) & inchikey==pc_inchikey ~ "full_key_match",
      TRUE ~ "no_unique_exact_key_resolution"
    ),
    identity_status=case_when(
      chembl_key_status=="full_key_match" & pubchem_key_status=="full_key_match" ~ "full_key_agreement_both_external_records",
      str_detect(chembl_key_status,"NOT_merged") ~ "external_structure_disagreement_requires_review",
      chembl_key_status=="full_key_match" ~ "full_key_agreement_ChEMBL_only",
      pubchem_key_status=="full_key_match" ~ "full_key_agreement_PubChem_only",
      TRUE ~ "unresolved_structure_identity"
    ),
    cid_recovered=!is.na(resolved_pubchem_cid) & is.na(pubchem_cid),
    cid_source_disagreement=!is.na(resolved_pubchem_cid) & !is.na(pubchem_cid) & resolved_pubchem_cid!=pubchem_cid,
    source_name_field_status=if_else(is.na(source_name_anomaly),"no_flag_from_09",
      "retain_source_name_discrepancy_separately_from_structure_check")
  )

# 3. Recover names by the checked PubChem records, in batches.
resolved_cids <- sort(unique(na.omit(joined$resolved_pubchem_cid)))
pc_syn <- list();groups <- chunks(resolved_cids,20L)
for (i in seq_along(groups)) {
  r <- fetch(paste0(pubchem_base,"cid/",paste(groups[[i]],collapse=","),"/synonyms/JSON"),
             paste0("pc_synonyms_",i))
  if (r$status!=200L) stop("Synonym batch unavailable; do not silently omit.")
  pc_syn[[i]] <- bind_rows(lapply(r$body$InformationList$Information,function(p)
    tibble(cid=rep(scalar(p$CID),length(p$Synonym)),alias=unlist(p$Synonym,use.names=FALSE))))
}
pc_synonyms <- bind_rows(pc_syn)
external_aliases <- bind_rows(
  joined |> select(seed_drugbank_id=drugbank_id,chembl_id,chembl_key_status) |>
    inner_join(chembl_synonyms,by="chembl_id",relationship="many-to-many") |>
    transmute(seed_drugbank_id,alias,alias_type,source="ChEMBL",source_record=chembl_id,
      record_match_status=chembl_key_status),
  joined |> select(seed_drugbank_id=drugbank_id,resolved_pubchem_cid,pubchem_key_status) |>
    inner_join(pc_synonyms,by=c("resolved_pubchem_cid"="cid"),relationship="many-to-many") |>
    transmute(seed_drugbank_id,alias,alias_type="PubChem_record_synonym",source="PubChem",
      source_record=resolved_pubchem_cid,record_match_status=pubchem_key_status)
) |> filter(!is.na(alias),nzchar(trimws(alias))) |>
  mutate(alias_validation="record-linked label; individual name not sufficient for entity merge")
aliases <- bind_rows(
  local_aliases |> transmute(seed_drugbank_id,alias,alias_type,source,
     source_record=record_drugbank_id,record_match_status=relationship,alias_validation),
  external_aliases
) |> distinct(seed_drugbank_id,source,source_record,alias,alias_type,.keep_all=TRUE) |>
  left_join(joined |> select(seed_drugbank_id=drugbank_id,identity_status),by="seed_drugbank_id")

# Alias names can nominate local records, but only full keys support equivalence.
local_index <- read_table("data/curated/drugbank_structure_identity_index.csv") |>
  mutate(alias_key=str_to_lower(str_squish(record_name)))
alias_candidates <- aliases |> mutate(alias_key=str_to_lower(str_squish(alias))) |>
  inner_join(local_index,by="alias_key",relationship="many-to-many") |>
  filter(seed_drugbank_id!=record_drugbank_id) |>
  left_join(registry |> select(seed_drugbank_id=drugbank_id,seed_inchikey=inchikey),by="seed_drugbank_id") |>
  mutate(match_disposition=if_else(!is.na(seed_inchikey)&!is.na(record_inchikey)&seed_inchikey==record_inchikey,
      "full_key_match_already_captured_by_10a","name_only_or_different_structure_NOT_merged"))

write_csv(joined,"data/curated/seed_identity_crosscheck.csv",na="")
write_csv(aliases,"data/curated/seed_aliases_enriched.csv",na="")
write_csv(chembl,"data/curated/chembl_seed_identity_records.csv",na="")
write_csv(pc,"data/curated/pubchem_seed_identity_records.csv",na="")
write_csv(alias_candidates,"outputs/10b_alias_candidate_record_matches.csv",na="")
write_csv(joined |> filter(identity_status!="full_key_agreement_both_external_records"),
          "outputs/10b_identity_followup_queue.csv",na="")
write_csv(bind_rows(request_log),"outputs/10b_api_provenance.csv")
summary <- joined |> count(identity_status,name="unique_drugs")
write_csv(summary,"outputs/10b_identity_summary.csv")
primary_ids <- read_table("data/curated/primary_compound_mapping_candidates.csv")$drugbank_id
primary_summary <- joined |> filter(drugbank_id %in% primary_ids) |>
  count(identity_status,name="unique_drugs")
write_csv(primary_summary,"outputs/10b_primary_identity_summary.csv")

# Register limited ChEMBL acquisition and supplementary PubChem provenance.
manifest <- read_table("config/database_manifest.csv")
manifest$status[manifest$source_id=="chembl"] <- "partial_identity_API_records_acquired"
manifest$retrieved_date[manifest$source_id=="chembl"] <- format(Sys.time(),tz="Europe/Rome",format="%Y-%m-%d")
manifest$download_url[manifest$source_id=="chembl"] <- chembl_base
manifest$local_file[manifest$source_id=="chembl"] <- "data/curated/chembl_seed_identity_records.csv"
write_csv(manifest,"config/database_manifest.csv",na="")
report <- c(
  "SCRIPT 10b - CROSS-SOURCE IDENTITY AND ALIAS CHECK",capture.output(print(summary)),
  "PRIMARY STRUCTURE QUEUE",capture.output(print(primary_summary)),
  paste("Recovered previously missing source CIDs:",sum(joined$cid_recovered)),
  paste("Resolved CIDs different from source CID:",sum(joined$cid_source_disagreement)),
  paste("Alias rows retained with provenance:",nrow(aliases)),
  "Result meaning: full-key agreement compares DrugBank, ChEMBL and PubChem",
  "record representations; it is not an independent laboratory identity assay.",
  "Database records may share upstream evidence. Source-name errors are retained.",
  "Aliases/brand labels can refer to related salts, formulations or products and",
  "are not individually validated molecular equivalences. No name-only merges occur.",
  "ChEMBL parents and connectivity-only matches remain explicit, not equivalent.",
  "Additional exact-key DrugBank record indications/mechanisms/xrefs are in the",
  "10a enrichment tables and retain their own source IDs; seed decisions are unchanged.",
  "No eEF1A binding, inhibition, allostery or activity data were inferred from names.",
  "Next: inspect primary identity disagreements, then use verified IDs and synonym",
  "leads for SIDER mapping with explicit confidence, coverage and missingness."
)
writeLines(report,"outputs/10b_identity_alias_summary.txt")
cat(paste0("\nSCRIPT 10b RUN - ",format(Sys.time(),tz="Europe/Rome",usetz=TRUE),"\n",
 "Script: scripts/10b_crosscheck_identities.R\n",
 "Inputs: source identity registry and 10a local alias/index tables.\n",
 "Outputs: seed_identity_crosscheck.csv; seed_aliases_enriched.csv;\n",
 "chembl_seed_identity_records.csv; pubchem_seed_identity_records.csv (data/curated/).\n",
 "API snapshots: data/raw/identity_api/*.rds; provenance and reports: outputs/10b_*.\n",
 paste(report,collapse="\n"),"\n"),file="Notes.txt",append=TRUE)
cat(paste(report,collapse="\n"),"\n")
