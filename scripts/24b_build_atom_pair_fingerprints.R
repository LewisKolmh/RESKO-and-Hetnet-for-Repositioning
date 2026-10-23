# All workflow analysis is R; R packages may use native compiled routines.
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr)})
source('scripts/helpers/prepare_R_chemistry_runtime.R')
source('scripts/helpers/RESKO_atom_pair_chemistry.R')
prepare_R_chemistry_runtime(install_missing=TRUE)
suppressPackageStartupMessages(library(ChemmineR))
if(as.character(packageVersion('ChemmineR'))!='3.64.0')
 stop('Runtime version changed; review and version the method before proceeding.')
inputs <- c('data/curated/RESKO24_ChEMBL_structure_records.rds',
 'data/curated/RESKO24_chemistry_entity_coverage.csv',
 'data/curated/RESKO24_entity_ChEMBL_structure_links.csv',
 'data/curated/RESKO24_proposal_cancer_cohort_context.csv',
 'config/RESKO_chemistry_structure_policy.csv')
stopifnot(all(file.exists(inputs)))
rec <- readRDS(inputs[1]);entities <- read_csv(inputs[2],show_col_types=FALSE)
links <- read_csv(inputs[3],show_col_types=FALSE)
sdf_paths <- paste0('data/raw/chemistry24/sdf/',rec$chembl_id,'.sdf')
stopifnot(all(file.exists(sdf_paths)))
protected <- c(inputs,sdf_paths);hashes <- tools::md5sum(protected)
settings <- tibble(parameter=c('R_runtime','fingerprint','vocabulary','uniquePairs','similarity',
 'addone','fragment_policy','AP_errors','empty_fingerprint','entity_projection','interpretation'),
 value=c('ChemmineR 3.64.0; Bioconductor 3.23; no Python or Java backend',
 'binary 2048-bit ChemmineR atom-pair vocabulary fingerprint',
 'first 2048 built-in apfp descriptors; frozen order, not selected from this cohort',
 'TRUE; preserve ChemmineR descriptor encoding including its multiplicity conventions',
 'shared set bits / union set bits','0; no smoothing',
 'hold disconnected graphs; do not strip salts/fragments or invent cross-fragment distances',
 'hold parse/graph/AP failures, sentinel returns and warnings',
 'hold all-zero bits; undefined empty/empty comparison is NA',
 'keep all 204 entities; missing/held similarities NA; same-key labels remain separate',
 '2D descriptors only; encoding choice not bit-for-bit original-paper reproduction'))
write_csv(settings,'config/RESKO_atom_pair_fingerprint_settings.csv')
vocab_env <- new.env();data('apfp',package='ChemmineR',envir=vocab_env)
vocabulary <- as.character(vocab_env$apfp$AP)[seq_len(2048L)]
stopifnot(length(vocabulary)==2048L,!anyDuplicated(vocabulary))
write_csv(tibble(bit=seq_along(vocabulary),descriptor=vocabulary),
 'data/curated/RESKO24_atom_pair_vocabulary.csv')
qc <- list();ap_list <- list()
for(i in seq_len(nrow(rec))) {
 id <- rec$chembl_id[i];message('24b parse/atom pairs ',i,'/',nrow(rec),': ',id)
 warnings <- character();stage <- 'parse';atom_count <- bond_count <- components <- NA_integer_
 result <- tryCatch(withCallingHandlers({
  s <- read.SDFset(sdf_paths[i])
  if(length(s)!=1L||!all(validSDF(s)))stop('Invalid parsed SDF record/counts.')
  atoms <- atomblock(s[[1]]);bonds <- bondblock(s[[1]])
  counts_line <- header(s[[1]])[['Counts_Line']]
  if(!grepl('V2000',counts_line,fixed=TRUE))stop('Undeclared SDF format requires review.')
  expected_atoms <- as.integer(substr(counts_line,1,3))
  expected_bonds <- as.integer(substr(counts_line,4,6))
  if(expected_bonds==0L && all(bonds[,1:2]==0))bonds <- bonds[FALSE,,drop=FALSE]
  atom_count <- nrow(atoms);bond_count <- nrow(bonds)
  if(atom_count!=expected_atoms||bond_count!=expected_bonds)stop('Parsed graph/count-line disagreement.')
  stage <- 'atom_count'
  if(atom_count<2L)stop('Fewer than two atoms for an atom-pair fingerprint.')
  stage <- 'graph';components <- RESKO_component_count(atom_count,bonds)
  if(components!=1L)stop('Disconnected structure: supplied fragments retained under hold.')
  stage <- 'atom_pair';apset <- sdf2ap(s,type='AP',uniquePairs=TRUE)
  pair <- as(apset,'list')[[1]]
  if(!length(pair)||any(!is.finite(pair))||
    (length(pair)==1L&&pair %in% c(0,1)))stop('Empty/sentinel/invalid AP descriptors.')
  if(length(warnings))stop('Parser/atom-pair warning requires review.')
  pair
 },warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
 error=function(e)e)
 if(inherits(result,'error')) {
  status <- paste0(stage,'_hold');reason <- conditionMessage(result);pair_count <- NA_integer_
 } else {
  ap_list[[id]] <- result;status <- 'atom_pairs_available';reason <- NA_character_;pair_count <- length(result)
 }
 qc[[i]] <- tibble(chembl_id=id,standard_inchikey=rec$standard_inchikey[i],
  atom_count=atom_count,bond_count=bond_count,components=components,AP_descriptor_count=pair_count,
  status=status,reason=reason,warnings=paste(unique(warnings),collapse='|'))
}
qc <- bind_rows(qc);stopifnot(length(ap_list)>0L)
bits <- desc2fp(ap_list,descnames=2048,type='matrix')
stopifnot(ncol(bits)==2048L,identical(colnames(bits),vocabulary),all(bits %in% c(0,1)))
features <- tibble(chembl_id=rownames(bits),active_bits=rowSums(bits),
 vocabulary_descriptor_overlap=vapply(ap_list,function(x)sum(vocabulary %in% as.character(x)),integer(1)),
 AP_descriptors_in_vocabulary=vapply(ap_list,function(x)sum(as.character(x) %in% vocabulary),integer(1)))
qc <- qc |> left_join(features,by='chembl_id') |>
 mutate(status=if_else(status=='atom_pairs_available' & active_bits==0,'zero_fingerprint_hold',status),
 status=if_else(status=='atom_pairs_available','fingerprint_available',status))
keep <- qc$chembl_id[qc$status=='fingerprint_available'];bits <- bits[keep,,drop=FALSE]
result <- RESKO_binary_Tanimoto(bits)
package_scores <- t(vapply(seq_len(nrow(bits)),function(i)
 as.numeric(fpSim(bits[i,],bits,sorted=FALSE,method='Tanimoto',addone=0)),numeric(nrow(bits))))
stopifnot(isTRUE(all.equal(unname(result$similarity),unname(package_scores),tolerance=1e-12)))
projection <- links |> select(entity_id,chembl_id,expected_inchikey) |> distinct() |>
 left_join(qc |> select(chembl_id,fingerprint_status=status),by='chembl_id')
entity_qc <- lapply(seq_len(nrow(entities)),function(i) {
 e <- entities[i,];p <- projection |> filter(entity_id==e$entity_id)
 ids <- p$chembl_id[p$fingerprint_status=='fingerprint_available'];representative <- NA_character_
 state <- if(nrow(p)==0L)e$structure_status else 'structure_fingerprint_hold'
 if(length(ids)) {
  all_ready <- length(ids)==nrow(p)
  agree <- all_ready && all(vapply(ids,function(id)identical(as.numeric(bits[id,]),
   as.numeric(bits[ids[1],])),logical(1)))
  if(agree) {representative <- sort(ids)[1];state <- 'fingerprint_available'} else
   state <- 'multiple_record_or_partial_fingerprint_hold'
 }
 tibble(entity_id=e$entity_id,scope=e$scope,entity_status=state,
 representative_chembl_id=representative,linked_structure_records=nrow(p))
}) |> bind_rows()
entity_matrix <- matrix(NA_real_,nrow(entities),nrow(entities),
 dimnames=list(entities$entity_id,entities$entity_id))
ready <- entity_qc$entity_status=='fingerprint_available'
entity_ids <- entity_qc$entity_id[ready];structure_ids <- entity_qc$representative_chembl_id[ready]
entity_matrix[entity_ids,entity_ids] <- result$similarity[structure_ids,structure_ids,drop=FALSE]
write_csv(qc,'data/curated/RESKO24_atom_pair_structure_QC.csv')
write_csv(entity_qc,'data/curated/RESKO24_atom_pair_entity_coverage.csv')
saveRDS(ap_list,'data/curated/RESKO24_atom_pair_descriptors.rds')
saveRDS(bits,'data/curated/RESKO24_atom_pair_2048_bits.rds')
saveRDS(result,'data/curated/RESKO24_structure_Tanimoto.rds')
saveRDS(entity_matrix,'data/curated/RESKO24_entity_Tanimoto.rds')
write_csv(projection,'data/curated/RESKO24_entity_fingerprint_projection.csv')
write_csv(qc |> count(status),'outputs/24b_structure_fingerprint_coverage.csv')
write_csv(entity_qc |> count(scope,entity_status),'outputs/24b_entity_fingerprint_coverage.csv')
write_csv(qc |> filter(status=='fingerprint_available') |>
 summarise(structures=n(),fingerprint_bits=2048L,min_active_bits=min(active_bits),
 max_active_bits=max(active_bits),median_active_bits=median(active_bits),
 min_AP_vocabulary_fraction=min(AP_descriptors_in_vocabulary/AP_descriptor_count),
 median_AP_vocabulary_fraction=median(AP_descriptors_in_vocabulary/AP_descriptor_count)),
 'outputs/24b_fingerprint_summary.csv')
write_csv(tibble(structure_matrix_rows=nrow(bits),entity_matrix_rows=nrow(entity_matrix),
 structure_pair_comparisons=nrow(bits)*(nrow(bits)-1L)/2L,
 entity_rows_with_features=sum(ready),entity_rows_missing=sum(!ready),
 symmetry_max_error=max(abs(result$similarity-t(result$similarity))),
 diagonal_min=min(diag(result$similarity)),diagonal_max=max(diag(result$similarity))),
 'outputs/24b_Tanimoto_validation_summary.csv')
ip <- installed.packages(lib.loc=.libPaths());used <- c('ChemmineR','Rcpp','rsvg','RCurl','rjson','BH','DT')
write_csv(as.data.frame(ip[match(used,ip[,'Package']),c('Package','Version','LibPath'),drop=FALSE]),
 'outputs/24b_runtime_version_summary.csv')
writeLines(capture.output(sessionInfo()),'data/raw/chemistry24/24b_session_info.txt')
stopifnot(identical(hashes,tools::md5sum(protected)))
protected <- c(protected,'config/RESKO_atom_pair_fingerprint_settings.csv',
 'data/curated/RESKO24_atom_pair_vocabulary.csv')
write_csv(tibble(path=protected,md5=unname(tools::md5sum(protected))),
 'outputs/24b_input_hashes_summary.csv')
context <- read_csv(inputs[4],show_col_types=FALSE) |>
 left_join(entity_qc |> select(entity_id,entity_status),by='entity_id')
write_csv(context |> count(historical_name,entity_status),
 'outputs/24b_proposal_fingerprint_coverage.csv')
# Identical 2D bit vectors do not establish chemical identity.
equal_pairs <- which(result$similarity==1 & upper.tri(result$similarity),arr.ind=TRUE)
write_csv(tibble(chembl_id_1=rownames(bits)[equal_pairs[,1]],
 chembl_id_2=rownames(bits)[equal_pairs[,2]]),
 'data/curated/RESKO24_identical_fingerprint_record_pairs.csv')
print(qc |> count(status));print(entity_qc |> count(scope,entity_status),width=Inf)
message('24b complete. Fingerprints and unsmoothed Tanimoto are 2D descriptors, not binding evidence.')
