# R graph checks and unsmoothed binary Tanimoto; no external chemistry process.
RESKO_component_count <- function(atom_count,bonds) {
 if(atom_count<1L) return(0L)
 adjacency <- vector('list',atom_count)
 if(nrow(bonds)) for(i in seq_len(nrow(bonds))) {
  a <- as.integer(bonds[i,1]);b <- as.integer(bonds[i,2])
  if(is.na(a)||is.na(b)||a<1L||b<1L||a>atom_count||b>atom_count||a==b)
   stop('Invalid bond endpoints in parsed graph.')
  adjacency[[a]] <- c(adjacency[[a]],b);adjacency[[b]] <- c(adjacency[[b]],a)
 }
 seen <- rep(FALSE,atom_count);components <- 0L
 for(start in seq_len(atom_count)) if(!seen[start]) {
  components <- components+1L;queue <- start;seen[start] <- TRUE
  while(length(queue)) {
   node <- queue[1];queue <- queue[-1]
   neighbours <- adjacency[[node]];next_nodes <- neighbours[!seen[neighbours]]
   if(length(next_nodes)) {seen[next_nodes] <- TRUE;queue <- c(queue,unique(next_nodes))}
  }
 }
 components
}
RESKO_binary_Tanimoto <- function(bits) {
 stopifnot(is.matrix(bits),all(bits %in% c(0L,1L)))
 intersection <- tcrossprod(bits)
 union <- outer(rowSums(bits),rowSums(bits),'+')-intersection
 similarity <- intersection/union;similarity[union==0] <- NA_real_
 dimnames(similarity) <- list(rownames(bits),rownames(bits))
 list(similarity=similarity,intersection=intersection,union=union)
}
