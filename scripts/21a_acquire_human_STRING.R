# Acquire a version-pinned complete human functional protein-association network.
library(readr);library(dplyr);library(tibble);library(httr2)
root<-'data/raw/string/v12.5';dir.create(root,recursive=TRUE,showWarnings=FALSE)
version<-'12.5';kinds<-c('protein.links.detailed','protein.info','protein.aliases');sources<-list()
for(kind in kinds){
 name<-paste0('9606.',kind,'.v',version,'.txt.gz');path<-file.path(root,name)
 url<-paste0('https://stringdb-downloads.org/download/',kind,'.v',version,'/',name)
 meta_path<-paste0(path,'.metadata.rds')
 if(!file.exists(path)){
  message('Downloading human STRING ',kind)
  response<-request(url)|>req_timeout(300)|>req_retry(max_tries=3,retry_on_failure=TRUE)|>req_perform(path=paste0(path,'.partial'))
  stopifnot(resp_status(response)==200L,file.rename(paste0(path,'.partial'),path))
  saveRDS(list(url=url,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),last_modified=resp_header(response,'last-modified'),
   file_bytes=file.info(path)$size,md5=unname(tools::md5sum(path))),meta_path)
 }
 stopifnot(file.exists(meta_path));meta<-readRDS(meta_path)
 stopifnot(identical(meta$url,url),identical(meta$md5,unname(tools::md5sum(path))))
 con<-file(path,'rb');magic<-readBin(con,'raw',n=2L);close(con);stopifnot(identical(magic,as.raw(c(31,139))))
 sources[[length(sources)+1L]]<-tibble(source_kind=kind,source_release=version,taxon=9606L,url,
 local_file=path,file_bytes=meta$file_bytes,md5=meta$md5,retrieved_utc=meta$retrieved_utc,last_modified=meta$last_modified)
}
sources<-bind_rows(sources)
write_csv(sources,file.path(root,'acquisition_manifest.csv'),na='')
network_path<-sources$local_file[sources$source_kind=='protein.links.detailed']
info_path<-sources$local_file[sources$source_kind=='protein.info']
alias_path<-sources$local_file[sources$source_kind=='protein.aliases']
# Preserve every published channel; do not request a capped API neighbourhood.
network<-read_table(network_path,col_types=cols(protein1=col_character(),protein2=col_character(),.default=col_integer()),progress=FALSE)
stopifnot(nrow(problems(network))==0L,all(c('protein1','protein2','combined_score') %in% names(network)))
score_cols<-setdiff(names(network),c('protein1','protein2'))
stopifnot(all(grepl('^9606[.]',network$protein1)),all(grepl('^9606[.]',network$protein2)),
 all(vapply(network[score_cols],function(v)all(!is.na(v)&v>=0L&v<=1000L),logical(1))),
 !any(network$protein1==network$protein2),!anyDuplicated(network[c('protein1','protein2')]))
# Symmetric functional-network edges are exported in both directions. Validate, then retain one copy.
forward<-network |> filter(protein1<protein2) |> arrange(protein1,protein2)
reverse<-network |> filter(protein1>protein2) |> transmute(protein1=protein2,protein2=network$protein1[network$protein1>network$protein2],across(all_of(score_cols))) |>
 arrange(protein1,protein2)
stopifnot(nrow(network)==2L*nrow(forward),identical(as.data.frame(forward),as.data.frame(reverse)))
info<-read_tsv(info_path,col_types=cols(.default=col_character()),quote='',progress=FALSE)
stopifnot(nrow(problems(info))==0L,names(info)[1]=='#string_protein_id')
names(info)[1]<-'protein';stopifnot(!anyDuplicated(info$protein),all(grepl('^9606[.]',info$protein)),
 all(network$protein1 %in% info$protein),all(network$protein2 %in% info$protein))
aliases<-read_tsv(alias_path,col_types=cols(.default=col_character()),quote='',progress=FALSE)
stopifnot(nrow(problems(aliases))==0L,ncol(aliases)==3L)
names(aliases)<-c('protein','alias','alias_sources')
stopifnot(all(aliases$protein %in% info$protein))
# RDS retains integer channels without enormous duplicate CSV exports.
saveRDS(list(source_release=version,taxon=9606L,network_type='functional_association_not_necessarily_physical_binding',
 source_manifest=sources,source_directed_rows=nrow(network),edges=forward),
 'data/curated/STRING_human_functional_network_v12.5.rds',compress='gzip')
write_csv(info,'data/curated/STRING_human_protein_registry_v12.5.csv',na='')
# Alias raw file remains local; later mapping reads it with explicit source restrictions.
channels<-bind_rows(lapply(score_cols,function(col)tibble(channel=col,
 source_undirected_edges=nrow(forward),edges_with_channel_support=sum(forward[[col]]>0L),
 minimum_score=min(forward[[col]]),maximum_score=max(forward[[col]]))))
write_csv(channels,'outputs/21a_channel_coverage.csv')
write_csv(sources |> select(-local_file),'outputs/21a_source_summary.csv',na='')
write_csv(tibble(source_release=version,registry_proteins=nrow(info),proteins_with_network_edges=n_distinct(c(forward$protein1,forward$protein2)),
 directed_source_rows=nrow(network),undirected_functional_edges=nrow(forward),alias_rows=nrow(aliases),
 minimum_published_combined_score=min(forward$combined_score)), 'outputs/21a_network_summary.csv')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
manifest<-r('config/database_manifest.csv');i<-manifest$source_id=='string';stopifnot(sum(i)==1L)
manifest$status[i]<-'complete_human_v12.5_functional_network_registry_aliases_acquired'
manifest$release_version[i]<-'12.5';manifest$retrieved_date[i]<-'2026-10-20'
manifest$download_url[i]<-'https://string-db.org/cgi/download';manifest$local_file[i]<-file.path(root,'acquisition_manifest.csv')
write_csv(manifest,'config/database_manifest.csv',na='')
report<-c('SCRIPT 21a - HUMAN STRING FUNCTIONAL NETWORK ACQUISITION',
 paste('Pinned source release:',version),paste('Registry proteins:',nrow(info)),
 paste('Directed source edge rows:',nrow(network)),paste('Unique undirected edges after reciprocal validation:',nrow(forward)),
 paste('Alias rows acquired:',nrow(aliases)),paste('Published minimum combined score:',min(network$combined_score)),
 'Complete human detailed network acquired; no API neighbour cap and no added confidence cutoff.',
 'Every published score channel retained on original 0..1000 scale. Symmetric reciprocal rows agree exactly on all channels.',
 'Functional protein associations are not necessarily physical interactions or evidence of compound binding.',
 'This is a version-pinned current-source adaptation; legacy STITCH/STRING10 identifiers are not silently treated as current nodes.',
 'Raw files, aliases and full source-rich network remain local; compact source, channel and network reports in outputs/21a_*.',
 'Next: exact-accession and legacy-to-current protein harmonisation with missing/ambiguous mappings explicit.')
writeLines(report,'outputs/21a_acquisition_summary.txt')
writeLines(c('PASS: official gzip downloads and source hashes verified; full network, registry and aliases parse without problems.',
 'PASS: human endpoints, score ranges, no self-edges/duplicate directed edges, exact reciprocal channel equality and registry referential integrity.'),
 'outputs/21a_acquisition_validation.txt')
cat(paste0('\nSCRIPT 21a RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/21a_acquire_human_STRING.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
