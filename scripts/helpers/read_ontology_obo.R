# OBO parser retaining primary/alternative IDs, obsolete flags, all edge types and xrefs.
read_ontology_obo<-function(path,prefix){
 lines<-readLines(path,warn=FALSE,encoding='UTF-8');starts<-grep('^\\[',lines)
 stopifnot(length(starts)>0L);ends<-c(starts[-1L]-1L,length(lines))
 terms<-list();edges<-list();aliases<-list();xrefs<-list();synonyms<-list()
 values<-function(block,key)sub(paste0('^',key,': '),'',grep(paste0('^',key,': '),block,value=TRUE))
 one<-function(block,key){v<-values(block,key);if(length(v))v[1]else NA_character_}
 for(i in seq_along(starts)){
  if(lines[starts[i]]!='[Term]')next
  b<-lines[starts[i]:ends[i]];id<-one(b,'id');if(is.na(id)||!startsWith(id,paste0(prefix,':')))next
  terms[[length(terms)+1L]]<-tibble::tibble(term_id=id,name=one(b,'name'),namespace=one(b,'namespace'),
   is_obsolete=identical(one(b,'is_obsolete'),'true'),replaced_by=paste(values(b,'replaced_by'),collapse='|'),
   consider=paste(values(b,'consider'),collapse='|'))
  for(v in values(b,'alt_id'))aliases[[length(aliases)+1L]]<-tibble::tibble(alt_id=v,term_id=id)
  for(v in values(b,'xref'))xrefs[[length(xrefs)+1L]]<-tibble::tibble(term_id=id,xref=sub(' .*$', '',v))
  for(v in values(b,'synonym'))synonyms[[length(synonyms)+1L]]<-tibble::tibble(term_id=id,synonym_record=v)
  for(v in values(b,'is_a'))edges[[length(edges)+1L]]<-tibble::tibble(child=id,parent=sub(' .*$', '',v),relation='is_a')
  for(v in values(b,'relationship')){
   bits<-strsplit(v,' +')[[1]]
   if(length(bits)>=2L)edges[[length(edges)+1L]]<-tibble::tibble(child=id,parent=bits[2],relation=bits[1])
  }
 }
 term_table<-dplyr::bind_rows(terms);stopifnot(!anyDuplicated(term_table$term_id))
 list(terms=term_table,edges=dplyr::bind_rows(edges),aliases=dplyr::bind_rows(aliases),
 xrefs=dplyr::bind_rows(xrefs),synonyms=dplyr::bind_rows(synonyms),
 header=lines[seq_len(starts[1]-1L)],data_version=one(lines[seq_len(starts[1]-1L)],'data-version'))
}
resolve_ontology_ids<-function(ids,ontology){
 aliases<-dplyr::distinct(ontology$aliases,alt_id,term_id)
 ambiguous<-unique(aliases$alt_id[duplicated(aliases$alt_id)|duplicated(aliases$alt_id,fromLast=TRUE)])
 aliases<-aliases[!aliases$alt_id %in% ambiguous,,drop=FALSE]
 known<-ids %in% ontology$terms$term_id;pos<-match(ids,aliases$alt_id)
 canonical<-ifelse(known,ids,aliases$term_id[pos])
 pos<-match(canonical,ontology$terms$term_id)
 tibble::tibble(source_term_id=ids,term_id=canonical,
 term_status=ifelse(!known & ids %in% ambiguous,'ambiguous_alternative_ID_hold',
 ifelse(is.na(pos),'unknown_term_hold',ifelse(ontology$terms$is_obsolete[pos],
 'obsolete_term_hold',ifelse(known,'active_primary_ID','active_documented_alternative_ID')))))
}
