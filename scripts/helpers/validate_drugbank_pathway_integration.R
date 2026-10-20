# Independent full-record spot check, bypassing the selective streamer.
library(readr);library(dplyr);library(tibble);library(xml2)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
manifest<-r('config/database_manifest.csv');path<-manifest$local_file[manifest$source_id=='drugbank']
con<-file(path,'rt',encoding='UTF-8');prefix<-readLines(con,n=110000L,warn=FALSE);close(con)
id_line<-which(prefix=='  <drugbank-id primary="true">DB00014</drugbank-id>')
stopifnot(length(id_line)==1L)
start<-max(which(grepl('^<drug ',prefix) & seq_along(prefix)<id_line))
end<-min(which(prefix=='</drug>' & seq_along(prefix)>id_line))
node<-read_xml(paste(prefix[start:end],collapse='\n'));xml_ns_strip(node)
role_source<-r('data/curated/DrugBank_selected_protein_role_annotations.csv') |> filter(drugbank_id=='DB00014')
refs<-r('data/curated/DrugBank_selected_protein_role_references.csv')
for(role in c('target','enzyme','carrier','transporter')){
 items<-xml_find_all(node,paste0('./',role,'s/',role));local<-role_source |> filter(protein_role==role)
 stopifnot(n_distinct(local$source_role_row_id)==length(items))
 for(j in seq_along(items)){
  item<-items[[j]];key<-paste('DB00014',role,j,sep=':');got<-local |> filter(source_role_row_id==key)
  expected_ids<-xml_attr(xml_find_all(item,'./polypeptide'),'id')
  got_ids<-na.omit(got$polypeptide_id)
  stopifnot(identical(sort(as.character(got_ids)),sort(as.character(expected_ids))),
   all(got$source_interaction_id==xml_text(xml_find_first(item,'./id'))))
  expected_refs<-xml_text(xml_find_all(item,'./references/articles/article/pubmed-id'))
  got_refs<-refs$pubmed_id[refs$source_role_row_id==key]
  stopifnot(identical(sort(as.character(na.omit(got_refs))),sort(expected_refs[nzchar(expected_refs)])))
 }
}
coverage<-r('data/curated/inventory_DrugBank_Reactome_coverage.csv')
human<-r('data/curated/inventory_DrugBank_human_protein_links.csv')
features<-r('data/curated/inventory_DrugBank_Reactome_role_features.csv')
stopifnot(nrow(coverage)==200L,!anyDuplicated(coverage$entity_id),
 all(human$protein_taxon=='9606'),!anyNA(human$uniprot),
 all(human$mapping_scope!='named_biologic_context_only'),
 !any(human$entity_id=='sider:CID100135113'),
 !anyDuplicated(features[c('entity_id','mapping_scope','protein_role','pathway_id','pathway_name')]))
validation<-c('PASS: independent complete XML record DB00014 reproduces every selected role, polypeptide accession and PMID list.',
 'PASS: all 200 inventory entities preserved; strict human accession and blocked-structure restrictions hold; role/pathway features unique.')
cat(paste(validation,collapse='\n'),'\n')
cat(paste0('\nSCRIPT 20b - INDEPENDENT VALIDATION\n',paste(validation,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste0(paste(validation,collapse='\n'),'\n'),file='outputs/20b_integration_validation.txt',append=TRUE)
