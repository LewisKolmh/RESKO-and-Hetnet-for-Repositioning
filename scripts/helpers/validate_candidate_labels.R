library(readr);library(dplyr);library(xml2)
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
d<-r('data/curated/candidate_label_identity_review.csv');s<-r('data/curated/candidate_label_selection.csv')
cov<-r('outputs/14b_candidate_label_coverage.csv');ident<-r('data/curated/candidate_identity_registry.csv')
db<-r('data/curated/candidate_DrugBank_record_enrichment.csv')
stopifnot(nrow(cov)==30,all(as.integer(cov$sets_examined)<=3),nrow(s)==sum(as.integer(cov$supported_AR_sets)),n_distinct(s$stitch_flat_id)==sum(as.integer(cov$supported_AR_sets)>0),
 !any(s$stitch_flat_id %in% ident$stitch_flat_id[ident$automatic_label_mapping_blocked=='TRUE']))
for(i in seq_len(nrow(s))){
 a<-s[i,];z<-readRDS(a$raw_cache_file);doc<-read_xml(z$body);xml_ns_strip(doc)
 ing<-xml_find_all(doc,'//ingredient[starts-with(@classCode,"ACTI")]/ingredientSubstance')
 ids<-unique(na.omit(xml_attr(xml_find_all(ing,'./code'),'code')))
 moieties<-unique(na.omit(xml_attr(xml_find_all(ing,'.//activeMoiety/code'),'code')))
 src<-db |> filter(stitch_flat_id==a$stitch_flat_id)
 stopifnot(length(ids)==1L,length(moieties)<=1L,any(src$unii %in% c(ids,moieties)),
  xml_attr(xml_find_first(doc,'/document/setId'),'root')==a$set_id,
  xml_attr(xml_find_first(doc,'/document/versionNumber'),'value')==a$label_version)
 normal<-readRDS(a$cache_file)$body$results[[1]]
 source_text<-unique(xml_text(xml_find_all(doc,'//section[code[@code="34084-4"]]//text')))
 stopifnot(identical(as.character(unlist(normal$adverse_reactions)),source_text),
  nchar(paste(source_text,collapse=' '))==as.integer(a$adverse_text_chars))
}
paths<-c('data/curated/candidate_identity_registry.csv','data/curated/candidate_DrugBank_record_enrichment.csv',
 'config/candidate_identity_exceptions.csv','scripts/14b_collect_candidate_labels.R')
write_csv(tibble(file=paths,md5=unname(tools::md5sum(paths))),'outputs/14b_input_summary.csv')
invisible(parse('scripts/14b_collect_candidate_labels.R'))
cat('PASS: independent reread of all selected raw XML confirms ingredient UNIIs, single-ingredient scope, versions and exact normalised AR text.\nPASS: maximum three examined sets per group and original conflict exclusion confirmed; R parsing passed.\n',file='outputs/14b_candidate_label_validation.txt',append=TRUE)
cat(paste0('\nSCRIPT 14b - INDEPENDENT VALIDATION\nAll ',nrow(s),' selected group/label rows were rechecked against raw XML: ingredient UNIIs, single-ingredient scope, document versions and normalised adverse-reaction text agree.\nAll groups retain at most three examined sets; the sacrosidase conflict remains excluded. R parsing passed.\nOutputs: outputs/14b_candidate_label_validation.txt and outputs/14b_input_summary.csv (input hashes).\n'),file='Notes.txt',append=TRUE)
