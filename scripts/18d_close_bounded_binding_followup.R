# Close the bounded follow-up; affinity is annotation, never a RESKO inclusion gate.
library(readr);library(dplyr);library(tibble);library(xml2)
source('scripts/helpers/literature_public_get.R')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
input<-'data/curated/association_literature_search_hits.csv';before<-tools::md5sum(input)
scope<-r(input) |> filter(source_name=='Testosterone enanthate')
stopifnot(nrow(scope)==5L)
rows<-list();passages<-list()
for(i in seq_len(nrow(scope))){
 id<-scope$pmcid[i];url<-paste0('https://www.ebi.ac.uk/europepmc/webservices/rest/',id,'/fullTextXML')
 cache<-paste0('data/raw/target_evidence/association_sources/fulltext_',id,'.rds')
 z<-literature_public_get(url,cache,FALSE)
 status<-paste0('full_text_endpoint_HTTP_',z$status)
 if(z$status==200L){
  h<-read_xml(z$body);source_pmid<-xml_text(xml_find_first(h,'/article/front/article-meta/article-id[@pub-id-type="pmid"]'))
  stopifnot(identical(source_pmid,scope$pmid[i]));status<-'targeted_full_text_context_reviewed'
  p<-xml_text(xml_find_all(h,'//body//p|//body//table-wrap'))
  selected<-p[grepl('testosterone enanthate|testosterone enantate|eef1a|elongation factor.?1.?alpha',p,ignore.case=TRUE)]
  passages[[length(passages)+1L]]<-tibble(literature_key=scope$literature_key[i],pmcid=id,passage=selected)
 }
 rows[[i]]<-tibble(literature_key=scope$literature_key[i],pmid=scope$pmid[i],pmcid=id,
 full_text_status=status,http_status=z$status,source_url=url,cache,retrieved_utc=z$retrieved_utc)
}
access<-bind_rows(rows)
decisions<-tribble(~pmcid,~context_conclusion,~review_basis,
 'PMC9996746','Male-contraception review discusses testosterone enanthate separately from gamendazole/eEF1A1; no transfer of binding claim','Cached XML targeted paragraphs/table',
 'PMC6563988','Human bed-rest/exercise/testosterone proteomics: EEF1A1 abundance/phosphorylation measurements are not affinity assays; significant control-group changes cannot be assigned to drug','Cached XML intervention paragraph, EEF1A1 result paragraph and tables 1-4',
 'PMC9184684','Human cancer skeletal-muscle proteomics: EEF1A2 S-nitrosylation change appears in cancer-treatment plus placebo context; not evidence of testosterone-enanthate binding','Cached XML treatment and S-nitrosylation paragraphs',
 'PMC9426764','Unresolved: full-text endpoint failure; retain pending context review','Acquisition status only; no mechanistic conclusion',
 'PMC4381878','Nonhormonal-contraception review: testosterone formulations in introduction; eEF1A1 claim concerns biotinylated gamendazole in separate section','PMC HTML manually checked 2026-10-20: introduction and section 2.5.2; XML endpoint failed')
write_csv(decisions,'config/18d_testosterone_context_decisions.csv')
review<-access |> left_join(decisions,by='pmcid',relationship='one-to-one') |>
 mutate(manual_source_url=if_else(pmcid=='PMC4381878','https://pmc.ncbi.nlm.nih.gov/articles/PMC4381878/',NA_character_),
 human_binding_affinity='not_established_in_reviewed_context',STITCH_original_study_verified=FALSE,
 qualifies_RESKO_inclusion_gate='Binding affinity is not an inclusion gate; retain all original leads')
write_csv(bind_rows(passages),'data/curated/testosterone_enanthate_full_text_contexts.csv',na='')
write_csv(review,'data/curated/testosterone_enanthate_context_review.csv',na='')
write_csv(review |> select(-cache),'outputs/18d_context_review_summary.csv',na='')
stopifnot(sum(access$http_status==200L)==3L,identical(before,tools::md5sum(input)))
policy<-tibble(setting=c('require_published_eef1a_affinity_for_RESKO_inclusion','target_species','binding_review_role',
 'missing_affinity_encoding','candidate_retention'),
 value=c('FALSE','Homo sapiens','Optional evidence annotation; not RESKO completion criterion',
 'NA/unreported/unknown; never zero or inactive','Retain original candidates regardless of affinity availability'),
 authority='User clarification 20 October 2026')
write_csv(policy,'config/RESKO_evidence_scope.csv')
report<-c('SCRIPT 18d - BOUNDED BINDING FOLLOW-UP CLOSED; RETURN TO CORE RESKO INTEGRATION',
 'User confirms that published affinity is not required for RESKO completion or candidate retention.',
 'Five testosterone-enanthate contexts: three XML sources acquired/reviewed, one additional PMC HTML context reviewed, one access hold.',
 'Reviewed contexts concern separate compounds, expression/proteomics or placebo-group changes, not a measured testosterone-enanthate/human eEF1A affinity.',
 'This is a bounded contextual review, not proof that no binding publication exists. Original STITCH provenance remains unresolved.',
 'All six leads remain in the inventory; missing affinity does not exclude a drug or imply inactivity.',
 'Outputs: context source/status review in outputs/18d_*; source passages local; user scope policy in config/RESKO_evidence_scope.csv.',
 'Next: audit the actual RESKO paper against completed code, then acquire/integrate the remaining original data sources; do not continue open-ended affinity searching.')
writeLines(report,'outputs/18d_followup_summary.txt')
writeLines('PASS: five contexts accounted for, three XML PMIDs validated, source hash unchanged, optional-affinity scope explicitly configured.',
 'outputs/18d_followup_validation.txt')
cat(paste0('\nSCRIPT 18d RUN - ',format(Sys.time(),tz='Europe/London',usetz=TRUE),'\nScript: scripts/18d_close_bounded_binding_followup.R\n',paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
