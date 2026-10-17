# Resolve platinum-drug identity coverage using supplementary ChEBI evidence.
# Keep original 10b/10c audit tables unchanged. Full keys alone are not geometry proof.
library(readr);library(dplyr);library(tibble);library(xml2);library(httr2)
if (!all(dir.exists(c("scripts", "config", "data")))) stop('Run from the project root.')
r<-function(p)read_csv(p,col_types=cols(.default=col_character()))
fetch<-function(url,path,json=FALSE){
  if(file.exists(path)) {x<-readRDS(path);stopifnot(identical(x$url,url));return(x)}
  resp<-request(url)|>req_timeout(45)|>req_retry(max_tries=2,retry_on_failure=TRUE)|>req_perform()
  x<-list(url=url,retrieved_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
    body=if(json)resp_body_json(resp) else resp_body_string(resp))
  saveRDS(x,path);x
}
x<-r('data/curated/seed_identity_crosscheck.csv')
chebi<-bind_rows(lapply(c('27899','31355'),function(id){
  z<-fetch(paste0('https://www.ebi.ac.uk/chebi/CHEBI:',id),paste0('data/raw/identity_api/10d_chebi_',id,'.rds'))
  rows<-xml_find_all(read_html(z$body),'//tr')
  val<-function(label){v<-lapply(rows,function(row)trimws(xml_text(xml_find_all(row,'./th|./td'))));
    hits<-v[vapply(v,function(a)length(a)>1L&&identical(a[1],label),logical(1))]
    if(length(hits)!=1L) stop('Missing/ambiguous ChEBI field: ',label)
    hits[[1]][2]}
  tibble(drugbank_id=if(id=='27899')'DB00515' else 'DB00958',chebi_id=paste0('CHEBI:',id),
    chebi_inchikey=val('InChIKey'),chebi_inchi=val('InChI'),chebi_smiles=val('SMILES'),
    chebi_formula=val('Formula'),chebi_url=z$url,chebi_retrieved_utc=z$retrieved_utc)
}))
u<-'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/84691,5460033/property/InChIKey,InChI,MolecularFormula,Title/JSON'
z<-fetch(u,'data/raw/identity_api/10d_cisplatin_candidate_CIDs.rds',TRUE)
candidates<-bind_rows(lapply(z$body$PropertyTable$Properties,function(p)tibble(cid=as.character(p$CID),
  title=p$Title,inchikey=p$InChIKey,inchi=p$InChI,formula=p$MolecularFormula)))
syn<-fetch('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/5460033/synonyms/JSON',
  'data/raw/identity_api/10d_cisplatin_synonyms.rds',TRUE)
names<-unlist(syn$body$InformationList$Information[[1]]$Synonym,use.names=FALSE)
cis<-x |> filter(drugbank_id=='DB00515')
selected<-candidates |> filter(cid=='5460033')
stopifnot(nrow(selected)==1L,tolower(selected$title)=='cisplatin',
 selected$inchikey==cis$inchikey,cis$cas_number %in% names)
y<-x |> left_join(chebi,by='drugbank_id') |>
 mutate(chebi_key_status=case_when(is.na(chebi_id)~'not_queried',
   inchikey==chebi_inchikey ~ 'full_key_match',TRUE~'conflict_requires_review'),
   reviewed_pubchem_cid=if_else(drugbank_id=='DB00515','5460033',resolved_pubchem_cid),
   reviewed_identity_status=case_when(drugbank_id=='DB00515'~'full_key_ChEBI_and_named_CAS_supported_PubChem_record',
     drugbank_id=='DB00958'~'full_key_ChEBI_and_PubChem',TRUE~identity_status),
   review_note=case_when(drugbank_id=='DB00515'~paste('10b could not choose uniquely between PubChem CIDs 84691 and 5460033:',
      'both share the full key. Select named cisplatin record 5460033 with source CAS synonym support;',
      'retain 84691 separately. Geometry not validated using a chemistry toolkit.'),
     drugbank_id=='DB00958'~'Source key matches ChEBI and PubChem 10339178; ChEMBL structure_type NONE.',
     TRUE~NA_character_))
stopifnot(nrow(y)==170,all(y$chebi_key_status[y$drugbank_id %in% chebi$drugbank_id]=='full_key_match'),
  all(y$inchi[y$drugbank_id %in% chebi$drugbank_id]==y$chebi_inchi[y$drugbank_id %in% chebi$drugbank_id]))
ids<-r('data/curated/seed_external_identifiers_enriched.csv') |>
 bind_rows(chebi |> transmute(seed_drugbank_id=drugbank_id,resource='ChEBI',identifier=chebi_id,
   source='ChEBI public record',source_record=chebi_id,relationship='external_full_key_match'),
   tibble(seed_drugbank_id='DB00515',resource='PubChem Compound',identifier='5460033',source='PubChem API',
     source_record='5460033',relationship='full_key_named_record_and_CAS_supported')) |> distinct()
a<-r('data/curated/seed_aliases_enriched.csv') |>
 bind_rows(tibble(seed_drugbank_id='DB00515',alias=names,alias_type='PubChem_record_synonym',source='PubChem',
   source_record='5460033',record_match_status='full_key_named_record_and_CAS_supported',
   alias_validation='record-linked label; individual name not sufficient for entity merge',
   identity_status='full_key_ChEBI_and_named_CAS_supported_PubChem_record')) |>
 distinct(seed_drugbank_id,source,source_record,alias,alias_type,.keep_all=TRUE)
write_csv(chebi,'data/curated/platinum_chebi_identity_records.csv',na='')
write_csv(candidates,'outputs/10d_cisplatin_PubChem_candidates.csv',na='')
write_csv(y,'data/curated/seed_identity_crosscheck_reviewed.csv',na='')
write_csv(ids,'data/curated/seed_external_identifiers_reviewed.csv',na='')
write_csv(a,'data/curated/seed_aliases_reviewed.csv',na='')
report<-c('SCRIPT 10d - PLATINUM IDENTITY REVIEW',
 'Cisplatin: CHEBI:27899; LXZZYRPGZAFOLE-UHFFFAOYSA-L; PubChem 5460033.',
 'Carboplatin: CHEBI:31355; OLESAACUTLOWQZ-UHFFFAOYSA-L; PubChem 10339178.',
 'Both ChEBI full InChI and InChIKey fields agree with DrugBank.',
 'Cisplatin was non-unique in PubChem key lookup, not absent: 84691 and 5460033',
 'share that key. Named record 5460033 is supported by CAS 15663-27-1 among synonyms.',
 'ChEMBL records CHEMBL11359 and CHEMBL1351 have structure_type NONE.',
 'Result: all 55 primary compounds now have an external full-key match, with',
 '53 matching ChEMBL/PubChem and two supported by ChEBI/PubChem review.',
 'Standard InChI/key agreement alone does not establish platinum coordination geometry;',
 'keep curated compound identities and source structure drawings for geometry-sensitive work.',
 'No different CID records were collapsed; original 10b/10c outputs remain audit snapshots.',
 'Reviewed outputs (data/curated/): seed_identity_crosscheck_reviewed.csv;',
 'seed_external_identifiers_reviewed.csv; seed_aliases_reviewed.csv;',
 'platinum_chebi_identity_records.csv. Candidate comparison: outputs/10d_*.',
 'These data establish mapping support, not eEF1A binding or inhibition.',
 'Next: SIDER mapping using the reviewed tables; retain source-specific representation caveats.')
writeLines(report,'outputs/10d_platinum_identity_summary.txt')
cat(paste0('\nSCRIPT 10d RUN - ',format(Sys.time(),tz='Europe/Rome',usetz=TRUE),'\n',
 'Script: scripts/10d_resolve_platinum_identities.R\n',
 'Inputs: 10b/10c outputs plus cached ChEBI pages and PubChem candidate/synonym API records.\n',
 paste(report,collapse='\n'),'\n'),file='Notes.txt',append=TRUE)
cat(paste(report,collapse='\n'),'\n')
