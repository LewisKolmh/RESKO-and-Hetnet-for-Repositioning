# Meaningful boundary/selection test independent of the actual extraction output.
validate_selected_drugbank_stream <- function() {
  path <- tempfile(fileext = '.xml'); on.exit(unlink(path), add = TRUE)
  writeLines(c('<drugbank>', '<drug type="small molecule">',
    '  <drugbank-id primary="true">DB00001</drugbank-id>', '  <name>Selected one</name>',
    '  <targets>', '    <target><id>BE1</id></target>', '  </targets>',
    '  <unwanted><x>Must not retain</x></unwanted>', '</drug>',
    '<drug type="small molecule">', '  <drugbank-id primary="true">DB00002</drugbank-id>',
    '  <name>Unselected</name>', '  <targets>', rep('    <ignored/>', 200), '  </targets>', '</drug>',
    '<drug type="biotech">', '  <drugbank-id primary="true">DB00003</drugbank-id>',
    '  <name>Selected three</name>', '  <targets/>', '</drug>', '</drugbank>'), path)
  outputs <- list()
  z <- stream_selected_drugbank_fields(path, c('DB00001', 'DB00003'),
    c('drugbank-id', 'name', 'targets'), function(node) {
      outputs[[length(outputs) + 1L]] <<- node
    }, chunk_lines = 3L)
  stopifnot(z$records_scanned == 3L, length(outputs) == 2L,
    xml2::xml_text(xml2::xml_find_first(outputs[[1]], './name')) == 'Selected one',
    xml2::xml_text(xml2::xml_find_first(outputs[[1]], './targets/target/id')) == 'BE1',
    length(xml2::xml_find_all(outputs[[1]], './unwanted')) == 0L,
    length(xml2::xml_find_all(outputs[[2]], './targets')) == 1L)
  TRUE
}
