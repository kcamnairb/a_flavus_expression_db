sra_id2df <- function(id) {
  xml_text = rentrez::entrez_fetch(db='sra', id=id, rettype = 'xml')
  sra_xml2 = xml2::as_xml_document(xml_text)
  run_df = xmlconvert::xml_to_df(text=xml_text, records.tags = '//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE') %>%
    pivot_wider(names_from = TAG, values_from = VALUE)
  bioproject = xml2::xml_find_first(sra_xml2, '//IDENTIFIERS/EXTERNAL_ID[@namespace="BioProject"]') %>% xml2::xml_text()
  biosample = xml2::xml_find_first(sra_xml2, '//IDENTIFIERS/EXTERNAL_ID[@namespace="BioSample"]') %>% xml2::xml_text()
  experiment = xml2::xml_find_first(sra_xml2, '//EXPERIMENT/IDENTIFIERS') %>% xml2::xml_text()
  run = xml2::xml_find_first(sra_xml2, '//RUN_SET/RUN/IDENTIFIERS/PRIMARY_ID') %>% xml2::xml_text()
  sample_name = xml2::xml_find_first(sra_xml2, '//Member') %>% xml2::xml_attr('sample_name')
  sample_title = xml2::xml_find_first(sra_xml2, '//Member') %>% xml2::xml_attr('sample_title')
  filename = xml2::xml_find_first(sra_xml2, '//SRAFile') %>% xml2::xml_attr('filename')
  email = xml2::xml_find_first(sra_xml2, '//Contact') %>% xml2::xml_attr('email')
  first_name = xml2::xml_find_first(sra_xml2, '//Contact/Name/First') %>% xml2::xml_text()
  last_name = xml2::xml_find_first(sra_xml2, '//Contact/Name/Last') %>% xml2::xml_text()
  alias = xml2::xml_find_first(sra_xml2, '//EXPERIMENT')  %>% xml2::xml_attr('alias')
  title2 = xml2::xml_find_first(sra_xml2, '//TITLE') %>% xml2::xml_text()
  design_description = xml2::xml_find_first(sra_xml2, '//DESIGN/DESIGN_DESCRIPTION') %>% xml2::xml_text()
  layout = xml2::xml_find_first(sra_xml2, '//LIBRARY_LAYOUT') %>% xml2::xml_child() %>% as.character() %>% str_remove_all('[</>]')
  pubmed = xml2::xml_find_first(sra_xml2, '//STUDY_LINKS/STUDY_LINK/XREF_LINK') %>% xml2::xml_text() %>% str_remove('pubmed')
  selection = xml2::xml_find_first(sra_xml2, '//LIBRARY_SELECTION') %>% xml2::xml_text()
  study_title = xml2::xml_find_first(sra_xml2, '//STUDY_TITLE') %>% xml2::xml_text()
  instrument = xml2::xml_find_first(sra_xml2, '//INSTRUMENT_MODEL') %>% xml2::xml_text()
  abstract = xml2::xml_find_first(sra_xml2, '//STUDY_ABSTRACT') %>% xml2::xml_text()
  date  = xml2::xml_find_first(sra_xml2, '//SRAFiles/SRAFile') %>% xml2::xml_attr('date')
  run_df = run_df %>%
    bind_cols(bioproject = bioproject,
              biosample = biosample,
              experiment = experiment,
              run = run,
              sample_name = sample_name,
              sample_title = sample_title,
              filename = filename,
              email = email,
              first_name = first_name,
              last_name = last_name,
              alias = alias,
              title2 = title2,
              design_description = design_description,
              layout = layout,
              pubmed = pubmed,
              selection = selection,
              study_title = study_title,
              instrument = instrument,
              abstract = abstract,
              date  = date)
  return(run_df)
} 
