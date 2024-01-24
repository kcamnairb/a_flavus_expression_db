library(data.table)
library(shiny)
library(bslib)
library(plotly)
library(tidyHeatmap)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(kableExtra)
library(igraph)
require(visNetwork)
library(openxlsx)
library(tidyverse)
source('bc3net_enrichment_modified.R')
#setwd('shiny_app')
shinyOptions(cache = cachem::cache_disk('cache'))
tpm_jcvi = read_csv('data/A_flavus_jcvi_tpm.csv', show_col_types = FALSE)
tpm_chrom_level = read_csv('data/A_flavus_chrom_level_tpm.csv', show_col_types = FALSE)
vst_jcvi = read_csv('data/vst_jcvi.csv', show_col_types = FALSE)
vst_chrom_level = read_csv('data/vst_chrom_level.csv', show_col_types = FALSE)

metadata = read_csv('data/sra_metadata_filtered.csv', show_col_types = FALSE) %>%
  select(-collection_date, run_title = title2) %>%
  dplyr::rename(growth_condition = `growth condition`) %>%
  mutate(across(c(strain, sample_type,  source_name, genotype), 
                ~str_remove(.x, regex('^A.* flavus ', ignore_case = TRUE)))) %>%
  mutate(date = as.character(date)) %>%
  
  group_by(bioproject) %>%
  mutate(title = case_when(n_distinct(sample_title) == 1 ~ run_title,
                           n_distinct(run_title) == 1 & n_distinct(sample_title) > 1 ~ sample_title,
                           TRUE ~ run_title)) %>%
  mutate(title = str_remove(title, 'GSM.*?: ') %>% str_remove('; Aspergillus flavus; RNA-Seq')) %>%
  ungroup() %>%
  unite('sample_description', title, strain, sample_type, source_name, genotype, treatment, timepoint, growth_condition, isolate,
        sep = ';', na.rm =TRUE, remove = FALSE) %>%
  filter(run %in% colnames(tpm_jcvi)) %>%
  mutate(submitter = str_c(first_name, last_name, collapse = ' ')) %>%
  select(-first_name, -last_name, -BioSampleModel, -pubmed) %>%
  mutate(doi = ifelse(is.na(doi), NA, str_glue('https://doi.org/{doi}'))) %>%
  mutate(doi = ifelse(is.na(doi), NA, str_glue("<a href='{doi}'>{doi}</a>"))) %>%
  select(-filename)
functional_annotation_jcvi = read_csv('data/functional_annotation_jcvi.csv') %>%
  dplyr::rename(gene_id = gene, `Gene Ontology` = `gene ontology`, `KEGG pathways` = `kegg pathways`, 
         `biosynthetic gene clusters` = smurf_and_known_metabolite, 
         `Subcellular localization (DeepLoc)` = deeploc_location, 
         `Interpro domains` = `interpro domains`) %>%
  mutate(`Gene Ontology` = str_remove_all(`Gene Ontology`, "'de novo'"))
functional_annotation_chrom_level = read_csv('data/functional_annotation_chrom_level.csv') %>% 
  filter(compound != 'aflatoxin G1')
annotation_col_to_df = function(column, annotation){
  df = annotation %>% 
    select(gene_id, column_to_select = {{ column }}) %>%
    filter(!is.na(column_to_select)) %>%
    separate_longer_delim(column_to_select, delim=';') %>%
    group_by(column_to_select) %>%
    summarise(gene_id = list(gene_id), num_genes = n()) %>%
    mutate(display_text = str_glue('{column_to_select} ({num_genes} genes)') %>% as.character())
  if (column %in% c('Gene Ontology', 'KEGG pathways', 
                    'Subcellular localization (DeepLoc)', 'Interpro domains')) {
    df = df %>% filter(num_genes > 10)
  }
  if (column %in% c('biosynthetic gene clusters')){
    df = df %>% filter(num_genes > 2)
  }
  return(df)
}
annotation_categories = c('Gene Ontology', 'KEGG pathways', 'biosynthetic gene clusters', 
                          'Subcellular localization (DeepLoc)', 'Interpro domains')
annotation_list = annotation_categories %>% purrr::set_names() %>% 
  map(~annotation_col_to_df(.x, annotation = functional_annotation_jcvi))
annotation_categories = c(annotation_categories, 'Gene list (Comma separated)')
annotation_list[['Gene list (Comma separated)']] = functional_annotation_jcvi %>%
  select(gene_id) %>%
  mutate(display_text = gene_id, column_to_select = gene_id) %>%
  rowwise() %>%
  mutate(gene_id = list(gene_id)) %>%
  ungroup()
bioproject_sizes = metadata %>% dplyr::count(bioproject, study_title, sort = TRUE)
mpn65 = c('#ff0029','#377eb8','#66a61e','#984ea3','#00d2d5','#ff7f00','#af8d00','#7f80cd','#b3e900','#c42e60','#a65628',
          '#f781bf','#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#fccde5','#bc80bd','#ffed6f','#c4eaff','#cf8c00',
          '#1b9e77','#d95f02','#e7298a','#e6ab02','#a6761d','#0097ff','#00d067','#000000','#252525','#525252','#737373',
          '#969696','#bdbdbd','#f43600','#4ba93b','#5779bb','#927acc','#97ee3f','#bf3947','#9f5b00','#f48758','#8caed6',
          '#f2b94f','#eff26e','#e43872','#d9b100','#9d7a00','#698cff','#d9d9d9','#00d27e','#d06800','#009f82','#c49200',
          '#cbe8ff','#fecddf','#c27eb6','#8cd2ce','#c4b8d9','#f883b0','#a49100','#f48800','#27d0df','#a04a9b')
AF_genes = c('AFLA_139150', 'AFLA_139160', 'AFLA_139170', 'AFLA_139180', 'AFLA_139190', 'AFLA_139200', 'AFLA_139210', 'AFLA_139220', 'AFLA_139230', 'AFLA_139240', 'AFLA_139250', 'AFLA_139260', 'AFLA_139270', 'AFLA_139280', 'AFLA_139290', 'AFLA_139300', 'AFLA_139310', 'AFLA_139320', 'AFLA_139330', 'AFLA_139340', 'AFLA_139360', 'AFLA_139370', 'AFLA_139380', 'AFLA_139390', 'AFLA_139400', 'AFLA_139410', 'AFLA_139420', 'AFLA_139430', 'AFLA_139440')
modules = readRDS('data/modules.rds')
network = readRDS('data/network.rds')
annotation_list_network = c('Gene Ontology', 'KEGG pathways', 'biosynthetic gene clusters', 
                            'Subcellular localization (DeepLoc)', 'Interpro domains') %>% purrr::set_names() %>% 
  map(~annotation_col_to_df(.x, annotation = functional_annotation_jcvi %>% 
                              filter(gene_id %in% V(network)$name)))
annotation_list_network[['Gene list (Comma separated)']] = functional_annotation_jcvi %>%
  select(gene_id) %>%
  filter(gene_id %in% V(network)$name) %>%
  mutate(display_text = gene_id, column_to_select = gene_id) %>%
  rowwise() %>%
  mutate(gene_id = list(gene_id)) %>%
  ungroup()
