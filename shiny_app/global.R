library(data.table)
library(shiny)
library(shinyjs)
library(bslib)
library(plotly)
library(tidyHeatmap)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(kableExtra)
library(igraph)
require(visNetwork)
library(openxlsx)
library(JBrowseR)
library(fst)
library(ggeasy)
library(shinyWidgets)
library(tidyverse)
source('bc3net_enrichment_modified.R')
pdf(file = NULL)
#setwd('shiny_app')
metadata = read_fst('data/sra_metadata_filtered_hand_edited.fst') %>% as_tibble() |>
  mutate(growth_medium = replace_na(growth_medium, "Unknown"))
functional_annotation_jcvi = read_fst('data/functional_annotation_jcvi.fst') %>%
  as_tibble() %>%
  dplyr::rename(gene_id = gene, `protein name` = `protein names`, `Gene Ontology` = `gene ontology`, 
                `KEGG pathways` = `kegg pathways`, `biosynthetic gene clusters` = smurf_and_known_metabolite, 
                `Subcellular localization` = deeploc_location, 
                `Interpro domains` = `interpro domains`) %>%
  mutate(`Gene Ontology` = str_remove_all(`Gene Ontology`, "'de novo'"), 
         genome = 'JCVI') %>% 
  select(-c(`best blast hit`, hgtector_potential_donor, hgtector_silhouette_score))
functional_annotation_chrom_level = read_fst('data/functional_annotation_chrom_level.fst') %>% 
  as_tibble() %>%
  filter(compound != 'aflatoxin G1' | is.na(compound)) %>% 
  mutate(`Gene Ontology (GO)`= str_remove_all(`Gene Ontology (GO)`, ' \\[GO:.*?\\]'),
         `Subcellular location [CC]` = str_remove(`Subcellular location [CC]`, 'SUBCELLULAR LOCATION: '), 
         genome = 'chrom_level') %>%
  select(gene_id, `protein name` = `Protein names`, 
         `Gene Ontology` = `Gene Ontology (GO)`, `biosynthetic gene clusters` = smurf_and_known_metabolite, 
         `KEGG pathways` = kegg_pathway,
         `Subcellular localization` = `Subcellular location [CC]`, `Interpro domains` = interpro, 
         SignalP:ApoplastP_probability, CAZyme, Cofactor, Pathway, genome)
functional_annotation = functional_annotation_jcvi %>%
  bind_rows(functional_annotation_chrom_level)
annotation_col_to_df = function(column, annotation){
  df = annotation %>% 
    select(gene_id, genome, column_to_select = {{ column }}) %>%
    filter(!is.na(column_to_select)) %>%
    separate_longer_delim(column_to_select, delim=';') %>%
    group_by(genome, column_to_select) %>%
    summarise(gene_id = list(gene_id), num_genes = n(), .groups = 'drop_last') %>%
    mutate(display_text = str_glue('{column_to_select} ({num_genes} genes)') %>% as.character())
  if (column %in% c('Gene Ontology', 'KEGG pathways', 
                    'Subcellular localization', 'Interpro domains')) {
    df = df %>% filter(num_genes > 10)
  }
  # if (column %in% c('biosynthetic gene clusters')){
  #   df = df %>% filter(num_genes > 2)
  # }
  return(df)
}
annotation_categories = c('Gene Ontology', 'KEGG pathways', 'biosynthetic gene clusters',
                         'Subcellular localization', 'Interpro domains')
annotation_list = annotation_categories %>% purrr::set_names() %>% 
  map(~annotation_col_to_df(.x, annotation = functional_annotation))
annotation_categories = c(annotation_categories, 'Gene list (Comma separated)')
annotation_list[['Gene list (Comma separated)']] = functional_annotation %>%
  select(gene_id, genome) %>%
  mutate(display_text = gene_id, column_to_select = gene_id) %>%
  rowwise() %>%
  mutate(gene_id = list(gene_id)) %>%
  ungroup()
bioproject_sizes = metadata %>% dplyr::count(bioproject, study_title, sort = TRUE)
colors202 <- c('#ff0029','#377eb8','#66a61e','#984ea3','#00d2d5','#ff7f00','#af8d00','#7f80cd','#b3e900','#c42e60','#a65628','#f781bf','#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#fccde5','#bc80bd','#ffed6f','#c4eaff','#cf8c00','#1b9e77','#d95f02','#e7298a','#e6ab02','#a6761d','#0097ff','#00d067','#000000','#252525','#525252','#737373','#969696','#bdbdbd','#f43600','#4ba93b','#5779bb','#927acc','#97ee3f','#bf3947','#9f5b00','#f48758','#8caed6','#f2b94f','#eff26e','#e43872','#d9b100','#9d7a00','#698cff','#d9d9d9','#00d27e','#d06800','#009f82','#c49200','#cbe8ff','#fecddf','#c27eb6','#8cd2ce','#c4b8d9','#f883b0','#a49100','#f48800','#27d0df','#a04a9b','#75d3e1','#ffb3b3','#3879c4','#a6d75b','#d673a3','#45c1c7','#ff9857','#d2a700','#8f8fe1','#caf025','#da4f82','#c07536','#ffbbdd','#67eed2','#aeb4ff','#f29292','#ffcc55','#e4b8ff','#ffdeb3','#8cddff','#ff8c42','#5ee592','#ccff00','#448aff','#ffdfc9','#d7de9a','#6e63ff','#00ffbb','#ff5e57','#ffccbb','#99cc00','#8f00ff','#eeaaee','#00aaee','#ff4444','#009988','#cc00ff','#ff7733','#448800','#0077cc','#ff99dd','#aa99ff','#cc3333','#33aa99','#cc8800','#9900cc','#00ddff','#ffcc00','#6600ff','#bb6600','#99ff00','#ff77aa','#55aaff','#ddee00','#00aa55','#ff5500','#3300ff','#ff1155','#55ffaa','#ffaa33','#0055aa','#99bb00','#ee7700','#660099','#00dd99','#ff88ff','#55cc77','#0099cc','#ff0066','#00ff44','#ffbb77','#7799ff','#9933ff','#77ff66','#ff9966','#66ddff','#ff7744','#1188ff','#33ffbb','#cc5500','#1166aa','#ff5599','#00aa77','#dd55ff','#77cc88','#ff3322','#2299ff','#99dd33','#ff7722','#3366ff','#00aa99','#ffaa11','#66ccff','#ff3366','#00dd66','#33ff77','#0099ff','#cc0033','#44ffcc','#0066ff','#ff33cc','#00ff99','#7700ff','#ff4400','#33ccff','#ff7700','#ffcc44','#22ff88','#ff0044','#0077ff','#33aa55','#ff8888','#bbff00','#ff0033','#00ffcc','#ff9922','#22ccff','#11ccff','#ffbb00','#66ff33','#22cc99','#aa00ff','#ff2200','#44ddff','#ff6655','#00ff55','#ffbb88','#44ff77','#dd00ff','#ff9900','#3388ff','#66ff88','#ff6600','#33cc44','#ff2233','#00ffdd','#dd9900')

#AF_genes = c('AFLA_139150', 'AFLA_139160', 'AFLA_139170', 'AFLA_139180', 'AFLA_139190', 'AFLA_139200', 'AFLA_139210', 'AFLA_139220', 'AFLA_139230', 'AFLA_139240', 'AFLA_139250', 'AFLA_139260', 'AFLA_139270', 'AFLA_139280', 'AFLA_139290', 'AFLA_139300', 'AFLA_139310', 'AFLA_139320', 'AFLA_139330', 'AFLA_139340', 'AFLA_139360', 'AFLA_139370', 'AFLA_139380', 'AFLA_139390', 'AFLA_139400', 'AFLA_139410', 'AFLA_139420', 'AFLA_139430', 'AFLA_139440')
network_genes = read_fst('data/gene_in_networks.fst')
annotation_list_network = c('Gene Ontology', 'KEGG pathways', 'biosynthetic gene clusters', 
                            'Subcellular localization', 'Interpro domains') %>% purrr::set_names() %>% 
  map(~annotation_col_to_df(.x, annotation = functional_annotation %>% 
                              semi_join(network_genes, by='gene_id')))
annotation_list_network[['Gene list (Comma separated)']] = functional_annotation %>%
  select(gene_id, genome) %>%
  semi_join(network_genes, by='gene_id') %>%
  mutate(display_text = gene_id, column_to_select = gene_id) %>%
  rowwise() %>%
  mutate(gene_id = list(gene_id)) %>%
  ungroup()
## Use local storage for JBrowse if running from host within ARS to reduce Google billing costs
jbrowse_resources_base = ifelse(str_starts(Sys.info()['nodename'], regex('ars', ignore_case = TRUE)),
                                'http://10.114.143.20/data/',
                                'https://storage.googleapis.com/flavus_jbrowse/'
                                )
#data_server = serve_data('data')
#data_server$stop_server()



