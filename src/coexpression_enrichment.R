library(shiny)
library(bslib)
library(plotly)
library(tidyHeatmap)
library(ComplexHeatmap)
library(igraph)
require(visNetwork)
library(openxlsx)
library(ggraph)
library(tidyverse)
source('bc3net_enrichment_modified.R')
#setwd('shiny_app')

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

## Get neighbors of each set of BGC genes
## Test each set of neighbors for enrichment of GO terms.
## See if any amino acid catabolism GO terms are among the enriched terms.
## Identify BGCs that have coexpression relationships
get_node_neighbors_adaptive_cutoff = function(gene_ids, num_nodes_for_cutoff=50000){
  node_neighbors = ego(network, order = 1, nodes = gene_ids)
  g = induced_subgraph(network, unlist(node_neighbors))
  edge_weight_cutoff = E(g)$abs_r %>% sort(decreasing = TRUE) %>% .[num_nodes_for_cutoff]
  network2 = delete.edges(network, which(E(network)$abs_r < edge_weight_cutoff))
  gene_ids = gene_ids[gene_ids %in% V(network2)$name]
  node_neighbors = ego(network2, order = 1, nodes = gene_ids)
  g = induced_subgraph(network2, unlist(node_neighbors))
  return(V(g)$name)
}
get_node_neighbors_adaptive_cutoff('AFLA_023020')
bgc_neighbors = annotation_list_network$`biosynthetic gene clusters` %>%
  select(column_to_select, gene_id) %>%
  rowwise() %>%
  mutate(neighbors = list(get_node_neighbors_adaptive_cutoff(gene_id)))
bgc_neighbors_enrich = bgc_neighbors %>% 
  ungroup() %>%
  mutate(n_neighbors = lengths(neighbors)) %>%
  rowwise() %>%
  mutate(enrich_res = list(enrichment_test(gene_list = list('network' = neighbors),
                                           columns_list = c('Gene Ontology', 'KEGG pathways', 'biosynthetic gene clusters', 
                                                            'Subcellular localization (DeepLoc)', 'Interpro domains'),
                                           functional_annotation = functional_annotation_jcvi %>% filter(gene_id %in% V(network)$name)) %>%
                             filter(column_to_select != TermID) %>%
                             mutate(padjust = p.adjust(pval, method = 'fdr'))%>%
                             filter(padjust < 0.05) %>% 
                             select(-pval, -gene_list_name) %>%
                             relocate(c(padjust, annotation_category), .after = all) %>% 
                             rename(n_genes_subnetwork = genes, n_genes_all_network = all, pvalue_adjusted = padjust) %>% 
                             arrange(pvalue_adjusted))) %>%
  select(bgc = column_to_select, n_neighbors, enrich_res, neighbors) %>%
  unnest(enrich_res)
bgc_neighbors_enrich %>% 
  filter(str_detect(TermID, 'catabol')) %>% View()
bgc_neighbors_enrich %>% 
  filter(str_detect(TermID, 'smurf')) %>% View()
openxlsx::write.xlsx(list('bgc_neighbors_enrich'=bgc_neighbors_enrich, 
                          'filtered_for_catabol' = bgc_neighbors_enrich %>% filter(str_detect(TermID, 'catabol'))),
                     file = 'Z:/flavus_sra_rnaseq/output/coexpression/bgc_neighbors_enrichment.xlsx',
                     firstRow = TRUE)
### Enrichment of neighbors for each gene
gene_neighbors = annotation_list_network$`Gene list (Comma separated)` %>%
  select(column_to_select, gene_id) %>%
  rowwise() %>%
  mutate(neighbors = list(get_node_neighbors_adaptive_cutoff(gene_id)))
gene_neighbors_enrich = gene_neighbors %>% 
  ungroup() %>%
  mutate(n_neighbors = lengths(neighbors)) %>%
  rowwise() %>%
  mutate(enrich_res = list(enrichment_test(gene_list = list('network' = neighbors),
                                           columns_list = c('Gene Ontology', 'KEGG pathways', 'biosynthetic gene clusters', 
                                                            'Subcellular localization (DeepLoc)', 'Interpro domains'),
                                           functional_annotation = functional_annotation_jcvi %>% filter(gene_id %in% V(network)$name)) %>%
                             filter(column_to_select != TermID) %>%
                             mutate(padjust = p.adjust(pval, method = 'fdr'))%>%
                             filter(padjust < 0.05) %>% 
                             select(-pval, -gene_list_name) %>%
                             relocate(c(padjust, annotation_category), .after = all) %>% 
                             rename(n_genes_subnetwork = genes, n_genes_all_network = all, pvalue_adjusted = padjust) %>% 
                             arrange(pvalue_adjusted))) %>%
  select(bgc = column_to_select, n_neighbors, enrich_res, neighbors) %>%
  unnest(enrich_res)
openxlsx::write.xlsx(gene_neighbors_enrich, file = 'Z:/flavus_sra_rnaseq/output/coexpression/gene_neighbors_enrichment.xlsx',
                     firstRow = TRUE)
## RIP genes
## Core gene: AFLA_063260
## TF: AFLA_063330
## Tyrosinase: AFLA_063220
## Cyclotransferase: AFLA_063250
## Chlorinating enzyme: AFLA_063290
#AFLA_063260,AFLA_063330,AFLA_063220,AFLA_063250,AFLA_063290
## AFLA_06260 not in network
ripp_enrich = tibble(bgc = 'ripp', gene_id = list(c('AFLA_063260', 'AFLA_063330', 'AFLA_063220', 'AFLA_063250', 'AFLA_063290'))) %>%
  rowwise() %>%
  mutate(neighbors = list(get_node_neighbors_adaptive_cutoff(gene_id))) %>%
  ungroup() %>%
  mutate(n_neighbors = lengths(neighbors)) %>%
  rowwise() %>%
  mutate(enrich_res = list(enrichment_test(gene_list = list('network' = neighbors),
                                           columns_list = c('Gene Ontology', 'KEGG pathways', 'biosynthetic gene clusters', 
                                                            'Subcellular localization (DeepLoc)', 'Interpro domains'),
                                           functional_annotation = functional_annotation_jcvi %>% filter(gene_id %in% V(network)$name)) %>%
                             #filter(column_to_select != TermID) %>%
                             mutate(padjust = p.adjust(pval, method = 'fdr'))%>%
                             filter(padjust < 0.05) %>% 
                             select(-pval, -gene_list_name) %>%
                             relocate(c(padjust, annotation_category), .after = all) %>% 
                             rename(n_genes_subnetwork = genes, n_genes_all_network = all, pvalue_adjusted = padjust) %>% 
                             arrange(pvalue_adjusted))) %>%
  select(bgc, n_neighbors, enrich_res, neighbors) %>%
  unnest(enrich_res)
openxlsx::write.xlsx(ripp_enrich, file = 'Z:/flavus_sra_rnaseq/output/coexpression/ripp_neighbors_enrichment.xlsx',
                     firstRow = TRUE)
bgc_network = bgc_neighbors_enrich %>% 
  filter(str_detect(TermID, 'smurf')) %>% #View()
  mutate(weight = -log10(pvalue_adjusted) %>% round()) %>%
  select(from = bgc, to = TermID, weight) %>%
  graph_from_data_frame(directed = FALSE)
E(bgc_network)$weight
class(bgc_network)
data = bgc_network %>% toVisNetworkData() 
visNetwork(nodes = data$nodes, edges = data$edges) %>%
  visEdges(value='weight')
bgc_network %>%
  ggraph(layout = 'fr') %>%
  geom_edge_link(aes(width = weight), alpha = .25) +
  geom_node_point(color = "blue", size = 2) + 
  geom_node_text(aes(label = name),  repel = TRUE)+
  theme_graph()+
  labs(title = 'Graph with weighted edges', 
       subtitle = 'No scaling')

test = data.frame(w1=rep('like', 5), 
                w2 = c('apple', 'orange', 'pear','peach', 'banana'),
                weight= c(2,3,5,8, 15)) %>%
  graph_from_data_frame() 
class(test)
test %>% ggraph(layout = "fr") +
  geom_edge_link(alpha = .25, 
                 aes(width = weight)) +
  geom_node_point(color = "blue", size = 2) + 
  geom_node_text(aes(label = name),  repel = TRUE)+
  theme_graph()+
  labs(title = 'Graph with weighted edges', 
       subtitle = 'No scaling')
