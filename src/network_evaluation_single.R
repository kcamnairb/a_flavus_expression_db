library(tidyverse)
setwd('~/flavus_sra_rnaseq')
source('shiny_app/bc3net_enrichment_modified.R')
args = commandArgs(trailingOnly = TRUE)
functional_annotation_jcvi_network = read_rds('~/flavus_sra_rnaseq/output/coexpression/network_evaluation/functional_annotation_jcvi_network')
top_100_edges_single = read_rds(args[1])
enrich_res = top_100_edges_single %>% 
  left_join(functional_annotation_jcvi_network %>% 
               filter(!is.na(`biosynthetic gene clusters`) | !is.na(`Gene Ontology`)) %>%
               select(gene_id, `biosynthetic gene clusters`, `Gene Ontology`), by='gene_id') %>%
  rowwise() %>%
  mutate(enrich_res = list(enrichment_test(gene_list = list('network' = top100),
                                           columns_list = c('Gene Ontology', 'biosynthetic gene clusters'),
                                           functional_annotation = functional_annotation_jcvi_network) %>%
                             mutate(padjust = p.adjust(pval, method = 'fdr'))%>%
                             filter(padjust < 0.05) %>% 
                             select(-pval, -gene_list_name) %>%
                             relocate(c(padjust, annotation_category), .after = all) %>% 
                             rename(n_genes_subnetwork = genes, n_genes_all_network = all, pvalue_adjusted = padjust) %>% 
                             arrange(pvalue_adjusted))) %>%
  unnest(enrich_res)
enrich_res %>% select(-top100) %>% write_csv(paste0(args[1],'_enrich.csv'))