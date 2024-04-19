## used network_evaluation_single.R to run in parallel on the server.
top_100_edges_enrich = list.files(here('output/coexpression/network_evaluation/'), 
                                  pattern = '.*-.*enrich.csv', full.names = TRUE) %>%
  map(~read_csv(.x) %>% mutate(network_method = str_replace(.x, '.*/(.*?)_top_100.*', '\\1'))) %>%
  bind_rows()
## From Vandenbon 2022 describing GO enrichment frequency:
## "We define EnrichmentMF, EnrichmentBP,and EnrichmentCC as the fraction of genes in a network for 
## which setX contained one or more significantly enriched GO terms"
network_enrich_freq = top_100_edges_enrich  %>%
  group_by(network_method, annotation_category) %>%
  summarize(enrichment_frequency = n() / nrow(genes_with_mean_tpm_above_1))
network_enrich_freq %>% write_csv(here('output/coexpression/network_enrichment_frequency.csv'))
enrich_frequency = top_100_edges_enrich %>%
  mutate(neg_log_pvalue = -log2(pvalue_adjusted),
         network_method = str_remove(network_method, 'get_')) %>%
  group_by(network_method, annotation_category) %>%
  mutate(med_neg_log_pvalue = median(neg_log_pvalue), 
         count = n()) %>%
  ungroup() %>%
  mutate(network_method = fct_reorder(network_method, med_neg_log_pvalue, .desc=FALSE, .fun = max),
         annotation_category = str_to_title(annotation_category))
make_tabular_workflow = . %>% separate_wider_delim(network_method, delim = '-', cols_remove = FALSE,
                                                   names = c('Downsampling', 'Batch Correction', 'Normalization', 'Scaling', 'Correlation')) %>% 
  mutate(Downsampling = if_else(Downsampling == 'downsampled_data', 'Random', 'None'),
         `Batch Correction` = if_else(`Batch Correction` == 'combat_seq_correct', 'ComBat-seq', 'None'),
         Normalization = str_remove(Normalization, '_normalize') %>% str_to_upper() %>% 
           str_replace('UPPER_QUARTILE', 'UQ'),
         Scaling = str_remove(Scaling, '_scale') %>% str_replace('inverse_hyperbolic_sine', 'asinh') %>%
           str_replace('no_scaling', 'None'),
         Correlation = str_replace_all(Correlation, c('pearson' = 'Pearson', 'wgcna' = 'WGNCA')))
p_frequency =  enrich_frequency %>%
  ggplot(aes(network_method, neg_log_pvalue)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.data = function(x) c(y = 2.5, label = length(x)),
               geom = "text", size=2.5) +
  facet_wrap(~annotation_category) +
  coord_flip(ylim = c(2.5, 60)) +
  theme_bw() +
  ggeasy::easy_remove_legend() +
  ggeasy::easy_remove_y_axis(what = c('title', 'text')) +
  ggtitle('Enrichment frequency') +
  ggeasy::easy_center_title()
df_frequency = enrich_frequency %>% 
  distinct(network_method) %>% 
  make_tabular_workflow() %>%
  arrange(desc(network_method)) %>%
  select(-network_method)
tbl_theme = ttheme_minimal(base_size = 6.8, base_colour = "black", base_family = "",
                           parse = FALSE, padding = unit(c(.90, .90), "mm"))
tbl_grob = gridExtra::tableGrob(df_frequency, theme= tbl_theme, rows = NULL)
p = gridExtra::grid.arrange(tbl_grob, p_frequency, ncol = 2, widths = c(2.5/10, 7.5/10))
ggsave(here('output/coexpression/enrichment_frequency.png'), p, width = 14, height = 8, units='in', dpi = 300)

## From Vandenbon 2022 describing GO enrichment accuracy:
## Where we found setX to have enriched GO terms, we checked if the enriched
## terms overlapped with the GO terms of gene X. AccuracyMF, AccuracyBP,and AccuracyCC were 
## defined as the fraction of genes in the network for which this was the case
top_100_edges_enrich %>% 
  filter(TermID == `biosynthetic gene clusters`) %>% 
  count(network_method, sort=TRUE) %>% View()
top_100_edges_enrich %>% filter(str_detect(`Gene Ontology`, 'ribosom')) %>%
  separate_longer_delim(`Gene Ontology`, ';') %>%
  filter(TermID == `Gene Ontology`) %>%
  count(network_method, sort=TRUE) %>% View()
enrich_accuracy = top_100_edges_enrich %>% 
  mutate(`Gene Ontology` = if_else(annotation_category == 'biosynthetic gene clusters', 
                                   NA_character_, `Gene Ontology`)) %>%
  separate_longer_delim(`Gene Ontology`, ';') %>%
  filter((TermID == `biosynthetic gene clusters`) | (TermID == `Gene Ontology`), 
         annotation_category %in% c('biosynthetic gene clusters', 'Gene Ontology')) %>% 
  mutate(neg_log_pvalue = -log2(pvalue_adjusted),
         network_method = str_remove(network_method, 'get_')) %>%
  group_by(network_method, annotation_category) %>%
  mutate(med_neg_log_pvalue = median(neg_log_pvalue),
         sum_neg_log_pvalue = sum(neg_log_pvalue),
         count = n()) %>%
  ungroup() %>%
  mutate(network_method = fct_reorder(network_method, med_neg_log_pvalue, .desc=FALSE, .fun = max),
         annotation_category = str_to_title(annotation_category)) 
p_accuracy = enrich_accuracy %>%
  filter(annotation_category == 'Biosynthetic Gene Clusters') %>%
  ggplot(aes(network_method, neg_log_pvalue)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.data = function(x) c(y = -5, label = length(x)), 
               geom = "text", size=2.5) +
  facet_wrap(~annotation_category, nrow = 1) +
  coord_flip(ylim = c(-5, 200)) +
  theme_bw() +
  ggeasy::easy_remove_legend() +
  ggeasy::easy_remove_y_axis(what = c('title', 'text')) +
  ggtitle('Enrichment accuracy') +
  ggeasy::easy_center_title() 
df_accuracy = enrich_accuracy %>% 
  distinct(network_method) %>% 
  make_tabular_workflow() %>%
  arrange(desc(network_method)) 
tbl_theme = ttheme_minimal(base_size = 6.8, base_colour = "black", base_family = "",
                           parse = FALSE, padding = unit(c(.90, .90), "mm"))
tbl_grob = gridExtra::tableGrob(df_accuracy %>% select(-network_method), theme= tbl_theme, rows = NULL)
p = gridExtra::grid.arrange(tbl_grob, p_accuracy, ncol = 2, widths = c(2.5/10, 7.5/10))
ggsave(here('output/coexpression/enrichment_accuracy.png'), p, width = 14, height = 8, units='in', dpi = 300)

p_total_genes_from_same_term = enrich_accuracy %>% 
  filter(n_genes_subnetwork >= 4) %>%
  group_by(annotation_category, network_method) %>%
  summarize(n_genes_from_same_term_total = sum(n_genes_subnetwork)) %>%
  mutate(network_method = fct_reorder(network_method, n_genes_from_same_term_total)) %>%
  ggplot(aes(network_method, n_genes_from_same_term_total)) +
  geom_col() +
  facet_wrap(~annotation_category, scales = 'free_x') +
  coord_flip() +
  theme_bw() +
  ggeasy::easy_remove_legend() +
  ggeasy::easy_remove_y_axis(what = c('title', 'text')) +
  ggtitle('Total number of genes from from same term among top 100 edges for each gene') +
  ggeasy::easy_center_title() 
df_total_genes_from_same_term = enrich_accuracy %>% 
  filter(n_genes_subnetwork >= 4) %>%
  group_by(annotation_category, network_method) %>%
  summarize(n_genes_from_same_term_total = sum(n_genes_subnetwork)) %>%
  mutate(network_method = fct_reorder(network_method, n_genes_from_same_term_total)) %>%
  distinct(network_method) %>% 
  make_tabular_workflow() %>%
  arrange(desc(network_method)) 
tbl_theme = ttheme_minimal(base_size = 6.8, base_colour = "black", base_family = "",
                           parse = FALSE, padding = unit(c(.90, .90), "mm"))
tbl_grob = gridExtra::tableGrob(df_accuracy %>% select(-network_method), theme= tbl_theme, rows = NULL)
p = gridExtra::grid.arrange(tbl_grob, p_total_genes_from_same_term, ncol = 2, widths = c(2.5/10, 7.5/10))
ggsave(here('output/coexpression/enrichment_accuracy_total_genes.png'), p, width = 14, height = 8, units='in', dpi = 300)
enrich_accuracy %>% 
  filter(n_genes_subnetwork >= 4,
         annotation_category == 'Biosynthetic Gene Clusters') %>%
  group_by(network_method, TermID) %>% 
  summarize(total_genes_same_term = sum(n_genes_subnetwork)) %>%
  ungroup() %>%
  group_by(TermID) %>%
  mutate(variance = var(total_genes_same_term), 
         num_methods_detected = n()) %>%
  write_csv(here('output/coexpression/total_genes_from_each_BGC_per_method.csv'))

## 'get_downsampled_data-combat_seq_correct-vst_normalize-no_scaling-wgcna'
## 'get_nondownsampled_data-no_correction-upper_quartile_normalize-no_scaling-pearson'
## Above two methods are the best when ranking by accuracy (how often the top 100 genes
## are enriched for an annotation term that the query gene itself shares).
## Since BGCs are known to be co-expressed, used the BGC enrichment accuracy for ranking the
## network construction methods. 

gene_name_2_jcvi = read_csv("Z:/genomes/af3357_annotation/uniprot_filtered.csv") %>%
  select(uniprot_id = Entry, gene_name = `Gene Names`) %>%
  mutate(gene_id = str_replace(gene_name, '.*(AFLA_\\d{6}).*', '\\1'),
         gene_name = str_remove(gene_name, 'AFLA_\\d{6}') %>% str_trim()) %>%
  separate_longer_delim(gene_name, delim = ' ')
pannzer = read_tsv('Z:/genomes/af3357_annotation/A_flavus_3357_pannzer2.tsv') %>%
  filter(str_length(desc) < 6) %>% 
  distinct(desc, qpid) %>%
  deframe()
string_2_jcvi = read_tsv('Z:/misc_db/string_db/a_flavus/332952.protein.info.v12.0.txt') %>%
  mutate(string_id = str_remove(`#string_protein_id`, '332952.')) %>%
  select(string_id, gene_id = preferred_name) %>%
  mutate(gene_id = str_replace_all(gene_id, 
                                   gene_name_2_jcvi %>% 
                                     filter(gene_name != '') %>%
                                     select(gene_name, gene_id) %>% deframe()),
         gene_id = str_replace_all(gene_id, pannzer)) %>% 
  write_csv('Z:/misc_db/string_db/a_flavus/string_id_2_jcvi.csv')
string_db_long = read_delim('Z:/misc_db/string_db/a_flavus/332952.protein.physical.links.full.v12.0.txt',
                            delim = ' ') %>%
  rowid_to_column('id') %>% 
  pivot_longer(cols = protein1:protein2, names_to = 'pair_num', values_to = 'protein_id') %>% 
  mutate(protein_id = str_remove(protein_id, '332952.')) %>%
  left_join(string_2_jcvi, by=c('protein_id' = 'string_id')) %>% 
  select(-protein_id) 
string_db = string_db_long  %>%
  pivot_wider(id_cols = id, names_from = pair_num, values_from = gene_id, unused_fn = first) %>%
  select(-id) %>%
  write_csv('Z:/misc_db/string_db/a_flavus/332952.protein.physical.links.full.v12.0_jcvi.csv')
string_db_filtered = string_db %>%
  filter(combined_score > 700) %>% 
  select(from = protein1, to = protein2, combined_score) %>%
  mutate(temp_from = if_else(from < to, from, to),
         temp_to = if_else(to > from, to, from)) %>%
  select(from = temp_from, to = temp_to, combined_score)
sort_from_to = function(df){
  df %>%
    mutate(temp_from = if_else(from < to, from, to),
           temp_to = if_else(to > from, to, from)) %>%
    select(-from, -to) %>%
    rename(from = temp_from, to = temp_to) %>%
    relocate(from, to)
}
string_db_filtered %>%
  left_join(test_net, by=c('from', 'to')) %>% 
  filter(weight > 0) %>%
  summarize(sum_weight = sum(weight, na.rm = TRUE))
method_combinations = method_combinations %>%
  rowwise() %>%
  mutate(net_filtered = list(filter_cor_and_create_igraph(workflow, cors) %>% 
                               igraph::as_data_frame() %>% as_tibble() %>% sort_from_to()))
method_combinations = method_combinations %>% 
  rowwise() %>%
  mutate(stringdb_weight_sum = string_db_filtered %>% 
           left_join(net_filtered, by=c('from', 'to')) %>% 
           filter(weight > 0) %>%
           {sum(.$weight, na.rm=TRUE)})
method_combinations %>% 
  select(workflow, stringdb_weight_sum) %>% 
  write_csv(here('output/coexpression/stringdb_network_weight_sums.csv'))

string_eval_lm = lm(data = method_combinations %>% filter(!str_detect(workflow, 'wgcna')), 
                    formula = stringdb_weight_sum ~ downsampling + batch_correction + library_size_normalization + scaling)
gtsummary::tbl_regression(string_eval_lm) %>% gtsummary::as_gt() %>%
  gt::gtsave(file = here('output/coexpression/lm_stringdb_network_weight.png'))
enrich_freq_bgc_lm = lm(data = network_enrich_freq %>%
                         filter(annotation_category == 'biosynthetic gene clusters') %>%
                         left_join(method_combinations %>% select(downsampling:workflow), 
                                   by=c('network_method' = 'workflow')), 
                       formula = enrichment_frequency ~ downsampling + batch_correction + library_size_normalization + scaling + correlation)
gtsummary::tbl_regression(enrich_freq_bgc_lm) %>% gtsummary::as_gt() %>%
  gt::gtsave(file = here('output/coexpression/lm_enrichment_frequency_bgc.png'))
## Best results: downsampling not much effect, no batch correction, TMM, no scaling, wgcna
enrich_freq_go_lm = lm(data = network_enrich_freq %>%
                             filter(annotation_category == 'Gene Ontology') %>%
                             left_join(method_combinations %>% select(downsampling:workflow), 
                                       by=c('network_method' = 'workflow')), 
                           formula = enrichment_frequency ~ downsampling + batch_correction + library_size_normalization + scaling + correlation)
gtsummary::tbl_regression(enrich_freq_go_lm) %>% gtsummary::as_gt() %>%
  gt::gtsave(file = here('output/coexpression/lm_enrichment_frequency_go.png'))
## Best results: no downsampling, no batch correction, TPM or TMM, no scaling, wgcna
enrich_accuracy_single = enrich_accuracy %>% 
  distinct(network_method, annotation_category, count) %>% 
  mutate(enrichment_accuracy = count / nrow(genes_with_mean_tpm_above_1)) %>%
  write_csv(here('output/coexpression/enrichment_accuracy.csv'))
enrich_accuracy_bgc_lm = lm(data = enrich_accuracy_single %>%
                              filter(annotation_category == 'Biosynthetic Gene Clusters') %>%
                              left_join(method_combinations %>% select(downsampling:workflow) %>%
                                          mutate(workflow = str_remove(workflow, 'get_')), 
                                        by=c('network_method' = 'workflow')), 
                            formula = count ~ downsampling + batch_correction + library_size_normalization + scaling + correlation)
gtsummary::tbl_regression(enrich_accuracy_bgc_lm) %>% gtsummary::as_gt() %>%
  gt::gtsave(file = here('output/coexpression/lm_enrichment_accuracy_bgc.png'))
## Best results: no downsampling, no batch correction, VST and TMM, no scaling, pearson
enrich_accuracy_go_lm = lm(data = enrich_accuracy_single %>%
                              filter(annotation_category == 'Gene Ontology') %>%
                              left_join(method_combinations %>% select(downsampling:workflow) %>%
                                          mutate(workflow = str_remove(workflow, 'get_')), 
                                        by=c('network_method' = 'workflow')), 
                            formula = count ~ downsampling + batch_correction + library_size_normalization + scaling + correlation)
gtsummary::tbl_regression(enrich_accuracy_go_lm) %>% gtsummary::as_gt() %>%
  gt::gtsave(file = here('output/coexpression/lm_enrichment_accuracy_go.png'))
## Best results: no downsampling, batch correction helped a lot, TPM, all scaling helped, pearson