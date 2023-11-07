library(igraph)
#counts_jcvi[[2]]
#counts_sums_jcvi
downsample_counts_col = function(counts_col, sample_name, sample_num){
  genes_rep = counts_col %>% 
    set_names(counts_jcvi$gene_id) %>%
    imap(~rep(.y, .x)) %>% 
    unlist() %>%
    sample(sample_num)
  tibble(gene_id = genes_rep) %>%
    dplyr::count(gene_id, name = sample_name)
}
samples_above_10m  = counts_sums_jcvi %>% 
  filter(total_counts >= 10e6, run %in% colnames(tpm_jcvi)) %>% pull(run)
downsampled_counts_jcvi = imap(counts_jcvi %>% select(all_of(samples_above_10m)), 
                              ~downsample_counts_col(.x, .y, 10e6))
downsampled_counts_jcvi = downsampled_counts_jcvi %>% purrr::reduce(full_join, by='genes')
downsampled_counts_jcvi[is.na(downsampled_counts_jcvi)] = 0
downsampled_counts_jcvi %>% write_csv(here('output/coexpression/downsampled_counts_jcvi.csv'))
downsampled_counts_jcvi %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='downsampled_count') %>%
  filter(sample_id == 'DRR452438') %>%
  left_join(counts_jcvi %>%
              pivot_longer(-gene_id, names_to = 'sample_id', values_to='original_count') %>%
              filter(sample_id == 'DRR452438')) %>%
  ggplot(aes(downsampled_count, original_count)) +
  geom_point(alpha=0.2) +
  geom_abline(slope = 1, intercept=0) +
  theme_bw()
downsampled_by_div_counts_jcvi = counts_jcvi %>%
  mutate(across(all_of(samples_above_10m),
                ~.x / (sum(.x)/10e6)))
downsampled_counts_jcvi %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='downsampled_count') %>%
  filter(sample_id == 'DRR452438') %>%
  left_join(downsampled_by_div_counts_jcvi %>%
              pivot_longer(-gene_id, names_to = 'sample_id', values_to='downsampled_by_div_count') %>%
              filter(sample_id == 'DRR452438')) %>%
  ggplot(aes(downsampled_count, downsampled_by_div_count)) +
  geom_point(alpha=0.2) +
  geom_abline(slope = 1, intercept=0) +
  theme_bw()
## The downsampleing method results in almost the same outcome as dividing each count and is much slower.
# UQ + ComBat + Pearson performed the best in Vandenbon, A. (2022). Evaluation of critical data processing steps for reliable prediction of gene co-expression from large collections of RNA-seq data. Plos one, 17(1), e0263344.
#https://bioinformatics.stackexchange.com/questions/2586/how-to-apply-upperquartile-normalization-on-rsem-expected-counts
downsample_uq_jcvi = downsampled_counts_jcvi %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='downsampled_count') %>%
  group_by(sample_id) %>%
  filter(downsampled_count > 0) %>% # get rid of genes that have 0 counts
  mutate(uq = quantile(downsampled_count, 0.75)) %>% 
  ungroup() %>% 
  mutate(uq_count = downsampled_count / uq) %>%
  select(-uq, -downsampled_count) %>%
  pivot_wider(names_from = sample_id, values_from = uq_count, values_fill = 0)
downsample_uq_jcvi %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='uq') %>%
  group_by(gene_id) %>%
  summarize(uq_mean = mean(uq), uq_sum = sum(uq), uq_max = max(uq),
         num_samples_above_1 = sum(uq > 1)) %>% 
  ungroup() %>%
  filter(uq_mean < 1) %>% View()
genes_with_mean_tpm_above_1 = tpm_jcvi %>%
  select(all_of(c('gene_id', samples_above_10m))) %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='tpm') %>%
  group_by(gene_id) %>%
  summarize(tpm_mean = mean(tpm), tpm_sum = sum(tpm), tpm_max = max(tpm),
            num_samples_above_1 = sum(tpm > 1)) %>% 
  ungroup() %>%
  filter(tpm_mean >= 1)
downsampled_counts_jcvi %>%
  group_by(gene_id) %>%
  summarise(var = var(mean.logTPM)) %>% 
  ungroup() %>% 
  filter(var > quantile(var, 0.667))
tpm_jcvi %>% select(gene_id, tpm=DRR452451) %>%
  full_join(downsample_uq_jcvi %>% select(gene_id, uq=DRR452451)) %>%
  ggplot(aes(tpm, uq)) +
  geom_point(alpha=0.1) +
  coord_cartesian(xlim=c(0, 500), ylim=c(0, 25)) +
  geom_smooth() +
  theme_bw()
downsample_uq_jcvi_sub = downsample_uq_jcvi %>% semi_join(genes_with_mean_tpm_above_1)
uq_cor = downsample_uq_jcvi_sub %>% as.data.frame() %>% column_to_rownames('gene_id') %>% t() %>% cor()
uq_cor_upper_tri = uq_cor
uq_cor_upper_tri[lower.tri(uq_cor_upper_tri)] = NA
n_samples = downsample_uq_jcvi_sub %>% ncol() - 1
edges = uq_cor_upper_tri %>% 
  as.data.frame() %>% 
  mutate(from = row.names(uq_cor)) %>% 
  pivot_longer(-from, names_to = "to", values_to = "r") %>% 
  filter(!is.na(r)) %>% 
  filter(from != to) %>% 
  mutate(t = r*sqrt((n_samples-2)/(1-r^2))) %>% 
  mutate(p.value = case_when(
    t > 0 ~ pt(t, df = n_samples-2, lower.tail = F),
    t <=0 ~ pt(t, df = n_samples-2, lower.tail = T)
  )) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr")) 

AF_genes = c('AFLA_139150', 'AFLA_139160', 'AFLA_139170', 'AFLA_139180', 'AFLA_139190', 'AFLA_139200', 'AFLA_139210', 'AFLA_139220', 'AFLA_139230', 'AFLA_139240', 'AFLA_139250', 'AFLA_139260', 'AFLA_139270', 'AFLA_139280', 'AFLA_139290', 'AFLA_139300', 'AFLA_139310', 'AFLA_139320', 'AFLA_139330', 'AFLA_139340', 'AFLA_139360', 'AFLA_139370', 'AFLA_139380', 'AFLA_139390', 'AFLA_139400', 'AFLA_139410', 'AFLA_139420', 'AFLA_139430', 'AFLA_139440')
edges %>%
  filter(r > 0) %>% 
  filter(FDR < 0.05) %>% 
  slice_min(order_by = abs(r), n = 10)
edges %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.01) %>% 
  slice_min(order_by = abs(r), n = 10)
edges %>% 
  filter((from %in% AF_genes) & (to %in% AF_genes))
edges %>% 
  slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.2, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )
quantiles = edges %>% mutate(abs_r = abs(r)) %>% pull(abs_r) %>% quantile(probs = 1:100 / 100)
edges_filtered = edges %>% mutate(abs_r = abs(r)) %>% filter(abs_r > 0.5) 
node_table = tibble(gene_id = c(edges_filtered$from, edges_filtered$to) %>% unique()) %>%
  left_join(functional_annotation_jcvi, by='gene_id')
                             