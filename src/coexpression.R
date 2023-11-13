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
# UQ also recommended in this paper evaluating different normalization methods for coexpression analysis: Johnson, Kayla A., and Arjun Krishnan. "Robust normalization and transformation techniques for constructing gene coexpression networks from RNA-seq data." Genome biology 23 (2022): 1-26.
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
tpm_jcvi %>%
  select(all_of(c('gene_id', samples_above_10m))) %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='tpm') %>%
  group_by(gene_id) %>%
  summarise(var = var(log1p(tpm)), 
            gene_id %in% genes_with_mean_tpm_above_1$gene_id) %>%  
  ungroup() %>%
  mutate(var > quantile(var, 0.667))  %>% View()
  ungroup() %>% 
  filter(var > quantile(var, 0.667))
tpm_jcvi %>% select(gene_id, tpm=DRR452451) %>%
  full_join(downsample_uq_jcvi %>% select(gene_id, uq=DRR452451)) %>%
  ggplot(aes(tpm, uq)) +
  geom_point(alpha=0.1) +
  coord_cartesian(xlim=c(0, 500), ylim=c(0, 25)) +
  geom_smooth() +
  theme_bw()
############## Gather inputs for constructing networks #############
downsample_uq_jcvi_sub = downsample_uq_jcvi %>% semi_join(genes_with_mean_tpm_above_1)
vst_jcvi_sub = vsd_jcvi %>% semi_join(genes_with_mean_tpm_above_1) %>%
  dplyr::select(all_of(c('gene_id', samples_above_10m))) 
downsample_vst_jcvi = DESeq2::vst(downsampled_counts_jcvi %>% column_to_rownames('gene_id') %>%
                                    as.matrix(), blind=TRUE) %>%
  as.data.frame() %>% rownames_to_column('gene_id') %>% 
  semi_join(genes_with_mean_tpm_above_1)
zscore_jcvi = tpm_jcvi %>%
  semi_join(genes_with_mean_tpm_above_1) %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to = 'tpm') %>%
  mutate(log_tpm = log10(tpm + 1)) %>%
  group_by(gene_id) %>%
  mutate(z_score = (log_tpm - mean(log_tpm))/ sd(log_tpm)) %>%
  ungroup() %>% 
  select(-log_tpm, -tpm) %>% 
  pivot_wider(names_from = 'sample_id', values_from = 'z_score') %>%
  dplyr::select(all_of(c('gene_id', samples_above_10m)))
downsample_tmm_jcvi = downsampled_counts_jcvi %>%
  dplyr::select(all_of(c('gene_id', samples_above_10m))) %>%
  column_to_rownames('gene_id') %>%
  DGEList()
downsample_tmm_jcvi = calcNormFactors(downsample_tmm_jcvi, method = 'TMM')
downsample_tmm_jcvi = cpm(downsample_tmm_jcvi, log=FALSE)
downsample_tmm_jcvi = downsample_tmm_jcvi %>% as.data.frame() %>% rownames_to_column('gene_id') %>%
  semi_join(genes_with_mean_tpm_above_1) %>% column_to_rownames('gene_id')
downsample_tmm_asinh_jcvi = downsample_tmm_jcvi %>% asinh() # hyperbolic arcsine recommended by Johnson, Kayla A., and Arjun Krishnan. "Robust normalization and transformation techniques for constructing gene coexpression networks from RNA-seq data." Genome biology 23 (2022): 1-26.
tmm_jcvi = counts_jcvi %>%
  dplyr::select(all_of(c('gene_id', samples_above_10m))) %>%
  column_to_rownames('gene_id') %>%
  DGEList()
tmm_jcvi = calcNormFactors(tmm_jcvi, method = 'TMM') # Only changes norm.factors column
log_tmm_jcvi = cpm(tmm_jcvi, log=TRUE) %>% as.data.frame() %>% rownames_to_column('gene_id') %>%
  semi_join(genes_with_mean_tpm_above_1)
tmm_jcvi = cpm(tmm_jcvi, log=FALSE) # Uses normalized library sizes
tmm_jcvi = tmm_jcvi %>% as.data.frame() %>% rownames_to_column('gene_id') %>%
  semi_join(genes_with_mean_tpm_above_1) %>% column_to_rownames('gene_id')
tmm_asinh_jcvi = tmm_jcvi %>% asinh()
normalized_inputs = list('uq_downsamp' = downsample_uq_jcvi_sub,
                         'vst' = vst_jcvi_sub,
                         'vst_downsamp' = downsample_vst_jcvi,
                         'zscore' = zscore_jcvi,
                         'tmm' = tmm_jcvi %>% rownames_to_column('gene_id'),
                         'tmm_asinh' = tmm_asinh_jcvi%>% rownames_to_column('gene_id'),
                         'downsample_tmm' = tmm_jcvi%>% rownames_to_column('gene_id'),
                         'downsample_tmm_asinh' = tmm_asinh_jcvi%>% rownames_to_column('gene_id'),
                         'tmm_log' = log_tmm_jcvi)
################### Make networks ###############################
cors = #normalized_inputs %>% 
  list('tmm_log' = log_tmm_jcvi) %>%
  map(~.x %>% as.data.frame() %>% column_to_rownames('gene_id') %>% t() %>% cor())

n_samples = downsample_uq_jcvi_sub %>% ncol() - 1
cor2edges = function(cor_matrix, n_samples){
  cor_upper_tri = cor_matrix
  cor_upper_tri[lower.tri(cor_upper_tri)] = NA
  cor_upper_tri %>% 
    as.data.frame() %>% 
    mutate(from = row.names(cor_matrix)) %>% 
    pivot_longer(-from, names_to = "to", values_to = "r") %>% 
    filter(!is.na(r)) %>% 
    filter(from != to) %>% 
    mutate(t = r*sqrt((n_samples-2)/(1-r^2))) %>% 
    mutate(p.value = case_when(
      t > 0 ~ pt(t, df = n_samples-2, lower.tail = F),
      t <=0 ~ pt(t, df = n_samples-2, lower.tail = T)
    )) %>% 
    mutate(FDR = p.adjust(p.value, method = "fdr")) 
}
edges = cors %>% map(~cor2edges(.x, n_samples))
edges[['tmm_log']] = cor2edges(cors[['tmm_log']], n_samples)
saveRDS(edges, here('output/coexpression/edges.rds'))
AF_genes = c('AFLA_139150', 'AFLA_139160', 'AFLA_139170', 'AFLA_139180', 'AFLA_139190', 'AFLA_139200', 'AFLA_139210', 'AFLA_139220', 'AFLA_139230', 'AFLA_139240', 'AFLA_139250', 'AFLA_139260', 'AFLA_139270', 'AFLA_139280', 'AFLA_139290', 'AFLA_139300', 'AFLA_139310', 'AFLA_139320', 'AFLA_139330', 'AFLA_139340', 'AFLA_139360', 'AFLA_139370', 'AFLA_139380', 'AFLA_139390', 'AFLA_139400', 'AFLA_139410', 'AFLA_139420', 'AFLA_139430', 'AFLA_139440')
quantiles = edges %>% map(~.x %>% mutate(abs_r = abs(r)) %>% pull(abs_r) %>% quantile(probs = 1:100 / 100))
edges_filtered = edges %>% map(~.x %>% mutate(abs_r = abs(r)) %>% filter(abs_r > 0.5))
node_tables = edges_filtered %>% 
  map(~.x %>% {tibble(gene_id = c(.$from, .$to) %>% unique())} %>%
        left_join(functional_annotation_jcvi, by='gene_id'))
networks = map2(edges_filtered, node_tables, 
                ~graph_from_data_frame(.x, vertices = .y, directed = FALSE))
modules = networks %>% 
  map(~cluster_leiden(.x, resolution_parameter = 2, objective_function = 'modularity') %>%
  {tibble(gene_id = names(membership(.)), module = as.vector(membership(.)))})
# https://github.com/cxli233/SimpleTidy_GeneCoEx#gene-wise-correlation
optimize_resolution = function(network, resolution, normalization_type){
  modules = network %>%
    cluster_leiden(resolution_parameter = resolution, objective_function = 'modularity') %>%
    {tibble(gene_id = names(membership(.)), module = as.vector(membership(.)))}
  num_module_5 = modules %>% 
    dplyr::count(module) %>% 
    arrange(-n) %>% 
    filter(n >= 5) %>% 
    nrow() 
  num_genes_contained = modules %>%  
    dplyr::count(module) %>% 
    arrange(-n) %>% 
    filter(n >= 5) %>% 
    ungroup() %>% 
    summarise(sum = sum(n)) 
  num_af_gene_modules = modules %>%
    filter(gene_id %in% AF_genes) %>%
    dplyr::count(module)
  bgc_overlap = functional_annotation_jcvi %>% 
    filter(!is.na(`biosynthetic gene clusters`)) %>% 
    group_by(`biosynthetic gene clusters`) %>%
    summarize(bgc_size = n(), 
           module_with_most_overlap = list(modules %>% 
                                             dplyr::rename(gene_id2 = gene_id) %>%
                                             group_by(module) %>% 
                                             mutate(module_size = n(),
                                                    overlap = sum(gene_id2 %in% gene_id),
                                                    frac_genes_in_mod = overlap / bgc_size) %>%
                                             ungroup() %>%
                                             slice_max(overlap, with_ties = FALSE) %>%
                                             dplyr::select(module_2 = module, overlap, module_size, 
                                                           frac_genes_in_mod))) %>% 
    unnest_wider(module_with_most_overlap) 
  tibble(num_module= num_module_5, num_contained_gene = num_genes_contained, num_edges = ecount(network),
         resolution = resolution, 
         most_af_genes_in_mod = num_af_gene_modules %>% slice_max(n) %>% pull(n),
         size_of_af_mod = modules %>% semi_join(num_af_gene_modules %>% slice_max(n)) %>% nrow(),
         bgc_overlap = list(bgc_overlap),
         normalization = normalization_type
  )
}
optimization_results = networks %>% imap(function(network, normalization_type){
  map(seq(from = 0.25, to = 5, by = 0.25), 
      ~optimize_resolution(network, .x, normalization_type))}) %>%
  bind_rows()

#uq_vst_comparison = modules %>% 
#  group_by(module) %>% 
#  summarise(module_1_size = n(), 
#            module_with_most_overlap = list(modules_vst %>% 
#                                              dplyr::rename(gene_id_vst = gene_id) %>%
#                                              group_by(module) %>% 
#                                              mutate(module_2_size = n(),
#                                                     overlap = sum(gene_id_vst %in% gene_id)) %>%
#                                              ungroup() %>%
#                                              slice_max(overlap, with_ties = FALSE) %>%
#                                              dplyr::select(module_2 = module, overlap, module_2_size))) %>% 
#  unnest_wider(module_with_most_overlap)

bgc_optimization_results = optimization_results %>% 
  dplyr::select(normalization, resolution, bgc_overlap) %>% unnest(bgc_overlap) %>%
  bind_rows() 
bgc_optimization_results %>% 
  group_by(normalization) %>% 
  filter(resolution >= 2, resolution <=3, bgc_size>=5) %>%
  #summarize(mean_frac_bgc_in_mod = mean(frac_genes_in_mod))
  ggplot(aes(resolution, frac_genes_in_mod, group=resolution, color=`biosynthetic gene clusters`)) +
  geom_boxplot() +
  geom_point(alpha=0.2) +
  facet_wrap(~normalization, nrow = 1) +
  theme_bw() +
  ggeasy::easy_remove_legend()
  list('uq' = optimization_results, 'zscore' = optimization_results_zscore, 
  'vst' = optimization_results_vst, 'tmm' = optimization_results_tmm, 
  'tmm_asinh' = optimization_results_tmmasinh) %>%
    imap(~.x %>% mutate(normalization = .y)) %>%
    bind_rows() %>% 
    ggplot(aes(resolution, num_module, group=normalization, fill=normalization)) +
    geom_col(position = 'dodge')
  optimization_results %>%
    filter(resolution >= 2, resolution <=3) %>%
    ggplot(aes(resolution, num_module, fill = normalization)) +
    geom_col(position = 'dodge')
 optimization_results %>% 
    pivot_longer()
    group_by(normalization) %>% 
    filter(resolution >= 2, resolution <=3, bgc_size>=5) %>%