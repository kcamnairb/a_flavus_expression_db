library(igraph)
library(tidyverse)
library(edgeR)
library(here)
counts_jcvi = read_csv(here('output/salmon_quant/jcvi_counts.csv'))
counts_sums_jcvi = counts_jcvi %>% 
  pivot_longer(-gene_id, names_to='run', values_to='counts') %>%
  group_by(run) %>%
  summarize(total_counts = sum(counts)) 
functional_annotation_jcvi = read_csv(here('shiny_app/data/functional_annotation_jcvi.csv')) %>%
  dplyr::rename(gene_id = gene, `Gene Ontology` = `gene ontology`, `KEGG pathways` = `kegg pathways`, 
                `biosynthetic gene clusters` = smurf_and_known_metabolite, 
                `Subcellular localization (DeepLoc)` = deeploc_location, 
                `Interpro domains` = `interpro domains`) %>%
  mutate(`Gene Ontology` = str_remove_all(`Gene Ontology`, "'de novo'"))
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
metadata = read_csv(here('shiny_app/data/sra_metadata_filtered.csv'))
tpm_jcvi= read_csv(here('output', 'salmon_quant', 'jcvi', 'A_flavus_jcvi_tpm_unfiltered.csv'))
tpm_jcvi = tpm_jcvi %>%
  select(gene_id, all_of(metadata$run))
#downsampled_counts_jcvi = imap(counts_jcvi %>% select(all_of(samples_above_10m)), 
#                              ~downsample_counts_col(.x, .y, 10e6))
#downsampled_counts_jcvi = downsampled_counts_jcvi %>% purrr::reduce(full_join, by='genes')
#downsampled_counts_jcvi[is.na(downsampled_counts_jcvi)] = 0
#downsampled_counts_jcvi %>% write_csv(here('output/coexpression/downsampled_counts_jcvi.csv'))
downsampled_counts_jcvi = read_csv(here('output/coexpression/downsampled_counts_jcvi.csv'))
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
genes_with_mean_tpm_above_1 = tpm_jcvi %>%
  select(all_of(c('gene_id', samples_above_10m))) %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='tpm') %>%
  group_by(gene_id) %>%
  summarize(tpm_mean = mean(tpm), tpm_sum = sum(tpm), tpm_max = max(tpm),
            num_samples_above_1 = sum(tpm > 1)) %>% 
  ungroup() %>%
  filter(tpm_mean >= 1)
downsample_uq_jcvi = downsampled_counts_jcvi %>%
  semi_join(genes_with_mean_tpm_above_1) %>%
  dplyr::select(all_of(c('gene_id', samples_above_10m))) %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='downsampled_count') %>%
  group_by(sample_id) %>%
  filter(downsampled_count > 0) %>% # get rid of genes that have 0 counts
  mutate(uq = quantile(downsampled_count, 0.75)) %>% 
  ungroup() %>% 
  mutate(uq_count = (downsampled_count / uq)*1e6) %>%
  select(-uq, -downsampled_count) %>%
  pivot_wider(names_from = sample_id, values_from = uq_count, values_fill = 0)
downsample_uq_jcvi %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='uq') %>%
  group_by(gene_id) %>%
  summarize(uq_mean = mean(uq), uq_sum = sum(uq), uq_max = max(uq),
            num_samples_above_1 = sum(uq > 1)) %>% 
  ungroup() %>%
  filter(uq_mean < 1) %>% View()
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
log_downsample_uq_jcvi= downsample_uq_jcvi %>%
  mutate(across(-gene_id, ~log2(.x + 1)))
asinh_downsample_uq_jcvi= downsample_uq_jcvi %>%
  mutate(across(-gene_id, asinh))
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
normalized_inputs = list('uq_downsamp' = downsample_uq_jcvi,
                         'uq_downsamp_log' = log_downsample_uq_jcvi,
                         'uq_downsamp_asinh' = asinh_downsample_uq_jcvi,
                         'vst_downsamp' = downsample_vst_jcvi,
                         'zscore' = zscore_jcvi,
                         'tmm' = tmm_jcvi %>% rownames_to_column('gene_id'),
                         'tmm_asinh' = tmm_asinh_jcvi%>% rownames_to_column('gene_id'),
                         'downsample_tmm' = tmm_jcvi%>% rownames_to_column('gene_id'),
                         'downsample_tmm_asinh' = tmm_asinh_jcvi%>% rownames_to_column('gene_id'),
                         'tmm_log' = log_tmm_jcvi)
################### Make networks ###############################
cors = #normalized_inputs %>% 
  map(~.x %>% as.data.frame() %>% column_to_rownames('gene_id') %>% t() %>% cor())
#cors[['uq_downsamp_log']] = normalized_inputs[['uq_downsamp_log']] %>% as.data.frame() %>% column_to_rownames('gene_id') %>% t() %>% cor()

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
#edges[['uq_downsamp_log']] = cor2edges(cors[['uq_downsamp_log']], n_samples)
saveRDS(edges, here('output/coexpression/edges.rds'))
edges = readRDS(here('output/coexpression/edges.rds'))
AF_genes = c('AFLA_139150', 'AFLA_139160', 'AFLA_139170', 'AFLA_139180', 'AFLA_139190', 'AFLA_139200', 'AFLA_139210', 'AFLA_139220', 'AFLA_139230', 'AFLA_139240', 'AFLA_139250', 'AFLA_139260', 'AFLA_139270', 'AFLA_139280', 'AFLA_139290', 'AFLA_139300', 'AFLA_139310', 'AFLA_139320', 'AFLA_139330', 'AFLA_139340', 'AFLA_139360', 'AFLA_139370', 'AFLA_139380', 'AFLA_139390', 'AFLA_139400', 'AFLA_139410', 'AFLA_139420', 'AFLA_139430', 'AFLA_139440')
quantiles = edges %>% map(~.x %>% mutate(abs_r = abs(r)) %>% pull(abs_r) %>% quantile(probs = 1:100 / 100))
edges_filtered = edges %>% map(~.x %>% mutate(abs_r = abs(r)) %>% filter(abs_r > 0.5))
node_tables = edges_filtered %>% 
  map(~.x %>% {tibble(gene_id = c(.$from, .$to) %>% unique())} %>%
        left_join(functional_annotation_jcvi, by='gene_id'))
networks = map2(edges_filtered, node_tables, 
                ~graph_from_data_frame(.x, vertices = .y, directed = FALSE))
saveRDS(networks, here('output/coexpression/networks.rds'))

# https://github.com/cxli233/SimpleTidy_GeneCoEx#gene-wise-correlation
create_communities_at_multiple_resolutions = function(network, resolutions){
  resolutions %>% set_names() %>%
    map(~network %>%
    cluster_leiden(resolution_parameter = .x, objective_function = 'modularity') %>%
      {tibble(gene_id = names(membership(.)), module = as.vector(membership(.)))}
  )
}
modules_resolutions = networks %>% map(~create_communities_at_multiple_resolutions(.x, seq(2, 4, by = 0.25)))
saveRDS(modules_resolutions, here('output/coexpression/modules_resolutions.rds'))
# https://github.com/krishnanlab/RNAseq_coexpression/tree/main/gold_standards
go_specific = read_tsv(here('data/raw/gobp_specific_terms.txt'), col_names = c('go_id', 'go_description'))
go_intermediate = read_tsv(here('data/raw/gobp_intermediate_terms.txt'), col_names = c('go_id', 'go_description'))  
go_terms = functional_annotation_jcvi %>% dplyr::select(gene_id, go_description = `Gene Ontology`) %>%
  separate_longer_delim(go_description, ';')
go_terms %>% 
  mutate(specific = go_description %in% c(go_specific$go_description, go_intermediate$go_description)) %>%
  filter(specific) %>% distinct(go_description)
# Only 6 terms are present in the GO specific terms from the coexpression methods comparison paper
# Most go terms between that have between 5-20 genes annotated with them seem fairly specific and would 
# likely be co-expressed.
go_terms_specific = go_terms %>% add_count(go_description, name = 'n_go') %>% filter(n_go >= 5, n_go <= 100)
af_edges = edges_filtered  %>%
  imap(~.x %>%
        filter(((from %in% AF_genes) & (to %in% AF_genes)) |
                 ((to %in% AF_genes) & (from %in% AF_genes))) %>%
        mutate(r_mean = mean(r), r_max = max(r), r_min = min(r)) %>%
        mutate(normalization_type = .y)
  ) %>%
  bind_rows()
af_edges %>% 
  mutate(normalization_type = fct_reorder(normalization_type, r_mean, .desc = TRUE)) %>%
  ggplot(aes(normalization_type, r)) +
  geom_boxplot() +
  theme_bw() +
  ggeasy::easy_rotate_x_labels(side='right')
ggsave(here('output/coexpression/AF_genes_edges_r_values.png'))
quantiles %>% imap(~tibble(abs_r = .x) %>% mutate(normalization_type = .y)) %>% bind_rows() %>%
  group_by(normalization_type) %>%
  mutate(median_abs_r = median(abs_r)) %>%
  ungroup() %>%
  mutate(normalization_type = fct_reorder(normalization_type, median_abs_r, .desc = TRUE)) %>%
  ggplot(aes(normalization_type, abs_r)) +
  geom_boxplot() + 
  theme_bw() +
  ggeasy::easy_rotate_x_labels(side='right') 
score_modules = function(modules, resolution, normalization_type){
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
    summarise(sum = sum(n)) %>%
    pull(sum)
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
  go_overlap = go_terms_specific %>% 
    group_by(go_description) %>%
    summarize(go_size = n(), 
              module_with_most_overlap = list(modules %>% 
                                                dplyr::rename(gene_id2 = gene_id) %>%
                                                group_by(module) %>% 
                                                mutate(module_size = n(),
                                                       overlap = sum(gene_id2 %in% gene_id),
                                                       frac_genes_in_mod = overlap / go_size) %>%
                                                ungroup() %>%
                                                slice_max(overlap, with_ties = FALSE) %>%
                                                dplyr::select(module_2 = module, overlap, module_size, 
                                                              frac_genes_in_mod))) %>% 
    unnest_wider(module_with_most_overlap) 
  tibble(num_module= num_module_5, num_contained_gene = num_genes_contained,
         resolution = resolution, 
         most_af_genes_in_mod = num_af_gene_modules %>% slice_max(n) %>% pull(n),
         size_of_af_mod = modules %>% semi_join(num_af_gene_modules %>% slice_max(n)) %>% nrow(),
         bgc_overlap = list(bgc_overlap),
         go_overlap = list(go_overlap),
         normalization_type = normalization_type
  )
}

optimization_results= modules_resolutions %>% 
  imap(function(modules_set, normalization_type){
    modules_set %>% imap(function(module, resolution){score_modules(module, resolution, normalization_type)})
}) %>%
  flatten() %>% bind_rows()
bgc_optimization_results = optimization_results %>% 
  dplyr::select(normalization_type, resolution, bgc_overlap) %>% unnest(bgc_overlap) %>%
  bind_rows() 
go_optimization_results = optimization_results %>% 
  dplyr::select(normalization_type, resolution, go_overlap) %>% unnest(go_overlap) %>%
  bind_rows() 
patchwork::wrap_plots(
  optimization_results %>%
    ggplot(aes(resolution, num_module, group=resolution, fill=normalization_type)) +
    geom_col() +
    facet_wrap(~normalization_type, nrow = 1) +
    theme_bw() +
    ggeasy::easy_remove_legend() +
    scale_fill_manual(values = mpn65) +
    ggeasy::easy_remove_x_axis(),
  optimization_results %>%
    ggplot(aes(resolution, most_af_genes_in_mod, group=resolution, fill=normalization_type)) +
    geom_col() +
    facet_wrap(~normalization_type, nrow = 1) +
    theme_bw() +
    ggeasy::easy_remove_legend() +
    scale_fill_manual(values = mpn65) +
    ggeasy::easy_remove_x_axis(),
  go_optimization_results %>% 
    dplyr::rename(`Fraction of BGC in module` = frac_genes_in_mod) %>%
    ggplot(aes(resolution, `Fraction of BGC in module`, group=resolution, fill=normalization_type)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~normalization_type, nrow = 1) +
    theme_bw() +
    ggeasy::easy_remove_legend() +
    scale_fill_manual(values = mpn65) +
    ggeasy::easy_remove_x_axis(),
  bgc_optimization_results %>% 
    dplyr::rename(`Fraction of GO category in module` = frac_genes_in_mod) %>%
    ggplot(aes(resolution, `Fraction of GO category in module`, group=resolution, fill=normalization_type)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~normalization_type, nrow = 1) +
    theme_bw() +
    ggeasy::easy_remove_legend() +
    scale_fill_manual(values = mpn65) +
    ggeasy::easy_rotate_x_labels(side = 'right'),
  ncol=1
)
ggsave(here('output/coexpression/network_optimization.png'), width=14, height=8, units='in')
go_optimization_results %>% 
  group_by(normalization_type, resolution) %>% 
  summarize(mean_frac_go_in_mod = mean(frac_genes_in_mod), median_frac_go_in_mod = median(frac_genes_in_mod)) %>% 
  arrange(desc(mean_frac_go_in_mod)) 
#normalization_type resolution mean_frac_go_in_mod
#<chr>              <chr>                    <dbl>
# 1 downsample_tmm     3.75                     0.387
#2 tmm                3.5                      0.387
#3 downsample_tmm     3.25                     0.385
#4 downsample_tmm     2.75                     0.385
#5 downsample_tmm     3.5                      0.383
#6 tmm                3.25                     0.382
#7 downsample_tmm     2                        0.381
#8 tmm                2.25                     0.380
#9 tmm                3.75                     0.370
#10 tmm                4                        0.370
bgc_optimization_results %>% 
  group_by(normalization_type, resolution) %>% 
  summarize(mean_frac_bgc_in_mod = mean(frac_genes_in_mod), median_frac_bgc_in_mod = median(frac_genes_in_mod)) %>% 
  arrange(desc(mean_frac_bgc_in_mod))%>%  View()#print(n=15)
#normalization_type resolution mean_frac_bgc_in_mod median_frac_bgc_in_mod
#<chr>              <chr>                     <dbl>                  <dbl>
#1 zscore             3                         0.498                  0.471
#2 zscore             3.25                      0.488                  0.455
#3 zscore             4                         0.486                  0.455
#4 zscore             2.25                      0.484                  0.444
#5 zscore             3.5                       0.484                  0.444
#6 zscore             2.5                       0.482                  0.444
#7 zscore             3.75                      0.482                  0.444
#8 zscore             2.75                      0.481                  0.444
#9 zscore             2                         0.478                  0.444
#10 downsample_tmm     2.75                      0.468                  0.5  
#11 downsample_tmm     3.25                      0.466                  0.5  
#12 downsample_tmm     3.75                      0.466                  0.5  
#13 tmm                3.5                       0.466                  0.5  
#14 tmm                3.25                      0.465                  0.5  
#15 tmm                2.75                      0.463                  0.5 
# Downsampled TMM with a resolution of 2.75 scored the best overall for BGC and GO categories being in the same module
# Zscore also scored very well for BGC. 
modules = modules_resolutions$downsample_tmm$`2.75` # or modules_resolutions$zscore$`3` 
saveRDS(modules, here('output/coexpression/modules.rds'))
saveRDS(modules, here('shiny_app/data/modules.rds'))
network = networks$downsample_tmm
saveRDS(network, here('output/coexpression/network.rds'))
saveRDS(network, here('shiny_app/data/network.rds'))
random_optimization_results = modules_resolutions %>% 
  imap(function(modules_set, normalization_type){
    modules_set %>% imap(function(module, resolution){
      score_modules(module %>% mutate(module = sample(module)), 
                    resolution, normalization_type)})
  }) %>%
  flatten() %>% bind_rows() 
patchwork::wrap_plots(
  random_optimization_results %>%
    ggplot(aes(resolution, num_module, group=resolution, fill=normalization_type)) +
    geom_col() +
    facet_wrap(~normalization_type, nrow = 1) +
    theme_bw() +
    ggeasy::easy_remove_legend() +
    scale_fill_manual(values = mpn65) +
    ggeasy::easy_remove_x_axis(),
  random_optimization_results %>%
    ggplot(aes(resolution, most_af_genes_in_mod, group=resolution, fill=normalization_type)) +
    geom_col() +
    facet_wrap(~normalization_type, nrow = 1) +
    theme_bw() +
    ggeasy::easy_remove_legend() +
    scale_fill_manual(values = mpn65) +
    ggeasy::easy_remove_x_axis(),
  random_optimization_results %>% 
    dplyr::select(normalization_type, resolution, bgc_overlap) %>% unnest(bgc_overlap) %>%
    bind_rows() %>%
    dplyr::rename(`Fraction of BGC in module` = frac_genes_in_mod) %>%
    ggplot(aes(resolution, `Fraction of BGC in module`, group=resolution, fill=normalization_type)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~normalization_type, nrow = 1) +
    theme_bw() +
    ggeasy::easy_remove_legend() +
    scale_fill_manual(values = mpn65) +
    ggeasy::easy_remove_x_axis(),
  random_optimization_results %>% 
    dplyr::select(normalization_type, resolution, go_overlap) %>% unnest(go_overlap) %>%
    bind_rows() %>%
    dplyr::rename(`Fraction of GO category in module` = frac_genes_in_mod) %>%
    ggplot(aes(resolution, `Fraction of GO category in module`, group=resolution, fill=normalization_type)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~normalization_type, nrow = 1) +
    theme_bw() +
    ggeasy::easy_remove_legend() +
    scale_fill_manual(values = mpn65) +
    ggeasy::easy_rotate_x_labels(side = 'right'),
  ncol=1
)
# Randomizing the modules reduced the scores by about one half.
ggsave(here('output/coexpression/network_optimization_randomized.png'), width=14, height=8, units='in')
modules_expression = modules %>% 
  left_join(normalized_inputs$downsample_tmm %>% as.tibble(), by='gene_id')
zscore_plot = modules %>% 
  left_join(normalized_inputs$zscore, by='gene_id') %>%
  pivot_longer(-c(gene_id, module), names_to = 'sample_id', values_to = 'expression') %>%
  ggplot(aes(sample_id, expression, group = gene_id)) +
  geom_line(alpha = 0.3, color = "grey70") +
  facet_wrap(~module, ncol = 1)
ggsave(here('output/coexpression/modules_zscore_plot.png'), plot = zscore_plot, height=16, width = 16, units='in', dpi = 300)  
zscore_plot = modules %>% 
  left_join(normalized_inputs$vst, by='gene_id') %>%
  pivot_longer(-c(gene_id, module), names_to = 'sample_id', values_to = 'expression') %>%
  ggplot(aes(sample_id, expression, group = gene_id)) +
  geom_line(alpha = 0.3, color = "grey70") +
  facet_wrap(~module, ncol = 1)
ggsave(here('output/coexpression/modules_vst_plot.png'), plot = zscore_plot, height=16, width = 16, units='in', dpi = 300)  


