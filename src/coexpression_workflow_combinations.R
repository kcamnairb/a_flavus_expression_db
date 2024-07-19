library(igraph)
library(edgeR)
library(WGCNA)
library(sva)
library(here)
library(mmtable2)
library(patchwork)
library(grid)
library(gridExtra)
library(WGCNA)
library(fst)
library(tidyverse)
mpn65 = c('#ff0029','#377eb8','#66a61e','#984ea3','#00d2d5','#ff7f00','#af8d00','#7f80cd','#b3e900','#c42e60','#a65628',
          '#f781bf','#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#fccde5','#bc80bd','#ffed6f','#c4eaff','#cf8c00',
          '#1b9e77','#d95f02','#e7298a','#e6ab02','#a6761d','#0097ff','#00d067','#000000','#252525','#525252','#737373',
          '#969696','#bdbdbd','#f43600','#4ba93b','#5779bb','#927acc','#97ee3f','#bf3947','#9f5b00','#f48758','#8caed6',
          '#f2b94f','#eff26e','#e43872','#d9b100','#9d7a00','#698cff','#d9d9d9','#00d27e','#d06800','#009f82','#c49200',
          '#cbe8ff','#fecddf','#c27eb6','#8cd2ce','#c4b8d9','#f883b0','#a49100','#f48800','#27d0df','#a04a9b')
metadata = read_csv('shiny_app/data/sra_metadata_filtered.csv')
counts_jcvi = read_csv(here('output/star/jcvi_counts.csv'))
counts_sums_jcvi = counts_jcvi %>% 
  pivot_longer(-gene_id, names_to='run', values_to='counts') %>%
  group_by(run) %>%
  summarize(total_counts = sum(counts)) 
downsample_counts_col = function(counts_col, sample_name, gene_ids, sample_num){
  genes_rep = counts_col %>% 
    set_names(gene_ids) %>%
    imap(~rep(.y, .x)) %>% 
    unlist() %>%
    sample(sample_num)
  tibble(gene_id = genes_rep) %>%
    dplyr::count(gene_id, name = sample_name)
}
samples_above_10m  = counts_sums_jcvi %>% 
  filter(total_counts >= 10e6) %>% pull(run)
length(samples_above_10m)
transcript_lengths_jcvi = read_csv('Z:/genomes/af3357_annotation/af3357_cds_lengths.csv', 
                                   col_names = c('gene_id', 'length'))
tpm_normalize = function(df, transcript_lengths_df = transcript_lengths_jcvi) {
  tpm_one_col = function(counts, gene_length)(counts / gene_length) / (sum(counts / gene_length) / 10^6)
  df %>%
    pivot_longer(-gene_id, names_to='sample_id', values_to='counts') %>%
    left_join(transcript_lengths_df, by='gene_id') %>%
    group_by(sample_id) %>%
    mutate(tpm = tpm_one_col(counts, length)) %>%
    ungroup() %>%
    select(-counts, -length) %>%
    pivot_wider(names_from=sample_id, values_from=tpm) 
}
tpm_jcvi = counts_jcvi %>% tpm_normalize(transcript_lengths_jcvi)
genes_with_mean_tpm_above_1 = tpm_jcvi %>%
  select(all_of(c('gene_id', samples_above_10m))) %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='tpm') %>%
  group_by(gene_id) %>%
  summarize(tpm_mean = mean(tpm), tpm_sum = sum(tpm), tpm_max = max(tpm),
            num_samples_above_1 = sum(tpm > 1)) %>% 
  ungroup() %>%
  filter(tpm_mean >= 1)
# downsampled_counts_jcvi = imap(counts_jcvi %>% select(all_of(samples_above_10m)),
#                              ~downsample_counts_col(.x, .y, counts_jcvi$gene_id, 10e6))
# downsampled_counts_jcvi = downsampled_counts_jcvi %>% purrr::reduce(full_join, by='gene_id')
# downsampled_counts_jcvi[is.na(downsampled_counts_jcvi)] = 0
# downsampled_counts_jcvi %>% write_csv(here('output/coexpression/downsampled_counts_jcvi.csv'))
downsampled_counts_jcvi = read_csv(here('output/coexpression/downsampled_counts_jcvi.csv'))
counts_chrom_level = read_csv(here('output/star_chrom_level/chrom_level_counts.csv'))
# downsampled_counts_chrom_level = imap(counts_chrom_level %>% select(all_of(samples_above_10m)),
#                              ~downsample_counts_col(.x, .y, counts_chrom_level$gene_id, 10e6))
# downsampled_counts_chrom_level = downsampled_counts_chrom_level %>% purrr::reduce(full_join, by='gene_id')
# downsampled_counts_chrom_level[is.na(downsampled_counts_chrom_level)] = 0
# downsampled_counts_chrom_level %>% write_csv(here('output/coexpression/downsampled_counts_chrom_level.csv'))
downsampled_counts_chrom_level = read_csv(here('output/coexpression/downsampled_counts_chrom_level.csv'))
transcript_lengths_chrom_level = read_csv(
  'Z:/genomes/A_flavus_3357_chromosome_level_JGI/download_from_ncbi/transcript_lengths.csv')
tpm_chrom_level = counts_chrom_level %>% tpm_normalize(transcript_lengths_chrom_level)
genes_with_mean_tpm_above_1 = tpm_jcvi %>%
  select(all_of(c('gene_id', samples_above_10m))) %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='tpm') %>%
  group_by(gene_id) %>%
  summarize(tpm_mean = mean(tpm), tpm_sum = sum(tpm), tpm_max = max(tpm),
            num_samples_above_1 = sum(tpm > 1)) %>% 
  ungroup() %>%
  filter(tpm_mean >= 1)
genes_with_mean_tpm_above_1_chrom_level = tpm_chrom_level %>%
  select(all_of(c('gene_id', samples_above_10m))) %>%
  pivot_longer(-gene_id, names_to = 'sample_id', values_to='tpm') %>%
  group_by(gene_id) %>%
  summarize(tpm_mean = mean(tpm), tpm_sum = sum(tpm), tpm_max = max(tpm),
            num_samples_above_1 = sum(tpm > 1)) %>% 
  ungroup() %>%
  filter(tpm_mean >= 1)
genes_with_mean_tpm_above_1 %>% select(gene_id) %>%
  bind_rows(genes_with_mean_tpm_above_1_chrom_level %>% select(gene_id)) %>%
  write_fst(here('shiny_app/data/gene_in_networks.fst'))
object_storage = list()
check_storage = function(inp, func){
  inp_func_hash = rlang::hash(list(inp, func))
  return(object_storage[[inp_func_hash]])
}
update_storage = function(inp, func, res){
  ## creates hash of the input and function, then updates object_storage with the return value 
  inp_func_hash = rlang::hash(list(inp, func))
  object_storage[[inp_func_hash]] <<- res
}
get_downsampled_data = function(){
  df = downsampled_counts_jcvi %>% dplyr::select(all_of(c('gene_id', samples_above_10m))) 
  runs_in_bioprojects_w_more_than_one_sample = metadata %>% 
    filter(run %in% colnames(df)) %>%
    add_count(bioproject) %>% 
    filter(n > 1) %>% 
    pull(run)
  df %>% dplyr::select(all_of(c('gene_id', runs_in_bioprojects_w_more_than_one_sample)))
}
get_nondownsampled_data = function(){
  df = counts_jcvi %>% dplyr::select(all_of(c('gene_id', samples_above_10m))) 
  runs_in_bioprojects_w_more_than_one_sample = metadata %>% 
    filter(run %in% colnames(df)) %>%
    add_count(bioproject) %>% 
    filter(n > 1) %>% 
    pull(run)
  df %>% dplyr::select(all_of(c('gene_id', runs_in_bioprojects_w_more_than_one_sample)))
}
get_downsampled_data_chrom_level = function(){
  df = downsampled_counts_chrom_level %>% dplyr::select(all_of(c('gene_id', samples_above_10m))) 
  runs_in_bioprojects_w_more_than_one_sample = metadata %>% 
    filter(run %in% colnames(df)) %>%
    add_count(bioproject) %>% 
    filter(n > 1) %>% 
    pull(run)
  df %>% dplyr::select(all_of(c('gene_id', runs_in_bioprojects_w_more_than_one_sample)))
}
get_nondownsampled_data_chrom_level = function(){
  df = counts_chrom_level %>% dplyr::select(all_of(c('gene_id', samples_above_10m))) 
  runs_in_bioprojects_w_more_than_one_sample = metadata %>% 
    filter(run %in% colnames(df)) %>%
    add_count(bioproject) %>% 
    filter(n > 1) %>% 
    pull(run)
  df %>% dplyr::select(all_of(c('gene_id', runs_in_bioprojects_w_more_than_one_sample)))
}
combat_seq_correct = function(df){
  storage_res = check_storage(df, deparse(sys.call()))
  if (!is.null(storage_res)) return(storage_res)
  res = ComBat_seq(counts = df %>% column_to_rownames('gene_id') %>% as.matrix(),
             batch = metadata %>% filter(run %in% colnames(df)) %>% pull(bioproject)) %>%
    as.data.frame() %>% rownames_to_column('gene_id')
  update_storage(df, deparse(sys.call()), res)
  return(res)
}
tmm_normalize = function(df){
  tmm = df %>%
    column_to_rownames('gene_id') %>%
    DGEList()
  tmm = calcNormFactors(tmm, method = 'TMM') # Only changes norm.factors column
  tmm = cpm(tmm, log=FALSE) %>% 
    as.data.frame() %>%
    rownames_to_column('gene_id') # Uses normalized library sizes
  return(tmm)
}
upper_quartile_normalize = function(df) {
  df %>%
    pivot_longer(-gene_id, names_to = 'sample_id', values_to='count') %>%
    group_by(sample_id) %>%
    filter(count > 0) %>% # get rid of genes that have 0 counts
    mutate(uq = quantile(count, 0.75)) %>% 
    ungroup() %>% 
    mutate(uq_count = (count / uq)*1e6) %>%
    select(-uq, -count) %>%
    pivot_wider(names_from = sample_id, values_from = uq_count, values_fill = 0)
}
vst_normalize = function(df){
  DESeq2::vst(df %>% column_to_rownames('gene_id') %>%
                as.matrix(), blind=TRUE) %>%
    as.data.frame() %>% rownames_to_column('gene_id') 
}
zscore_scale = function(df){
  df %>%
    pivot_longer(-gene_id, names_to = 'sample_id', values_to = 'tpm') %>%
    mutate(log_tpm = log10(tpm + 1)) %>%
    group_by(gene_id) %>%
    mutate(z_score = (log_tpm - mean(log_tpm))/ sd(log_tpm)) %>%
    ungroup() %>% 
    select(-log_tpm, -tpm) %>% 
    pivot_wider(names_from = 'sample_id', values_from = 'z_score')
}
no_scaling = no_correction = function(df) df
log_plus_one = function(df) df %>% mutate(across(-gene_id, log1p))
inverse_hyperbolic_sine = function(df) df %>% mutate(across(-gene_id, asinh))
pearson = function(df) df%>% as.data.frame() %>% column_to_rownames('gene_id') %>% t() %>% cor()
get_top_100_edges = function(cor_matrix, gene_id){
  rownames(cor_matrix)[order(cor_matrix[,gene_id], decreasing = TRUE)[1:100]]
}
get_top_100_edges_for_all_genes = function(cor_matrix){
  tibble(gene_id = colnames(cor_matrix)) %>%
    rowwise() %>%
    mutate(top100 = list(get_top_100_edges(cor_matrix, gene_id)))
}
wgcna = function(df){
  options(stringsAsFactors = FALSE)
  datExpr = as.data.frame(t(df %>%column_to_rownames('gene_id')))
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed hybrid")
  corOptions = list(use = 'p', maxPOutliers = 0.1)
  adjacency = adjacency(datExpr, power = sft$powerEstimate,
                        type = "signed hybrid", corFnc = "bicor")
  TOM = TOMsimilarity(adjacency)
  rownames(TOM) = df$gene_id
  colnames(TOM) = df$gene_id
  return(TOM)
}
method_combinations = expand_grid(downsampling = c('get_nondownsampled_data', 'get_downsampled_data'),
            batch_correction = c('combat_seq_correct', 'no_correction'),
            library_size_normalization = c('upper_quartile_normalize', 'tmm_normalize', 'tpm_normalize', 'vst_normalize'), 
            scaling =  c('no_scaling', 'log_plus_one', 'inverse_hyperbolic_sine', 'zscore_scale'),
            correlation = 'pearson') %>%
  bind_rows(
    expand_grid(downsampling = c('get_nondownsampled_data', 'get_downsampled_data'),
                batch_correction = c('combat_seq_correct', 'no_correction'),
                library_size_normalization =  'vst_normalize', 
                scaling =  'no_scaling',
                correlation = 'wgcna')
  ) %>%
  mutate(workflow = str_c(downsampling, batch_correction, library_size_normalization, scaling, correlation, sep='-')) %>%
  rowwise() %>%
  mutate(res = list(get(downsampling)() %>% 
                      get(batch_correction)() %>% 
                      get(library_size_normalization)() %>% 
                      get(scaling)())) %>%
  rowwise() %>%
  mutate(genes_filtered = list(res %>% semi_join(genes_with_mean_tpm_above_1)),
         cors = list(genes_filtered %>% get(correlation)())) %>% 
  rowwise() %>% 
  mutate(top_100 = list(get_top_100_edges_for_all_genes(cors)))
method_combinations = method_combinations %>%
  mutate(batch_correction = factor(batch_correction, levels = c('no_correction', 'combat_seq_correct')),
         scaling =  factor(scaling, levels = c('no_scaling', 'log_plus_one', 'inverse_hyperbolic_sine', 'zscore_scale')),
         downsampling = factor(downsampling, levels = c('get_nondownsampled_data', 'get_downsampled_data')))
functional_annotation_jcvi_network = functional_annotation_jcvi %>% semi_join(genes_with_mean_tpm_above_1)
## 1844 genes removed
functional_annotation_jcvi_network %>% write_rds('Z:/flavus_sra_rnaseq/output/coexpression/network_evaluation/functional_annotation_jcvi_network')
method_combinations %>%
  pwalk(function(workflow, top_100, ...){
    write_rds(top_100, here(paste0('output/coexpression/network_evaluation/', workflow, '_top_100.rds')))
  })
method_combinations %>%
  pwalk(function(workflow, res, ...){
    write_rds(res, here(paste0('output/coexpression/network_evaluation/', workflow, '.rds')))
  })
top_methods  = c('get_nondownsampled_data-no_correction-upper_quartile_normalize-no_scaling-pearson')
method_combinations %>%
  filter(workflow == 'get_nondownsampled_data-no_correction-upper_quartile_normalize-no_scaling-pearson') %>%
  pull(cors) %>%
  pluck(1) %>% 
  {ecdf(abs(.))(0.5)}
## [1] 0.9805056 
## An R of 0.5 is the 0.9805056 percentile in the Pearson correlation based network.
filter_cor_and_create_igraph = function(workflow, cor_matrix){
  threshold = ifelse(str_detect(workflow, 'wgcna'), 
                     quantile(abs(cor_matrix), 0.9805056), 
                     0.5)
  ## Used the percentile of a 0.5 pearson coefficient cutoff in the pearson network,
  ## which is the 0.9805056 percentile, as a heuristic to determine equivalent cut off in WGCNA's TOM
  cor_matrix[ abs(cor_matrix) <= threshold ] = 0
  return(graph.adjacency(cor_matrix, diag = FALSE, weighted=TRUE, mode="lower"))
}
method_combinations %>% 
  filter(workflow %in% top_methods) %>%
  pwalk(function(workflow, cors, ...){
    filter_cor_and_create_igraph(workflow, cors) %>%
      write_rds(here(paste0('shiny_app/data/', workflow, '.rds')))
  })
method_combinations_top_chrom_level = expand_grid(downsampling = c('get_nondownsampled_data_chrom_level', 'get_downsampled_data_chrom_level'),
                                  batch_correction = c('combat_seq_correct', 'no_correction'),
                                  library_size_normalization = c('upper_quartile_normalize', 'tmm_normalize', 'tpm_normalize', 'vst_normalize'), 
                                  scaling =  c('no_scaling', 'log_plus_one', 'inverse_hyperbolic_sine', 'zscore_scale'),
                                  correlation = 'pearson') %>%
  bind_rows(
    expand_grid(downsampling = c('get_nondownsampled_data_chrom_level', 'get_downsampled_data_chrom_level'),
                batch_correction = c('combat_seq_correct', 'no_correction'),
                library_size_normalization =  'vst_normalize', 
                scaling =  'no_scaling',
                correlation = 'wgcna')
  ) %>%
  mutate(workflow = str_c(downsampling, batch_correction, library_size_normalization, scaling, correlation, sep='-')) %>%
  filter(workflow %in% (top_methods %>% str_replace('data', 'data_chrom_level'))) %>%
  rowwise() %>%
  mutate(res = list(get(downsampling)() %>% 
                      get(batch_correction)() %>% 
                      get(library_size_normalization)() %>% 
                      get(scaling)())) %>% 
  rowwise() %>%
  mutate(genes_filtered = list(res %>% semi_join(genes_with_mean_tpm_above_1_chrom_level)),
         cors = list(genes_filtered %>% get(correlation)()))
method_combinations_top_chrom_level %>%
  pwalk(function(workflow, cors, ...){
    filter_cor_and_create_igraph(workflow, cors) %>%
      write_rds(here(paste0('shiny_app/data/', workflow, '.rds')))
  })
read_rds(here(paste0('shiny_app/data/', top_methods, '.rds'))) %>%
  upgrade_graph() %>%
  write_rds(here(paste0('shiny_app/data/', top_methods, '.rds')))
