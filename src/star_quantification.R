library(tidyverse)
library(here)
library(fst)
counts_star_jcvi = list.files(here('output/star/'), pattern = 'ReadsPerGene.out.tab', full.names = TRUE) %>%
  set_names() %>%
  imap(~read_tsv(.x, skip = 4, col_names = c('gene_id', 'count_unstranded', 'count_strand1', 'count_strand2')) %>%
         mutate(file = str_replace(.y, '.*/(.*)_Reads.*', '\\1')) %>%
         rename(run = file))
sra_df = read_csv(here('output/sra_metadata.csv'), col_types = cols(.default = col_character()))
count_strand_totals = counts_star_jcvi %>% map(~.x %>%
  summarize(count_unstranded= sum(count_unstranded),
            count_strand1 = sum(count_strand1), 
            count_strand2 = sum(count_strand2),
            run = first(run)
  )) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(frac_strand1 = count_strand1 / sum(count_strand1, count_strand2),
         frac_strand2 = count_strand2 / sum(count_strand1, count_strand2)) %>% 
  mutate(strand_specificity = case_when(frac_strand1 > 0.7 ~ 'strand1',
                                        frac_strand1 < 0.3 ~ 'strand2',
                                        frac_strand1 >= 0.3 & frac_strand1 <= 0.7 ~ 'unstranded')) %>%
  left_join(sra_df, by='run')
count_strand_totals %>% write_csv(here('output/star/strand_totals.csv'))
count_strand_totals = read_csv(here('output/star/strand_totals.csv'))
p = count_strand_totals %>%
  ggplot(aes(frac_strand1, frac_strand2, key=bioproject)) +
  geom_point(alpha=0.1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
plotly::ggplotly(p)

select_star_count_col = function(df){
  current_run = df[1,]$run
  current_strand_specificity = count_strand_totals %>% 
    filter(run == current_run) %>%
    pull(strand_specificity)
  df %>% select(gene_id, !!quo_name(current_run) := all_of(paste0('count_', current_strand_specificity)))
}
counts_star_jcvi = counts_star_jcvi %>% 
  map(select_star_count_col) %>% 
  reduce(full_join)
counts_star_jcvi = counts_star_jcvi %>% 
  filter(!str_detect(gene_id, 'tRNA'))
counts_sums_jcvi = counts_star_jcvi %>% 
  pivot_longer(-gene_id, names_to='run', values_to='counts') %>%
  group_by(run) %>%
  summarize(total_counts = sum(counts)) 
low_coding_samples = read_lines(here('output/picard/low_coding_samples.txt'))
q30_below_75 = read_lines(here('output/fastp_reports/q30_below_75_samples.txt'))
samples_to_keep = counts_sums_jcvi %>% 
  filter(total_counts >= 1e6, 
         !run %in% low_coding_samples,
         !run %in% q30_below_75) %>%
  pull(run)
length(samples_to_keep)
counts_star_jcvi  = counts_star_jcvi %>% select(all_of(c('gene_id', samples_to_keep)))
counts_star_jcvi %>% write_csv(here('output/star/jcvi_counts.csv'))
counts_star_jcvi = read_csv(here('output/star/jcvi_counts.csv'))
metadata = sra_df %>% 
  filter(run %in% samples_to_keep)
metadata %>% 
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
  mutate(submitter = str_c(first_name, last_name, collapse = ' ')) %>%
  select(-first_name, -last_name, -BioSampleModel, -pubmed) %>%
  mutate(doi = ifelse(is.na(doi), NA, str_glue('https://doi.org/{doi}'))) %>%
  mutate(doi = ifelse(is.na(doi), NA, str_glue("<a href='{doi}'>{doi}</a>"))) %>%
  select(-filename) %>%
  mutate(strain = if_else(bioproject == 'PRJNA1085269', 'NRRL 3357', strain)) %>%
  write_csv(here('shiny_app/data/sra_metadata_filtered.csv'))
counts_star_chrom_level = list.files(here('output/star_chrom_level/'), pattern = 'ReadsPerGene.out.tab', full.names = TRUE) %>%
  set_names() %>%
  imap(~read_tsv(.x, skip = 4, col_names = c('gene_id', 'count_unstranded', 'count_strand1', 'count_strand2')) %>%
         mutate(file = str_replace(.y, '.*/(.*)_Reads.*', '\\1')) %>%
         rename(run = file))
count_strand_totals_chrom_level = counts_star_chrom_level %>% map(~.x %>%
                                                 summarize(count_unstranded= sum(count_unstranded),
                                                           count_strand1 = sum(count_strand1), 
                                                           count_strand2 = sum(count_strand2),
                                                           run = first(run)
                                                 )) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(frac_strand1 = count_strand1 / sum(count_strand1, count_strand2),
         frac_strand2 = count_strand2 / sum(count_strand1, count_strand2)) %>% 
  mutate(strand_specificity = case_when(frac_strand1 > 0.7 ~ 'strand1',
                                        frac_strand1 < 0.3 ~ 'strand2',
                                        frac_strand1 >= 0.3 & frac_strand1 <= 0.7 ~ 'unstranded')) %>%
  left_join(sra_df, by='run')
counts_star_chrom_level  = counts_star_chrom_level %>% select(all_of(c('gene_id', samples_to_keep)))
count_strand_totals_chrom_level %>% write_csv(here('output/star_chrom_level/strand_totals.csv'))
count_strand_totals %>%
  left_join(count_strand_totals_chrom_level, by='run', suffix = c('_jcvi', '_chrom')) %>% 
  summarize(sum(strand_specificity_jcvi != strand_specificity_chrom))
counts_star_chrom_level = counts_star_chrom_level %>% 
  map(select_star_count_col) %>% 
  reduce(full_join)
counts_star_chrom_level %>% 
  mutate(gene_id = str_remove(gene_id, 'gene-')) %>%
  write_csv(here('output/star_chrom_level/chrom_level_counts.csv'))
counts_star_chrom_level = read_csv(here('output/star_chrom_level/chrom_level_counts.csv'))
counts_sums_chrom_level = counts_star_chrom_level %>% 
  pivot_longer(-gene_id, names_to='run', values_to='counts') %>%
  group_by(run) %>%
  summarize(total_counts = sum(counts)) 
counts_sums_jcvi %>%
  left_join(counts_sums_chrom_level, by='run', suffix = c('_jcvi', '_chrom'))
### PCA ###
metadata_pca = metadata %>% 
  mutate(run = factor(run, levels = colnames(counts_star_jcvi)[2:ncol(counts_star_jcvi)])) %>%
  arrange(run) %>%
  left_join(qc, by='run') %>%
  column_to_rownames('run') %>%
  mutate(run = rownames(.))
  
all(rownames(metadata_pca) == colnames(counts_star_jcvi)[2:ncol(counts_star_jcvi)])
dds_jcvi = DESeq2::DESeqDataSetFromMatrix(counts_star_jcvi %>% column_to_rownames('gene_id'), 
                                          colData = metadata_pca, design = ~ 1)
dds_jcvi = DESeq2::DESeq(dds_jcvi)
vst_jcvi = DESeq2::vst(dds_jcvi, blind=TRUE)
SummarizedExperiment::assay(vst_jcvi) %>% as_tibble(rownames = 'gene_id') %>%
 write_csv(here('shiny_app/data', 'vst_jcvi.csv'))

pca_data_jcvi = DESeq2::plotPCA(vst_jcvi, intgroup='bioproject', returnData=T, ntop=4000)
percentVar_jcvi = round(100 * attr(pca_data_jcvi, "percentVar"))
ggplotly(
  pca_data_jcvi %>% as.data.frame() %>% rownames_to_column('run') %>%
    left_join(metadata_pca %>%
                select(-bioproject) %>%
                mutate(across(c(strain, sample_type,  source_name, genotype), 
                              ~str_remove(.x, regex('^A.* flavus ', ignore_case = TRUE)))) %>%
                #mutate(Date = as.character(Date)) %>%
                mutate(across(c(strain, sample_type,  source_name, genotype, treatment, isolate), 
                              ~na_if(.x, 'missing'))) %>%
                unite('sample_description', strain, sample_type,  source_name, genotype, treatment, isolate,
                      sep = ';', na.rm =TRUE, remove = FALSE), by='run') %>% 
    ggplot(aes(PC1, PC2, color=pct_mrna_bases,  label=sample_description)) +
    geom_point(alpha=0.5) +
    xlab(paste0("PC1: ",percentVar_jcvi[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar_jcvi[2],"% variance")) +
    coord_fixed(sqrt(percentVar_jcvi[2] / percentVar_jcvi[1])) +
    theme_bw() +
    #scale_fill_manual(values = mpn65) +
    #scale_color_gradient(low='blue', high='red') +
    scale_color_gradientn(colours = rainbow(5)) #+
    #ggeasy::easy_remove_legend()
)

all(rownames(metadata_pca) == colnames(counts_star_chrom_level)[2:ncol(counts_star_chrom_level)])
dds_chrom_level = DESeq2::DESeqDataSetFromMatrix(counts_star_chrom_level %>% column_to_rownames('gene_id'), 
                                          colData = metadata_pca, design = ~ 1)
dds_chrom_level = DESeq2::DESeq(dds_chrom_level)
vst_chrom_level = DESeq2::vst(dds_chrom_level, blind=TRUE)
SummarizedExperiment::assay(vst_chrom_level) %>% as_tibble(rownames = 'gene_id') %>%
  mutate(gene_id = str_remove(gene_id, 'gene-')) %>%
  write_csv(here('shiny_app/data', 'vst_chrom_level.csv'))

transcript_lengths_jcvi = read_csv('Z:/genomes/af3357_annotation/af3357_cds_lengths.csv', 
                                   col_names = c('gene_id', 'length'))
tpm_normalize = function(df, transcript_lengths_df) {
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

tpm_jcvi = counts_star_jcvi %>% tpm_normalize(transcript_lengths_jcvi)
tpm_jcvi %>% write_csv(here('shiny_app/data/A_flavus_jcvi_tpm.csv'))
transcript_lengths_chrom_level = read_csv(
  'Z:/genomes/A_flavus_3357_chromosome_level_JGI/download_from_ncbi/transcript_lengths.csv')
tpm_chrom_level = counts_star_chrom_level %>% tpm_normalize(transcript_lengths_chrom_level)
tpm_chrom_level %>% write_csv(here('shiny_app/data/A_flavus_chrom_level_tpm.csv'))
convert_file_to_fst = function(file) read_csv(file) %>% write_fst(str_replace(file, '\\.csv', '.fst'))
c('data/A_flavus_jcvi_tpm.csv', 'data/A_flavus_chrom_level_tpm.csv',
  'data/vst_jcvi.csv', 'data/vst_chrom_level.csv') %>% walk(convert_file_to_fst)
'data/sra_metadata_filtered.csv' %>% walk(convert_file_to_fst)
'data/sra_metadata_filtered_hand_edited.csv' %>% walk(convert_file_to_fst)
'data/functional_annotation_jcvi.csv' %>% walk(convert_file_to_fst)
c('data/functional_annotation_chrom_level.csv') %>% walk(convert_file_to_fst)
