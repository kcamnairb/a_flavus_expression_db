library(here)
library(tidyverse)

fastp_reports = read_tsv(here('output/fastp_reports/filtered/multiqc_data/multiqc_general_stats.txt')) %>%
  rename_with(~str_remove(.x, 'fastp_mqc-generalstats-fastp-')) %>%
  dplyr::rename(run = Sample) %>%
  mutate(run = str_remove(run, '_1'))
counts_jcvi = read_csv(here('output/star/jcvi_counts.csv'))
sra_df = read_csv(here('output/sra_metadata.csv'))

counts_sums_jcvi = counts_jcvi %>% 
  pivot_longer(-gene_id, names_to='run', values_to='counts') %>%
  group_by(run) %>%
  summarize(total_counts = sum(counts))
counts_sums_jcvi %>% dplyr::count(total_counts > 10e6)
fastp = fastp_reports %>%
  left_join(counts_sums_jcvi, by='run') %>% 
  left_join(sra_df, by='run')
mpn65 = c('#ff0029','#377eb8','#66a61e','#984ea3','#00d2d5','#ff7f00','#af8d00','#7f80cd','#b3e900','#c42e60','#a65628',
          '#f781bf','#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#fccde5','#bc80bd','#ffed6f','#c4eaff','#cf8c00',
          '#1b9e77','#d95f02','#e7298a','#e6ab02','#a6761d','#0097ff','#00d067','#000000','#252525','#525252','#737373',
          '#969696','#bdbdbd','#f43600','#4ba93b','#5779bb','#927acc','#97ee3f','#bf3947','#9f5b00','#f48758','#8caed6',
          '#f2b94f','#eff26e','#e43872','#d9b100','#9d7a00','#698cff','#d9d9d9','#00d27e','#d06800','#009f82','#c49200',
          '#cbe8ff','#fecddf','#c27eb6','#8cd2ce','#c4b8d9','#f883b0','#a49100','#f48800','#27d0df','#a04a9b')

plotly::ggplotly(fastp %>%
           filter(total_counts >= 1e6, run != 'SRR8526599') %>%
           pivot_longer(cols=c(pct_duplication, after_filtering_gc_content, pct_surviving, total_counts),
                        names_to = 'metric', values_to = 'value') %>%
  ggplot(aes(x=run, y=value, color=bioproject, label=title2)) +
  geom_col() +
  facet_wrap(~metric, ncol = 1, scales = 'free') +
  ggeasy::easy_remove_legend() +
  ggeasy::easy_remove_x_axis(what = c('ticks', 'text')) +
  ggplot2::scale_color_manual(values = mpn65)
)
picard = read_tsv(here('output/picard/multiqc_data/multiqc_picard_RnaSeqMetrics.txt')) %>% 
  janitor::clean_names() %>%
  rename(run = sample)
picard %>% select(run, pct_mrna_bases:median_5prime_to_3prime_bias) %>% 
  left_join(counts_sums_jcvi, by='run') %>%
  filter(run != 'SRR8526599') %>%
  pivot_longer(cols = pct_mrna_bases:total_counts, 
               names_to = 'metric', values_to='value') %>%
  mutate(run = fct_reorder(run, frac_strand2)) %>%
  ggplot(aes(run, value, fill = strand_specificity)) +
  geom_col() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  facet_wrap(~metric, scales = 'free_y')
ggsave(here('output/picard/picard_qc.png'), unit='in', width = 12, height=8, dpi = 300)
picard %>% count(pct_mrna_bases > 30)
## SRR10160904 has a high median 3' bias of 3.783471. 
## It is a clear outlier on the gene coverage visualization, going to remove.
gene_coverage = read_tsv('output/picard/multiqc_data/mqc_picard_rna_coverage_1.txt')
gene_coverage %>%
  filter(Sample %in% c('SRR10160904', 'SRR283858')) %>%
  pivot_longer(-Sample, names_to = 'percent through gene', values_to = 'coverage') %>%
  ggplot(aes(`percent through gene`, coverage, color=Sample, group=Sample)) +
  geom_smooth(se=F, alpha = 0.2) +
  theme_bw() +
  ggeasy::easy_remove_legend() 
## https://bioinformatics.ccr.cancer.gov/btep/wp-content/uploads/RNA-seq_BETP_2019rev.pdf#page=21.00  
picard %>% left_join(metadata, by='run') %>% 
  count(low_coding = pct_coding_bases < 40, bioproject, sort=T) %>% 
  filter(low_coding == TRUE, 
         #bioproject %in% c('PRJNA573552'),
         )
## Filtering out samples with coding bases less than 40% results in 15 samples being removed. 
## Only one of these samples (SRR10160904) is in the co-expression network.
picard %>% 
  left_join(metadata, by='run') %>%
  left_join(counts_sums_jcvi, by='run') %>%
  filter(pct_coding_bases < 40) %>%
  pull(run) %>%
  write_lines(here('output/picard/low_coding_samples.txt'))

