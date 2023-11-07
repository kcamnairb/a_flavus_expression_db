library(here)
library(tidyverse)

fastp_reports = read_tsv(here('output/fastp_reports/multiqc_data/multiqc_general_stats.txt')) %>%
  rename_with(~str_remove(.x, 'fastp_mqc-generalstats-fastp-')) %>%
  dplyr::rename(run = Sample) %>%
  mutate(run = str_remove(run, '_1'))
counts_jcvi = read_csv(here('output', 'salmon_quant', 'jcvi', 'A_flavus_jcvi_counts_unfiltered.csv'))
sra_df = read_csv(here('output/sra_metadata.csv'))

counts_sums_jcvi = counts_jcvi %>% 
  pivot_longer(-gene_id, names_to='run', values_to='counts') %>%
  group_by(run) %>%
  summarize(total_counts = sum(counts))
counts_sums_jcvi %>% dplyr::count(total_counts > 10e6)
qc = fastp_reports %>%
  left_join(counts_sums_jcvi, by='run') %>% 
  left_join(sra_df, by='run')
qc %>%
  mutate(low_count = total_counts < 1e6) %>%
  dplyr::count(bioproject, low_count, sort=T) %>%
  pivot_wider(names_from = low_count, values_from = n, values_fill = 0) %>% 
  select(-`NA`) %>% 
  arrange(desc(`TRUE`))
  left_join(sra_df %>% distinct(bioproject, abstract)) %>% View()
## Almost all of the bioprojects that have a lot of samples with low counts consist of mixed organisms, 
## so the reads are probably coming mostly from the other organism.
counts_sums_jcvi %>%
  anti_join(fastp_reports) 

ggplotly(qc %>%
           filter(total_counts >= 1e6, run != 'SRR8526599') %>%
           pivot_longer(cols=c(pct_duplication, after_filtering_gc_content, pct_surviving, total_counts),
                        names_to = 'metric', values_to = 'value') %>%
  ggplot(aes(x=run, y=value, color=bioproject, label=title2)) +
  geom_col() +
  facet_wrap(~metric, ncol = 1, scales = 'free') +
  ggeasy::easy_remove_legend()
)
