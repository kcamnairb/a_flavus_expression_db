#!/bin/bash
## QC filtering cutoffs from https://github.com/yongchao/motrpac_rnaseq/blob/master/MOP_details.md
## Low RIN scores, e.g. RIN < 6
## Low number of sequenced reads after trimming adaptors, e.g. samples with < 20M reads
## Abnormal GC content, e.g. GC% > 80% or < 20%
## High percentage of rRNA, e.g. rRNA > 20%
## Low number of mapped reads, e.g. < 60%
## Low number of % exonic read, % exonic read < 50%
## Mapped read count < 50% of average mapped read count per sample
.libPaths('/project/aflavus_field_isolates/R_packages/4.3')
library(tidyverse)
count_strand_totals = read_csv('output/star/strand_totals.csv')
count_strand_totals %>%
  mutate(strand_specificity = str_replace_all(strand_specificity,
    c('unstranded' = 'NONE', 
    'strand1' = 'FIRST_READ_TRANSCRIPTION_STRAND', 
    'strand2'= 'SECOND_READ_TRANSCRIPTION_STRAND'))) %>%
  pmap(function(run, strand_specificity, ...){
    system(command = paste('sbatch src/picard_strand_aware.sh', run, strand_specificity),
           wait = FALSE)
  })
