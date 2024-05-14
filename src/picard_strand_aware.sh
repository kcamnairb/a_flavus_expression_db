#!/bin/bash
#SBATCH --job-name="picard"
#SBATCH -n 2
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"

module load picard
run=$1
strand=$2
java -jar /software/el9/apps/picard/3.0.0/picard.jar CollectRnaSeqMetrics \
  I=output/star/${run}_Aligned.sortedByCoord.out.bam \
  STRAND=${strand} \
  O=output/picard/$run \
  REF_FLAT=../FungiDB-3.2_Aflavus_NRRL3357_cufflinks_NW2EQ_replaced.refflat
