#!/bin/bash
#SBATCH --job-name="picard"
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -o "logs/picard_%A_%a.out"
#SBATCH -e "logs/picard_%A_%a.err"
#SBATCH --partition="ceres"
#SBATCH --time=1-00:00:00

module load picard
run=$1
strand=$2
if ! ls output/picard/$run 1> /dev/null 2>&1 ; then
  java -jar /software/el9/apps/picard/3.0.0/picard.jar CollectRnaSeqMetrics \
    I=output/star/${run}_Aligned.sortedByCoord.out.bam \
    STRAND=${strand} \
    O=output/picard/$run \
    REF_FLAT=../FungiDB-3.2_Aflavus_NRRL3357_cufflinks_NW2EQ_replaced.refflat
fi