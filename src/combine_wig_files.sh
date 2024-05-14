#!/bin/bash
#SBATCH --job-name="wig"
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -p mem
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"
module load miniconda
source activate wiggletools_env
wiggletools mean output/star/*.unstranded.out.bw > output/star/jcvi_mean_unstranded.wig
wigToBigWig output/star/jcvi_mean_unstranded.wig ../A_flavus_NRRL3357_chrom.sizes output/star/jcvi_mean_unstranded.bw
wiggletools mean output/star_chrom_level/*_unstrandedSignal.Unique.str1.out.bw > \
  output/star_chrom_level/chrom_level_mean_unstranded.wig
wigToBigWig output/star_chrom_level/chrom_level_mean_unstranded.wig \
  ../Af3357_chrom_level_chrom.sizes output/star_chrom_level/chrom_level_mean_unstranded.bw

