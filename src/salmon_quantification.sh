#!/bin/bash
#SBATCH --job-name="salmon"
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"

set -x
set -e
set -u
module load salmon/1.10.0

base=$(basename $1 | cut -d _ -f1 | cut -d . -f1)
if (( $# == 1 )); then
  salmon quant -i data/salmon_db/salmon_index_jcvi -l A -o output/salmon_quant/jcvi/$base --gcBias -p 12 -r $1 
  salmon quant -i data/salmon_db/salmon_chrom_level -l A -o output/salmon_quant/chrom_level/$base --gcBias -p 12 -r $1
elif (( $# == 2 )); then
  salmon quant -i data/salmon_db/salmon_index_jcvi -l A -o output/salmon_quant/jcvi/$base --gcBias -p 12 -1 $1 -2 $2
  salmon quant -i data/salmon_db/salmon_chrom_level -l A -o output/salmon_quant/chrom_level/$base --gcBias -p 12 -1 $1 -2 $2
else 
  >&2 echo "Wrong number of parameters"
fi
