#!/bin/bash
#SBATCH --job-name="fastq_trim"
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"

set -x
set -e
set -u
 
module load miniconda
source activate fastp_env

base=$(basename $1 | cut -d _ -f1 | cut -d . -f1)
if ! ls output/tmp_data/fastq_trimmed/$(basename $1)* 1> /dev/null 2>&1 ; then
  if (( $# == 1 )); then
    fastp --thread 12 -l 40 -i $1 -o output/tmp_data/fastq_trimmed/$(basename $1) \
      --report_title $base --html ${base}.html --json ${base}.json
  elif (( $# == 2 )); then
    fastp --thread 12 -l 40 -i $1 -I $2 -o output/tmp_data/fastq_trimmed/$(basename $1) \
      -O output/tmp_data/fastq_trimmed/$(basename $2) --report_title $base \
      --html ${base}.html --json ${base}.json
  else 
    >&2 echo "Wrong number of parameters"
  fi
fi
