#!/bin/bash
#SBATCH --job-name="fastq_repair"
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"

set -x
set -e
set -u
 
module load bbtools/38.79

base=$(basename $1 | cut -d _ -f1 | cut -d . -f1)
if (( $# == 1 )); then
  repair.sh -Xmx35g in=output/tmp_data/fastq/$(basename $1) \
    out=output/tmp_data/fastq_repaired/$(basename $1) \
    repair=t tossbrokenreads=t
elif (( $# == 2 )); then
  repair.sh -Xmx35g in=output/tmp_data/fastq/$(basename $1) in2=output/tmp_data/fastq/$(basename $2) \
    out=output/tmp_data/fastq_repaired/$(basename $1) out2=output/tmp_data/fastq_repaired/$(basename $2) \
    repair=t tossbrokenreads=t
else 
    >&2 echo "Wrong number of parameters"
fi
