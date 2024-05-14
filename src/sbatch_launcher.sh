#!/bin/bash
#SBATCH --job-name="multiqc"
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"
set -x
set -e
set -u

module load parallel
sra_runs_file=output/sra_runs.txt
parallel -a $sra_runs_file sbatch src/download_trim_count.sh {}
parallel -a $sra_runs_file sbatch src/align_count_chrom_level.sh {}
#for f in $(grep -l error stderr*) ; do  grep ' --accession' $f | perl -p -e 's/.* (.*)/\1/' >> output/runs_w_errors.txt ; done



