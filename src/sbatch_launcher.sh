#!/bin/bash
set -x
set -e
set -u

sra_runs_file=output/sra_runs.txt
# sra_runs_file=output/sra_runs_additional.txt
# sra_runs_file=output/runs_w_errors.txt
first_job_id=$(sbatch --array=1-$(wc -l < $sra_runs_file)%8 \
  --export=ALL,SRA_FILE=$sra_runs_file \
  src/download_trim.sh | awk '{print $4}')
second_job_id=$(sbatch --depend=afterok:$first_job_id \
  --array=1-$(wc -l < $sra_runs_file) \
  --export=ALL,SRA_FILE=$sra_runs_file \
  src/align_count.sh | awk '{print $4}')
# second_job_id=$(sbatch \
#   --array=1-$(wc -l < $sra_runs_file) \
#   --export=ALL,SRA_FILE=$sra_runs_file \
#   src/align_count.sh | awk '{print $4}')
#sbatch --time=2-00:00:00  --wrap="module load miniconda ; source activate multiqc_env ; multiqc -o output/fastp_reports output/fastp_reports"
#for f in $(grep -l error stderr*) ; do  grep ' --accession' $f | perl -p -e 's/.* (.*)/\1/' >> output/runs_w_errors.txt ; done

# find logs/ -type f -mtime -4 -exec grep -li "error" {} + | xargs -I {} sed -n '4p' {} | grep "SRA_RUN=" | cut -d "=" -f2 | sort -u

