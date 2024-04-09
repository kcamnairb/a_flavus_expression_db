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
#sra_runs_file=output/sra_runs.txt
#sra_runs_file=output/sra_runs_additional.txt
#sra_runs_file=output/runs_w_errors.txt
sra_runs_file=output/sra_runs_wig_and_timecourse.txt
parallel -a $sra_runs_file sbatch src/download_trim_count.sh {}
#for f in $(grep -l error stderr*) ; do  grep ' --accession' $f | perl -p -e 's/.* (.*)/\1/' >> output/runs_w_errors.txt ; done
## remove empty files because "--force" doesn't seem to work with fastq-dl
#find  output/tmp_data/fastq -type f -empty -delete 
#job_ids=()
#for run in $(cat $sra_runs_file) ; do
#  job_id=$(sbatch --parsable src/trim_reads.sh data/tmp_data/fastq/${run}*)
#  job_ids+=($job_id)
#  sbatch src/trim_reads.sh data/tmp_data/fastq/${run}*
#done
#dependency_list=$(IFS=,; echo "${job_ids[*]}")

#for run in $(cat output/runs_with_trim_errors.txt) ; do
#  sbatch src/repair_fastq.sh output/tmp_data/fastq_trimmed/${run}*
#done
#mv output/tmp_data/fastq_repaired/* output/tmp_data/fastq/
#for run in $(cat output/runs_with_trim_errors.txt) ; do
#  sbatch src/trim_reads.sh output/tmp_data/fastq/${run}*
#done
#mkdir output/fastp_reports
#mv *{html,json} output/fastp_reports/
#module load miniconda
#source activate multiqc_env
#multiqc -o output/fastp_reports output/fastp_reports

#for run in $(cat $sra_runs_file) ; do
#  if ! ls output/salmon_quant/jcvi/${run} 1> /dev/null 2>&1 ; then
##    sbatch --dependency=afterok:$dependency_list src/salmon_quantification.sh output/tmp_data/fastq_trimmed/${run}*
#    sbatch src/salmon_quantification.sh output/tmp_data/fastq_trimmed/${run}*
#  fi
#done
