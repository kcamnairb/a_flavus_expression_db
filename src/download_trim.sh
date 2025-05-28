#!/bin/bash
#SBATCH --job-name="sra_download"
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -o "logs/download_trim_%A_%a.out"
#SBATCH -e "logs/download_trim_%A_%a.err"
#SBATCH --time=1-00:00:00
#SBATCH --partition="ceres"

set -x
set -e
set -u
SRA_RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SRA_FILE})
module load miniconda
source activate kingfisher_env
if ! ls data/fastq/${SRA_RUN}*fastq.gz 1> /dev/null 2>&1 ; then 
  kingfisher get -r $SRA_RUN -m ena-ascp aws-http ena-ftp \
    --download-threads 6 --output_directory data/fastq/ -f fastq.gz
fi
source activate fastp_env
if ! ls output/fastq_trimmed/${SRA_RUN}* 1> /dev/null 2>&1 ; then
  if [[ -f data/fastq/${SRA_RUN}.fastq.gz ]]; then
    fastp --thread 12 -l 40 -i data/fastq/${SRA_RUN}.fastq.gz -o output/fastq_trimmed/${SRA_RUN}.fastq.gz \
      --report_title $SRA_RUN --html output/fastp_reports/${SRA_RUN}.html --json output/fastp_reports/${SRA_RUN}.json
  elif [[ -f data/fastq/${SRA_RUN}_1.fastq.gz && -f data/fastq/${SRA_RUN}_2.fastq.gz ]]; then
    fastp --thread 12 -l 40 -i data/fastq/${SRA_RUN}_1.fastq.gz -I data/fastq/${SRA_RUN}_2.fastq.gz -o output/fastq_trimmed/${SRA_RUN}_1.fastq.gz \
      -O output/fastq_trimmed/${SRA_RUN}_2.fastq.gz --report_title $SRA_RUN \
      --html output/fastp_reports/${SRA_RUN}.html --json output/fastp_reports/${SRA_RUN}.json
  else 
    >&2 echo "Wrong number of parameters"
  fi
fi

