#!/bin/bash
#SBATCH --job-name="sra_download"
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"

set -x
set -e
set -u
 
module load miniconda
#source activate fastq-dl_env
source activate kingfisher
if ! ls data/tmp_data/fastq/${1}* 1> /dev/null 2>&1 ; then 
  #fastq-dl --outdir data/tmp_data/fastq --cpus 6 --accession $1 --provider sra
  kingfisher get -r $1 -m ena-ftp  --download-threads 6 
  mv ${1}.fastq.gz data/tmp_data/fastq
fi
