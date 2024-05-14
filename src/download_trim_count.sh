#!/bin/bash
#SBATCH --job-name="sra_download"
#SBATCH -n 72
#SBATCH -N 1
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"

set -x
set -e
set -u
 
module load miniconda
source activate kingfisher_env
if ! ls data/fastq/${1}* 1> /dev/null 2>&1 ; then 
  kingfisher get -r $1 -m ena-ascp aws-http ena-ftp \
    --download-threads 6 --output_directory data/fastq/ -f fastq.gz
fi
source activate fastp_env
base=$1
if ! ls output/fastq_trimmed/${base}* 1> /dev/null 2>&1 ; then
  if [[ -f data/fastq/${base}.fastq.gz ]]; then
    fastp --thread 12 -l 40 -i data/fastq/${base}.fastq.gz -o output/fastq_trimmed/${base}.fastq.gz \
      --report_title $base --html ${base}.html --json ${base}.json
  elif [[ -f data/fastq/${base}_1.fastq.gz && -f data/fastq/${base}_2.fastq.gz ]]; then
    fastp --thread 12 -l 40 -i data/fastq/${base}_1.fastq.gz -I data/fastq/${base}_2.fastq.gz -o output/fastq_trimmed/${base}_1.fastq.gz \
      -O output/fastq_trimmed/${base}_2.fastq.gz --report_title $base \
      --html ${base}.html --json ${base}.json
  else 
    >&2 echo "Wrong number of parameters"
  fi
fi
module load star/2.7.10b
if [[ -f output/fastq_trimmed/${base}.fastq.gz ]]; then
  STAR --readFilesIn output/fastq_trimmed/${base}.fastq.gz \
    --genomeDir ../star_db/Af3357/ \
    --outFileNamePrefix output/star/${base}_ \
    --runThreadN 72 --alignIntronMax 1000  --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic --quantMode GeneCounts --readFilesCommand zcat \
    --outWigType wiggle --outBAMsortingBinsN 200 --limitBAMsortRAM 17179861000
elif [[ -f output/fastq_trimmed/${base}_1.fastq.gz && -f output/fastq_trimmed/${base}_2.fastq.gz ]]; then
  STAR --readFilesIn output/fastq_trimmed/${base}_1.fastq.gz \
    output/fastq_trimmed/${base}_2.fastq.gz \
    --genomeDir ../star_db/Af3357/ \
    --outFileNamePrefix output/star/${base}_ \
    --runThreadN 72 --alignIntronMax 1000  --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic --quantMode GeneCounts --readFilesCommand zcat \
    --outWigType wiggle --outBAMsortingBinsN 200 --limitBAMsortRAM 17179861000
else 
  >&2 echo "Wrong number of parameters"
fi
 STAR --runMode inputAlignmentsFromBAM \
   --inputBAMfile output/star/${base}_Aligned.sortedByCoord.out.bam \
   --outWigType wiggle \
   --outWigStrand Unstranded \
   --outFileNamePrefix output/star/${base}_
## rename  str1 unstranded output/star/*_Signal.Unique.str1.out.bw
