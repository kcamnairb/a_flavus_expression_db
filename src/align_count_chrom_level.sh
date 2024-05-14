#!/bin/bash
#SBATCH --job-name="align_count"
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"

set -x
set -e
set -u

module load star/2.7.10b
base=$1
if [[ -f output/star_chrom_level/${base}_Log.final.out ]] ; then
  echo $base already finished
  exit
fi
if [[ -f output/fastq_trimmed/${base}.fastq.gz ]]; then
  STAR --readFilesIn output/fastq_trimmed/${base}.fastq.gz \
  --genomeDir ../star_db/Af3357_chrom_level/ \
  --outFileNamePrefix output/star_chrom_level/${base}_ \
  --runThreadN 12 --alignIntronMax 1000  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic --quantMode GeneCounts --readFilesCommand zcat \
  --outWigType wiggle --outBAMsortingBinsN 200 --limitBAMsortRAM 3000000000
elif [[ -f output/fastq_trimmed/${base}_1.fastq.gz && -f output/fastq_trimmed/${base}_2.fastq.gz ]]; then
  STAR --readFilesIn output/fastq_trimmed/${base}_1.fastq.gz \
  output/fastq_trimmed/${base}_2.fastq.gz \
  --genomeDir ../star_db/Af3357_chrom_level/ \
  --outFileNamePrefix output/star_chrom_level/${base}_ \
  --runThreadN 12 --alignIntronMax 1000  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic --quantMode GeneCounts --readFilesCommand zcat \
  --outWigType wiggle --outBAMsortingBinsN 200 --limitBAMsortRAM 3000000000
else
  >&2 echo "Wrong number of parameters"
fi
 STAR --runMode inputAlignmentsFromBAM \
   --inputBAMfile output/star_chrom_level/${base}_Aligned.sortedByCoord.out.bam \
   --outWigType wiggle \
   --outWigStrand Unstranded \
   --outFileNamePrefix output/star_chrom_level/${base}_unstranded
## parallel sbatch --wrap \"wigToBigWig {} ../Af3357_chrom_level_chrom.sizes {.}.bw\" ::: output/star_chrom_level/*_unstrandedSignal.Unique.str1.out.wig
