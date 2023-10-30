#!/bin/bash
#SBATCH --job-name="salmon_db"
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"

set -x
set -e
set -u

module load salmon/1.9.0

wget -P data/salmon_db https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/275/GCA_000006275.2_JCVI-afl1-v2.0/GCA_000006275.2_JCVI-afl1-v2.0_rna_from_genomic.fna.gz
wget -P data/salmon_db https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/275/GCA_000006275.2_JCVI-afl1-v2.0/GCA_000006275.2_JCVI-afl1-v2.0_genomic.fna.gz
wget -P data/salmon_db https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/017/415/GCA_009017415.1_ASM901741v1/GCA_009017415.1_ASM901741v1_rna_from_genomic.fna.gz
wget -P data/salmon_db https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/017/415/GCA_009017415.1_ASM901741v1/GCA_009017415.1_ASM901741v1_genomic.fna.gz
gunzip data/salmon_db/*


perl -i -pe "s/>.*(AFLA_\d{6}).*/>\1/" data/salmon_db/GCA_000006275.2_JCVI-afl1-v2.0_rna_from_genomic.fna 
grep "^>" data/salmon_db/GCA_000006275.2_JCVI-afl1-v2.0_genomic.fna | cut -d " " -f 1 \
  | cut -d ">" -f2 > data/salmon_db/decoys_jcvi.txt
cat data/salmon_db/GCA_000006275.2_JCVI-afl1-v2.0_rna_from_genomic.fna \
  data/salmon_db/GCA_000006275.2_JCVI-afl1-v2.0_genomic.fna \
   > data/salmon_db/gentrome_jcvi.fa
salmon index -t data/salmon_db/gentrome_jcvi.fa -d data/salmon_db/decoys_jcvi.txt -p 2 \
  -i data/salmon_db/salmon_index_jcvi

perl -i -pe "s/>.*(F9C07_\d.*?)\].*/>\1/" data/salmon_db/GCA_009017415.1_ASM901741v1_rna_from_genomic.fna
grep "^>" data/salmon_db/GCA_009017415.1_ASM901741v1_genomic.fna | cut -d " " -f 1 \
  | cut -d ">" -f2 > data/salmon_db/decoys_chrom_level.txt
cat data/salmon_db/GCA_009017415.1_ASM901741v1_rna_from_genomic.fna \
  data/salmon_db/GCA_009017415.1_ASM901741v1_genomic.fna \
   > data/salmon_db/gentrome_chrom_level.fa
salmon index -t data/salmon_db/gentrome_chrom_level.fa -d data/salmon_db/decoys_chrom_level.txt -p 2 \
  -i data/salmon_db/salmon_chrom_level