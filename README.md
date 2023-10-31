# Aspergillus flavus expression database

This repository contains code to download all Aspergillus flavus RNA-Seq data from the NCBI SRA, process it to generate gene expression values, and visualize the results in a Shiny app.

## Pipeline

The pipeline code is in the `pipeline` directory. It performs the following steps:

1. Identifies A. flavus RNA-Seq experiments available on NCBI and downloads metadata using rentrez. Metadata about the samples is stored in `output/sra_metadata.csv`.
2. Downloads FASTQ files from the NCBI SRA using fastq-dl. 
3. Removes adapters and low quality sequence from the reads using fastp.
4. Generates gene counts and TPM using Salmon.  

## Shiny App
The Shiny app code is in the `shiny_app` directory. It displays the expression data along with sample metadata in an interactive web interface.
The app has two tabs:

1. Barplot: Displays expression for a selected gene across samples as a barplot. Samples are grouped and colored by their bioproject. Clicking on any bar shows the available metadata for that sample.
2. Heatmap: Shows a heatmap of genes and bioprojects selected by the user. The selection of the genes displayed in the heatmap can be made using various functional annotations available or a comma separated list of genes can be provided.

