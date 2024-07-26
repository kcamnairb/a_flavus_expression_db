# Aspergillus flavus expression database

This repository contains code to download all Aspergillus flavus RNA-Seq data from the NCBI SRA, process it to generate gene expression values, and visualize the results in a Shiny app. The app can be viewed at https://a-flavus-expression-db-jyqnpeuvta-uc.a.run.app/.

## Pipeline

The code to download and process the RNA-Seq data is in the `src` directory. It performs the following steps:

1. `get_sra_metadata.Rmd` identifies A. flavus RNA-Seq experiments available on NCBI and downloads metadata using rentrez.
2. `sbatch_launcher.sh`, `download_trim_count.sh`, and `align_count_chrom_level.sh` download FASTQ files from the NCBI SRA using Kingfisher, removes adapters and low quality sequence from the reads using fastp, aligns reads and generates read counts using STAR.
3. `star_quantification.R` combines count data, examines quality control results, and creates TPM and variance stabilized transformed counts.
4. `coexpression_workflow_combinations.R`  creates combinations of different workflow steps for making different co-expression networks. `shiny_app/network_evaluation_single.R` runs functional enrichment on the top 100 edges for each gene in each network. `coexpression_network_evaluations.R` combines and examines the enrichment results.

## Shiny App
The Shiny app code is in the `shiny_app` directory. It displays the expression values and co-expression network along with sample metadata in an interactive web interface.
The app has five tools, each on a different tab:

1. Barplot: Displays expression for a selected gene across all samples as a barplot. Samples are grouped and colored by their bioproject. Clicking on any bar shows the available metadata for that sample.
2. Heatmap: Shows a heatmap of genes and bioprojects selected by the user. The selection of the genes displayed in the heatmap can be made using various functional annotations available or a comma separated list of genes can be provided.
3. Co-expression network: Shows co-expression network creating using upper quartile normalized counts and pearsons correlation. Query genes can be selected using same method as for the heatmap and their first order neihbors will be displayed. Functional enrichment of the subnetwork can also be performed.
4. PCA: Shows PCA of variance stabilized transformed counts and colors samples by various metadata.
5. JBrowse: Shows genome browser with genes and RNA-Seq coverage track representing the average across all samples.

