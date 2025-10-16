# Master's Thesis on Transcriptomic Analysis in Chromohalobacter-salexigens
Repository containing Bash and R scripts for transcriptomic analysis of Chromohalobacter salexigens.

## Analysis Description
This repository contains the Bash and R scripts used to perform a transcriptomic (RNA-seq) analysis in the halophilic bacterium Chromohalobacter salexigens. The purpose of this analysis was to determine the molecular mechanisms involving a key regulatory molecule in this bacterium, EupR. The repository provides a reproducible pipeline for data preprocessing, differential gene expression analysis, and functional analysis to facilitate future studies on halophilic gene regulation.

Initially, the analysis was performed on a remote Linux environment accessed from a Windows system through an SSH connection using MobaXterm as the terminal client. Subsequently, the analysis was carried out in R (v4.3.1) within the RStudio environment.

## RNA-seq data

A total of 14 samples of RNAseq data were analyzed. Two datasets were used:
 - The first comprised six samples from the Wt strain, three of which were cultivated at 0.6 M NaCl (Wt06A, Wt06B, Wt06C) and the remaining three at 2.5 M NaCl (Wt25A, Wt25B, Wt25C). The files used were .bam files. These data were obtained using the SOLiD sequencing platform.
 - The second dataset comprised six samples from the eupR mutant strain, three when it was cultivated at 0.6 M NaCl (EupR06A, EupR06B, EupR06C) and three at 2.5 M NaCl (EupR25A, EupR25B, EupR25C), and two samples from the Wt strain (Wt06I, Wt25I), one corresponding to each salinity. The files used were the FASTQ files generados por la plataforma de secuenciación Illumina.

## Analyses performed in a remote Linux environment
1) Sample processing using the cat command to obtain the complete R1 and R2 (paired-end) files for each sample in .fastq format.
2) Quality assessment of the reads using the FastQC software.
3) Indexing of the reference genome with the Burrows-Wheeler Aligner (BWA), a required step prior to mapping the reads to the reference genome.
4) Mapping of reads to the reference genome using BWA-MEM. The SAMTOOLS software was used to convert the mapping .sam files into .bam format and to analyse mapping statistics.
5) Quantification of reads mapped to the reference genome using the FeatureCounts software.

## Analyses performed in R (RStudio)
1) Normalization of read counts using the DESeq2 package.
2) Principal Component Analysis (PCA) visualization of samples grouped according to experimental conditions.
3) Differential expression analysis for each salinity condition using the DESeq2 package. Visualization of results with Venn diagrams and Volcano plots.
4) Functional analyses: Clusters of Orthologous Genes (COG), Kyoto Encyclopedia of Genes and Genomes (KEGG) and transcription factors. R was used to compare tables and quantify which categories from each analysis appeared in the set of differentially expressed genes.

## Scripts information
- bash_script.sh: includes all Bash code used for sample processing, written as a command-line workflow. Two additional scripts were generated for specific processes, which are also included.
- analysis_fastqc.sh: includes only the quality analysis of the reads performed with the FASTQC software.
- map_bwa.sh: includes only the mapping of reads to the reference genome using the Burrows-Wheeler Aligner (BWA) with the BWA-MEM algorithm.
- r_analysis.sh: includes all code used in R.
