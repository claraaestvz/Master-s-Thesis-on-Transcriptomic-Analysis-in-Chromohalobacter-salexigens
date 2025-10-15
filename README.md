# Master's Thesis on Transcriptomic Analysis in Chromohalobacter-salexigens
Repository containing Bash and R scripts for transcriptomic analysis of Chromohalobacter salexigens.

## Analysis Description
This repository contains the Bash and R scripts used to perform a transcriptomic (RNA-seq) analysis in the halophilic bacterium Chromohalobacter salexigens. The purpose of this analysis was to determine the molecular mechanisms involving a key regulatory molecule in this bacterium, EupR. The repository provides a reproducible pipeline for data preprocessing, differential gene expression analysis, and functional analysis to facilitate future studies on halophilic gene regulation.

Initially, the analysis was performed on a remote Linux environment accessed from a Windows system through an SSH connection using MobaXterm as the terminal client. Subsequently, the analysis was carried out in R (v4.3.1) within the RStudio environment.

## RNA-seq data

A total of 14 samples of RNAseq data were analyzed. Two datasets were used:
 - The first comprised six samples from the Wt strain, three of which were cultivated at 0.6 M NaCl (Wt06A, Wt06B, Wt06C) and the remaining three at 2.5 M NaCl (Wt25A, Wt25B, Wt25C). The files used were .bam files. These data were obtained using the SOLiD sequencing platform.
 - The second dataset comprised six samples from the eupR mutant strain, three when it was cultivated at 0.6 M NaCl (EupR06A, EupR06B, EupR06C) and three at 2.5 M NaCl (EupR25A, EupR25B, EupR25C), and two samples from the Wt strain (Wt06I, Wt25I), one corresponding to each salinity. The files used were the FASTQ files generados por la plataforma de secuenciación Illumina.

## Análisis llevados a cabo en remote Linux environment
1) 

## Scripts information
- bash_script.sh: includes all Bash code used for sample processing, written as a command-line workflow. Two additional scripts were generated for specific processes, which are also included.
- analysis_fastqc.sh: includes only the quality analysis of the reads performed with the FASTQC software.
- map_bwa.sh: includes only the mapping of reads to the reference genome using the Burrows-Wheeler Aligner (BWA) with the BWA-MEM algorithm.
- r_analysis.sh: includes all code used in R.
