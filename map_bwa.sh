#!/bin/bash

# Arguments
DIR="$1"   # folder path
NAME="$2"  # base sample name

# Fixed paths to genome and annotation
GENOME="../../genome/index_bwa_mem/genome_reference_israelensis.fa"
ANNOTATION="../../genome/annotation_israelensis.gtf"

echo "[INFO] Mapping reads with BWA..."
bwa mem -t 4 -M -T 30 "$GENOME" \
  "${NAME}_R1.fastq.gz" "${NAME}_R2.fastq.gz" > "${NAME}.sam"

echo "[INFO] Converting SAM to BAM..."
samtools view -bS "${NAME}.sam" > "${NAME}.bam"

echo "[INFO] Generating mapping statistics..."
samtools flagstat "${NAME}.bam" > "${NAME}_flagstat.txt"
samtools stats "${NAME}.bam" > "${NAME}_stats.txt"

echo "[INFO] Sorting BAM by coordinates..."
samtools sort "${NAME}.bam" -o "${NAME}_sorted.bam"

echo "[COMPLETED] Files generated in: $DIR"

