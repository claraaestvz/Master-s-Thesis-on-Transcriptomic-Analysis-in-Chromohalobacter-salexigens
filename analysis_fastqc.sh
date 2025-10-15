#!/bin/bash

# Relative paths from scripts/
SAMPLES_DIR="../samples"
RESULTS_DIR="../results/fastqc"

# Iterate over all FASTQ files in the sample folders
for sample in "$SAMPLES_DIR"/*; do
    for file in "$sample"/*.fastq.gz; do
        echo "Analyzing $file..."
        fastqc "$file" -o "$RESULTS_DIR"
    done
done

echo "FastQC analysis completed. Results saved in $RESULTS_DIR"

