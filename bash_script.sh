# !/bin/bash

# First, merge all R1 files together and all R2 files separately (paired-end). Use the cat command to combine them into new R1 and R2 files.

# Wt06I
cat Wt06_ATGTCA_L008_R1_*.fastq.gz > Wt06_illumina_R1.fastq.gz && cat Wt06_ATGTCA_L008_R2_*.fastq.gz > Wt06_illumina_R2.fastq.gz

# Wt25 I
cat Wt25_CCGTCC_L008_R1_*.fastq.gz > Wt25_illumina_R1.fastq.gz && cat Wt25_CCGTCC_L008_R2_*.fastq.gz > Wt25_illumina_R2.fastq.gz

# EupR06A
cat EupR06A_CAGATC_L008_R1_*.fastq.gz > EupR06A_illumina_R1.fastq.gz && cat EupR06A_CAGATC_L008_R2_*.fastq.gz > EupR06A_illumina_R2.fastq.gz

# EupR06B
cat EupR06B_ACTTGA_L008_R1_*.fastq.gz > EupR06B_illumina_R1.fastq.gz && cat EupR06B_ACTTGA_L008_R2_*.fastq.gz > EupR06B_illumina_R2.fastq.gz

# EupR06C
cat EupR06C_GATCAG_L008_R1_*.fastq.gz > EupR06C_illumina_R1.fastq.gz && cat EupR06C_GATCAG_L008_R2_*.fastq.gz > EupR06C_illumina_R2.fastq.gz

# EupR25A
cat EupR25A_TAGCTT_L008_R1_*.fastq.gz > EupR25A_illumina_R1.fastq.gz && cat EupR25A_TAGCTT_L008_R2_*.fastq.gz > EupR25A_illumina_R2.fastq.gz

# EupR25B
cat EupR25B_GGCTAC_L008_R1_*.fastq.gz > EupR25B_illumina_R1.fastq.gz && cat EupR25B_GGCTAC_L008_R2_*.fastq.gz > EupR25B_illumina_R2.fastq.gz

# EupR25C
cat EupR25C_CTTGTA_L008_R1_*.fastq.gz > EupR25C_illumina_R1.fastq.gz && cat EupR25C_CTTGTA_L008_R2_*.fastq.gz > EupR25C_illumina_R2.fastq.gz


# Each pair of samples was stored in folders named sample01 through sample08, following the previous order. Other folders were generated for the execution of the analysis.
# Different scripts were generated for the various steps of the analysis and are indicated as follows: === script: .sh ===

# Quality analysis of FASTQ samples with FASTQC
# === script: analysis_fastq.sh ===

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
###########################################################################

# Reference genome index using Burrow-Wheeler Aligner (BWA)
bwa index -a is genome_reference_israelensis.fa

# Mapping of fastq.gz files to the reference genome (BWA-MEM)
# === script: map_bwa.sh ===

#!/bin/bash

# Arguments
DIR="$1"   # folder name
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
###############################################################################

# The .bam files from the samples sequenced with the SOLiD platform were added to a folder named sample09 (Wt06A.bam, Wt06B.bam, Wt06C.bam, Wt25A.bam, Wt25B.bam, y Wt25C.bam)
for bam in ./sample09/*.bam; do
    flagstat="${bam%.bam}_flagstat.txt"
    stats="${bam%.bam}_stats.txt"
    samtools flagstat "$bam" > "$flagstat" &
    samtools stats "$bam" > "$stats" &
done

# Quantification of reads with featureCounts
featureCounts -T 12 -a ../genome/annotation_israelensis.gtf -o gene_matrix_counts.txt -t CDS -g gene_id -p ./sample0*/*.bam
