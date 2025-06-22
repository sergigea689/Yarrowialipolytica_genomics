#!/bin/bash
#SBATCH -p long
#SBATCH -n 64
#SBATCH -J cel
#SBATCH --mem=131072
#SBATCH -t 9-23:30:15   
#SBATCH -o slurm/cel.%j.out
#SBATCH -e slurm/cel.%j.err

#Pangenome preparation pipeline
#we start from trimmed reads already for quality.
#Since we detected a contamination specific to Cellulosimicrobium cellulans in 70 samples, we are going to eliminate it to keep it from altering our pangenome analysis.
# Define variables
GENOME="celcellulans.fasta"
INPUT_DIR="/mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads"
OUTPUT_DIR="/mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/cleanedcel"
STATS_FILE="/mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/cleanedcel/mapping_summary.csv"

# Create output directories
mkdir -p "$OUTPUT_DIR"

# Create CSV file with a header
echo "Sample,Total_Reads,Mapped_Reads,Unmapped_Reads,Percentage_Mapped" > "$STATS_FILE"

# Iterate over each paired-end file fq.gz (R1 y R2)
for R1_FILE in "$INPUT_DIR"/*1.fq.gz; do
    # Verify initial string 'YL'
    BASENAME=$(basename "$R1_FILE" "1.fq.gz")
    if [[ ! "$BASENAME" =~ ^YL ]]; then
        continue  
    fi

    R2_FILE="$INPUT_DIR/${BASENAME}2.fq.gz"

    if [ ! -f "$R2_FILE" ]; then
        echo "Error: R2 file for $BASENAME not found!"
        continue
    fi

    echo "Processing $BASENAME..."

    # Alignment and extraction of unmapped reads
    MAPPED_READS=$(bwa mem "$GENOME" "$R1_FILE" "$R2_FILE" | samtools view -b -f 4 - | samtools fastq - > "$OUTPUT_DIR/${BASENAME}_clean.fq")

    # Calculate mapping strategies
    TOTAL_READS=$(zcat "$R1_FILE" "$R2_FILE" | echo $((`wc -l` / 4)))
    MAPPED_READS=$(samtools view -c -F 4 "$GENOME.bam")
    UNMAPPED_READS=$((TOTAL_READS - MAPPED_READS))
    PERCENTAGE_MAPPED=$(echo "scale=2; ($MAPPED_READS / $TOTAL_READS) * 100" | bc)

    # Save statistics in CSV
    echo "$BASENAME,$TOTAL_READS,$MAPPED_READS,$UNMAPPED_READS,$PERCENTAGE_MAPPED" >> "$STATS_FILE"

    echo "Processed $BASENAME: Total Reads: $TOTAL_READS, Mapped Reads: $MAPPED_READS, Unmapped Reads: $UNMAPPED_READS, Percentage Mapped: $PERCENTAGE_MAPPED%"
done

echo "Mapping summary saved to $STATS_FILE"


