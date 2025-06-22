#!/bin/bash

# Script to align and quantify reads mapping to each MAT locus
# This script classifies strains based on the relative number of reads mapping to MATA and MATB reference sequences

# Requirements: bwa, samtools, paste, sed
# Ensure required programs are available
for cmd in bwa samtools paste sed; do
  if ! command -v $cmd &> /dev/null; then
    echo "Error: Required program '$cmd' is not installed or not in PATH."
    exit 1
  fi
done

# Create output directory if it does not exist
mkdir -p mat

# Create results file with header
echo -e "strain\tMATA_reads\tMATB_reads\tclassification" > mat/ploidy_counts.tsv

# Read file pairs from the trimming list and process each pair
# 'trimlist.txt' is expected to contain alternating lines of forward and reverse read filenames
paste - - < trimlist.txt | while read R1 R2; do

  # Extract the strain name by removing the final character ('1' or '2') from the filename
  strain=$(echo "$R1" | sed 's/1$//')

  echo "Processing $strain..."

  # Define full filenames for read pairs
  R1_file="${strain}1.fq.gz"
  R2_file="${strain}2.fq.gz"

  # Align reads to MATA reference and count mapped reads
  count_MATA=$(bwa mem MATA_W29.fasta "$R1_file" "$R2_file" | samtools view -F 4 -c -)

  # Align reads to MATB reference and count mapped reads
  count_MATB=$(bwa mem MATB_E150.fasta "$R1_file" "$R2_file" | samtools view -F 4 -c -)

  # Classify based on which locus had more mapped reads
  if [[ "$count_MATA" -gt "$count_MATB" ]]; then
    classification="MatA"
  elif [[ "$count_MATA" -lt "$count_MATB" ]]; then
    classification="MatB"
  else
    classification="Mixed"
  fi

  # Append results to output file
  echo -e "$strain\t$count_MATA\t$count_MATB\t$classification" >> mat/ploidy_counts.tsv

done
