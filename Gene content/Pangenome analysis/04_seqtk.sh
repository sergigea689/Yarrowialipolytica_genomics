#!/bin/bash
#SBATCH -p short
#SBATCH -n 32
#SBATCH -J compgenom2
#SBATCH --mem 32768
#SBATCH -t 0-20:30:15   
#SBATCH -o slurm/fastp.%j.out
#SBATCH -e slurm/fastp.%j.err

#Script to filter out contigs of less than 1000 bp out of genome assemblies. It will make the pangenome analysis more restrictive.
# Folders
input_dir="/mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/assemblies/assembly"
output_dir="/mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/assemblies/filtassemblies"

# Create folder
#mkdir -p "$output_dir"

# Iterate over each fasta (assemblies)
for file in "$input_dir"/*.fasta; do
    # Obtain basenames
    filename=$(basename "$file")
    
    # Filter according to contig length
    seqtk seq -L 1000 "$file" > "$output_dir/$filename"
    
    echo "Procesado: $filename"
done