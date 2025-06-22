#!/bin/bash
#SBATCH -p short
#SBATCH -n 64
#SBATCH -J spades
#SBATCH --mem 64768
#SBATCH -t 0-20:30:15   
#SBATCH -o slurm/spades.%j.out
#SBATCH -e slurm/spades.%j.err
cd /mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/cleanedcel
# Vreate output directory
mkdir -p "./assemblies"

# Iterate over each sample in list
for Variable in $(< listgenome1.txt); do
    # Avoid re-doing assemblies
    if [ ! -d "./assemblies/${Variable}" ]; then
        echo "Processing ${Variable}..."
        spades.py -s "./${Variable}" -o "./assemblies/${Variable}" --isolate
    else
        echo "Directory ./assemblies/${Variable} already exists. Skipping..."
    fi
done

echo "Assembly process completed."