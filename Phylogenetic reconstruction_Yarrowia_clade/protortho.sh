#!/bin/bash
#SBATCH -p long
#SBATCH -n 64
#SBATCH -J braker
#SBATCH --mem=64072
#SBATCH -t 0-23:30:15   
#SBATCH -o slurm/braker.%j.out
#SBATCH -e slurm/braker.%j.err

# Cambiar al directorio de trabajo
cd /mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/assemblies/annotation/protfasta/weird 

#Run proteinOrtho
proteinortho ./format/*.aa