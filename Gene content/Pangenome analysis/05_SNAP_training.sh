#!/bin/bash
#Script for training a model for SNAP gene prediction
#SBATCH -p short
#SBATCH -n 32
#SBATCH -J extractgenomes
#SBATCH --mem 32768
#SBATCH -t 0-12:30:15   
#SBATCH -o slurm/extractgenomes.%j.out
#SBATCH -e slurm/extractgenomes.%j.err



cd  /mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/assemblies/train_yl
#Create hmm model
../SNAP/hmm-assembler.pl Y.lipolytica . > at2.hmm