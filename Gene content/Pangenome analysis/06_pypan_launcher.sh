#!/bin/bash
#SBATCH -p long
#SBATCH -n 64
#SBATCH -J pypan2
#SBATCH --mem=131072
#SBATCH -t 9-23:30:15   
#SBATCH -o slurm/pypan2.%j.out
#SBATCH -e slurm/pypan2.%j.err

#Install Pypan dependencies accordingly to https://github.com/jsgounot/PyPan
#It is a complex pipeline involving Blast+, Augustus and SNAP, troubleshooting may be required to ensure proper linkage between all the steps of the program.
#Eventually, dividing the pipeline in sequential steps and running each separately can fix some of the errors related to segmentation fault.

# Comando para ejecutar el script de Python
python /mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/assemblies/PyPan-master/basic_launcher.py /mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/assemblies/yali_config.json yarr --ncore 64