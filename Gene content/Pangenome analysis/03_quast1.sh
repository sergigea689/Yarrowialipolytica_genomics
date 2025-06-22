#!/bin/bash
#SBATCH -p short
#SBATCH -n 32
#SBATCH -J compgenom2
#SBATCH --mem 32768
#SBATCH -t 0-22:30:15   
#SBATCH -o slurm/fastp.%j.out
#SBATCH -e slurm/fastp.%j.err

#Script de quast para la evaluación de la calidad de los ensamblajes.
#Se puede hacer una evaluación con referencia o sin ella.
#La evaluación con referencia permite encontrar mismatches, ya sean artificiales o variaciones estructurales verdaderas.
#Además da una serie de parámetros comparativos con respecto al ensamblaje de referencia.
#Necesitamos iterar por cada carpeta del directorio y trabajar sobre el archivo scaffolds output de spades.
#La evaluación de la calidad de los ensamblajes se realizará con quast, que se encuentra instalado en el env spades de conda.

cd /mnt/lustre/home/sergiz01/genomes/Ylipgenomes/fastq/rawreads/trimmed_reads/assemblies
for Variable in $(< ../listagenomes1.txt)
do
	mkdir -p ./"$Variable"/quast
	quast.py ./"$Variable"/scaffolds.fasta -r ./YALI0-1.fasta -o ./"$Variable"/quast
done
