#!/bin/bash
#SBATCH -p short
#SBATCH -n 32
#SBATCH -J iqtree
#SBATCH --mem=32072
#SBATCH -t 0-23:30:15   
#SBATCH -o slurm/iqtree.%j.out
#SBATCH -e slurm/iqtree.%j.err
#Tree building from the alignment
iqtree -s final_aligned.fasta -m Q.yeast+F+I+G4 -bb 1000 -nt 8
