#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH -p education

module load IQ-TREE

iqtree2 -s mito_abl.fasta -bb 1000 -nt AUTO