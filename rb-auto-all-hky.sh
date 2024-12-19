#!/bin/bash
#SBATCH --job-name=hky
#SBATCH --time=7:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=12
#SBATCH -p education

module load RevBayes/1.2.2-GCC-12.2.0

rb rb_base_hky.rev