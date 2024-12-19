#!/bin/bash
#SBATCH --job-name=gtr
#SBATCH --time=7:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=12
#SBATCH -p scavenge

module load RevBayes/1.2.2-GCC-12.2.0

rb rb_nooutgroup_gtr.rev