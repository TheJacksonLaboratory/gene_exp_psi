#!/usr/bin/env bash

#SBATCH -J gene_psi
#SBATCH -t 1:00
#SBATCH -n 1
#SBATCH --mem=100M
#SBATCH -o %u-%x-%j

# NOTE THE FOLLOWING ENVIROMENT MUST BE CREATED
# BEFORE RUNNING THIS SCRIPT -- SEE SLURM.md
conda activate esvenv

sbatch run_snake.sh