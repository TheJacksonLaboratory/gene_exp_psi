#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --job-name="pps"
#SBATCH -o ppr-%j.out
#SBATCH -e ppr-%j.err 
#SBATCH --mail-user=guy.karlebach@jax.org
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -N 1
#SBATCH -n 31
#SBATCH --mem-per-cpu=24G

cd $SLURM_SUBMIT_DIR

module load singularity

singularity exec sing.sif bash run_snakemake.sh 
















	







