#!/bin/bash
#SBATCH --job-name=submit_snakemake
#SBATCH --partition=shared,intermediate,hsph,sapphire
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=full_0.5_snakemake_fixed_%j.out
#SBATCH --error=full_0.5_snakemake_fixed_%j.err

# Load conda environment
source ~/.bashrc
conda activate phylo

# Change to the right directory
cd /n/netscratch/grad_lab/Lab/mkline/bridging_project/code/snakefiles

# Run the pipeline
snakemake -s basic_sim_run_0.5_rest.smk \
    --executor slurm \
    --default-resources slurm_account=grad_lab \
    --jobs 5 \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds \
    --verbose 
