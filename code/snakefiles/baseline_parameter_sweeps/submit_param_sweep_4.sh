#!/bin/bash
#SBATCH --job-name=param_sweep_4
#SBATCH --partition=hsph,shared,intermediate,sapphire
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=param_sweep_250_logs/param_sweep_baseline_4_%j.out
#SBATCH --error=param_sweep_250_logs/param_sweep_baseline_4_%j.err

source ~/.bashrc
conda activate phylo
cd /n/netscratch/grad_lab/Lab/mkline/bridging_project/code/snakefiles

echo "Starting parameter sweep pipeline 4 (250 sims per value, seed 3000)..."

snakemake -s param_sweep_baseline_250sims_4.smk --unlock

snakemake -s param_sweep_baseline_250sims_4.smk --executor slurm --default-resources slurm_account=grad_lab --retries 3 --restart-times 2 --jobs 25 --latency-wait 300 --rerun-incomplete --printshellcmds --verbose --keep-going

echo "Parameter sweep pipeline 4 completed!"
