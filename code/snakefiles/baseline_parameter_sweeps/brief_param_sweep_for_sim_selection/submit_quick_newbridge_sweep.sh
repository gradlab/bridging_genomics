#!/bin/bash
#SBATCH --job-name=quick_param_sweep
#SBATCH --partition=hsph,shared,intermediate,sapphire
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=param_sweep_10_28_bridge_%j.out
#SBATCH --error=param_sweep_10_28_bridge_%j.err

source ~/.bashrc
conda activate phylo
cd /n/netscratch/grad_lab/Lab/mkline/bridging_project/code/snakefiles

echo "Starting QUICK parameter sweep pipeline (100 sims per value)..."

snakemake -s param_sweep_batched_100_new_bridge.smk --unlock

snakemake -s param_sweep_batched_100_new_bridge.smk --executor slurm --default-resources slurm_account=grad_lab --retries 3 --restart-times 2 --jobs 50 --latency-wait 300 --rerun-incomplete --printshellcmds --verbose --keep-going

echo "Quick parameter sweep pipeline completed!"
