#!/bin/bash
#SBATCH --job-name=phylo_0.95
#SBATCH --partition=hsph,shared,intermediate,sapphire
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --output=finish_fullseq_logs/phylo_0.95_%j.out
#SBATCH --error=finish_fullseq_logs/phylo_0.95_%j.err

# Create logs directory
mkdir -p finish_fullseq_logs

echo "Starting phylo finish pipeline at $(date)"
echo "Job ID: $SLURM_JOB_ID"

# Unlock
snakemake -s phylo_finish_baseilne_0.95.smk --unlock

# Run with detailed logging
snakemake -s phylo_finish_baseline_0.95.smk --executor slurm --default-resources slurm_account=grad_lab --retries 3 --restart-times 2 --jobs 25 --latency-wait 300 --rerun-incomplete --printshellcmds --verbose --keep-going

echo "Pipeline completed at $(date)"
