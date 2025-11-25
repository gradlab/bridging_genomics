#!/bin/bash
#SBATCH --job-name=launch_remaining_47
#SBATCH --partition=hsph,shared,intermediate,sapphire
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=test_logs/launch_remaining_47_%j.out
#SBATCH --error=test_logs/launch_remaining_47_%j.err

source ~/.bashrc
conda activate phylo
cd /n/netscratch/grad_lab/Lab/mkline/bridging_project/code/snakefiles

echo "🚀 Launching remaining 47 simulations (PROVEN WORKING APPROACH) at $(date)"

# Get all sims, exclude already running/submitted
SIM_BASE="../../output/baseline_phylo_bridging_complete_fresh_v2/phylo"
ALL_SIMS=($(ls -1 $SIM_BASE | sort))
EXCLUDE=("0.5_060" "0.05_014" "0.05_016")  # Running + 2 tests

REMAINING_SIMS=()
for sim in "${ALL_SIMS[@]}"; do
    if [[ ! " ${EXCLUDE[*]} " =~ " $sim " ]]; then
        REMAINING_SIMS+=("$sim")
    fi
done

echo "🎯 Will submit ${#REMAINING_SIMS[@]} remaining simulations"
echo "📋 Simulations to submit:"
printf '  %s\n' "${REMAINING_SIMS[@]}" | head -10
echo "  ... (and $((${#REMAINING_SIMS[@]} - 10)) more)"

SUBMITTED=0
FAILED=0

for sim_id in "${REMAINING_SIMS[@]}"; do
    echo ""
    echo "📤 Submitting $((SUBMITTED + 1))/${#REMAINING_SIMS[@]}: $sim_id"
    
    # Create individual submission script
    cat > "temp_submit_${sim_id}.sh" << EOF
#!/bin/bash
#SBATCH --job-name=pipeline_${sim_id}
#SBATCH --partition=hsph,shared,intermediate,sapphire
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1  
#SBATCH --output=test_logs/pipeline_${sim_id}_%j.out
#SBATCH --error=test_logs/pipeline_${sim_id}_%j.err

source ~/.bashrc
conda activate phylo
cd /n/netscratch/grad_lab/Lab/mkline/bridging_project/code/snakefiles

echo "🚀 Starting pipeline for $sim_id at \$(date)"

snakemake --snakefile single_sim_phylo_existing.smk \\
          --config sim_id="$sim_id" \\
          --cores 16 --jobs 4 --executor slurm \\
          --default-resources slurm_account=grad_lab \\
          --retries 2 --restart-times 1 \\
          --latency-wait 300 --rerun-incomplete \\
          --keep-going \\
          --target-files "../../output/baseline_phylo_bridging_complete_fresh_v2/phylo/$sim_id/aggregated_sampled_bridging_30_results.csv" \\
                         "../../output/baseline_phylo_bridging_complete_fresh_v2/phylo/$sim_id/aggregated_phylogenetic_bridging_30_results.csv"

echo "✅ Pipeline completed for $sim_id at \$(date)"
EOF

    # Submit and capture job ID
    chmod +x "temp_submit_${sim_id}.sh"
    JOBID=$(sbatch "temp_submit_${sim_id}.sh" 2>&1)
    
    if echo "$JOBID" | grep -q "Submitted batch job"; then
        JOB_NUM=$(echo "$JOBID" | grep -o '[0-9]\+')
        echo "✅ Submitted $sim_id as job $JOB_NUM"
        ((SUBMITTED++))
    else
        echo "❌ Failed to submit $sim_id: $JOBID"
        ((FAILED++))
    fi
    
    # Clean up temp script
    rm "temp_submit_${sim_id}.sh"
    
    # Small delay to be gentle on scheduler
    sleep 5
    
    # Longer pause every 10 submissions
    if [ $((SUBMITTED % 10)) -eq 0 ]; then
        echo "😴 Brief pause after $SUBMITTED submissions..."
        sleep 20
    fi
done

echo ""
echo "📈 FINAL SUBMISSION SUMMARY:"
echo "   Total to submit: ${#REMAINING_SIMS[@]}"
echo "   Successfully submitted: $SUBMITTED"
echo "   Failed submissions: $FAILED"
echo "   Already running: 3 (0.5_060 + 2 tests)"
echo "   TOTAL ACTIVE PIPELINES: $((SUBMITTED + 3))"

echo ""
echo "🎉 LAUNCH COMPLETED at $(date)"
echo "📊 Monitor with: watch 'squeue -u mkline | wc -l'"
