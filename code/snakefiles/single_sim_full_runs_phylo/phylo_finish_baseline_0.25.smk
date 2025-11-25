# ============================================================================
# Phylogenetic Bridging Analysis Pipeline - Full Sequences Only
# Picks up from completed IQTree analysis and finishes downstream steps
# ============================================================================

# CONFIGURATION - CUSTOMIZE FOR EACH DIRECTORY
SIM_DIR = "../../output/baseline_params_0.25_rep_357"  # CHANGE THIS FOR EACH SNAKEFILE
PIPELINE_SEED = 42

# Software paths
LSD2_PATH = "/n/holylfs05/LABS/grad_lab/Users/mkline/bridging_sims/software/lsd2/src/lsd2"

# Resource scaling and partition functions
def get_compute_partition(wildcards, attempt):
    partitions = ["hsph", "sapphire", "intermediate", "shared"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

def get_long_partition(wildcards, attempt):
    partitions = ["intermediate", "hsph", "sapphire", "shared"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

# ============================================================================
# Rules
# ============================================================================

rule all:
    input:
        f"{SIM_DIR}/phylo_analysis/filtering/full_sequences_bridging_analysis.csv",
        f"{SIM_DIR}/phylo_analysis/filtering/full_sequences_bridging_analysis_counts.csv",
        f"{SIM_DIR}/phylo_analysis/filtering/full_sequences_bridging_analysis_totals.csv"

rule extract_dates_and_states:
    input:
        treefile = f"{SIM_DIR}/phylo_analysis/filtering/full_sequences.fasta.treefile"
    output:
        dates_file = f"{SIM_DIR}/phylo_analysis/filtering/dates.txt",
        msmw_msm_states = f"{SIM_DIR}/phylo_analysis/filtering/behaviors_msmw_msm_states.csv",
        msmw_msw_states = f"{SIM_DIR}/phylo_analysis/filtering/behaviors_msmw_msw_states.csv",
        extraction_flag = f"{SIM_DIR}/phylo_analysis/filtering/extraction_complete.flag"
    resources:
        mem_mb = 8000,
        runtime = 60,
        cpus_per_task = 1,
        slurm_partition = get_compute_partition
    conda: "phylo"
    shell:
        """
        echo "🕐 Starting date and state extraction at $(date)"
        
        # Create filtered transmission data for full sequences (post-burnin, non-superseded, sampled)
        python -c "
import pandas as pd
import json
import os

# Load full transmission data
transmission_df = pd.read_csv('{SIM_DIR}/transmission_df.csv')
print(f'Loaded {{len(transmission_df)}} total transmissions')

# Get burn-in cutoff from parameters
try:
    with open('{SIM_DIR}/parameters_used.json', 'r') as f:
        params = json.load(f)
    cutoff_day = (params['simulation']['partnership_burnin_days'] + 
                 params['simulation']['transmission_burnin_days'])
    print(f'Using burn-in cutoff: {{cutoff_day}}')
except Exception as e:
    print(f'Warning: Could not read parameters, using default cutoff 12000: {{e}}')
    cutoff_day = 12000

# Apply filtering (same as tree building)
filtered_transmissions = transmission_df[
    (transmission_df['superseded_simultaneous'] == False) &
    (transmission_df['day_of_transmission'] >= cutoff_day) &
    (pd.notna(transmission_df['day_of_sampling']))
]

print(f'After filtering: {{len(filtered_transmissions)}} transmissions')

# Save filtered data
filtered_transmissions.to_csv('{SIM_DIR}/phylo_analysis/filtering/temp_filtered_transmissions.csv', index=False)
"

        # Extract dates
        echo "Extracting tip dates..."
        python ../scripts/10_27_extract_tip_dates.py \\
            "{SIM_DIR}/phylo_analysis/filtering/temp_filtered_transmissions.csv" \\
            "{SIM_DIR}/phylo_analysis/filtering/full_sequences.fasta" \\
            "{SIM_DIR}" \\
            "{output.dates_file}"

        # Extract states for MSMW+MSM
        echo "Extracting MSMW+MSM states..."
        python ../scripts/10_27_extract_tip_states.py \\
            "{SIM_DIR}/phylo_analysis/filtering/temp_filtered_transmissions.csv" \\
            "{SIM_DIR}/phylo_analysis/filtering/full_sequences.fasta" \\
            "{SIM_DIR}" \\
            "{SIM_DIR}/phylo_analysis/filtering/behaviors"

        # Clean up temporary file
        rm -f "{SIM_DIR}/phylo_analysis/filtering/temp_filtered_transmissions.csv"

        # Verify outputs were created
        if [ ! -f "{output.dates_file}" ]; then
            echo "❌ ERROR: Dates file not created"
            exit 1
        fi

        if [ ! -f "{output.msmw_msm_states}" ]; then
            echo "❌ ERROR: MSMW+MSM states file not created" 
            exit 1
        fi

        if [ ! -f "{output.msmw_msw_states}" ]; then
            echo "❌ ERROR: MSMW+MSW states file not created"
            exit 1
        fi

        echo "EXTRACTION_COMPLETE=True" > {output.extraction_flag}
        echo "✅ Date and state extraction completed at $(date)"
        """

rule run_lsd2:
    input:
        treefile = f"{SIM_DIR}/phylo_analysis/filtering/full_sequences.fasta.treefile",
        dates_file = f"{SIM_DIR}/phylo_analysis/filtering/dates.txt",
        extraction_flag = f"{SIM_DIR}/phylo_analysis/filtering/extraction_complete.flag"
    output:
        lsd2_result_nwk = f"{SIM_DIR}/phylo_analysis/filtering/full_sequences.fasta.treefile.result.nwk",
        lsd2_result = f"{SIM_DIR}/phylo_analysis/filtering/full_sequences.fasta.treefile.result",
        lsd2_flag = f"{SIM_DIR}/phylo_analysis/filtering/lsd2_complete.flag"
    resources:
        mem_mb = 200000,  # 200GB for full sequences
        runtime = 2880,   # 48 hours max
        cpus_per_task = 8,
        slurm_partition = get_long_partition
    shell:
        """
        echo "🕐 Starting LSD2 analysis at $(date)"
        
        cd {SIM_DIR}/phylo_analysis/filtering
        
        {LSD2_PATH} \\
            -i full_sequences.fasta.treefile \\
            -d dates.txt \\
            -s 10000 -r a -l 0
        
        if [ ! -f "full_sequences.fasta.treefile.result.nwk" ]; then
            echo "❌ ERROR: LSD2 did not create result file"
            exit 1
        fi
        
        echo "LSD2_COMPLETE=True" > {output.lsd2_flag}
        echo "✅ LSD2 completed at $(date)"
        """

rule run_treetime_msmw_msm:
    input:
        lsd2_result_nwk = f"{SIM_DIR}/phylo_analysis/filtering/full_sequences.fasta.treefile.result.nwk",
        msmw_msm_states = f"{SIM_DIR}/phylo_analysis/filtering/behaviors_msmw_msm_states.csv",
        lsd2_flag = f"{SIM_DIR}/phylo_analysis/filtering/lsd2_complete.flag"
    output:
        msmw_msm_nexus = f"{SIM_DIR}/phylo_analysis/filtering/msmw_msm/annotated_tree.nexus",
        treetime_msmw_msm_flag = f"{SIM_DIR}/phylo_analysis/filtering/treetime_msmw_msm_complete.flag"
    resources:
        mem_mb = 16000,  
        runtime = 360,   
        cpus_per_task = 4,
        slurm_partition = get_long_partition
    conda: "phylo"
    shell:
        """
        echo "🌳 Starting TreeTime MSMW+MSM analysis at $(date)"
        
        mkdir -p {SIM_DIR}/phylo_analysis/filtering/msmw_msm
        
        treetime mugration \\
            --tree "{input.lsd2_result_nwk}" \\
            --states "{input.msmw_msm_states}" \\
            --attribute state \\
            --outdir "{SIM_DIR}/phylo_analysis/filtering/msmw_msm" \\
            --verbose 2 --confidence
        
        if [ ! -f "{output.msmw_msm_nexus}" ]; then
            echo "❌ ERROR: TreeTime MSMW+MSM did not create annotated tree"
            exit 1
        fi
        
        echo "TREETIME_MSMW_MSM_COMPLETE=True" > {output.treetime_msmw_msm_flag}
        echo "✅ TreeTime MSMW+MSM completed at $(date)"
        """

rule run_treetime_msmw_msw:
    input:
        lsd2_result_nwk = f"{SIM_DIR}/phylo_analysis/filtering/full_sequences.fasta.treefile.result.nwk",
        msmw_msw_states = f"{SIM_DIR}/phylo_analysis/filtering/behaviors_msmw_msw_states.csv",
        lsd2_flag = f"{SIM_DIR}/phylo_analysis/filtering/lsd2_complete.flag"
    output:
        msmw_msw_nexus = f"{SIM_DIR}/phylo_analysis/filtering/msmw_msw/annotated_tree.nexus",
        treetime_msmw_msw_flag = f"{SIM_DIR}/phylo_analysis/filtering/treetime_msmw_msw_complete.flag"
    resources:
        mem_mb = 16000,    
        runtime = 360,  
        cpus_per_task = 4,
        slurm_partition = get_long_partition
    conda: "phylo"
    shell:
        """
        echo "🌳 Starting TreeTime MSMW+MSW analysis at $(date)"
        
        mkdir -p {SIM_DIR}/phylo_analysis/filtering/msmw_msw
        
        treetime mugration \\
            --tree "{input.lsd2_result_nwk}" \\
            --states "{input.msmw_msw_states}" \\
            --attribute state \\
            --outdir "{SIM_DIR}/phylo_analysis/filtering/msmw_msw" \\
            --verbose 2 --confidence
        
        if [ ! -f "{output.msmw_msw_nexus}" ]; then
            echo "❌ ERROR: TreeTime MSMW+MSW did not create annotated tree"
            exit 1
        fi
        
        echo "TREETIME_MSMW_MSW_COMPLETE=True" > {output.treetime_msmw_msw_flag}
        echo "✅ TreeTime MSMW+MSW completed at $(date)"
        """

rule run_phylogenetic_bridging:
    input:
        msmw_msm_nexus = f"{SIM_DIR}/phylo_analysis/filtering/msmw_msm/annotated_tree.nexus",
        msmw_msw_nexus = f"{SIM_DIR}/phylo_analysis/filtering/msmw_msw/annotated_tree.nexus",
        lsd2_result_file = f"{SIM_DIR}/phylo_analysis/filtering/full_sequences.fasta.treefile.result",
        treetime_msmw_msm_flag = f"{SIM_DIR}/phylo_analysis/filtering/treetime_msmw_msm_complete.flag",
        treetime_msmw_msw_flag = f"{SIM_DIR}/phylo_analysis/filtering/treetime_msmw_msw_complete.flag"
    output:
        bridging_analysis = f"{SIM_DIR}/phylo_analysis/filtering/full_sequences_bridging_analysis.csv",
        bridging_counts = f"{SIM_DIR}/phylo_analysis/filtering/full_sequences_bridging_analysis_counts.csv", 
        bridging_totals = f"{SIM_DIR}/phylo_analysis/filtering/full_sequences_bridging_analysis_totals.csv",
        bridging_flag = f"{SIM_DIR}/phylo_analysis/filtering/bridging_complete.flag"
    resources:
        mem_mb = 32000,   # 32GB for bridging analysis
        runtime = 720,    # 12 hours
        cpus_per_task = 2,
        slurm_partition = get_compute_partition
    conda: "phylo"
    shell:
        """
        echo "🌲 Starting phylogenetic bridging analysis at $(date)"
        
        python ../scripts/10_30_treetime_bridging.py \\
            "{input.msmw_msm_nexus}" \\
            "{input.msmw_msw_nexus}" \\
            "{input.lsd2_result_file}" \\
            --output "{output.bridging_analysis}" \\
            --detailed-output "{SIM_DIR}/phylo_analysis/filtering/full_sequences_bridging_detailed.csv" \\
            --counts-output "{output.bridging_counts}" \\
            --totals-output "{output.bridging_totals}"
        
        if [ ! -f "{output.bridging_analysis}" ]; then
            echo "❌ ERROR: Bridging analysis did not create main output file"
            exit 1
        fi
        
        echo "BRIDGING_COMPLETE=True" > {output.bridging_flag}
        echo "✅ Phylogenetic bridging analysis completed at $(date)"
        """
