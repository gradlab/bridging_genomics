#!/usr/bin/env python3
"""
Complete Baseline Phylogenetic Bridging Pipeline
Simulation → Ground Truth → Phylogenetic Analysis → Comparison
"""

import json
import os
import glob
from pathlib import Path
import math
import pandas as pd

# Configuration
PARAM_DIR = "../parameters/baseline_parm_selected_seeds"
OUTPUT_BASE = "../../output/baseline_phylo_bridging_complete_fresh_v2"
PIPELINE_SEED = 42

# Phylo pipeline parameters
DOWNSAMPLE_RATES = [0.5, 0.25, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01]
N_REPLICATES = 10
HIGH_ACTIVITY_SUBSAMPLE = 0.25

# Software paths
SEQGEN_PATH = "/n/holylfs05/LABS/grad_lab/Users/mkline/bridging_sims/software/Seq-Gen/source/seq-gen"
LSD2_PATH = "/n/holylfs05/LABS/grad_lab/Users/mkline/bridging_sims/software/lsd2/src/lsd2"

# Resource limits
MAX_MEM_MB = 64000
MAX_RUNTIME = 2880
MAX_CPUS = 16

# Find all parameter files and create simulation list
param_files = glob.glob(f"{PARAM_DIR}/params_*.json")
print(f"Found {len(param_files)} parameter files")

SIMULATIONS = []
for param_file in param_files:
    filename = Path(param_file).stem
    # Extract p_msmw_w and replicate from filename: params_0.25_042.json
    parts = filename.replace("params_", "").split("_")
    if len(parts) == 2:
        p_val, rep = parts[0], parts[1]
        sim_id = f"{p_val}_{rep}"
        SIMULATIONS.append({
            'id': sim_id,
            'p_msmw_w': float(p_val),
            'replicate': int(rep),
            'param_file': param_file
        })

# Create mapping for rules
SIM_TO_PARAM = {sim['id']: sim['param_file'] for sim in SIMULATIONS}
SIM_IDS = [sim['id'] for sim in SIMULATIONS]

print(f"Configured {len(SIMULATIONS)} simulations:")
p_counts = {}
for sim in SIMULATIONS:
    p_val = sim['p_msmw_w']
    p_counts[p_val] = p_counts.get(p_val, 0) + 1
for p_val, count in sorted(p_counts.items()):
    print(f"  p_msmw_w = {p_val}: {count} simulations")

# Resource scaling functions
def get_mem_mb_sim(wildcards, attempt):
    return min(attempt * 16000, 32000)

def get_mem_mb_tree_building(wildcards, attempt):
    return min(attempt * 24000, 64000)  

def get_mem_mb_phylo(wildcards, attempt):
    return min(attempt * 20000, 40000)

def get_runtime_sim(wildcards, attempt):
    return min(attempt * 180, 480)  # 3-8 hours

def get_runtime_phylo(wildcards, attempt):
    return min(attempt * 240, 720)  # 4-12 hours

def get_partition(wildcards, attempt):
    partitions = ["hsph", "sapphire", "intermediate", "shared"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

def get_quick_partition(wildcards, attempt):
    partitions = ["hsph", "sapphire", "shared", "intermediate"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

def get_long_partition(wildcards, attempt):
    partitions = ["intermediate", "hsph", "sapphire", "shared"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

# Phylo batch configuration functions (from your phylo pipeline)
def get_batch_config(dataset_type, downsample_rate=None):
    """Return batch configuration based on dataset characteristics"""
    if dataset_type == "full_sequences":
        return {'batch_size': 1, 'mem_mb': 128000, 'runtime': 2880, 'cpus_per_task': 8}
    elif dataset_type == "complete":
        return {'batch_size': 1, 'mem_mb': 32000, 'runtime': 720, 'cpus_per_task': 4}
    elif downsample_rate in [50, 25]:
        return {'batch_size': 3, 'mem_mb': 24000, 'runtime': 360, 'cpus_per_task': 2}
    elif downsample_rate in [10, 5]:
        return {'batch_size': 8, 'mem_mb': 12000, 'runtime': 180, 'cpus_per_task': 1}
    else:
        return {'batch_size': 20, 'mem_mb': 4000, 'runtime': 120, 'cpus_per_task': 1}

def get_tier_name(dataset, downsample_rate):
    """Get tier name for batching"""
    if dataset == 'complete':
        return 'tier2_complete'
    elif dataset == 'full':
        return 'tier1_full'
    elif downsample_rate in [50, 25]:
        return 'tier3_large'
    elif downsample_rate in [10, 5]:
        return 'tier4_medium'
    else:
        return 'tier5_small'

def get_batch_resources(wildcards, input, resource_type):
    """Extract resource requirements for this batch"""
    batch_df = pd.read_csv(input.batch_assignments)
    batch_data = batch_df[batch_df['batch_id'] == int(wildcards.batch_id)]
    if len(batch_data) > 0:
        value = batch_data.iloc[0][resource_type]
        return int(value)
    return 4000 if resource_type == 'mem_mb' else (60 if resource_type == 'runtime' else 1)

def get_actual_batch_ids(sim_id):
    """Get actual batch IDs for a simulation"""
    try:
        batch_file = f"{OUTPUT_BASE}/phylo/{sim_id}/batch_assignments.csv"
        if os.path.exists(batch_file):
            batch_df = pd.read_csv(batch_file)
            return sorted(batch_df['batch_id'].unique().tolist())
        return []
    except:
        return []

def get_batch_flags_for_aggregation(wildcards):
    """Get batch flags after checkpoint completes - for aggregation"""
    checkpoints.create_phylo_batches.get(sim_id=wildcards.sim_id)
    batch_file = f"{OUTPUT_BASE}/phylo/{wildcards.sim_id}/batch_assignments.csv"
    batch_df = pd.read_csv(batch_file)
    batch_ids = sorted(batch_df['batch_id'].unique())
    return [f"{OUTPUT_BASE}/phylo/{wildcards.sim_id}/flags/batch_{bid}_tree_bridging_30_complete.flag" 
            for bid in batch_ids]

rule all:
    input:
        # Final comprehensive comparison and summary
        f"{OUTPUT_BASE}/final_analysis/comprehensive_bridging_comparison.csv",
        f"{OUTPUT_BASE}/final_analysis/summary_statistics.csv"

rule run_simulation:
    input:
        param_file = lambda wildcards: SIM_TO_PARAM[wildcards.sim_id]
    output:
        simulation_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/simulation_complete.flag"
    conda: "phylo"
    resources:
        mem_mb=get_mem_mb_sim,
        runtime=get_runtime_sim,
        cpus_per_task=2,
        slurm_partition=get_partition
    shell:
        """
        echo "🚀 Running simulation for {wildcards.sim_id}"
        echo "Using parameter file: {input.param_file}"
        
        # Create directory and modify parameter file
        mkdir -p {OUTPUT_BASE}/simulations/{wildcards.sim_id}
        
        python -c "
import json
import os
import glob

# Load original parameters
with open('{input.param_file}', 'r') as f:
    params = json.load(f)

# Update output directory to our pipeline structure
params['output']['output_dir'] = '{OUTPUT_BASE}/simulations/{wildcards.sim_id}'
params['output']['run_name'] = 'sim'

# Save modified parameters
temp_param_file = '{OUTPUT_BASE}/simulations/{wildcards.sim_id}/modified_params.json'
with open(temp_param_file, 'w') as f:
    json.dump(params, f, indent=2)

print(f'Modified parameter file: {{temp_param_file}}')
"
        
        # Run simulation with modified parameters
        python ../scripts/10_16_core_sim.py {OUTPUT_BASE}/simulations/{wildcards.sim_id}/modified_params.json
        
        # Find actual simulation output directory - FIXED VERSION
        python -c "
import glob
import os

# Try both patterns: timestamped and non-timestamped
patterns = [
    '{OUTPUT_BASE}/simulations/{wildcards.sim_id}/sim_*',
    '{OUTPUT_BASE}/simulations/{wildcards.sim_id}/sim'
]

actual_dir = None
for pattern in patterns:
    dirs = glob.glob(pattern)
    if dirs:
        if len(dirs) == 1 and os.path.isdir(dirs[0]):
            actual_dir = dirs[0]
            break
        elif len(dirs) > 1:
            actual_dir = max(dirs, key=os.path.getctime)  # Most recent
            break

if actual_dir:
    print(f'Found actual simulation directory: {{actual_dir}}')
    
    # Check required outputs exist
    required_files = [
        'transmission_df.csv', 'parameters_used.json', 'nodes_df.csv',
        'initial_infectors.csv', 'edge_df.csv'
    ]
    
    missing_files = []
    for file in required_files:
        if not os.path.exists(os.path.join(actual_dir, file)):
            missing_files.append(file)
    
    if missing_files:
        print(f'ERROR: Missing files: {{missing_files}}')
        exit(1)
    else:
        print('✅ All required simulation outputs found')
        
        # Write actual directory to flag file
        with open('{output.simulation_flag}', 'w') as f:
            f.write(f'ACTUAL_SIM_DIR={{actual_dir}}\\n')
        print('✅ Simulation {wildcards.sim_id} completed successfully')
else:
    print('ERROR: No simulation directory found')
    print('Checked patterns:')
    for pattern in patterns:
        print(f'  {{pattern}}')
    exit(1)
"
        """

rule check_lineage:
    input:
        simulation_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/simulation_complete.flag"
    output:
        lineage_results = f"{OUTPUT_BASE}/simulations/{{sim_id}}/lineage_test_results.json",
        lineage_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/lineage_check_complete.flag"
    conda: "phylo"
    resources:
        mem_mb=4000,
        runtime=30,
        cpus_per_task=1,
        slurm_partition="shared"
    shell:
        """
        echo "🧬 Checking lineage for {wildcards.sim_id}"
        
        # Get actual simulation directory from flag
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.simulation_flag} | cut -d'=' -f2)
        echo "Using simulation directory: $ACTUAL_SIM_DIR"
        
        # Run lineage test
        python ../scripts/10_17_single_lineage_test_script.py "$ACTUAL_SIM_DIR"
        
        if [ -f "$ACTUAL_SIM_DIR/lineage_test_results.json" ]; then
            # Copy results to our pipeline location
            cp "$ACTUAL_SIM_DIR/lineage_test_results.json" {output.lineage_results}
            
            # Check results and create flag
            SINGLE_LINEAGE=$(python -c "import json; print(json.load(open('{output.lineage_results}'))['single_lineage'])")
            echo "SINGLE_LINEAGE=$SINGLE_LINEAGE" > {output.lineage_flag}
            echo "ACTUAL_SIM_DIR=$ACTUAL_SIM_DIR" >> {output.lineage_flag}
            
            echo "✅ Lineage check {wildcards.sim_id} completed: $SINGLE_LINEAGE"
        else
            echo "❌ Lineage check {wildcards.sim_id} failed"
            exit 1
        fi
        """

rule ground_truth_analysis:
    input:
        simulation_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/simulation_complete.flag"
    output:
        behavior_classifications = f"{OUTPUT_BASE}/simulations/{{sim_id}}/infection_behavior_classifications.csv",
        bridging_values = f"{OUTPUT_BASE}/simulations/{{sim_id}}/ground_truth_bridging.csv",
        bridging_counts = f"{OUTPUT_BASE}/simulations/{{sim_id}}/ground_truth_bridging_counts.csv",
        bridging_totals = f"{OUTPUT_BASE}/simulations/{{sim_id}}/ground_truth_bridging_totals.csv",
        gt_analysis_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/ground_truth_analysis_complete.flag"
    conda: "phylo"
    resources:
        mem_mb=12000,
        runtime=90,
        cpus_per_task=2,
        slurm_partition=get_partition
    shell:
        """
        echo "🏷️ Running ground truth analysis for {wildcards.sim_id}"
        
        # Get actual simulation directory from flag
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.simulation_flag} | cut -d'=' -f2)
        echo "Using simulation directory: $ACTUAL_SIM_DIR"
        
        # Run behavior classification
        echo "Running behavior classification..."
        python ../scripts/10_17_infection_behavior_classification_analysis.py "$ACTUAL_SIM_DIR" --both
        
        # Copy behavior classifications to pipeline location
        if [ -f "$ACTUAL_SIM_DIR/infection_behavior_classifications.csv" ]; then
            cp "$ACTUAL_SIM_DIR/infection_behavior_classifications.csv" {output.behavior_classifications}
            echo "✅ Behavior classification completed"
        else
            echo "❌ Behavior classification failed"
            exit 1
        fi
        
        # Run ground truth bridging analysis
        echo "Running ground truth bridging analysis..."
        python ../scripts/10_30_ground_truth_bridging.py \\
            "$ACTUAL_SIM_DIR" \\
            "{output.behavior_classifications}" \\
            --output "{output.bridging_values}" \\
            --detailed-output "{OUTPUT_BASE}/simulations/{wildcards.sim_id}/ground_truth_bridging_detailed.csv"
        
        # Check bridging outputs
        if [ -f {output.bridging_values} ] && [ -f {output.bridging_counts} ] && [ -f {output.bridging_totals} ]; then
            echo "ACTUAL_SIM_DIR=$ACTUAL_SIM_DIR" > {output.gt_analysis_flag}
            echo "✅ Ground truth analysis {wildcards.sim_id} completed"
        else
            echo "❌ Ground truth bridging analysis failed"
            exit 1
        fi
        """

rule build_transmission_tree:
    input:
        lineage_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/lineage_check_complete.flag"
    output:
        final_trimmed_tree = f"{OUTPUT_BASE}/simulations/{{sim_id}}/final_trimmed_tree.nwk",
        tree_building_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/tree_building_complete.flag"
    conda: "phylo"
    resources:
        mem_mb=get_mem_mb_tree_building,  # 24GB → 48GB → 64GB on retries
        runtime=180,  # 3 hours (more than basic sim's 2 hours for analysis)
        cpus_per_task=4,  # More CPUs too
        slurm_partition=get_partition
    shell:
        """
        echo "🌳 Building transmission tree for {wildcards.sim_id}"
        
        # Get actual simulation directory and check lineage
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.lineage_flag} | cut -d'=' -f2)
        SINGLE_LINEAGE=$(grep SINGLE_LINEAGE {input.lineage_flag} | cut -d'=' -f2)
        
        echo "Using simulation directory: $ACTUAL_SIM_DIR"
        echo "Single lineage status: $SINGLE_LINEAGE"
        
        if [ "$SINGLE_LINEAGE" = "True" ]; then
            echo "✅ Proceeding with tree building (single lineage confirmed)"
            
            # Step 1: Build hybrid transmission tree
            echo "Building hybrid transmission tree..."
            python ../scripts/10_21_build_hybrid_transmission_tree.py "$ACTUAL_SIM_DIR"
            
            # Check if tree was created
            TREE_FILE=$(ls "$ACTUAL_SIM_DIR"/tree_*_hybrid_episodes*.newick 2>/dev/null | head -1)
            if [ -z "$TREE_FILE" ]; then
                echo "❌ No hybrid tree file found"
                exit 1
            fi
            echo "Found tree file: $TREE_FILE"
            
            # Step 2: First trimming
            echo "Running first tree trimming..."
            Rscript ../scripts/10_20_tree_trimmer.R "$TREE_FILE" "$ACTUAL_SIM_DIR" "$ACTUAL_SIM_DIR/trimmed_tree.nwk"
            
            if [ ! -f "$ACTUAL_SIM_DIR/trimmed_tree.nwk" ]; then
                echo "❌ First trimming failed"
                exit 1
            fi
            
            # Step 3: Additional trimming
            echo "Running additional tree trimming..."
            Rscript ../scripts/10_21_tree_additional_trimming.R \\
                "$ACTUAL_SIM_DIR/trimmed_tree.nwk" \\
                "$ACTUAL_SIM_DIR/transmission_df.csv" \\
                "$ACTUAL_SIM_DIR/final_trimmed_tree.nwk"
            
            if [ ! -f "$ACTUAL_SIM_DIR/final_trimmed_tree.nwk" ]; then
                echo "❌ Additional trimming failed"
                exit 1
            fi
            
            # Copy final tree to pipeline location
            cp "$ACTUAL_SIM_DIR/final_trimmed_tree.nwk" {output.final_trimmed_tree}
            
            # Create completion flag
            echo "ACTUAL_SIM_DIR=$ACTUAL_SIM_DIR" > {output.tree_building_flag}
            echo "TREE_READY=True" >> {output.tree_building_flag}
            
            echo "✅ Tree building and trimming completed for {wildcards.sim_id}"
        else
            echo "❌ Skipping tree building - multiple lineages detected"
            exit 1
        fi
        """

rule run_seqgen:
    input:
        tree_building_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/tree_building_complete.flag"
    output:
        sequences = f"{OUTPUT_BASE}/phylo/{{sim_id}}/sequences/sequences.fasta",
        seqgen_log = f"{OUTPUT_BASE}/phylo/{{sim_id}}/sequences/sequences.log",
        seqgen_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/seqgen_complete.flag"
    resources:
        mem_mb = 2000,
        runtime = 30,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    shell:
        """
        echo "🧬 Starting seq-gen for {wildcards.sim_id} at $(date)"
        
        # Get actual simulation directory and tree file
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.tree_building_flag} | cut -d'=' -f2)
        TREE_READY=$(grep TREE_READY {input.tree_building_flag} | cut -d'=' -f2)
        
        if [ "$TREE_READY" != "True" ]; then
            echo "❌ Tree not ready for {wildcards.sim_id}"
            exit 1
        fi
        
        echo "Using final trimmed tree from: $ACTUAL_SIM_DIR/final_trimmed_tree.nwk"
        echo "Starting seq-gen at $(date)" > {output.seqgen_log}
        
        mkdir -p $(dirname {output.sequences})
        
        {SEQGEN_PATH} \\
            -m GTR -r 1.00578 4.59475 0.74535 0.46171 4.55284 1.00000 \\
            -f 0.154 0.346 0.343 0.157 -a 0.45 -l 10000 -n 1 -s 2.254247e-06 -z {PIPELINE_SEED} \\
            < "$ACTUAL_SIM_DIR/final_trimmed_tree.nwk" > {output.sequences} 2>> {output.seqgen_log}
        
        if [ ! -s {output.sequences} ]; then
            echo "ERROR: seq-gen output is empty!" >> {output.seqgen_log}
            exit 1
        fi
        
        SEQ_COUNT=$(grep -c '^>' {output.sequences})
        echo "Generated $SEQ_COUNT sequences" >> {output.seqgen_log}
        echo "Completed seq-gen at $(date)" >> {output.seqgen_log}
        
        # Create completion flag
        echo "SEQUENCES_READY=True" > {output.seqgen_flag}
        echo "ACTUAL_SIM_DIR=$ACTUAL_SIM_DIR" >> {output.seqgen_flag}
        
        echo "✅ SeqGen completed for {wildcards.sim_id}"
        """
rule run_filtering_downsampling:
    input:
        seqgen_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/seqgen_complete.flag"
    output:
        filtering_summary = f"{OUTPUT_BASE}/phylo/{{sim_id}}/filtering/filtering_summary.csv",
        full_sequences = f"{OUTPUT_BASE}/phylo/{{sim_id}}/filtering/full_sequences.fasta",
        filtering_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/filtering_complete.flag"
    params:
        output_dir = f"{OUTPUT_BASE}/phylo/{{sim_id}}/filtering",
        rates_str = " ".join(map(str, DOWNSAMPLE_RATES))
    resources:
        mem_mb = 16000,
        runtime = 60,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    conda: "phylo"
    shell:
        """
        echo "🔍 Starting filtering and downsampling for {wildcards.sim_id}"
        
        # Get paths from flags
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.seqgen_flag} | cut -d'=' -f2)
        SEQUENCES_READY=$(grep SEQUENCES_READY {input.seqgen_flag} | cut -d'=' -f2)
        
        if [ "$SEQUENCES_READY" != "True" ]; then
            echo "❌ Sequences not ready for {wildcards.sim_id}"
            exit 1
        fi
        
        echo "Using sequences: {OUTPUT_BASE}/phylo/{wildcards.sim_id}/sequences/sequences.fasta"
        echo "Using simulation data from: $ACTUAL_SIM_DIR"
        echo "Using seed: {PIPELINE_SEED}"
        
        mkdir -p {params.output_dir}
        
        # Copy full sequences
        cp {OUTPUT_BASE}/phylo/{wildcards.sim_id}/sequences/sequences.fasta {output.full_sequences}
        
        # Run filtering and downsampling
        python ../scripts/10_26_infection_filtering_downsampling.py \\
            "$ACTUAL_SIM_DIR" \\
            {OUTPUT_BASE}/phylo/{wildcards.sim_id}/sequences/sequences.fasta \\
            {params.output_dir} \\
            --downsample-rates {params.rates_str} \\
            --n-replicates {N_REPLICATES} \\
            --high-activity-subsample {HIGH_ACTIVITY_SUBSAMPLE} \\
            --seed {PIPELINE_SEED}
        
        if [ ! -f {output.filtering_summary} ]; then
            echo "❌ Filtering summary not created for {wildcards.sim_id}"
            exit 1
        fi
        
        DATASET_COUNT=$(tail -n +2 {output.filtering_summary} | wc -l)
        echo "✅ Created $DATASET_COUNT filtered datasets for {wildcards.sim_id}"
        
        # Create completion flag
        echo "FILTERING_READY=True" > {output.filtering_flag}
        echo "ACTUAL_SIM_DIR=$ACTUAL_SIM_DIR" >> {output.filtering_flag}
        echo "DATASET_COUNT=$DATASET_COUNT" >> {output.filtering_flag}
        
        echo "✅ Filtering and downsampling completed for {wildcards.sim_id}"
        """
rule run_sampled_bridging_analysis:
    input:
        filtering_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/filtering_complete.flag"
    output:
        sampled_bridging_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/sampled_bridging_complete.flag"
    resources:
        mem_mb = 12000,
        runtime = 90,
        cpus_per_task = 2,
        slurm_partition = get_partition
    conda: "phylo"
    shell:
        """
        echo "📊 Running sampled bridging analysis for {wildcards.sim_id}"
        
        # Get paths from flag
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.filtering_flag} | cut -d'=' -f2)
        FILTERING_READY=$(grep FILTERING_READY {input.filtering_flag} | cut -d'=' -f2)
        
        if [ "$FILTERING_READY" != "True" ]; then
            echo "❌ Filtering not ready for {wildcards.sim_id}"
            exit 1
        fi
        
        echo "Using simulation directory: $ACTUAL_SIM_DIR"
        
        # Read filtering summary to find all datasets (excluding full sequences)
        python -c "
import pandas as pd
import os

summary_file = '{OUTPUT_BASE}/phylo/{wildcards.sim_id}/filtering/filtering_summary.csv'
summary_df = pd.read_csv(summary_file)

print(f'Processing {{len(summary_df)}} filtered datasets for sampled bridging...')

success_count = 0
total_count = 0

for _, row in summary_df.iterrows():
    dataset_id = f\"{{row['population']}}_{{row['dataset']}}\"
    infections_file = row['infection_file']
    dataset_dir = os.path.dirname(row['fasta_file'])
    
    total_count += 1
    print(f'Running sampled bridging for {{dataset_id}}...', flush=True)
    
    # Run 10_30_sampled_bridging.py
    bridging_output = os.path.join(dataset_dir, 'sampled_bridging_30.csv')
    bridging_detailed = os.path.join(dataset_dir, 'sampled_bridging_30_detailed.csv')
    
    cmd = 'python ../scripts/10_30_sampled_bridging.py "' + infections_file + '" "$ACTUAL_SIM_DIR" --output "' + bridging_output + '" --detailed-output "' + bridging_detailed + '"'

    result = os.system(cmd)
    
    if result == 0:
        print(f'SUCCESS: Sampled bridging completed for {{dataset_id}}')
        success_count += 1
    else:
        print(f'ERROR: Sampled bridging failed for {{dataset_id}}')

print(f'Sampled bridging completed: {{success_count}}/{{total_count}} datasets successful')
"
        
        echo "SAMPLED_BRIDGING_READY=True" > {output.sampled_bridging_flag}
        echo "ACTUAL_SIM_DIR=$ACTUAL_SIM_DIR" >> {output.sampled_bridging_flag}
        
        echo "✅ Sampled bridging analysis completed for {wildcards.sim_id}"
        """

checkpoint create_phylo_batches:
    input:
        filtering_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/filtering_complete.flag"
    output:
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv",
        batching_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batching_complete.flag"
    run:
        import pandas as pd
        import os

        # Read ACTUAL_SIM_DIR from filtering flag
        actual_sim_dir = None
        with open(input.filtering_flag, 'r') as f:
            for line in f:
                if line.startswith('ACTUAL_SIM_DIR='):
                    actual_sim_dir = line.strip().split('=', 1)[1]
                    break
        
        if not actual_sim_dir:
            raise ValueError(f"ACTUAL_SIM_DIR not found in {input.filtering_flag}")

        # Load filtering summary
        summary_file = f"{OUTPUT_BASE}/phylo/{wildcards.sim_id}/filtering/filtering_summary.csv"
        summary_df = pd.read_csv(summary_file)

        batch_assignments = []

        # Process filtered datasets
        for _, row in summary_df.iterrows():
            population = row['population']
            dataset = row['dataset']
            fasta_file = row['fasta_file']

            # Determine dataset type and batch config
            if dataset == 'complete':
                dataset_type = 'complete'
                downsample_rate = None
                batch_config = get_batch_config('complete')
            else:
                dataset_type = 'downsampled'
                downsample_rate = int(dataset.split('pct')[0])
                batch_config = get_batch_config('downsampled', downsample_rate)

            batch_assignments.append({
                'dataset_id': f"{population}_{dataset}",
                'population': population,
                'dataset': dataset,
                'fasta_file': fasta_file,
                'infections_file': row['infection_file'],
                'batch_tier': get_tier_name(dataset, downsample_rate),
                'batch_size': batch_config['batch_size'],
                'mem_mb': batch_config['mem_mb'],
                'runtime': batch_config['runtime'],
                'cpus_per_task': batch_config['cpus_per_task']
            })

        # Add full sequences as special case
        full_seq_row = {
            'dataset_id': 'full_sequences',
            'population': 'all',
            'dataset': 'full',
            'fasta_file': f"{OUTPUT_BASE}/phylo/{wildcards.sim_id}/filtering/full_sequences.fasta",
            'infections_file': None,
            'batch_tier': 'tier1_full',
            'batch_size': 1,
            'mem_mb': 128000,
            'runtime': 2880,
            'cpus_per_task': 8
        }

        batch_df = pd.DataFrame(batch_assignments)
        batch_df = pd.concat([pd.DataFrame([full_seq_row]), batch_df], ignore_index=True)

        # Assign actual batch numbers within tiers (DYNAMIC)
        batch_df['batch_id'] = 0
        current_batch = 1

        for tier in batch_df['batch_tier'].unique():
            tier_data = batch_df[batch_df['batch_tier'] == tier]
            batch_size = tier_data.iloc[0]['batch_size']

            for i in range(0, len(tier_data), batch_size):
                batch_df.loc[tier_data.index[i:i+batch_size], 'batch_id'] = current_batch
                current_batch += 1

        # CRITICAL: Calculate actual max batch ID dynamically
        max_batch_id = batch_df['batch_id'].max()
        
        batch_df.to_csv(output.batch_assignments, index=False)
        
        # Create flag with DYNAMIC batch count AND ACTUAL_SIM_DIR
        with open(output.batching_flag, 'w') as f:
            f.write("BATCHING_READY=True\n")
            f.write(f"MAX_BATCH_ID={max_batch_id}\n")
            f.write(f"NUM_BATCHES={max_batch_id}\n")
            f.write(f"NUM_DATASETS={len(batch_df)}\n")
            f.write(f"ACTUAL_SIM_DIR={actual_sim_dir}\n")  # ← ADDED THIS LINE
        
        print(f"✅ Created {max_batch_id} batches for {wildcards.sim_id} phylogenetic analysis")
        print(f"   Using simulation directory: {actual_sim_dir}")
        print(f"   Datasets: {len(batch_df)} total")
        print(f"   Batch distribution:")
        for tier in sorted(batch_df['batch_tier'].unique()):
            tier_batches = batch_df[batch_df['batch_tier'] == tier]['batch_id'].nunique()
            tier_datasets = len(batch_df[batch_df['batch_tier'] == tier])
            print(f"     {tier}: {tier_datasets} datasets in {tier_batches} batches")

rule run_parallel_batch:
    input:
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv",
        batching_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batching_complete.flag"
    output:
        iqtree_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_iqtree_complete.flag",
        dates_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_dates_complete.flag",
        behaviors_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_behaviors_complete.flag"
    resources:
        mem_mb = lambda wildcards, input: max(get_batch_resources(wildcards, input, "mem_mb"), 8000),
        runtime = lambda wildcards, input: max(get_batch_resources(wildcards, input, "runtime"), 120),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_partition
    conda: "phylo"
    shell:
        """
        echo "🔄 Starting parallel batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"
        
        # Get actual simulation directory
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.batching_flag} | cut -d'=' -f2)
        
        mkdir -p {OUTPUT_BASE}/phylo/{wildcards.sim_id}/flags

        python -c "
import pandas as pd
import os
import json
from concurrent.futures import ThreadPoolExecutor

def process_dataset(row, sim_dir, sim_id):
    dataset_id = row['dataset_id']
    fasta_file = row['fasta_file']
    infections_file = row['infections_file']

    print(f'Processing {{dataset_id}} for {{sim_id}}...', flush=True)

    # Determine infections file to use
    if dataset_id == 'full_sequences':
        # FIXED: Filter transmission_df to post-burn-in, non-superseded, sampled infections
        print(f'  Applying post-burn-in filtering for full_sequences...', flush=True)
        
        # Load full transmission_df
        full_transmission_df = pd.read_csv(os.path.join(sim_dir, 'transmission_df.csv'))
        print(f'  Loaded {{len(full_transmission_df)}} total transmissions', flush=True)
        
        # Get burn-in cutoff from parameters
        try:
            with open(os.path.join(sim_dir, 'parameters_used.json'), 'r') as f:
                params = json.load(f)
            cutoff_day = (params['simulation']['partnership_burnin_days'] + 
                         params['simulation']['transmission_burnin_days'])
            print(f'  Using burn-in cutoff: {{cutoff_day}}', flush=True)
        except Exception as e:
            print(f'  Warning: Could not read parameters, using default cutoff 12000: {{e}}', flush=True)
            cutoff_day = 12000  # fallback value
        
        # Apply same filtering as tree building
        filtered_transmissions = full_transmission_df[
            (full_transmission_df['superseded_simultaneous'] == False) &
            (full_transmission_df['day_of_transmission'] >= cutoff_day) &
            (pd.notna(full_transmission_df['day_of_sampling']))
        ]
        
        print(f'  After filtering: {{len(filtered_transmissions)}} post-burn-in, non-superseded, sampled transmissions', flush=True)
        
        # Save filtered data to temporary file for extraction scripts
        dataset_dir = os.path.dirname(fasta_file)
        temp_infections_file = os.path.join(dataset_dir, 'temp_filtered_transmissions.csv')
        filtered_transmissions.to_csv(temp_infections_file, index=False)
        
        infections_to_use = temp_infections_file
    else:
        infections_to_use = infections_file
        
    dataset_dir = os.path.dirname(fasta_file)

    jobs = []

    # 1. IQ-TREE
    if os.path.exists(fasta_file):
        cmd1 = 'iqtree -s ' + fasta_file + ' -m GTR+G --fast -seed {PIPELINE_SEED}'
        jobs.append(('iqtree', cmd1))

    # 2-3. Dates and behaviors (for TreeTime - needed for both full and filtered)
    if os.path.exists(infections_to_use):
        dates_output = os.path.join(dataset_dir, 'dates.txt')
        cmd2 = 'python ../scripts/10_27_extract_tip_dates.py "' + infections_to_use + '" "' + fasta_file + '" "' + sim_dir + '" "' + dates_output + '"'
        jobs.append(('dates', cmd2))

        behaviors_prefix = os.path.join(dataset_dir, 'behaviors')
        cmd3 = 'python ../scripts/10_27_extract_tip_states.py "' + infections_to_use + '" "' + fasta_file + '" "' + sim_dir + '" "' + behaviors_prefix + '"'
        jobs.append(('behaviors', cmd3))
    else:
        print(f'  Warning: Infections file not found: {{infections_to_use}}', flush=True)

    # Run jobs in parallel
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [(job_type, executor.submit(os.system, cmd)) for job_type, cmd in jobs]

    for job_type, future in futures:
        result = future.result()
        if result == 0:
            print(f'SUCCESS: {{job_type}} completed for {{dataset_id}}')
        else:
            print(f'ERROR: {{job_type}} failed for {{dataset_id}}')
    
    # Clean up temporary file if created
    if dataset_id == 'full_sequences' and 'temp_filtered_transmissions.csv' in infections_to_use:
        try:
            os.remove(infections_to_use)
            print(f'  Cleaned up temporary file: {{infections_to_use}}', flush=True)
        except:
            pass

batch_df = pd.read_csv('{input.batch_assignments}')
batch_data = batch_df[batch_df['batch_id'] == {wildcards.batch_id}]

for _, row in batch_data.iterrows():
    process_dataset(row, '$ACTUAL_SIM_DIR', '{wildcards.sim_id}')
"

        echo "✅ Parallel batch {wildcards.batch_id} completed at $(date)"
        touch {output.iqtree_flag} {output.dates_flag} {output.behaviors_flag}
        """

rule run_lsd2_batch:
    input:
        iqtree_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_iqtree_complete.flag",
        dates_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_dates_complete.flag",
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv"
    output:
        lsd2_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_lsd2_complete.flag"
    resources:
        mem_mb = lambda wildcards, input: max(get_batch_resources(wildcards, input, "mem_mb"), 8000),
        runtime = lambda wildcards, input: max(get_batch_resources(wildcards, input, "runtime"), 120),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_partition
    conda: "phylo"
    shell:
        """
        echo "🕐 Starting LSD2 batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"

        python -c "
import pandas as pd
import os

batch_df = pd.read_csv('{input.batch_assignments}')
batch_data = batch_df[batch_df['batch_id'] == {wildcards.batch_id}]

for _, row in batch_data.iterrows():
    dataset_id = row['dataset_id']
    fasta_file = row['fasta_file']

    dataset_dir = os.path.dirname(fasta_file)
    tree_file = fasta_file + '.treefile'
    dates_file = os.path.join(dataset_dir, 'dates.txt')

    if os.path.exists(tree_file) and os.path.exists(dates_file):
        print(f'Running LSD2 for {{dataset_id}}...', flush=True)

        cmd = '{LSD2_PATH} -i "' + tree_file + '" -d "' + dates_file + '" -s 10000 -r a -l 0'
        result = os.system(cmd)

        if result == 0:
            print(f'SUCCESS: LSD2 completed for {{dataset_id}}')
        else:
            print(f'ERROR: LSD2 failed for {{dataset_id}}')
    else:
        print(f'Missing files for {{dataset_id}} - tree: {{os.path.exists(tree_file)}}, dates: {{os.path.exists(dates_file)}}')
"

        echo "✅ LSD2 batch {wildcards.batch_id} completed at $(date)"
        touch {output.lsd2_flag}
        """

rule run_treetime_msmw_msm_batch:
    input:
        lsd2_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_lsd2_complete.flag",
        behaviors_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_behaviors_complete.flag",
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv"
    output:
        treetime_msmw_msm_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_treetime_msmw_msm_complete.flag"
    resources:
        mem_mb = lambda wildcards, input: max(get_batch_resources(wildcards, input, "mem_mb"), 16000),
        runtime = lambda wildcards, input: max(get_batch_resources(wildcards, input, "runtime"), 240),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_long_partition
    conda: "phylo"
    shell:
        """
        echo "🌳 Starting TreeTime MSMW+MSM batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"

        python -c "
import pandas as pd
import os

batch_df = pd.read_csv('{input.batch_assignments}')
batch_data = batch_df[batch_df['batch_id'] == {wildcards.batch_id}]

for _, row in batch_data.iterrows():
    dataset_id = row['dataset_id']
    fasta_file = row['fasta_file']

    dataset_dir = os.path.dirname(fasta_file)
    lsd2_tree = fasta_file + '.treefile.result.nwk'
    states_file = os.path.join(dataset_dir, 'behaviors_msmw_msm_states.csv')
    outdir = os.path.join(dataset_dir, 'msmw_msm')

    if os.path.exists(lsd2_tree) and os.path.exists(states_file):
        print(f'Running TreeTime MSMW+MSM for {{dataset_id}}...', flush=True)

        os.makedirs(outdir, exist_ok=True)
        cmd = 'treetime mugration --tree "' + lsd2_tree + '" --states "' + states_file + '" --attribute state --outdir "' + outdir + '" --verbose 2 --confidence'
        result = os.system(cmd)

        if result == 0:
            print(f'SUCCESS: TreeTime MSMW+MSM completed for {{dataset_id}}')
        else:
            print(f'ERROR: TreeTime MSMW+MSM failed for {{dataset_id}}')
    else:
        print(f'Missing files for {{dataset_id}} - lsd2_tree: {{os.path.exists(lsd2_tree)}}, states: {{os.path.exists(states_file)}}')
"

        echo "✅ TreeTime MSMW+MSM batch {wildcards.batch_id} completed at $(date)"
        touch {output.treetime_msmw_msm_flag}
        """

rule run_treetime_msmw_msw_batch:
    input:
        lsd2_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_lsd2_complete.flag",
        behaviors_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_behaviors_complete.flag",
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv"
    output:
        treetime_msmw_msw_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_treetime_msmw_msw_complete.flag"
    resources:
        mem_mb = lambda wildcards, input: max(get_batch_resources(wildcards, input, "mem_mb"), 16000),
        runtime = lambda wildcards, input: max(get_batch_resources(wildcards, input, "runtime"), 240),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_long_partition
    conda: "phylo"
    shell:
        """
        echo "🌳 Starting TreeTime MSMW+MSW batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"

        python -c "
import pandas as pd
import os

batch_df = pd.read_csv('{input.batch_assignments}')
batch_data = batch_df[batch_df['batch_id'] == {wildcards.batch_id}]

for _, row in batch_data.iterrows():
    dataset_id = row['dataset_id']
    fasta_file = row['fasta_file']

    dataset_dir = os.path.dirname(fasta_file)
    lsd2_tree = fasta_file + '.treefile.result.nwk'
    states_file = os.path.join(dataset_dir, 'behaviors_msmw_msw_states.csv')
    outdir = os.path.join(dataset_dir, 'msmw_msw')

    if os.path.exists(lsd2_tree) and os.path.exists(states_file):
        print(f'Running TreeTime MSMW+MSW for {{dataset_id}}...', flush=True)

        os.makedirs(outdir, exist_ok=True)
        cmd = 'treetime mugration --tree "' + lsd2_tree + '" --states "' + states_file + '" --attribute state --outdir "' + outdir + '" --verbose 2 --confidence'
        result = os.system(cmd)

        if result == 0:
            print(f'SUCCESS: TreeTime MSMW+MSW completed for {{dataset_id}}')
        else:
            print(f'ERROR: TreeTime MSMW+MSW failed for {{dataset_id}}')
    else:
        print(f'Missing files for {{dataset_id}} - lsd2_tree: {{os.path.exists(lsd2_tree)}}, states: {{os.path.exists(states_file)}}')
"

        echo "✅ TreeTime MSMW+MSW batch {wildcards.batch_id} completed at $(date)"
        touch {output.treetime_msmw_msw_flag}
        """

rule run_tree_bridging_30_batch:
    input:
        treetime_msmw_msm_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_treetime_msmw_msm_complete.flag",
        treetime_msmw_msw_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_treetime_msmw_msw_complete.flag",
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv"
    output:
        tree_bridging_30_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_tree_bridging_30_complete.flag"
    resources:
        mem_mb = lambda wildcards, input: max(get_batch_resources(wildcards, input, "mem_mb"), 16000),
        runtime = lambda wildcards, input: max(get_batch_resources(wildcards, input, "runtime"), 120),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_partition
    conda: "phylo"
    shell:
        """
        echo "🌲 Starting TreeTime bridging analysis 30 batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"

        python -c "
import pandas as pd
import os

batch_df = pd.read_csv('{input.batch_assignments}')
batch_data = batch_df[batch_df['batch_id'] == {wildcards.batch_id}]

for _, row in batch_data.iterrows():
    dataset_id = row['dataset_id']
    fasta_file = row['fasta_file']

    dataset_dir = os.path.dirname(fasta_file)
    msmw_msm_nexus = os.path.join(dataset_dir, 'msmw_msm', 'annotated_tree.nexus')
    msmw_msw_nexus = os.path.join(dataset_dir, 'msmw_msw', 'annotated_tree.nexus')
    lsd2_result = fasta_file + '.treefile.result'
    bridging_prefix = os.path.join(dataset_dir, 'tree_bridging_30')

    if os.path.exists(msmw_msm_nexus) and os.path.exists(msmw_msw_nexus) and os.path.exists(lsd2_result):
        print(f'Running TreeTime bridging analysis 30 for {{dataset_id}}...', flush=True)

        cmd = 'python ../scripts/10_30_treetime_bridging.py "' + msmw_msm_nexus + '" "' + msmw_msw_nexus + '" "' + lsd2_result + '" "' + bridging_prefix + '"'
        result = os.system(cmd)

        if result == 0:
            print(f'SUCCESS: TreeTime bridging analysis 30 completed for {{dataset_id}}')
            
            # Verify output files were created
            expected_files = [
                bridging_prefix + '_treetime_bridging.csv',
                bridging_prefix + '_treetime_bridging_counts.csv',
                bridging_prefix + '_treetime_bridging_totals.csv'
            ]
            
            missing_files = [f for f in expected_files if not os.path.exists(f)]
            if missing_files:
                print(f'WARNING: Missing output files for {{dataset_id}}: {{missing_files}}')
            else:
                print(f'All TreeTime bridging outputs verified for {{dataset_id}}')
        else:
            print(f'ERROR: TreeTime bridging analysis 30 failed for {{dataset_id}}')
    else:
        missing = []
        if not os.path.exists(msmw_msm_nexus): missing.append('msmw_msm_nexus')
        if not os.path.exists(msmw_msw_nexus): missing.append('msmw_msw_nexus') 
        if not os.path.exists(lsd2_result): missing.append('lsd2_result')
        print(f'Missing files for {{dataset_id}}: {{missing}}')
"

        echo "✅ TreeTime bridging analysis 30 batch {wildcards.batch_id} completed at $(date)"
        touch {output.tree_bridging_30_flag}
        """

rule aggregate_phylogenetic_bridging_30:
    input:
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv",
        batching_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batching_complete.flag",
        batch_flags = get_batch_flags_for_aggregation
    output:
        aggregated_phylo_bridging = f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_phylogenetic_bridging_30_results.csv"
    resources:
        mem_mb = 8000,
        runtime = 60,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    run:
        import pandas as pd
        import os

        print(f"🌳 Aggregating phylogenetic bridging 30 results for {wildcards.sim_id}...")

        # Read batch assignments to get actual batch IDs
        batch_df = pd.read_csv(input.batch_assignments)
        actual_batch_ids = sorted(batch_df['batch_id'].unique())
        
        print(f"Expected batches: {actual_batch_ids}")
        
        # Check that all batch flags exist
        missing_flags = []
        for batch_id in actual_batch_ids:
            flag_file = f"{OUTPUT_BASE}/phylo/{wildcards.sim_id}/flags/batch_{batch_id}_tree_bridging_30_complete.flag"
            if not os.path.exists(flag_file):
                missing_flags.append(f"batch_{batch_id}")
        
        if missing_flags:
            print(f"❌ ERROR: Missing batch completion flags: {missing_flags}")
            # Create empty output
            pd.DataFrame().to_csv(output.aggregated_phylo_bridging, index=False)
            return

        # Proceed with aggregation
        all_results = []

        for _, row in batch_df.iterrows():
            dataset_id = row['dataset_id']
            population = row['population']
            dataset = row['dataset']
            fasta_file = row['fasta_file']

            dataset_dir = os.path.dirname(fasta_file)
            bridging_summary_file = os.path.join(dataset_dir, 'tree_bridging_30_treetime_bridging.csv')

            if os.path.exists(bridging_summary_file):
                try:
                    bridging_df = pd.read_csv(bridging_summary_file, index_col=0)
                    bridging_df = bridging_df.reset_index()  # Convert pivot to long format

                    for _, bridging_row in bridging_df.iterrows():
                        result_row = {
                            'sim_id': wildcards.sim_id,
                            'dataset_id': dataset_id,
                            'population': population,
                            'dataset_type': 'complete' if dataset == 'complete' else 'downsampled',
                            'dataset': dataset,
                            'measure': bridging_row['measure']
                        }

                        if dataset == 'complete':
                            result_row['downsample_rate'] = '100'
                            result_row['replicate'] = 'all'
                        elif dataset == 'full':
                            result_row['downsample_rate'] = '100'
                            result_row['replicate'] = 'full_sequences'
                        else:
                            parts = dataset.split('_')
                            result_row['downsample_rate'] = parts[0].replace('pct', '')
                            result_row['replicate'] = parts[1]

                        for col in bridging_df.columns:
                            if col != 'measure':
                                result_row[f'phylo_30_{col}'] = bridging_row[col]

                        all_results.append(result_row)

                    print(f"✅ Collected phylo bridging results from {dataset_id}")
                except Exception as e:
                    print(f"❌ Error processing {bridging_summary_file}: {e}")

        if all_results:
            aggregated_df = pd.DataFrame(all_results)
            column_order = [
                'sim_id', 'dataset_id', 'population', 'dataset_type', 'dataset',
                'downsample_rate', 'replicate', 'measure'
            ] + [col for col in aggregated_df.columns if col.startswith('phylo_30_')]

            aggregated_df = aggregated_df[column_order]
            aggregated_df = aggregated_df.sort_values(['population', 'downsample_rate', 'replicate', 'measure'])
            aggregated_df.to_csv(output.aggregated_phylo_bridging, index=False)

            print(f"🎉 Phylogenetic bridging 30 aggregation complete! ({len(aggregated_df)} records)")
        else:
            pd.DataFrame().to_csv(output.aggregated_phylo_bridging, index=False)
            print("❌ No phylogenetic bridging results found")

rule aggregate_sampled_bridging_30:
    input:
        sampled_bridging_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/sampled_bridging_complete.flag",
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv"
    output:
        aggregated_sampled_bridging = f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_sampled_bridging_30_results.csv"
    resources:
        mem_mb = 4000,
        runtime = 30,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    run:
        import pandas as pd
        import os

        print(f"📊 Aggregating sampled bridging 30 results for {wildcards.sim_id}...")

        batch_df = pd.read_csv(input.batch_assignments)
        all_results = []

        for _, row in batch_df.iterrows():
            dataset_id = row['dataset_id']
            population = row['population']
            dataset = row['dataset']
            fasta_file = row['fasta_file']

            # Skip full sequences (we don't run sampled bridging on them)
            if dataset == 'full':
                continue

            dataset_dir = os.path.dirname(fasta_file)
            sampled_bridging_file = os.path.join(dataset_dir, 'sampled_bridging_30.csv')

            if os.path.exists(sampled_bridging_file):
                try:
                    sampled_df = pd.read_csv(sampled_bridging_file, index_col=0)
                    sampled_df = sampled_df.reset_index()  # Convert pivot back to long format

                    for _, sampled_row in sampled_df.iterrows():
                        result_row = {
                            'sim_id': wildcards.sim_id,
                            'dataset_id': dataset_id,
                            'population': population,
                            'dataset_type': 'complete' if dataset == 'complete' else 'downsampled',
                            'dataset': dataset,
                            'measure': sampled_row['measure']
                        }

                        if dataset == 'complete':
                            result_row['downsample_rate'] = '100'
                            result_row['replicate'] = 'all'
                        else:
                            parts = dataset.split('_')
                            result_row['downsample_rate'] = parts[0].replace('pct', '')
                            result_row['replicate'] = parts[1]

                        for col in sampled_df.columns:
                            if col != 'measure':
                                result_row[f'sampled_30_{col}'] = sampled_row[col]

                        all_results.append(result_row)

                    print(f"✅ Collected sampled bridging from {dataset_id}")
                except Exception as e:
                    print(f"❌ Error processing {sampled_bridging_file}: {e}")
            else:
                print(f"⚠️ Sampled bridging file not found for {dataset_id}: {sampled_bridging_file}")

        if all_results:
            aggregated_df = pd.DataFrame(all_results)

            column_order = [
                'sim_id', 'dataset_id', 'population', 'dataset_type', 'dataset',
                'downsample_rate', 'replicate', 'measure'
            ] + [col for col in aggregated_df.columns if col.startswith('sampled_30_')]

            aggregated_df = aggregated_df[column_order]
            aggregated_df = aggregated_df.sort_values(['population', 'downsample_rate', 'replicate'])
            aggregated_df.to_csv(output.aggregated_sampled_bridging, index=False)

            print(f"🎉 Sampled bridging 30 aggregation complete! ({len(aggregated_df)} records)")
        else:
            pd.DataFrame().to_csv(output.aggregated_sampled_bridging, index=False)
            print("❌ No sampled bridging results found")

rule comprehensive_aggregation:
    input:
        # Ground truth results from all simulations
        ground_truth_results = expand(f"{OUTPUT_BASE}/simulations/{{sim_id}}/ground_truth_bridging.csv", sim_id=SIM_IDS),
        # Phylogenetic results from all simulations  
        phylo_results = expand(f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_phylogenetic_bridging_30_results.csv", sim_id=SIM_IDS),
        # Sampled results from all simulations
        sampled_results = expand(f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_sampled_bridging_30_results.csv", sim_id=SIM_IDS)
    output:
        comprehensive_comparison = f"{OUTPUT_BASE}/final_analysis/comprehensive_bridging_comparison.csv",
        summary_statistics = f"{OUTPUT_BASE}/final_analysis/summary_statistics.csv"
    resources:
        mem_mb = 16000,
        runtime = 120,
        cpus_per_task = 2,
        slurm_partition = get_partition
    run:
        import pandas as pd
        import numpy as np
        import os
        from pathlib import Path

        print("🎯 Starting comprehensive cross-simulation aggregation...")
        
        # Create output directory
        os.makedirs(f"{OUTPUT_BASE}/final_analysis", exist_ok=True)
        
        all_comparisons = []
        sim_success_stats = {'total': len(SIM_IDS), 'ground_truth': 0, 'phylo': 0, 'sampled': 0, 'complete': 0}
        
        # Process each simulation
        for sim in SIMULATIONS:
            sim_id = sim['id']
            p_val = sim['p_msmw_w']  
            rep = sim['replicate']
            
            print(f"Processing {sim_id} (p_msmw_w={p_val}, rep={rep})...")
            
            # Load ground truth results
            gt_file = f"{OUTPUT_BASE}/simulations/{sim_id}/ground_truth_bridging.csv"
            gt_data = None
            if os.path.exists(gt_file) and os.path.getsize(gt_file) > 0:
                try:
                    gt_df = pd.read_csv(gt_file, index_col=0)
                    # Convert pivot to long format
                    gt_long = gt_df.stack().reset_index()
                    gt_long.columns = ['measure', 'grouping', 'ground_truth_value']
                    gt_data = gt_long
                    sim_success_stats['ground_truth'] += 1
                    print(f"  ✅ Ground truth: {len(gt_long)} records")
                except Exception as e:
                    print(f"  ❌ Ground truth error: {e}")
            
            # Load phylogenetic results  
            phylo_file = f"{OUTPUT_BASE}/phylo/{sim_id}/aggregated_phylogenetic_bridging_30_results.csv"
            phylo_data = None
            if os.path.exists(phylo_file) and os.path.getsize(phylo_file) > 0:
                try:
                    phylo_df = pd.read_csv(phylo_file)
                    sim_success_stats['phylo'] += 1
                    print(f"  ✅ Phylogenetic: {len(phylo_df)} records")
                    phylo_data = phylo_df
                except Exception as e:
                    print(f"  ❌ Phylogenetic error: {e}")
            
            # Load sampled results
            sampled_file = f"{OUTPUT_BASE}/phylo/{sim_id}/aggregated_sampled_bridging_30_results.csv"  
            sampled_data = None
            if os.path.exists(sampled_file) and os.path.getsize(sampled_file) > 0:
                try:
                    sampled_df = pd.read_csv(sampled_file)
                    sim_success_stats['sampled'] += 1
                    print(f"  ✅ Sampled: {len(sampled_df)} records")
                    sampled_data = sampled_df
                except Exception as e:
                    print(f"  ❌ Sampled error: {e}")
            
            # Create comprehensive comparison for this simulation
            if gt_data is not None:
                # Start with ground truth as base
                for _, gt_row in gt_data.iterrows():
                    measure = gt_row['measure']
                    grouping = gt_row['grouping']
                    
                    comparison_row = {
                        'sim_id': sim_id,
                        'p_msmw_w': p_val,
                        'replicate': rep,
                        'measure': measure,
                        'grouping': grouping,
                        'ground_truth_value': gt_row['ground_truth_value']
                    }
                    
                    # Add phylogenetic values (multiple datasets per measure/grouping)
                    if phylo_data is not None:
                        phylo_matches = phylo_data[
                            (phylo_data['measure'] == measure) 
                        ]
                        
                        # Group by dataset characteristics
                        for _, phylo_row in phylo_matches.iterrows():
                            row_copy = comparison_row.copy()
                            row_copy.update({
                                'dataset_type': phylo_row['dataset_type'],
                                'dataset': phylo_row['dataset'],
                                'downsample_rate': phylo_row['downsample_rate'],
                                'replicate_phylo': phylo_row['replicate'],
                                'population': phylo_row['population']
                            })
                            
                            # Add phylo values (find matching grouping column)
                            phylo_col = f"phylo_30_{grouping}"
                            if phylo_col in phylo_row:
                                row_copy['phylo_value'] = phylo_row[phylo_col]
                            
                            all_comparisons.append(row_copy)
                    
                    # Add sampled values (multiple datasets per measure/grouping)  
                    if sampled_data is not None:
                        sampled_matches = sampled_data[
                            (sampled_data['measure'] == measure)
                        ]
                        
                        for _, sampled_row in sampled_matches.iterrows():
                            # Find existing comparison row or create new one
                            matching_rows = [r for r in all_comparisons if 
                                           r['sim_id'] == sim_id and 
                                           r['measure'] == measure and
                                           r['grouping'] == grouping and
                                           r.get('dataset') == sampled_row['dataset']]
                            
                            if matching_rows:
                                # Add to existing row
                                sampled_col = f"sampled_30_{grouping}"
                                if sampled_col in sampled_row:
                                    matching_rows[0]['sampled_value'] = sampled_row[sampled_col]
                            else:
                                # Create new row for sampled-only data
                                row_copy = comparison_row.copy()
                                row_copy.update({
                                    'dataset_type': sampled_row['dataset_type'],
                                    'dataset': sampled_row['dataset'], 
                                    'downsample_rate': sampled_row['downsample_rate'],
                                    'replicate_sampled': sampled_row['replicate'],
                                    'population': sampled_row['population']
                                })
                                
                                sampled_col = f"sampled_30_{grouping}"
                                if sampled_col in sampled_row:
                                    row_copy['sampled_value'] = sampled_row[sampled_col]
                                
                                all_comparisons.append(row_copy)
                
                if gt_data is not None and phylo_data is not None and sampled_data is not None:
                    sim_success_stats['complete'] += 1
        
        # Save comprehensive comparison
        if all_comparisons:
            comparison_df = pd.DataFrame(all_comparisons)
            comparison_df = comparison_df.sort_values(['p_msmw_w', 'replicate', 'measure', 'grouping', 'dataset_type', 'downsample_rate'])
            comparison_df.to_csv(output.comprehensive_comparison, index=False)
            
            print(f"🎉 Comprehensive comparison saved: {len(comparison_df)} records")
        else:
            pd.DataFrame().to_csv(output.comprehensive_comparison, index=False)
            print("❌ No comparison data to save")
        
        # Create summary statistics
        summary_stats = {
            'total_simulations': sim_success_stats['total'],
            'successful_ground_truth': sim_success_stats['ground_truth'],
            'successful_phylogenetic': sim_success_stats['phylo'],
            'successful_sampled': sim_success_stats['sampled'],
            'complete_analyses': sim_success_stats['complete'],
            'ground_truth_success_rate': sim_success_stats['ground_truth'] / sim_success_stats['total'],
            'phylogenetic_success_rate': sim_success_stats['phylo'] / sim_success_stats['total'],
            'sampled_success_rate': sim_success_stats['sampled'] / sim_success_stats['total'],
            'complete_success_rate': sim_success_stats['complete'] / sim_success_stats['total']
        }
        
        # Add parameter value breakdown
        p_value_stats = {}
        for p_val in sorted(set(sim['p_msmw_w'] for sim in SIMULATIONS)):
            p_sims = [sim for sim in SIMULATIONS if sim['p_msmw_w'] == p_val]
            p_value_stats[f'p_msmw_w_{p_val}_total'] = len(p_sims)
        
        summary_stats.update(p_value_stats)
        
        summary_df = pd.DataFrame([summary_stats])
        summary_df.to_csv(output.summary_statistics, index=False)
        
        print(f"📊 Summary Statistics:")
        print(f"  Total simulations: {sim_success_stats['total']}")
        print(f"  Complete analyses: {sim_success_stats['complete']} ({sim_success_stats['complete']/sim_success_stats['total']:.1%})")
        print(f"  Ground truth success: {sim_success_stats['ground_truth']} ({sim_success_stats['ground_truth']/sim_success_stats['total']:.1%})")
        print(f"  Phylogenetic success: {sim_success_stats['phylo']} ({sim_success_stats['phylo']/sim_success_stats['total']:.1%})")
        print(f"  Sampled success: {sim_success_stats['sampled']} ({sim_success_stats['sampled']/sim_success_stats['total']:.1%})")
        
        print("🎉 Comprehensive aggregation complete!")




