#!/usr/bin/env python3
"""
Complete Single Simulation Pipeline: Simulation → Ground Truth → Phylogenetic Analysis
Uses working f-string approach for downstream steps + processes full_sequences too
"""

import pandas as pd
import json
import os

wildcard_constraints:
    sim_id=r"[0-9]+\.[0-9]+_[0-9]+",
    batch_id=r"[0-9]+"

# Configuration
OUTPUT_BASE = "/n/netscratch/grad_lab/Lab/mkline/bridging_project/output/baseline_phylo_bridging_complete_50each"
PARAM_DIR = "/n/netscratch/grad_lab/Lab/mkline/bridging_project/code/parameters/11_13_selected_baseline_params_50/params"
PIPELINE_SEED = 42

# Phylo pipeline parameters
DOWNSAMPLE_RATES = [0.5, 0.25, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01]
N_REPLICATES = 10
HIGH_ACTIVITY_SUBSAMPLE = 0.25

# Software paths
SEQGEN_PATH = "/n/holylfs05/LABS/grad_lab/Users/mkline/bridging_sims/software/Seq-Gen/source/seq-gen"
LSD2_PATH = "/n/holylfs05/LABS/grad_lab/Users/mkline/bridging_sims/software/lsd2/src/lsd2"

# Resource limits
MAX_MEM_MB = 128000
MAX_RUNTIME = 2880
MAX_CPUS = 16


# Partition functions
def get_partition(wildcards, attempt):
    partitions = ["hsph", "sapphire", "intermediate", "shared"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

def get_quick_partition(wildcards, attempt):
    partitions = ["hsph", "sapphire", "shared", "intermediate"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

def get_long_partition(wildcards, attempt):
    partitions = ["intermediate", "hsph", "sapphire", "shared"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

# Batch configuration functions
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
    batch_df = pd.read_csv(input.batch_assignments)
    batch_data = batch_df[batch_df['batch_id'] == int(wildcards.batch_id)]
    if len(batch_data) > 0:
        value = batch_data.iloc[0][resource_type]
        return int(value)
    return 4000 if resource_type == 'mem_mb' else (60 if resource_type == 'runtime' else 1)

# ============================================================================
# COMPLETE PIPELINE RULES
# ============================================================================

rule all:
    input:
        f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_sampled_bridging_30_results.csv".format(sim_id=config["sim_id"]),
        f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_phylogenetic_bridging_30_results.csv".format(sim_id=config["sim_id"])

rule run_simulation:
    input:
        param_file = f"{PARAM_DIR}/params_{{sim_id}}.json"
    output:
        simulation_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/simulation_complete.flag"
    conda: "phylo"
    resources:
        mem_mb = 32000,
        runtime = 480,
        cpus_per_task = 2,
        slurm_partition = get_partition
    shell:
        """
        echo "🚀 Running simulation for {wildcards.sim_id}"
        mkdir -p {OUTPUT_BASE}/simulations/{wildcards.sim_id}

        python -c "
import json

# Load and modify parameters
with open('{input.param_file}', 'r') as f:
    params = json.load(f)

params['output']['output_dir'] = '{OUTPUT_BASE}/simulations/{wildcards.sim_id}'
params['output']['run_name'] = 'sim'

# Save modified parameters
temp_param_file = '{OUTPUT_BASE}/simulations/{wildcards.sim_id}/modified_params.json'
with open(temp_param_file, 'w') as f:
    json.dump(params, f, indent=2)
"

        # Run simulation
        python ../scripts/10_16_core_sim.py {OUTPUT_BASE}/simulations/{wildcards.sim_id}/modified_params.json

        # Find actual simulation directory
        python -c "
import glob
import os

patterns = [
    '{OUTPUT_BASE}/simulations/{wildcards.sim_id}/sim_*',
    '{OUTPUT_BASE}/simulations/{wildcards.sim_id}/sim'
]

actual_dir = None
for pattern in patterns:
    dirs = glob.glob(pattern)
    if dirs:
        actual_dir = dirs[0] if len(dirs) == 1 else max(dirs, key=os.path.getctime)
        break

if actual_dir and os.path.isdir(actual_dir):
    required_files = ['transmission_df.csv', 'parameters_used.json', 'nodes_df.csv', 'initial_infectors.csv', 'edge_df.csv']
    missing = [f for f in required_files if not os.path.exists(os.path.join(actual_dir, f))]

    if not missing:
        with open('{output.simulation_flag}', 'w') as f:
            f.write(f'ACTUAL_SIM_DIR={{actual_dir}}\\n')
        print(f'✅ Simulation completed: {{actual_dir}}')
    else:
        print(f'ERROR: Missing files: {{missing}}')
        exit(1)
else:
    print('ERROR: No simulation directory found')
    exit(1)
"
        """

rule check_lineage:
    input:
        simulation_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/simulation_complete.flag"
    output:
        lineage_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/lineage_check_complete.flag"
    conda: "phylo"
    resources:
        mem_mb = 4000,
        runtime = 60,
        cpus_per_task = 1,
        slurm_partition = "shared"
    shell:
        """
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.simulation_flag} | cut -d'=' -f2)

        python ../scripts/10_17_single_lineage_test_script.py "$ACTUAL_SIM_DIR"

        if [ -f "$ACTUAL_SIM_DIR/lineage_test_results.json" ]; then
            SINGLE_LINEAGE=$(python -c "import json; print(json.load(open('$ACTUAL_SIM_DIR/lineage_test_results.json'))['single_lineage'])")

            echo "SINGLE_LINEAGE=$SINGLE_LINEAGE" > {output.lineage_flag}
            echo "ACTUAL_SIM_DIR=$ACTUAL_SIM_DIR" >> {output.lineage_flag}

            if [ "$SINGLE_LINEAGE" = "True" ]; then
                echo "✅ Single lineage confirmed - continuing pipeline"
            else
                echo "❌ Multiple lineages detected - stopping pipeline"
                exit 1
            fi
        else
            echo "❌ Lineage check failed"
            exit 1
        fi
        """

rule ground_truth_analysis:
    input:
        simulation_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/simulation_complete.flag"
    output:
        bridging_values = f"{OUTPUT_BASE}/simulations/{{sim_id}}/ground_truth_bridging.csv",
        gt_analysis_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/ground_truth_analysis_complete.flag"
    conda: "phylo"
    resources:
        mem_mb = 12000,
        runtime = 90,
        cpus_per_task = 2,
        slurm_partition = get_partition
    shell:
        """
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.simulation_flag} | cut -d'=' -f2)

        # Run behavior classification
        python ../scripts/10_17_infection_behavior_classification_analysis.py "$ACTUAL_SIM_DIR" --both

        # Run ground truth bridging analysis
        python ../scripts/10_30_ground_truth_bridging.py \\
            "$ACTUAL_SIM_DIR" \\
            "$ACTUAL_SIM_DIR/infection_behavior_classifications.csv" \\
            --output "{output.bridging_values}" \\
            --detailed-output "{OUTPUT_BASE}/simulations/{wildcards.sim_id}/ground_truth_bridging_detailed.csv"

        echo "ACTUAL_SIM_DIR=$ACTUAL_SIM_DIR" > {output.gt_analysis_flag}
        """

rule build_transmission_tree:
    input:
        lineage_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/lineage_check_complete.flag"
    output:
        tree_building_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/tree_building_complete.flag"
    conda: "phylo"
    resources:
        mem_mb = 64000,
        runtime = 180,
        cpus_per_task = 4,
        slurm_partition = get_partition
    shell:
        """
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.lineage_flag} | cut -d'=' -f2)

        # Build hybrid transmission tree
        python ../scripts/10_21_build_hybrid_transmission_tree.py "$ACTUAL_SIM_DIR"

        # First trimming
        TREE_FILE=$(ls "$ACTUAL_SIM_DIR"/tree_*_hybrid_episodes*.newick | head -1)
        Rscript ../scripts/10_20_tree_trimmer.R "$TREE_FILE" "$ACTUAL_SIM_DIR" "$ACTUAL_SIM_DIR/trimmed_tree.nwk"

        # Additional trimming
        Rscript ../scripts/10_21_tree_additional_trimming.R \\
            "$ACTUAL_SIM_DIR/trimmed_tree.nwk" \\
            "$ACTUAL_SIM_DIR/transmission_df.csv" \\
            "$ACTUAL_SIM_DIR/final_trimmed_tree.nwk"

        echo "ACTUAL_SIM_DIR=$ACTUAL_SIM_DIR" > {output.tree_building_flag}
        echo "TREE_READY=True" >> {output.tree_building_flag}
        """

rule run_seqgen:
    input:
        tree_building_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/tree_building_complete.flag"
    output:
        seqgen_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/seqgen_complete.flag"
    resources:
        mem_mb = 2000,
        runtime = 30,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    shell:
        """
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.tree_building_flag} | cut -d'=' -f2)

        mkdir -p {OUTPUT_BASE}/phylo/{wildcards.sim_id}/sequences

        {SEQGEN_PATH} \\
            -m GTR -r 1.00578 4.59475 0.74535 0.46171 4.55284 1.00000 \\
            -f 0.154 0.346 0.343 0.157 -a 0.45 -l 10000 -n 1 -s 2.254247e-06 -z {PIPELINE_SEED} \\
            < "$ACTUAL_SIM_DIR/final_trimmed_tree.nwk" > {OUTPUT_BASE}/phylo/{wildcards.sim_id}/sequences/sequences.fasta

        echo "SEQUENCES_READY=True" > {output.seqgen_flag}
        echo "ACTUAL_SIM_DIR=$ACTUAL_SIM_DIR" >> {output.seqgen_flag}
        """

rule run_filtering_downsampling:
    input:
        seqgen_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/seqgen_complete.flag"
    output:
        filtering_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/filtering_complete.flag"
    params:
        rates_str = " ".join(map(str, DOWNSAMPLE_RATES))
    resources:
        mem_mb = 16000,
        runtime = 60,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    conda: "phylo"
    shell:
        """
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.seqgen_flag} | cut -d'=' -f2)

        mkdir -p {OUTPUT_BASE}/phylo/{wildcards.sim_id}/filtering

        # Copy full sequences
        cp {OUTPUT_BASE}/phylo/{wildcards.sim_id}/sequences/sequences.fasta {OUTPUT_BASE}/phylo/{wildcards.sim_id}/filtering/full_sequences.fasta

        # Run filtering and downsampling
        python ../scripts/10_26_infection_filtering_downsampling.py \\
            "$ACTUAL_SIM_DIR" \\
            {OUTPUT_BASE}/phylo/{wildcards.sim_id}/sequences/sequences.fasta \\
            {OUTPUT_BASE}/phylo/{wildcards.sim_id}/filtering \\
            --downsample-rates {params.rates_str} \\
            --n-replicates {N_REPLICATES} \\
            --high-activity-subsample {HIGH_ACTIVITY_SUBSAMPLE} \\
            --seed {PIPELINE_SEED}

        echo "FILTERING_READY=True" > {output.filtering_flag}
        echo "ACTUAL_SIM_DIR=$ACTUAL_SIM_DIR" >> {output.filtering_flag}
        """

rule create_full_sequences_downsamples:
    input:
        filtering_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/filtering_complete.flag"
    output:
        full_sequences_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/full_sequences_complete.flag"
    resources:
        mem_mb = 8000,
        runtime = 60,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    conda: "phylo"
    shell:
        """
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {input.filtering_flag} | cut -d'=' -f2)

        echo "📊 Creating full_sequences dataset + downsamples using external script..."

        python ../scripts/create_full_sequences_downsamples.py \\
            "$ACTUAL_SIM_DIR" \\
            "{wildcards.sim_id}" \\
            "{OUTPUT_BASE}" \\
            "{PIPELINE_SEED}" \\
            "{output.full_sequences_flag}"
        """

checkpoint create_phylo_batches:
    input:
        filtering_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/filtering_complete.flag",
        full_sequences_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/full_sequences_complete.flag"
    output:
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv",
        batching_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batching_complete.flag"
    run:
        import pandas as pd
        import os
        import re

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
                dataset_id = f"{population}_{dataset}"
            elif dataset == 'full':
                dataset_type = 'full_sequences'
                downsample_rate = None
                batch_config = get_batch_config('full_sequences')
                
                # FIX: Extract downsample info from the directory path
                # Pattern: .../full_sequences/50pct/rep01/sequences.fasta or .../full_sequences/complete/sequences.fasta
                if '/complete/' in fasta_file:
                    # Original full sequences (no downsampling)
                    dataset_id = 'full_sequences'
                else:
                    # Downsampled full sequences - extract rate and replicate from path
                    match = re.search(r'/full_sequences/(\d+pct)/(rep\d+)/', fasta_file)
                    if match:
                        rate = match.group(1)  # e.g., "50pct"
                        replicate = match.group(2)  # e.g., "rep01"
                        dataset_id = f"full_sequences_{rate}_{replicate}"
                    else:
                        # Fallback if pattern doesn't match
                        print(f"Warning: Could not parse downsample info from {fasta_file}")
                        dataset_id = 'full_sequences'
            else:
                dataset_type = 'downsampled'
                downsample_rate = int(dataset.split('pct')[0])
                batch_config = get_batch_config('downsampled', downsample_rate)
                dataset_id = f"{population}_{dataset}"

            batch_assignments.append({
                'dataset_id': dataset_id,  # ← NOW PROPERLY SET
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

        # Rest stays the same...
        batch_df = pd.DataFrame(batch_assignments)

        # Assign actual batch numbers within tiers (DYNAMIC)
        batch_df['batch_id'] = 0
        current_batch = 1

        for tier in batch_df['batch_tier'].unique():
            tier_data = batch_df[batch_df['batch_tier'] == tier]
            batch_size = tier_data.iloc[0]['batch_size']

            for i in range(0, len(tier_data), batch_size):
                batch_df.loc[tier_data.index[i:i+batch_size], 'batch_id'] = current_batch
                current_batch += 1

        # Calculate actual max batch ID dynamically
        max_batch_id = batch_df['batch_id'].max()

        batch_df.to_csv(output.batch_assignments, index=False)

        # Create flag with DYNAMIC batch count AND ACTUAL_SIM_DIR
        with open(output.batching_flag, 'w') as f:
            f.write("BATCHING_READY=True\n")
            f.write(f"MAX_BATCH_ID={max_batch_id}\n")
            f.write(f"NUM_BATCHES={max_batch_id}\n")
            f.write(f"NUM_DATASETS={len(batch_df)}\n")
            f.write(f"ACTUAL_SIM_DIR={actual_sim_dir}\n")

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
        gt_analysis_flag = f"{OUTPUT_BASE}/simulations/{{sim_id}}/ground_truth_analysis_complete.flag"
    output:
        iqtree_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_iqtree_complete.flag",
        sampled_bridging_30_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_sampled_bridging_30_complete.flag",
        dates_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_dates_complete.flag",
        behaviors_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_behaviors_complete.flag"
    resources:
        mem_mb = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "mem_mb"), 4000), MAX_MEM_MB),
        runtime = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "runtime"), 60), MAX_RUNTIME),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_partition
    conda: "phylo"
    shell:
        """
        echo "🔄 Starting parallel batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"
        mkdir -p {OUTPUT_BASE}/phylo/{wildcards.sim_id}/flags

        # Get actual simulation directory from batching flag
        ACTUAL_SIM_DIR=$(grep ACTUAL_SIM_DIR {OUTPUT_BASE}/phylo/{wildcards.sim_id}/batching_complete.flag | cut -d'=' -f2)

        python -c "
import pandas as pd
import os

def process_dataset(row, sim_dir, sim_id):
    dataset_id = row['dataset_id']
    fasta_file = row['fasta_file']
    infections_file = row['infections_file']

    print(f'Processing {{dataset_id}} for {{sim_id}}...', flush=True)

    # Determine infections file to use
    if dataset_id == 'full_sequences':
        # For full sequences, apply post-burn-in filtering
        print(f'  Applying post-burn-in filtering for full_sequences...', flush=True)

        full_transmission_df = pd.read_csv(os.path.join(sim_dir, 'transmission_df.csv'))
        print(f'  Loaded {{len(full_transmission_df)}} total transmissions', flush=True)

        try:
            import json
            with open(os.path.join(sim_dir, 'parameters_used.json'), 'r') as f:
                params = json.load(f)
            cutoff_day = (params['simulation']['partnership_burnin_days'] +
                         params['simulation']['transmission_burnin_days'])
            print(f'  Using burn-in cutoff: {{cutoff_day}}', flush=True)
        except Exception as e:
            print(f'  Warning: Could not read parameters, using default cutoff 12000: {{e}}', flush=True)
            cutoff_day = 12000

        filtered_transmissions = full_transmission_df[
            (full_transmission_df['superseded_simultaneous'] == False) &
            (full_transmission_df['day_of_transmission'] >= cutoff_day) &
            (pd.notna(full_transmission_df['day_of_sampling']))
        ]

        print(f'  After filtering: {{len(filtered_transmissions)}} transmissions', flush=True)

        dataset_dir = os.path.dirname(fasta_file)
        temp_infections_file = os.path.join(dataset_dir, 'temp_filtered_transmissions.csv')
        filtered_transmissions.to_csv(temp_infections_file, index=False)
        infections_to_use = temp_infections_file
    else:
        infections_to_use = infections_file

    dataset_dir = os.path.dirname(fasta_file)

    # 1. IQ-TREE (WORKING F-STRING APPROACH)
    if os.path.exists(fasta_file):
        treefile_expected = fasta_file + '.treefile'
        if not os.path.exists(treefile_expected) or os.path.getsize(treefile_expected) == 0:
            print(f'Running IQ-TREE for {{dataset_id}}...', flush=True)
            cmd1 = f'iqtree -s {{fasta_file}} -m GTR+G --fast -seed {PIPELINE_SEED}'
            result1 = os.system(cmd1)
            print(f'IQTree result: {{result1}}')
        else:
            print(f'SKIP: iqtree already completed for {{dataset_id}}')

    # 2. Dates extraction (WORKING APPROACH)
    if os.path.exists(infections_to_use):
        dates_output = os.path.join(dataset_dir, 'dates.txt')
        if not os.path.exists(dates_output) or os.path.getsize(dates_output) == 0:
            print(f'Running dates extraction for {{dataset_id}}...', flush=True)
            cmd2 = f'python ../scripts/10_27_extract_tip_dates.py {{infections_to_use}} {{fasta_file}} {{sim_dir}} {{dates_output}}'
            result2 = os.system(cmd2)
            print(f'Dates extraction result: {{result2}}')
        else:
            print(f'SKIP: dates already extracted for {{dataset_id}}')

        # 3. Behaviors extraction (WORKING APPROACH)
        behaviors_msmw_msm = os.path.join(dataset_dir, 'behaviors_msmw_msm_states.csv')
        behaviors_msmw_msw = os.path.join(dataset_dir, 'behaviors_msmw_msw_states.csv')

        if (not os.path.exists(behaviors_msmw_msm) or os.path.getsize(behaviors_msmw_msm) == 0 or
            not os.path.exists(behaviors_msmw_msw) or os.path.getsize(behaviors_msmw_msw) == 0):
            print(f'Running behaviors extraction for {{dataset_id}}...', flush=True)
            behaviors_prefix = os.path.join(dataset_dir, 'behaviors')
            cmd3 = f'python ../scripts/10_27_extract_tip_states.py {{infections_to_use}} {{fasta_file}} {{sim_dir}} {{behaviors_prefix}}'
            result3 = os.system(cmd3)
            print(f'Behaviors extraction result: {{result3}}')
        else:
            print(f'SKIP: behaviors already extracted for {{dataset_id}}')

        # 4. Sampled bridging analysis (skip for full_sequences)
        if dataset_id != 'full_sequences':
            bridging_30_output = os.path.join(dataset_dir, 'sampled_bridging_30.csv')
            bridging_30_detailed = os.path.join(dataset_dir, 'sampled_bridging_30_detailed.csv')

            if (not os.path.exists(bridging_30_output) or os.path.getsize(bridging_30_output) == 0 or
                not os.path.exists(bridging_30_detailed) or os.path.getsize(bridging_30_detailed) == 0):
                print(f'Running sampled bridging 30 for {{dataset_id}}...', flush=True)
                cmd4 = f'python ../scripts/10_30_sampled_bridging.py {{infections_to_use}} {{sim_dir}} --output {{bridging_30_output}} --detailed-output {{bridging_30_detailed}}'
                result4 = os.system(cmd4)
                print(f'Sampled bridging 30 result: {{result4}}')
            else:
                print(f'SKIP: sampled bridging already completed for {{dataset_id}}')
        else:
            print(f'SKIP: sampled bridging not applicable for full_sequences')

    # Clean up temporary file
    if dataset_id == 'full_sequences' and 'temp_filtered_transmissions.csv' in infections_to_use:
        try:
            os.remove(infections_to_use)
        except:
            pass

batch_df = pd.read_csv('{input.batch_assignments}')
batch_data = batch_df[batch_df['batch_id'] == {wildcards.batch_id}]

for _, row in batch_data.iterrows():
    process_dataset(row, '$ACTUAL_SIM_DIR', '{wildcards.sim_id}')
"

        echo "✅ Parallel batch {wildcards.batch_id} completed at $(date)"
        touch {output.iqtree_flag} {output.sampled_bridging_30_flag} {output.dates_flag} {output.behaviors_flag}
        """

rule run_lsd2_batch:
    input:
        iqtree_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_iqtree_complete.flag",
        dates_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_dates_complete.flag",
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv"
    output:
        lsd2_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_lsd2_complete.flag"
    resources:
        mem_mb = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "mem_mb"), 64000), MAX_MEM_MB),
        runtime = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "runtime"), 720), MAX_RUNTIME),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_long_partition
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

        # Change to dataset directory for LSD2 (WORKING APPROACH)
        original_dir = os.getcwd()
        os.chdir(dataset_dir)

        tree_basename = os.path.basename(tree_file)
        dates_basename = os.path.basename(dates_file)

        cmd = f'{LSD2_PATH} -i {{tree_basename}} -d {{dates_basename}} -s 10000 -r a -l 0'
        result = os.system(cmd)

        os.chdir(original_dir)
        print(f'LSD2 result: {{result}}')
    else:
        print(f'Missing prerequisites for {{dataset_id}}')
"

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
        mem_mb = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "mem_mb"), 16000), MAX_MEM_MB),
        runtime = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "runtime"), 240), MAX_RUNTIME),
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
        cmd = f'treetime mugration --tree {{lsd2_tree}} --states {{states_file}} --attribute state --outdir {{outdir}} --verbose 2 --confidence'
        result = os.system(cmd)
        print(f'TreeTime MSMW+MSM result: {{result}}')
    else:
        print(f'Missing prerequisites for {{dataset_id}}')
"

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
        mem_mb = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "mem_mb"), 16000), MAX_MEM_MB),
        runtime = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "runtime"), 240), MAX_RUNTIME),
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
        cmd = f'treetime mugration --tree {{lsd2_tree}} --states {{states_file}} --attribute state --outdir {{outdir}} --verbose 2 --confidence'
        result = os.system(cmd)
        print(f'TreeTime MSMW+MSW result: {{result}}')
    else:
        print(f'Missing prerequisites for {{dataset_id}}')
"

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
        mem_mb = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "mem_mb"), 16000), MAX_MEM_MB),
        runtime = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "runtime"), 120), MAX_RUNTIME),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_partition
    conda: "phylo"
    shell:
        """
        echo "🌲 Starting Tree Bridging 30 batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"

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

    # Output files from 10_30_treetime_bridging.py
    bridging_output = os.path.join(dataset_dir, 'tree_bridging_30_treetime_bridging.csv')
    bridging_detailed = os.path.join(dataset_dir, 'tree_bridging_30_treetime_bridging_detailed.csv')

    if os.path.exists(msmw_msm_nexus) and os.path.exists(msmw_msw_nexus) and os.path.exists(lsd2_result):
        print(f'Running Tree Bridging 30 for {{dataset_id}}...', flush=True)

        bridging_prefix = os.path.join(dataset_dir, 'tree_bridging_30')
        cmd = f'python ../scripts/10_30_treetime_bridging.py {{msmw_msm_nexus}} {{msmw_msw_nexus}} {{lsd2_result}} {{bridging_prefix}}'
        result = os.system(cmd)
        print(f'Tree Bridging 30 result: {{result}}')
    else:
        print(f'Missing prerequisites for {{dataset_id}}')
"

        touch {output.tree_bridging_30_flag}
        """

# Remove the old get_batch_flags_for_aggregation function entirely
# Replace the aggregation rules with these:

def get_all_batch_flags(wildcards, flag_type):
    """Get batch flags dynamically from batch assignments"""
    batch_file = f"{OUTPUT_BASE}/phylo/{wildcards.sim_id}/batch_assignments.csv"
    if os.path.exists(batch_file):
        import pandas as pd
        batch_df = pd.read_csv(batch_file)
        batch_ids = sorted(batch_df['batch_id'].unique())
        return [f"{OUTPUT_BASE}/phylo/{wildcards.sim_id}/flags/batch_{bid}_{flag_type}_complete.flag" 
                for bid in batch_ids]
    return []

rule aggregate_sampled_bridging_30:
    input:
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv",
        # Use a function that gets called at runtime, not during DAG building
        sampled_bridging_flags = lambda wildcards: get_all_batch_flags(wildcards, "sampled_bridging_30") 
    output:
        aggregated_sampled_bridging_30 = f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_sampled_bridging_30_results.csv"
    resources:
        mem_mb = 4000,
        runtime = 60,
        cpus_per_task = 1,
        slurm_partition = get_partition
    run:
        import pandas as pd
        import os

        print("📊 Aggregating sampled bridging 30 results...")

        batch_df = pd.read_csv(input.batch_assignments)
        all_results = []

        for _, row in batch_df.iterrows():
            dataset_id = row['dataset_id']
            population = row['population']
            dataset = row['dataset']
            fasta_file = row['fasta_file']

            # Skip full sequences for sampled bridging
            if dataset == 'full' or dataset_id == 'full_sequences':
                continue

            dataset_dir = os.path.dirname(fasta_file)
            sampled_bridging_file = os.path.join(dataset_dir, 'sampled_bridging_30.csv')

            if os.path.exists(sampled_bridging_file) and os.path.getsize(sampled_bridging_file) > 0:
                try:
                    sampled_df = pd.read_csv(sampled_bridging_file, index_col=0)
                    sampled_df = sampled_df.reset_index()

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
                            result_row['replicate'] = parts[1] if len(parts) > 1 else 'unknown'

                        for col in sampled_df.columns:
                            if col != 'measure':
                                result_row[f'sampled_30_{col}'] = sampled_row[col]

                        all_results.append(result_row)

                    print(f"✅ Collected sampled bridging 30 from {dataset_id}")
                except Exception as e:
                    print(f"❌ Error processing {sampled_bridging_file}: {e}")

        if all_results:
            aggregated_df = pd.DataFrame(all_results)
            column_order = [
                'sim_id', 'dataset_id', 'population', 'dataset_type', 'dataset',
                'downsample_rate', 'replicate', 'measure'
            ] + [col for col in aggregated_df.columns if col.startswith('sampled_30_')]

            aggregated_df = aggregated_df[column_order]
            aggregated_df = aggregated_df.sort_values(['population', 'downsample_rate', 'replicate', 'measure'])
            aggregated_df.to_csv(output.aggregated_sampled_bridging_30, index=False)

            print(f"🎉 Sampled bridging 30 aggregation complete! ({len(aggregated_df)} records)")
        else:
            pd.DataFrame().to_csv(output.aggregated_sampled_bridging_30, index=False)
            print("❌ No sampled bridging results found")

rule aggregate_phylogenetic_bridging_30:
    input:
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv",
        # Use a function that gets called at runtime, not during DAG building  
        bridging_flags = lambda wildcards: get_all_batch_flags(wildcards, "tree_bridging_30")
    output:
        aggregated_phylo_bridging_30 = f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_phylogenetic_bridging_30_results.csv"
    resources:
        mem_mb = 4000,
        runtime = 60,
        cpus_per_task = 1,
        slurm_partition = get_partition
    run:
        import pandas as pd
        import os
        import re

        print("🌳 Aggregating phylogenetic bridging 30 results...")

        batch_df = pd.read_csv(input.batch_assignments)
        all_results = []

        for _, row in batch_df.iterrows():
            dataset_id = row['dataset_id']
            population = row['population']
            dataset = row['dataset']
            fasta_file = row['fasta_file']

            dataset_dir = os.path.dirname(fasta_file)
            bridging_summary_file = os.path.join(dataset_dir, 'tree_bridging_30_treetime_bridging.csv')

            if os.path.exists(bridging_summary_file) and os.path.getsize(bridging_summary_file) > 0:
                try:
                    bridging_df = pd.read_csv(bridging_summary_file, index_col=0)
                    bridging_df = bridging_df.reset_index()

                    for _, bridging_row in bridging_df.iterrows():
                        result_row = {
                            'sim_id': wildcards.sim_id,
                            'dataset_id': dataset_id,
                            'population': population,
                            'dataset_type': 'complete' if dataset == 'complete' else ('full_sequences' if dataset == 'full' else 'downsampled'),
                            'dataset': dataset,
                            'measure': bridging_row['measure']
                        }

                        # FIX: Properly handle full_sequences downsamples
                        if dataset == 'complete':
                            result_row['downsample_rate'] = '100'
                            result_row['replicate'] = 'all'
                        elif dataset == 'full':
                            # Check if this is a downsampled full_sequences dataset
                            if '_' in dataset_id and 'pct' in dataset_id:
                                # Parse: full_sequences_50pct_rep01 -> downsample_rate=50, replicate=rep01
                                parts = dataset_id.split('_')
                                for i, part in enumerate(parts):
                                    if 'pct' in part:
                                        result_row['downsample_rate'] = part.replace('pct', '')
                                        # Get replicate from next part if it exists
                                        if i + 1 < len(parts):
                                            result_row['replicate'] = parts[i + 1]
                                        else:
                                            result_row['replicate'] = 'unknown'
                                        break
                                else:
                                    # Fallback if pattern not found
                                    result_row['downsample_rate'] = '100'
                                    result_row['replicate'] = 'full_sequences'
                            else:
                                # Original full sequences (no downsampling)
                                result_row['downsample_rate'] = '100'
                                result_row['replicate'] = 'full_sequences'
                        else:
                            # Regular downsampled datasets
                            parts = dataset.split('_')
                            result_row['downsample_rate'] = parts[0].replace('pct', '')
                            result_row['replicate'] = parts[1] if len(parts) > 1 else 'unknown'

                        for col in bridging_df.columns:
                            if col != 'measure':
                                result_row[f'phylo_30_{col}'] = bridging_row[col]

                        all_results.append(result_row)

                    print(f"✅ Collected phylogenetic bridging 30 from {dataset_id}")
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
            aggregated_df.to_csv(output.aggregated_phylo_bridging_30, index=False)

            print(f"🎉 Phylogenetic bridging 30 aggregation complete! ({len(aggregated_df)} records)")
        else:
            pd.DataFrame().to_csv(output.aggregated_phylo_bridging_30, index=False)
            print("❌ No phylogenetic bridging results found")
