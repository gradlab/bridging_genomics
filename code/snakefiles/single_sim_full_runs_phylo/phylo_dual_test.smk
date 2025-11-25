# ============================================================================
# Phylogenetic Bridging Analysis Pipeline (Dual Methods: 10_27 + 10_28)
# ============================================================================
#
# USAGE: This pipeline should be run on a SIMULATION DIRECTORY that contains:
#   - final_trimmed_tree.nwk
#   - transmission_df.csv  
#   - nodes_df.csv
#   - parameters_used.json
#   - infection_behavior_classifications.csv (for 10_28 analysis)
#
# Example: snakemake --snakefile phylo_pipeline_dual_bridge.smk \
#            --config sim_dir="/path/to/baseline_params_0.75_rep_454" \
#            --cores 16 --cluster "sbatch --partition {resources.slurm_partition} ..."
#
# OUTPUT: Creates both 10_27 (original) and 10_28 (directional) bridging results
# ============================================================================

# Pipeline parameters
PIPELINE_SEED = 42
DOWNSAMPLE_RATES = [0.5, 0.25, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01]
N_REPLICATES = 10
HIGH_ACTIVITY_SUBSAMPLE = 0.25

# Software paths
SEQGEN_PATH = "/n/holylfs05/LABS/grad_lab/Users/mkline/bridging_sims/software/Seq-Gen/source/seq-gen"
LSD2_PATH = "/n/holylfs05/LABS/grad_lab/Users/mkline/bridging_sims/software/lsd2/src/lsd2"

# Estimate maximum batch ID (will be determined dynamically)
MAX_BATCH_ID = 70  # THIS CAUSED PROBLEMS, SHOULD MAKE DYNAMIC

# Resource scaling functions
MAX_MEM_MB = 64000    # 64GB max
MAX_RUNTIME = 2880    # 48 hours max
MAX_CPUS = 16         # 16 CPUs max

def get_mem_mb_scaling(base_mem_func):
    """Dynamic memory allocation with scaling"""
    def mem_func(wildcards, input):
        base_value = base_mem_func(wildcards, input)  # Call the function first
        return min(base_value, MAX_MEM_MB)
    return mem_func

def get_runtime_scaling(base_runtime_func):
    """Dynamic runtime allocation with scaling"""  
    def runtime_func(wildcards, input):
        base_value = base_runtime_func(wildcards, input)  # Call the function first
        return min(base_value, MAX_RUNTIME)
    return runtime_func

def get_cpus_scaling(base_cpus_func):
    """Dynamic CPU allocation with scaling"""
    def cpu_func(wildcards, input):
        base_value = base_cpus_func(wildcards, input)  # Call the function first
        return min(base_value, MAX_CPUS)
    return cpu_func

# Partition selection functions (fallback on retry)
def get_compute_partition(wildcards, attempt):
    partitions = ["hsph", "sapphire", "intermediate", "shared"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

def get_quick_partition(wildcards, attempt):
    partitions = ["hsph", "sapphire", "shared", "intermediate"] 
    return partitions[min(attempt - 1, len(partitions) - 1)]

def get_long_partition(wildcards, attempt):
    partitions = ["intermediate", "hsph", "sapphire", "shared"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

# ============================================================================
# Helper functions
# ============================================================================

def get_batch_config(dataset_type, downsample_rate=None):
    """Return batch configuration based on dataset characteristics"""
    
    if dataset_type == "full_sequences":
        return {
            'batch_size': 1,
            'mem_mb': 64000,
            'runtime': 2880,     # 48 hours for full sequences
            'cpus_per_task': 8
        }
    elif dataset_type == "complete":
        return {
            'batch_size': 1, 
            'mem_mb': 32000,
            'runtime': 720,      # 12 hours for complete datasets
            'cpus_per_task': 4
        }
    elif downsample_rate in [50, 25]:
        return {
            'batch_size': 3,
            'mem_mb': 24000,
            'runtime': 360,      # 6 hours for large downsamples
            'cpus_per_task': 2
        }
    elif downsample_rate in [10, 5]:
        return {
            'batch_size': 8,
            'mem_mb': 12000,
            'runtime': 180,      # 3 hours for medium downsamples
            'cpus_per_task': 1
        }
    else:
        return {
            'batch_size': 20,
            'mem_mb': 4000,
            'runtime': 120,      # 2 hours for small downsamples
            'cpus_per_task': 1
        }

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
    import pandas as pd
    batch_df = pd.read_csv(input.batch_assignments)
    batch_data = batch_df[batch_df['batch_id'] == int(wildcards.batch_id)]
    if len(batch_data) > 0:
        value = batch_data.iloc[0][resource_type]
        return int(value)  # Convert to int!
    return 4000 if resource_type == 'mem_mb' else (60 if resource_type == 'runtime' else 1)

# ============================================================================
# Pipeline Rules
# ============================================================================

rule all:
    input:
        # Original 10_27 bridging results
        config["sim_dir"] + "/phylo_analysis/aggregated_phylogenetic_bridging_27_results.csv",
        config["sim_dir"] + "/phylo_analysis/aggregated_sampled_bridging_27_results.csv",
        # New 10_28 directional bridging results  
        config["sim_dir"] + "/phylo_analysis/aggregated_phylogenetic_bridging_28_results.csv",
        config["sim_dir"] + "/phylo_analysis/aggregated_sampled_bridging_28_results.csv"

rule run_seqgen:
    input:
        tree = "{sim_dir}/final_trimmed_tree.nwk"
    output:
        sequences = "{sim_dir}/phylo_analysis/sequences/sequences.fasta",
        log = "{sim_dir}/phylo_analysis/sequences/sequences.log"
    resources:
        mem_mb = 2000,
        runtime = 30,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    shell:
        """
        echo "Starting seq-gen at $(date)" > {output.log}
        echo "Input tree: {input.tree}" >> {output.log}
        
        mkdir -p $(dirname {output.sequences})
        
        {SEQGEN_PATH} \
            -m GTR -r 1.00578 4.59475 0.74535 0.46171 4.55284 1.00000 \
            -f 0.154 0.346 0.343 0.157 -a 0.45 -l 10000 -n 1 -s 2.254247e-06 -z {PIPELINE_SEED} \
            < {input.tree} > {output.sequences} 2>> {output.log}
        
        if [ ! -s {output.sequences} ]; then
            echo "ERROR: seq-gen output is empty!" >> {output.log}
            exit 1
        fi
        
        SEQ_COUNT=$(grep -c '^>' {output.sequences})
        echo "Generated $SEQ_COUNT sequences" >> {output.log}
        echo "Completed seq-gen at $(date)" >> {output.log}
        """

rule run_filtering_downsampling:
    input:
        sequences = "{sim_dir}/phylo_analysis/sequences/sequences.fasta",
        transmission_df = "{sim_dir}/transmission_df.csv",
        nodes_df = "{sim_dir}/nodes_df.csv", 
        params = "{sim_dir}/parameters_used.json"
    output:
        summary = "{sim_dir}/phylo_analysis/filtering/filtering_summary.csv",
        full_sequences = "{sim_dir}/phylo_analysis/filtering/full_sequences.fasta"
    params:
        output_dir = "{sim_dir}/phylo_analysis/filtering",
        rates_str = " ".join(map(str, DOWNSAMPLE_RATES))
    resources:
        mem_mb = 16000,
        runtime = 60,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    shell:
        """
        echo "Starting filtering and downsampling at $(date)"
        echo "Using seed: {PIPELINE_SEED}"
        
        mkdir -p {params.output_dir}
        cp {input.sequences} {output.full_sequences}
        
        python ../scripts/10_26_infection_filtering_downsampling.py \
            {wildcards.sim_dir} \
            {input.sequences} \
            {params.output_dir} \
            --downsample-rates {params.rates_str} \
            --n-replicates {N_REPLICATES} \
            --high-activity-subsample {HIGH_ACTIVITY_SUBSAMPLE} \
            --seed {PIPELINE_SEED}
        
        if [ ! -f {output.summary} ]; then
            echo "ERROR: Filtering summary not created!"
            exit 1
        fi
        
        DATASET_COUNT=$(tail -n +2 {output.summary} | wc -l)
        echo "Created $DATASET_COUNT filtered datasets"
        echo "Filtering and downsampling completed at $(date)"
        """

rule create_phylo_batches:
    input:
        summary = "{sim_dir}/phylo_analysis/filtering/filtering_summary.csv"
    output:
        batch_assignments = "{sim_dir}/phylo_analysis/batch_assignments.csv"
    run:
        import pandas as pd
        import os
        
        # Load filtering summary
        summary_df = pd.read_csv(input.summary)
        
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
            'fasta_file': f"{wildcards.sim_dir}/phylo_analysis/filtering/full_sequences.fasta",
            'infections_file': None,
            'batch_tier': 'tier1_full',
            'batch_size': 1,
            'mem_mb': 32000,
            'runtime': 480,
            'cpus_per_task': 8
        }
        
        batch_df = pd.DataFrame(batch_assignments)
        batch_df = pd.concat([pd.DataFrame([full_seq_row]), batch_df], ignore_index=True)
        
        # Assign actual batch numbers within tiers
        batch_df['batch_id'] = 0
        current_batch = 1
        
        for tier in batch_df['batch_tier'].unique():
            tier_data = batch_df[batch_df['batch_tier'] == tier]
            batch_size = tier_data.iloc[0]['batch_size']
            
            for i in range(0, len(tier_data), batch_size):
                batch_df.loc[tier_data.index[i:i+batch_size], 'batch_id'] = current_batch
                current_batch += 1
        
        batch_df.to_csv(output.batch_assignments, index=False)
        print(f"Created {current_batch-1} batches for phylogenetic analysis")

rule run_parallel_batch:
    input:
        batch_assignments = "{sim_dir}/phylo_analysis/batch_assignments.csv"
    output:
        iqtree_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_iqtree_complete.flag",
        bridging_27_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_sampled_bridging_27_complete.flag", 
        bridging_28_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_sampled_bridging_28_complete.flag",
        dates_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_dates_complete.flag",
        behaviors_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_behaviors_complete.flag"
    resources:
        mem_mb = get_mem_mb_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "mem_mb"), 64000)),
        runtime = get_runtime_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "runtime"), 2880)),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_compute_partition
    conda: "phylo"
    shell:
        """
        echo "Starting parallel batch {wildcards.batch_id} at $(date)"
        mkdir -p {wildcards.sim_dir}/phylo_analysis/flags
        
        python -c "
import pandas as pd
import os
from concurrent.futures import ThreadPoolExecutor

def process_dataset(row, sim_dir):
    dataset_id = row['dataset_id']
    fasta_file = row['fasta_file'] 
    infections_file = row['infections_file']
    
    print(f'Processing {{dataset_id}}...', flush=True)
    
    # Determine infections file to use
    if dataset_id == 'full_sequences':
        infections_to_use = os.path.join(sim_dir, 'transmission_df.csv')
        dataset_dir = os.path.dirname(fasta_file)
    else:
        infections_to_use = infections_file
        dataset_dir = os.path.dirname(fasta_file)
    
    jobs = []
    
    # 1. IQ-TREE
    if os.path.exists(fasta_file):
        cmd1 = f'iqtree -s {{fasta_file}} -m GTR+G --fast -seed {PIPELINE_SEED}'
        jobs.append(('iqtree', cmd1))
    
    # 2-6. Other analyses (now includes both 27 and 28 bridging methods)
    if os.path.exists(infections_to_use):
        dates_output = os.path.join(dataset_dir, 'dates.txt')
        cmd2 = f'python ../scripts/10_27_extract_tip_dates.py {{infections_to_use}} {{fasta_file}} {{sim_dir}} {{dates_output}}'
        jobs.append(('dates', cmd2))
        
        behaviors_prefix = os.path.join(dataset_dir, 'behaviors')
        cmd3 = f'python ../scripts/10_27_extract_tip_states.py {{infections_to_use}} {{fasta_file}} {{sim_dir}} {{behaviors_prefix}}'
        jobs.append(('behaviors', cmd3))
        
        # 10_27 bridging analysis (original)
        bridging_27_prefix = os.path.join(dataset_dir, 'sampled_bridging_27')
        cmd4 = f'python ../scripts/10_27_calculate_sampled_bridging.py {{infections_to_use}} {{sim_dir}} {{bridging_27_prefix}}'
        jobs.append(('bridging_27', cmd4))
        
        # 10_28 bridging analysis (directional)
        bridging_28_output = os.path.join(dataset_dir, 'sampled_bridging_28.csv')
        bridging_28_detailed = os.path.join(dataset_dir, 'sampled_bridging_28_detailed.csv')
        cmd5 = f'python ../scripts/10_28_calculate_sampled_bridging.py {{infections_to_use}} {{sim_dir}} --output {{bridging_28_output}} --detailed-output {{bridging_28_detailed}}'
        jobs.append(('bridging_28', cmd5))
    
    # Run jobs in parallel
    with ThreadPoolExecutor(max_workers=6) as executor:
        futures = [(job_type, executor.submit(os.system, cmd)) for job_type, cmd in jobs]
    
    for job_type, future in futures:
        result = future.result()
        if result == 0:
            print(f'SUCCESS: {{job_type}} completed for {{dataset_id}}')
        else:
            print(f'ERROR: {{job_type}} failed for {{dataset_id}}')

batch_df = pd.read_csv('{input.batch_assignments}')
batch_data = batch_df[batch_df['batch_id'] == {wildcards.batch_id}]

for _, row in batch_data.iterrows():
    process_dataset(row, '{wildcards.sim_dir}')
"

        echo "Parallel batch {wildcards.batch_id} completed at $(date)"
        touch {output.iqtree_flag} {output.bridging_27_flag} {output.bridging_28_flag} {output.dates_flag} {output.behaviors_flag}
        """

rule run_lsd2_batch:
    input:
        iqtree_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_iqtree_complete.flag",
        dates_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_dates_complete.flag",
        batch_assignments = "{sim_dir}/phylo_analysis/batch_assignments.csv"
    output:
        lsd2_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_lsd2_complete.flag"
    resources:
        mem_mb = get_mem_mb_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "mem_mb"), 64000)),
        runtime = get_runtime_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "runtime"), 2880)),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_compute_partition
    shell:
        """
        echo "Starting LSD2 batch {wildcards.batch_id} at $(date)"
        
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
        
        cmd = f'{LSD2_PATH} -i {{tree_file}} -d {{dates_file}} -s 10000 -r a -l 0'
        result = os.system(cmd)
        
        if result == 0:
            print(f'SUCCESS: LSD2 completed for {{dataset_id}}')
        else:
            print(f'ERROR: LSD2 failed for {{dataset_id}}')
    else:
        print(f'Missing files for {{dataset_id}}')
"
        
        echo "LSD2 batch {wildcards.batch_id} completed at $(date)"
        touch {output.lsd2_flag}
        """

rule run_treetime_msmw_msm_batch:
    input:
        lsd2_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_lsd2_complete.flag",
        behaviors_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_behaviors_complete.flag",
        batch_assignments = "{sim_dir}/phylo_analysis/batch_assignments.csv"
    output:
        treetime_msmw_msm_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_treetime_msmw_msm_complete.flag"
    resources:
        mem_mb = get_mem_mb_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "mem_mb"), 24000)),
        runtime = get_runtime_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "runtime"), 720)),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_long_partition
    conda: "phylo"
    shell:
        """
        echo "Starting TreeTime MSMW+MSM batch {wildcards.batch_id} at $(date)"
        
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
        
        if result == 0:
            print(f'SUCCESS: TreeTime MSMW+MSM completed for {{dataset_id}}')
        else:
            print(f'ERROR: TreeTime MSMW+MSM failed for {{dataset_id}}')
"
        
        echo "TreeTime MSMW+MSM batch {wildcards.batch_id} completed at $(date)"
        touch {output.treetime_msmw_msm_flag}
        """

rule run_treetime_msmw_msw_batch:
    input:
        lsd2_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_lsd2_complete.flag",
        behaviors_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_behaviors_complete.flag",
        batch_assignments = "{sim_dir}/phylo_analysis/batch_assignments.csv"
    output:
        treetime_msmw_msw_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_treetime_msmw_msw_complete.flag"
    resources:
        mem_mb = get_mem_mb_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "mem_mb"), 24000)),
        runtime = get_runtime_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "runtime"), 720)),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_long_partition
    conda: "phylo"
    shell:
        """
        echo "Starting TreeTime MSMW+MSW batch {wildcards.batch_id} at $(date)"
        
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
        
        if result == 0:
            print(f'SUCCESS: TreeTime MSMW+MSW completed for {{dataset_id}}')
        else:
            print(f'ERROR: TreeTime MSMW+MSW failed for {{dataset_id}}')
"
        
        echo "TreeTime MSMW+MSW batch {wildcards.batch_id} completed at $(date)"
        touch {output.treetime_msmw_msw_flag}
        """

rule run_tree_bridging_27_batch:
    input:
        treetime_msmw_msm_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_treetime_msmw_msm_complete.flag",
        treetime_msmw_msw_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_treetime_msmw_msw_complete.flag",
        batch_assignments = "{sim_dir}/phylo_analysis/batch_assignments.csv"
    output:
        tree_bridging_27_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_tree_bridging_27_complete.flag"
    resources:
        mem_mb = get_mem_mb_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "mem_mb"), 24000)),
        runtime = get_runtime_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "runtime"), 720)),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_compute_partition
    conda: "phylo"
    shell:
        """
        echo "Starting Tree Bridging Analysis 27 batch {wildcards.batch_id} at $(date)"
        
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
    bridging_prefix = os.path.join(dataset_dir, 'tree_bridging_27')
    
    if os.path.exists(msmw_msm_nexus) and os.path.exists(msmw_msw_nexus) and os.path.exists(lsd2_result):
        print(f'Running Tree Bridging Analysis 27 for {{dataset_id}}...', flush=True)
        
        cmd = f'python ../scripts/10_27_treetime_bridging.py {{msmw_msm_nexus}} {{msmw_msw_nexus}} {{lsd2_result}} {{bridging_prefix}}'
        result = os.system(cmd)
        
        if result == 0:
            print(f'SUCCESS: Tree Bridging Analysis 27 completed for {{dataset_id}}')
        else:
            print(f'ERROR: Tree Bridging Analysis 27 failed for {{dataset_id}}')
"
        
        echo "Tree Bridging Analysis 27 batch {wildcards.batch_id} completed at $(date)"
        touch {output.tree_bridging_27_flag}
        """

rule run_tree_bridging_28_batch:
    input:
        treetime_msmw_msm_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_treetime_msmw_msm_complete.flag",
        treetime_msmw_msw_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_treetime_msmw_msw_complete.flag",
        batch_assignments = "{sim_dir}/phylo_analysis/batch_assignments.csv"
    output:
        tree_bridging_28_flag = "{sim_dir}/phylo_analysis/flags/batch_{batch_id}_tree_bridging_28_complete.flag"
    resources:
        mem_mb = get_mem_mb_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "mem_mb"), 24000)),
        runtime = get_runtime_scaling(lambda wildcards, input: max(get_batch_resources(wildcards, input, "runtime"), 720)),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_compute_partition
    conda: "phylo"
    shell:
        """
        echo "Starting Tree Bridging Analysis 28 batch {wildcards.batch_id} at $(date)"
        
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
    bridging_output = os.path.join(dataset_dir, 'tree_bridging_28.csv')
    bridging_detailed = os.path.join(dataset_dir, 'tree_bridging_28_detailed.csv')
    
    if os.path.exists(msmw_msm_nexus) and os.path.exists(msmw_msw_nexus) and os.path.exists(lsd2_result):
        print(f'Running Tree Bridging Analysis 28 for {{dataset_id}}...', flush=True)
        
        cmd = f'python ../scripts/10_28_treetime_bridging.py {{msmw_msm_nexus}} {{msmw_msw_nexus}} {{lsd2_result}} --output {{bridging_output}} --detailed-output {{bridging_detailed}}'
        result = os.system(cmd)
        
        if result == 0:
            print(f'SUCCESS: Tree Bridging Analysis 28 completed for {{dataset_id}}')
        else:
            print(f'ERROR: Tree Bridging Analysis 28 failed for {{dataset_id}}')
"
        
        echo "Tree Bridging Analysis 28 batch {wildcards.batch_id} completed at $(date)"
        touch {output.tree_bridging_28_flag}
        """

rule aggregate_phylogenetic_bridging_27:
    input:
        lambda wildcards: expand("{sim_dir}/phylo_analysis/flags/batch_{batch_id}_tree_bridging_27_complete.flag",
                                sim_dir=wildcards.sim_dir, 
                                batch_id=[i for i in range(1, MAX_BATCH_ID + 1)]),
        batch_assignments = "{sim_dir}/phylo_analysis/batch_assignments.csv"
    output:
        aggregated_bridging_27 = "{sim_dir}/phylo_analysis/aggregated_phylogenetic_bridging_27_results.csv"
    resources:
        mem_mb = 4000,
        runtime = 30,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    run:
        import pandas as pd
        import os
        
        print("🌳 Aggregating phylogenetic bridging 27 results...")
        
        batch_df = pd.read_csv(input.batch_assignments)
        all_results = []
        
        for _, row in batch_df.iterrows():
            dataset_id = row['dataset_id']
            population = row['population']
            dataset = row['dataset']
            fasta_file = row['fasta_file']
            
            dataset_dir = os.path.dirname(fasta_file)
            bridging_summary_file = os.path.join(dataset_dir, 'tree_bridging_27_bridging_comparison.csv')
            
            if os.path.exists(bridging_summary_file):
                try:
                    bridging_df = pd.read_csv(bridging_summary_file)
                    
                    for _, bridging_row in bridging_df.iterrows():
                        result_row = {
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
                                result_row[f'phylo_27_{col}'] = bridging_row[col]
                        
                        all_results.append(result_row)
                    
                    print(f"✅ Collected bridging 27 results from {dataset_id}")
                except Exception as e:
                    print(f"❌ Error processing {bridging_summary_file}: {e}")
        
        if all_results:
            aggregated_df = pd.DataFrame(all_results)
            column_order = [
                'dataset_id', 'population', 'dataset_type', 'dataset', 
                'downsample_rate', 'replicate', 'measure'
            ] + [col for col in aggregated_df.columns if col.startswith('phylo_27_')]
            
            aggregated_df = aggregated_df[column_order]
            aggregated_df = aggregated_df.sort_values(['population', 'downsample_rate', 'replicate', 'measure'])
            aggregated_df.to_csv(output.aggregated_bridging_27, index=False)
            
            print(f"🎉 Phylogenetic bridging 27 aggregation complete! ({len(aggregated_df)} records)")
        else:
            pd.DataFrame().to_csv(output.aggregated_bridging_27, index=False)

rule aggregate_phylogenetic_bridging_28:
    input:
        lambda wildcards: expand("{sim_dir}/phylo_analysis/flags/batch_{batch_id}_tree_bridging_28_complete.flag",
                                sim_dir=wildcards.sim_dir, 
                                batch_id=[i for i in range(1, MAX_BATCH_ID + 1)]),
        batch_assignments = "{sim_dir}/phylo_analysis/batch_assignments.csv"
    output:
        aggregated_bridging_28 = "{sim_dir}/phylo_analysis/aggregated_phylogenetic_bridging_28_results.csv"
    resources:
        mem_mb = 4000,
        runtime = 30,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    run:
        import pandas as pd
        import os
        
        print("🌳 Aggregating phylogenetic bridging 28 results...")
        
        batch_df = pd.read_csv(input.batch_assignments)
        all_results = []
        
        for _, row in batch_df.iterrows():
            dataset_id = row['dataset_id']
            population = row['population']
            dataset = row['dataset']
            fasta_file = row['fasta_file']
            
            dataset_dir = os.path.dirname(fasta_file)
            bridging_summary_file = os.path.join(dataset_dir, 'tree_bridging_28.csv')
            
            if os.path.exists(bridging_summary_file):
                try:
                    bridging_df = pd.read_csv(bridging_summary_file, index_col=0)
                    bridging_df = bridging_df.reset_index()  # Convert pivot back to long format
                    
                    for _, bridging_row in bridging_df.iterrows():
                        result_row = {
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
                                result_row[f'phylo_28_{col}'] = bridging_row[col]
                        
                        all_results.append(result_row)
                    
                    print(f"✅ Collected bridging 28 results from {dataset_id}")
                except Exception as e:
                    print(f"❌ Error processing {bridging_summary_file}: {e}")
        
        if all_results:
            aggregated_df = pd.DataFrame(all_results)
            column_order = [
                'dataset_id', 'population', 'dataset_type', 'dataset', 
                'downsample_rate', 'replicate', 'measure'
            ] + [col for col in aggregated_df.columns if col.startswith('phylo_28_')]
            
            aggregated_df = aggregated_df[column_order]
            aggregated_df = aggregated_df.sort_values(['population', 'downsample_rate', 'replicate', 'measure'])
            aggregated_df.to_csv(output.aggregated_bridging_28, index=False)
            
            print(f"🎉 Phylogenetic bridging 28 aggregation complete! ({len(aggregated_df)} records)")
        else:
            pd.DataFrame().to_csv(output.aggregated_bridging_28, index=False)

rule aggregate_sampled_bridging_27:
    input:
        lambda wildcards: expand("{sim_dir}/phylo_analysis/flags/batch_{batch_id}_sampled_bridging_27_complete.flag",
                                sim_dir=wildcards.sim_dir,
                                batch_id=[i for i in range(1, MAX_BATCH_ID + 1)]),
        batch_assignments = "{sim_dir}/phylo_analysis/batch_assignments.csv"
    output:
        aggregated_sampled_bridging_27 = "{sim_dir}/phylo_analysis/aggregated_sampled_bridging_27_results.csv"
    resources:
        mem_mb = 4000,
        runtime = 30,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    run:
        import pandas as pd
        import os
        
        print("📊 Aggregating sampled bridging 27 results...")
        
        batch_df = pd.read_csv(input.batch_assignments)
        all_results = []
        
        for _, row in batch_df.iterrows():
            dataset_id = row['dataset_id']
            population = row['population']
            dataset = row['dataset']
            fasta_file = row['fasta_file']
            
            dataset_dir = os.path.dirname(fasta_file)
            sampled_bridging_file = os.path.join(dataset_dir, 'sampled_bridging_27_sampled_bridging_comparison.csv')
            
            if os.path.exists(sampled_bridging_file):
                try:
                    sampled_df = pd.read_csv(sampled_bridging_file)
                    
                    for _, sampled_row in sampled_df.iterrows():
                        result_row = {
                            'dataset_id': dataset_id,
                            'population': population,
                            'dataset_type': 'complete' if dataset == 'complete' else 'downsampled',
                            'dataset': dataset,
                            'measure': sampled_row['measure']
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
                        
                        for col in sampled_df.columns:
                            if col != 'measure':
                                result_row[f'sampled_27_{col}'] = sampled_row[col]
                        
                        all_results.append(result_row)
                    
                    print(f"✅ Collected sampled bridging 27 from {dataset_id}")
                except Exception as e:
                    print(f"❌ Error processing {sampled_bridging_file}: {e}")
        
        if all_results:
            aggregated_df = pd.DataFrame(all_results)
            
            column_order = [
                'dataset_id', 'population', 'dataset_type', 'dataset', 
                'downsample_rate', 'replicate', 'measure'
            ] + [col for col in aggregated_df.columns if col.startswith('sampled_27_')]
            
            aggregated_df = aggregated_df[column_order]
            aggregated_df = aggregated_df.sort_values(['population', 'downsample_rate', 'replicate'])
            aggregated_df.to_csv(output.aggregated_sampled_bridging_27, index=False)
            
            print(f"🎉 Sampled bridging 27 aggregation complete! ({len(aggregated_df)} records)")
        else:
            pd.DataFrame().to_csv(output.aggregated_sampled_bridging_27, index=False)
            print("❌ No sampled bridging 27 results found!")

rule aggregate_sampled_bridging_28:
    input:
        lambda wildcards: expand("{sim_dir}/phylo_analysis/flags/batch_{batch_id}_sampled_bridging_28_complete.flag",
                                sim_dir=wildcards.sim_dir,
                                batch_id=[i for i in range(1, MAX_BATCH_ID + 1)]),
        batch_assignments = "{sim_dir}/phylo_analysis/batch_assignments.csv"
    output:
        aggregated_sampled_bridging_28 = "{sim_dir}/phylo_analysis/aggregated_sampled_bridging_28_results.csv"
    resources:
        mem_mb = 4000,
        runtime = 30,
        cpus_per_task = 1,
        slurm_partition = get_quick_partition
    run:
        import pandas as pd
        import os
        
        print("📊 Aggregating sampled bridging 28 results...")
        
        batch_df = pd.read_csv(input.batch_assignments)
        all_results = []
        
        for _, row in batch_df.iterrows():
            dataset_id = row['dataset_id']
            population = row['population']
            dataset = row['dataset']
            fasta_file = row['fasta_file']
            
            dataset_dir = os.path.dirname(fasta_file)
            sampled_bridging_file = os.path.join(dataset_dir, 'sampled_bridging_28.csv')
            
            if os.path.exists(sampled_bridging_file):
                try:
                    sampled_df = pd.read_csv(sampled_bridging_file, index_col=0)
                    sampled_df = sampled_df.reset_index()  # Convert pivot back to long format
                    
                    for _, sampled_row in sampled_df.iterrows():
                        result_row = {
                            'dataset_id': dataset_id,
                            'population': population,
                            'dataset_type': 'complete' if dataset == 'complete' else 'downsampled',
                            'dataset': dataset,
                            'measure': sampled_row['measure']
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
                        
                        for col in sampled_df.columns:
                            if col != 'measure':
                                result_row[f'sampled_28_{col}'] = sampled_row[col]
                        
                        all_results.append(result_row)
                    
                    print(f"✅ Collected sampled bridging 28 from {dataset_id}")
                except Exception as e:
                    print(f"❌ Error processing {sampled_bridging_file}: {e}")
        
        if all_results:
            aggregated_df = pd.DataFrame(all_results)
            
            column_order = [
                'dataset_id', 'population', 'dataset_type', 'dataset', 
                'downsample_rate', 'replicate', 'measure'
            ] + [col for col in aggregated_df.columns if col.startswith('sampled_28_')]
            
            aggregated_df = aggregated_df[column_order]
            aggregated_df = aggregated_df.sort_values(['population', 'downsample_rate', 'replicate'])
            aggregated_df.to_csv(output.aggregated_sampled_bridging_28, index=False)
            
            print(f"🎉 Sampled bridging 28 aggregation complete! ({len(aggregated_df)} records)")
        else:
            pd.DataFrame().to_csv(output.aggregated_sampled_bridging_28, index=False)
            print("❌ No sampled bridging 28 results found!")
