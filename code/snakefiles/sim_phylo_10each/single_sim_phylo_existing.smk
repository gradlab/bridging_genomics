# ============================================================================
# Single Simulation Phylogenetic Bridging Pipeline (Uses Existing Structure)
# ============================================================================
#
# USAGE: Run on ONE simulation using existing phylo structure
# snakemake --snakefile single_sim_phylo_existing.smk \
#           --config sim_id="0.5_060" \
#           --cores 16 --jobs 10 --executor slurm --default-resources slurm_account=grad_lab
#
# ============================================================================

# Pipeline parameters
PIPELINE_SEED = 42
OUTPUT_BASE = "../../output/baseline_phylo_bridging_complete_fresh_v2"

# Software paths
LSD2_PATH = "/n/holylfs05/LABS/grad_lab/Users/mkline/bridging_sims/software/lsd2/src/lsd2"

# Resource limits
MAX_MEM_MB = 128000
MAX_RUNTIME = 2880
MAX_CPUS = 16
MAX_BATCH_ID = 70

# Partition functions
def get_compute_partition(wildcards, attempt):
    partitions = ["hsph", "sapphire", "intermediate", "shared"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

def get_long_partition(wildcards, attempt):
    partitions = ["intermediate", "hsph", "sapphire", "shared"]
    return partitions[min(attempt - 1, len(partitions) - 1)]

def get_batch_resources(wildcards, input, resource_type):
    import pandas as pd
    batch_df = pd.read_csv(input.batch_assignments)
    batch_data = batch_df[batch_df['batch_id'] == int(wildcards.batch_id)]
    if len(batch_data) > 0:
        value = batch_data.iloc[0][resource_type]
        return int(value)
    return 4000 if resource_type == 'mem_mb' else (60 if resource_type == 'runtime' else 1)

# ============================================================================
# Pipeline Rules (Using EXISTING directory structure)
# ============================================================================

rule all:
    input:
        f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_sampled_bridging_30_results.csv".format(sim_id=config["sim_id"]),
        f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_phylogenetic_bridging_30_results.csv".format(sim_id=config["sim_id"])

# Updated parallel batch rule - uses EXISTING paths
rule run_parallel_batch:
    input:
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv"
    output:
        iqtree_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_iqtree_complete.flag",
        sampled_bridging_30_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_sampled_bridging_30_complete.flag",
        dates_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_dates_complete.flag",
        behaviors_flag = f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_behaviors_complete.flag"
    resources:
        mem_mb = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "mem_mb"), 4000), MAX_MEM_MB),
        runtime = lambda wildcards, input: min(max(get_batch_resources(wildcards, input, "runtime"), 60), MAX_RUNTIME),
        cpus_per_task = lambda wildcards, input: min(get_batch_resources(wildcards, input, "cpus_per_task"), MAX_CPUS),
        slurm_partition = get_compute_partition
    conda: "phylo"
    shell:
        """
        echo "Starting parallel batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"
        mkdir -p {OUTPUT_BASE}/phylo/{wildcards.sim_id}/flags

        # Get actual simulation directory
        SIM_DIR="{OUTPUT_BASE}/simulations/{wildcards.sim_id}/sim"

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

        # 4. UPDATED sampled bridging analysis (10_30 script)
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

    # Clean up temporary file
    if dataset_id == 'full_sequences' and 'temp_filtered_transmissions.csv' in infections_to_use:
        try:
            os.remove(infections_to_use)
        except:
            pass

batch_df = pd.read_csv('{input.batch_assignments}')
batch_data = batch_df[batch_df['batch_id'] == {wildcards.batch_id}]

for _, row in batch_data.iterrows():
    process_dataset(row, '$SIM_DIR', '{wildcards.sim_id}')
"

        echo "Parallel batch {wildcards.batch_id} completed at $(date)"
        touch {output.iqtree_flag} {output.sampled_bridging_30_flag} {output.dates_flag} {output.behaviors_flag}
        """

# Add these rules to your single_sim_phylo_existing.smk file:

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
        echo "Starting LSD2 batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"

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
        echo "Starting TreeTime MSMW+MSM batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"

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
        echo "Starting TreeTime MSMW+MSW batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"

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
        slurm_partition = get_compute_partition
    conda: "phylo"
    shell:
        """
        echo "Starting Tree Bridging 30 batch {wildcards.batch_id} for {wildcards.sim_id} at $(date)"

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

rule aggregate_sampled_bridging_30:
    input:
        lambda wildcards: expand(f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_sampled_bridging_30_complete.flag",
                                sim_id=wildcards.sim_id,
                                batch_id=[i for i in range(1, MAX_BATCH_ID + 1)]),
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv"
    output:
        aggregated_sampled_bridging_30 = f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_sampled_bridging_30_results.csv"
    resources:
        mem_mb = 4000,
        runtime = 30,
        cpus_per_task = 1,
        slurm_partition = get_compute_partition
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
            if dataset == 'full':
                continue

            dataset_dir = os.path.dirname(fasta_file)
            sampled_bridging_file = os.path.join(dataset_dir, 'sampled_bridging_30.csv')

            if os.path.exists(sampled_bridging_file) and os.path.getsize(sampled_bridging_file) > 0:
                try:
                    sampled_df = pd.read_csv(sampled_bridging_file, index_col=0)
                    sampled_df = sampled_df.reset_index()

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
                'dataset_id', 'population', 'dataset_type', 'dataset',
                'downsample_rate', 'replicate', 'measure'
            ] + [col for col in aggregated_df.columns if col.startswith('sampled_30_')]

            aggregated_df = aggregated_df[column_order]
            aggregated_df = aggregated_df.sort_values(['population', 'downsample_rate', 'replicate', 'measure'])
            aggregated_df.to_csv(output.aggregated_sampled_bridging_30, index=False)

            print(f"🎉 Sampled bridging 30 aggregation complete! ({len(aggregated_df)} records)")
        else:
            pd.DataFrame().to_csv(output.aggregated_sampled_bridging_30, index=False)

rule aggregate_phylogenetic_bridging_30:
    input:
        lambda wildcards: expand(f"{OUTPUT_BASE}/phylo/{{sim_id}}/flags/batch_{{batch_id}}_tree_bridging_30_complete.flag",
                                sim_id=wildcards.sim_id,
                                batch_id=[i for i in range(1, MAX_BATCH_ID + 1)]),
        batch_assignments = f"{OUTPUT_BASE}/phylo/{{sim_id}}/batch_assignments.csv"
    output:
        aggregated_phylo_bridging_30 = f"{OUTPUT_BASE}/phylo/{{sim_id}}/aggregated_phylogenetic_bridging_30_results.csv"
    resources:
        mem_mb = 4000,
        runtime = 30,
        cpus_per_task = 1,
        slurm_partition = get_compute_partition
    run:
        import pandas as pd
        import os

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
                'dataset_id', 'population', 'dataset_type', 'dataset',
                'downsample_rate', 'replicate', 'measure'
            ] + [col for col in aggregated_df.columns if col.startswith('phylo_30_')]

            aggregated_df = aggregated_df[column_order]
            aggregated_df = aggregated_df.sort_values(['population', 'downsample_rate', 'replicate', 'measure'])
            aggregated_df.to_csv(output.aggregated_phylo_bridging_30, index=False)

            print(f"🎉 Phylogenetic bridging 30 aggregation complete! ({len(aggregated_df)} records)")
        else:
            pd.DataFrame().to_csv(output.aggregated_phylo_bridging_30, index=False)
