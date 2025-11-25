import json
import os
import glob
from pathlib import Path
import datetime
import math

# CONFIGURABLE PARAMETERS
STARTING_SEED = 1000  # CHANGE THIS FOR EACH OF YOUR 4 JOBS (1, 1000, 2000, 3000)
P_MSMW_W_VALUES = [0.05, 0.25, 0.5, 0.75, 0.95]
N_REPS = 250  # UPDATED: 250 sims per p_msmw_w value
BASE_PARAM_FILE = "../parameters/11_3_base_params.json"  # UPDATED: Use baseline params

# Batching configuration
SIMULATION_BATCH_SIZE = 10
ANALYSIS_BATCH_SIZE = 10

# Updated timestamp with starting seed
TIMESTAMP = f"baseline_250each_seed_{STARTING_SEED}"
SWEEP_OUTPUT_DIR = f"../../output/param_sweep_{TIMESTAMP}"

# Generate all parameter combinations
PARAM_COMBINATIONS = []
for i, (p_val, rep) in enumerate([(p, rep) for p in P_MSMW_W_VALUES for rep in range(1, N_REPS+1)]):
    PARAM_COMBINATIONS.append({
        'id': i,
        'p_msmw_w': p_val,
        'rep': f"{rep:03d}",
        'seed_offset': rep - 1
    })

# Calculate batches
NUM_SIM_BATCHES = math.ceil(len(PARAM_COMBINATIONS) / SIMULATION_BATCH_SIZE)
NUM_ANALYSIS_BATCHES = math.ceil(len(PARAM_COMBINATIONS) / ANALYSIS_BATCH_SIZE)

# DYNAMIC RESOURCE FUNCTIONS with multi-partition support
def get_mem_mb_simulation(wildcards, attempt):
    base_mem = 15000
    return min(attempt * base_mem, 32000)

def get_mem_mb_analysis(wildcards, attempt):
    base_mem = 15000
    return min(attempt * base_mem, 24000)

def get_runtime_scaling(base_minutes):
    return lambda wildcards, attempt: min(attempt * base_minutes, 720)

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

rule all:
    input:
        f"{SWEEP_OUTPUT_DIR}/aggregated/bridging_values.csv",
        f"{SWEEP_OUTPUT_DIR}/aggregated/bridging_counts.csv",
        f"{SWEEP_OUTPUT_DIR}/aggregated/bridging_totals.csv",
        f"{SWEEP_OUTPUT_DIR}/aggregated/prevalence_timeseries.csv",
        f"{SWEEP_OUTPUT_DIR}/aggregated/prevalence_averages.csv"

# 1. Generate all parameter files with baseline params and starting seed
rule generate_all_params:
    input:
        BASE_PARAM_FILE
    output:
        expand(f"{SWEEP_OUTPUT_DIR}/params/params_{{p_msmw_w}}_{{rep}}.json",
               p_msmw_w=P_MSMW_W_VALUES, rep=[f"{r:03d}" for r in range(1, N_REPS+1)]),
        f"{SWEEP_OUTPUT_DIR}/batch_mapping.txt"
    run:
        import copy

        # Load baseline parameters
        with open(input[0], 'r') as f:
            base_params = json.load(f)

        base_seed = base_params["simulation"]["seed"]
        os.makedirs(f"{SWEEP_OUTPUT_DIR}/params", exist_ok=True)

        # Generate parameter files
        for combo in PARAM_COMBINATIONS:
            params = copy.deepcopy(base_params)
            params["partnerships"]["p_MSMW_W"] = combo['p_msmw_w']
            # UPDATED: Use starting seed + base seed + offset for unique seeds across jobs
            params["simulation"]["seed"] = base_seed + STARTING_SEED + combo['seed_offset']
            params["output"]["output_dir"] = f"{SWEEP_OUTPUT_DIR}/{combo['p_msmw_w']}"
            params["output"]["run_name"] = f"rep_{combo['rep']}"

            param_file = f"{SWEEP_OUTPUT_DIR}/params/params_{combo['p_msmw_w']}_{combo['rep']}.json"
            with open(param_file, 'w') as f:
                json.dump(params, f, indent=2)

        # Create batch mapping file
        with open(f"{SWEEP_OUTPUT_DIR}/batch_mapping.txt", 'w') as f:
            for combo in PARAM_COMBINATIONS:
                f.write(f"{combo['id']},{combo['p_msmw_w']},{combo['rep']}\n")

# 2. Batched simulation jobs (unchanged from calibrated version)
def get_simulation_batch_params(wildcards):
    """Get parameter files for this simulation batch"""
    batch_id = int(wildcards.batch_id)
    start_idx = (batch_id - 1) * SIMULATION_BATCH_SIZE
    end_idx = min(start_idx + SIMULATION_BATCH_SIZE, len(PARAM_COMBINATIONS))

    param_files = []
    for i in range(start_idx, end_idx):
        combo = PARAM_COMBINATIONS[i]
        param_files.append(f"{SWEEP_OUTPUT_DIR}/params/params_{combo['p_msmw_w']}_{combo['rep']}.json")

    return param_files

rule run_simulation_batch:
    input:
        param_files = get_simulation_batch_params,
        batch_mapping = f"{SWEEP_OUTPUT_DIR}/batch_mapping.txt"
    output:
        f"{SWEEP_OUTPUT_DIR}/batches/sim_batch_{{batch_id}}_complete.flag"
    conda: "phylo"
    resources:
        mem_mb=get_mem_mb_simulation,
        runtime=get_runtime_scaling(600),
        cpus_per_task=2,
        slurm_partition=get_compute_partition
    shell:
        """
        echo "Simulation batch {wildcards.batch_id}: Running {SIMULATION_BATCH_SIZE} simulations"

        # Create batch output directory
        mkdir -p {SWEEP_OUTPUT_DIR}/batches

        # Run each simulation in this batch
        for param_file in {input.param_files}; do
            echo "Running simulation with: $param_file"
            python ../scripts/10_16_core_sim.py "$param_file"
        done

        echo "Batch {wildcards.batch_id} simulations complete"
        touch {output}
        """

# 3. Batched lineage tests (unchanged)
def get_lineage_batch_inputs(wildcards):
    """Get simulation outputs for this lineage batch"""
    batch_id = int(wildcards.batch_id)
    start_idx = (batch_id - 1) * ANALYSIS_BATCH_SIZE
    end_idx = min(start_idx + ANALYSIS_BATCH_SIZE, len(PARAM_COMBINATIONS))

    # Find which sim batches we depend on
    sim_batch_flags = set()
    for i in range(start_idx, end_idx):
        sim_batch_id = (i // SIMULATION_BATCH_SIZE) + 1
        sim_batch_flags.add(f"{SWEEP_OUTPUT_DIR}/batches/sim_batch_{sim_batch_id}_complete.flag")

    return list(sim_batch_flags)

rule run_lineage_batch:
    input:
        sim_batches = get_lineage_batch_inputs
    output:
        f"{SWEEP_OUTPUT_DIR}/batches/lineage_batch_{{batch_id}}_complete.flag"
    conda: "phylo"
    resources:
        mem_mb=6000,
        runtime=60,
        cpus_per_task=1,
        slurm_partition=get_quick_partition
    shell:
        """
        echo "Lineage batch {wildcards.batch_id}: Processing lineage tests"

        # Calculate which parameter combinations this batch should process
        BATCH_ID={wildcards.batch_id}
        BATCH_SIZE={ANALYSIS_BATCH_SIZE}
        START_IDX=$(( (BATCH_ID - 1) * BATCH_SIZE ))

        # Read batch mapping and process the relevant combinations
        tail -n +$((START_IDX + 1)) {SWEEP_OUTPUT_DIR}/batch_mapping.txt | head -n $BATCH_SIZE | while IFS=, read -r combo_id p_val rep; do
            OUTPUT_DIR="{SWEEP_OUTPUT_DIR}/$p_val/rep_$rep"
            echo "Testing lineage for: $OUTPUT_DIR"
            python ../scripts/10_17_single_lineage_test_script.py "$OUTPUT_DIR"
        done

        echo "Lineage batch {wildcards.batch_id} complete"
        touch {output}
        """

# 4. Analysis batch with proven 10_30 bridging script
rule run_analysis_batch:
    input:
        lineage_batch = f"{SWEEP_OUTPUT_DIR}/batches/lineage_batch_{{batch_id}}_complete.flag"
    output:
        f"{SWEEP_OUTPUT_DIR}/batches/analysis_batch_{{batch_id}}_complete.flag"
    conda: "phylo"
    resources:
        mem_mb=get_mem_mb_analysis,
        runtime=get_runtime_scaling(120),
        cpus_per_task=1,
        slurm_partition=get_long_partition
    shell:
        """
        echo "Analysis batch {wildcards.batch_id}: Processing {ANALYSIS_BATCH_SIZE} analyses"

        # Calculate which parameter combinations this batch should process
        BATCH_ID={wildcards.batch_id}
        BATCH_SIZE={ANALYSIS_BATCH_SIZE}
        START_IDX=$(( (BATCH_ID - 1) * BATCH_SIZE ))

        # Read batch mapping and process the relevant combinations
        tail -n +$((START_IDX + 1)) {SWEEP_OUTPUT_DIR}/batch_mapping.txt | head -n $BATCH_SIZE | while IFS=, read -r combo_id p_val rep; do
            OUTPUT_DIR="{SWEEP_OUTPUT_DIR}/$p_val/rep_$rep"
            LINEAGE_FILE="$OUTPUT_DIR/lineage_test_results.json"

            echo "Processing analysis for: $OUTPUT_DIR"

            # Check if single lineage
            SINGLE_LINEAGE=$(python -c "import json; data=json.load(open('$LINEAGE_FILE')); print(data.get('success', False) and data.get('single_lineage', False))")

            if [ "$SINGLE_LINEAGE" = "True" ]; then
                echo "Running analysis for single lineage simulation: $OUTPUT_DIR"

                # Behavior classification
                python ../scripts/10_17_infection_behavior_classification_analysis.py "$OUTPUT_DIR" --both

                # Proven 10_30 bridging analysis
                python ../scripts/10_30_ground_truth_bridging.py \
                    "$OUTPUT_DIR" \
                    "$OUTPUT_DIR/infection_behavior_classifications.csv" \
                    --output "$OUTPUT_DIR/bridging_analysis.csv" \
                    --detailed-output "$OUTPUT_DIR/bridging_detailed.csv"
            else
                echo "Skipping analysis - epidemic died out or multiple lineages: $OUTPUT_DIR"
                touch "$OUTPUT_DIR/bridging_analysis.csv"
            fi
        done

        echo "Analysis batch {wildcards.batch_id} complete"
        touch {output}
        """

# 5. Aggregation for 10_30 bridging output structure (values, counts, totals)
rule aggregate_bridging:
    input:
        analysis_batches = expand(f"{SWEEP_OUTPUT_DIR}/batches/analysis_batch_{{batch_id}}_complete.flag",
                                 batch_id=range(1, NUM_ANALYSIS_BATCHES + 1))
    output:
        values = f"{SWEEP_OUTPUT_DIR}/aggregated/bridging_values.csv",
        counts = f"{SWEEP_OUTPUT_DIR}/aggregated/bridging_counts.csv",
        totals = f"{SWEEP_OUTPUT_DIR}/aggregated/bridging_totals.csv"
    conda: "phylo"
    resources:
        mem_mb=12000,
        runtime=120,
        cpus_per_task=1,
        slurm_partition=get_compute_partition
    run:
        import pandas as pd

        os.makedirs(f"{SWEEP_OUTPUT_DIR}/aggregated", exist_ok=True)

        all_values = []
        all_counts = []
        all_totals = []

        for combo in PARAM_COMBINATIONS:
            p_val = combo['p_msmw_w']
            rep = int(combo['rep'])

            # Check for 10_30 script outputs
            values_file = f"{SWEEP_OUTPUT_DIR}/{p_val}/rep_{combo['rep']}/bridging_analysis.csv"
            counts_file = f"{SWEEP_OUTPUT_DIR}/{p_val}/rep_{combo['rep']}/bridging_analysis_counts.csv"
            totals_file = f"{SWEEP_OUTPUT_DIR}/{p_val}/rep_{combo['rep']}/bridging_analysis_totals.csv"

            if os.path.exists(values_file) and os.path.getsize(values_file) > 0:
                try:
                    # Load values (pivot table format - has measure as index)
                    values_df = pd.read_csv(values_file, index_col=0)
                    # Reset index to make measure a column, then add identifiers
                    values_df = values_df.reset_index()
                    values_df['p_msmw_w'] = p_val
                    values_df['replicate'] = rep
                    all_values.append(values_df)

                    # Load counts if exists
                    if os.path.exists(counts_file):
                        counts_df = pd.read_csv(counts_file, index_col=0)
                        counts_df = counts_df.reset_index()
                        counts_df['p_msmw_w'] = p_val
                        counts_df['replicate'] = rep
                        all_counts.append(counts_df)

                    # Load totals if exists
                    if os.path.exists(totals_file):
                        totals_df = pd.read_csv(totals_file, index_col=0)
                        totals_df = totals_df.reset_index()
                        totals_df['p_msmw_w'] = p_val
                        totals_df['replicate'] = rep
                        all_totals.append(totals_df)

                    print(f"Added results for p_msmw_w={p_val}, rep={rep}")
                except Exception as e:
                    print(f"Failed to read {values_file}: {e}")
                    continue

        # Save aggregated results
        if all_values:
            combined_values = pd.concat(all_values, ignore_index=True)
            combined_values.to_csv(output.values, index=False)
            print(f"Combined values from {len(all_values)} successful simulations")
        else:
            pd.DataFrame().to_csv(output.values, index=False)
            print("No successful simulations found for values")

        if all_counts:
            combined_counts = pd.concat(all_counts, ignore_index=True)
            combined_counts.to_csv(output.counts, index=False)
            print(f"Combined counts from {len(all_counts)} successful simulations")
        else:
            pd.DataFrame().to_csv(output.counts, index=False)
            print("No counts files found")

        if all_totals:
            combined_totals = pd.concat(all_totals, ignore_index=True)
            combined_totals.to_csv(output.totals, index=False)
            print(f"Combined totals from {len(all_totals)} successful simulations")
        else:
            pd.DataFrame().to_csv(output.totals, index=False)
            print("No totals files found")

# 6. Prevalence aggregation rules (unchanged)
rule aggregate_prevalence_timeseries:
    input:
        expand(f"{SWEEP_OUTPUT_DIR}/batches/analysis_batch_{{batch_id}}_complete.flag",
               batch_id=range(1, NUM_ANALYSIS_BATCHES + 1))
    output:
        f"{SWEEP_OUTPUT_DIR}/aggregated/prevalence_timeseries.csv"
    conda: "phylo"
    resources:
        mem_mb=12000,
        runtime=30,
        cpus_per_task=1,
        slurm_partition=get_quick_partition
    run:
        import pandas as pd

        all_prevalence = []
        for combo in PARAM_COMBINATIONS:
            p_val = combo['p_msmw_w']
            rep = int(combo['rep'])

            # Check if analysis was successful
            bridging_file = f"{SWEEP_OUTPUT_DIR}/{p_val}/rep_{combo['rep']}/bridging_analysis.csv"
            if os.path.exists(bridging_file) and os.path.getsize(bridging_file) > 0:
                # Look for prevalence time series file
                prev_file = f"{SWEEP_OUTPUT_DIR}/{p_val}/rep_{combo['rep']}/prevalence_time_series.csv"
                if os.path.exists(prev_file):
                    try:
                        df = pd.read_csv(prev_file)
                        df['p_msmw_w'] = p_val
                        df['replicate'] = rep
                        all_prevalence.append(df)
                        print(f"Added prevalence timeseries for p_msmw_w={p_val}, rep={rep}")
                    except Exception as e:
                        print(f"Failed to read {prev_file}: {e}")
                        continue

        if all_prevalence:
            combined_df = pd.concat(all_prevalence, ignore_index=True)
            combined_df.to_csv(output[0], index=False)
            print(f"Combined prevalence timeseries from {len(all_prevalence)} successful simulations")
        else:
            pd.DataFrame().to_csv(output[0], index=False)
            print("No prevalence timeseries found - created empty file")

rule aggregate_prevalence_averages:
    input:
        expand(f"{SWEEP_OUTPUT_DIR}/batches/analysis_batch_{{batch_id}}_complete.flag",
               batch_id=range(1, NUM_ANALYSIS_BATCHES + 1))
    output:
        f"{SWEEP_OUTPUT_DIR}/aggregated/prevalence_averages.csv"
    conda: "phylo"
    resources:
        mem_mb=12000,
        runtime=30,
        cpus_per_task=1,
        slurm_partition=get_quick_partition
    run:
        import pandas as pd

        all_averages = []
        for combo in PARAM_COMBINATIONS:
            p_val = combo['p_msmw_w']
            rep = int(combo['rep'])

            # Check if analysis was successful
            bridging_file = f"{SWEEP_OUTPUT_DIR}/{p_val}/rep_{combo['rep']}/bridging_analysis.csv"
            if os.path.exists(bridging_file) and os.path.getsize(bridging_file) > 0:
                # Look for average prevalence file
                avg_file = f"{SWEEP_OUTPUT_DIR}/{p_val}/rep_{combo['rep']}/average_prevalences.csv"
                if os.path.exists(avg_file):
                    try:
                        df = pd.read_csv(avg_file)
                        df['p_msmw_w'] = p_val
                        df['replicate'] = rep
                        all_averages.append(df)
                        print(f"Added prevalence averages for p_msmw_w={p_val}, rep={rep}")
                    except Exception as e:
                        print(f"Failed to read {avg_file}: {e}")
                        continue

        if all_averages:
            combined_df = pd.concat(all_averages, ignore_index=True)
            combined_df.to_csv(output[0], index=False)
            print(f"Combined prevalence averages from {len(all_averages)} successful simulations")
        else:
            pd.DataFrame().to_csv(output[0], index=False)
            print("No prevalence averages found - created empty file")
