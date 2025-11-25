# local_snakefile.smk
import json
import os
import glob
from pathlib import Path

# Configuration
PARAM_FILE = "../parameters/10_27_0.75_rep454_params.json"

# Load parameters to get base output directory
with open(PARAM_FILE, 'r') as f:
    params = json.load(f)

BASE_OUTPUT_DIR = params["output"]["output_dir"].rstrip("/")
RUN_NAME = params["output"].get("run_name", "simulation")
PIPELINE_OUTPUTS_DIR = f"{BASE_OUTPUT_DIR}/{RUN_NAME}/pipeline_outputs"

# Function to find the actual timestamped directory
def get_actual_output_dir():
    """Find the timestamped output directory created by simulation"""
    pattern = f"{BASE_OUTPUT_DIR}/{RUN_NAME}_*"
    dirs = glob.glob(pattern)
    if dirs:
        # Return the most recent one
        return max(dirs, key=os.path.getctime)
    else:
        # If no timestamped dir exists yet, return the pattern for Snakemake
        return f"{BASE_OUTPUT_DIR}/{{timestamp}}"

# Main target rule
rule all:
    input:
        f"{PIPELINE_OUTPUTS_DIR}/plots_complete.flag",
        f"{PIPELINE_OUTPUTS_DIR}/analysis_complete.flag"

# Step 1: Run simulation and capture actual output directory
rule run_simulation:
    input:
        params_file = PARAM_FILE
    output:
        flag = f"{PIPELINE_OUTPUTS_DIR}/simulation_complete.flag"
    params:
        base_dir = BASE_OUTPUT_DIR,
        run_name = RUN_NAME,
        pipeline_dir = PIPELINE_OUTPUTS_DIR
    conda:
        "phylo"
    resources:
        mem_mb=8000,
        runtime=300,  # 5 hours
        cpus_per_task=2,
        slurm_partition="shared"
    shell:
        """
        mkdir -p {params.pipeline_dir}
        python ../scripts/10_16_core_sim.py {input.params_file}
        
        # Point directly to the expected directory (no timestamp searching)
        echo "ACTUAL_OUTPUT_DIR={params.base_dir}/{params.run_name}" > {output.flag}
        """

# Step 2: Run plotting using actual directory
rule generate_plots:
    input:
        f"{PIPELINE_OUTPUTS_DIR}/simulation_complete.flag"
    output:
        f"{PIPELINE_OUTPUTS_DIR}/plots_complete.flag"
    conda:
        "phylo"   
    resources:
        mem_mb=12000,
        runtime=60,  # 1 hour
        cpus_per_task=1,
        slurm_partition="shared"
    shell:
        """
        # Read actual output directory
        ACTUAL_DIR=$(grep ACTUAL_OUTPUT_DIR {input} | cut -d'=' -f2)
        
        echo "Running plotting on: $ACTUAL_DIR"
        python ../scripts/10_20_plotting_script.py "$ACTUAL_DIR"
        echo "PLOTS_COMPLETE=true" > {output}
        """

# Step 3: Test for single lineage
rule test_single_lineage:
    input:
        f"{PIPELINE_OUTPUTS_DIR}/simulation_complete.flag"
    output:
        f"{PIPELINE_OUTPUTS_DIR}/lineage_complete.flag"
    conda:
        "phylo"
    resources:
        mem_mb=4000,
        runtime=30,  # 30 minutes
        cpus_per_task=1,
        slurm_partition="shared"
    shell:
        """
        # Read actual output directory
        ACTUAL_DIR=$(grep ACTUAL_OUTPUT_DIR {input} | cut -d'=' -f2)
        echo "Running lineage test on: $ACTUAL_DIR"

        python ../scripts/10_17_single_lineage_test_script.py "$ACTUAL_DIR"

        # Check if single lineage and store result
        if [ -f "$ACTUAL_DIR/lineage_test_results.json" ]; then
            SINGLE_LINEAGE=$(python -c "import json; print(json.load(open('$ACTUAL_DIR/lineage_test_results.json'))['single_lineage'])")
            echo "SINGLE_LINEAGE=$SINGLE_LINEAGE" > {output}
            echo "ACTUAL_OUTPUT_DIR=$ACTUAL_DIR" >> {output}
        else
            echo "SINGLE_LINEAGE=False" > {output}
            echo "ACTUAL_OUTPUT_DIR=$ACTUAL_DIR" >> {output}
        fi
        """

# Step 4: Conditional analysis based on lineage results
rule run_analysis:
    input:
        f"{PIPELINE_OUTPUTS_DIR}/lineage_complete.flag"
    output:
        f"{PIPELINE_OUTPUTS_DIR}/analysis_complete.flag"
    conda:
        "phylo"
    resources:
        mem_mb=12000,  # More memory for tree building and analysis
        runtime=240,   # 4 hours
        cpus_per_task=2,
        slurm_partition="shared"
    shell:
        """
        # Read lineage results
        SINGLE_LINEAGE=$(grep SINGLE_LINEAGE {input} | cut -d'=' -f2)
        ACTUAL_DIR=$(grep ACTUAL_OUTPUT_DIR {input} | cut -d'=' -f2)

        echo "Single lineage: $SINGLE_LINEAGE"
        echo "Output directory: $ACTUAL_DIR"

        if [ "$SINGLE_LINEAGE" = "True" ]; then
            echo "✅ Running full analysis pipeline..."

            # Tree building
            echo "Building transmission tree..."
            python ../scripts/10_21_build_hybrid_transmission_tree.py "$ACTUAL_DIR"

            # Tree trimming
            echo "Trimming tree..."
            TREE_FILE=$(ls "$ACTUAL_DIR"/tree_*_hybrid_episodes*.newick | head -1)
            Rscript ../scripts/10_20_tree_trimmer.R "$TREE_FILE" "$ACTUAL_DIR" "$ACTUAL_DIR/trimmed_tree.nwk"

            Rscript ../scripts/10_21_tree_additional_trimming.R "$ACTUAL_DIR/trimmed_tree.nwk" "$ACTUAL_DIR/transmission_df.csv" "$ACTUAL_DIR/final_trimmed_tree.nwk"

            # Behavior classification
            echo "Running behavior classification..."
            python ../scripts/10_17_infection_behavior_classification_analysis.py "$ACTUAL_DIR" --both

            # Bridging analysis
            echo "Calculating bridging measures..."
            python ../scripts/10_23_calculate_bridging.py "$ACTUAL_DIR" "$ACTUAL_DIR/infection_behavior_classifications.csv" --output "$ACTUAL_DIR/bridging_analysis.csv" --detailed-output "$ACTUAL_DIR/bridging_detailed.csv"

            echo "ANALYSIS_COMPLETE=true" > {output}
        else
            echo "⚠️  Skipping analysis - multiple lineages detected"
            echo "ANALYSIS_SKIPPED=multiple_lineages" > {output}
        fi
        """

# Utility rules
rule clean:
    shell:
        """
        rm -rf {BASE_OUTPUT_DIR}/{RUN_NAME}_*
        rm -rf {PIPELINE_OUTPUTS_DIR}
        """
