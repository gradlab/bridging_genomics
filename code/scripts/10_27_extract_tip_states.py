#!/usr/bin/env python3
"""
Extract behavior states for tree tips using dynamic behavior classifications
Creates TreeTime-compatible behavior files with different state groupings
"""

import pandas as pd
import argparse
import os
from Bio import SeqIO

def load_behavior_classifications(sim_dir):
    """Load the infection behavior classifications from the original simulation"""
    behavior_file = os.path.join(sim_dir, "infection_behavior_classifications.csv")
    
    if not os.path.exists(behavior_file):
        raise FileNotFoundError(f"Behavior classifications not found: {behavior_file}")
    
    behavior_df = pd.read_csv(behavior_file)
    print(f"Loaded {len(behavior_df)} behavior classifications from simulation")
    
    # Create mapping from transmission_id to classified behavior
    behavior_mapping = dict(zip(behavior_df['transmission_id'], behavior_df['classified_behavior']))
    
    return behavior_mapping

def extract_tip_behaviors(infections_csv, fasta_file, sim_dir, output_prefix):
    """Extract behavior states for tree tips using dynamic classifications"""
    
    print(f"Reading filtered infections from: {infections_csv}")
    infections_df = pd.read_csv(infections_csv)
    
    print(f"Reading sequences from: {fasta_file}")
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    tip_names = [seq.id for seq in sequences]
    
    print(f"Found {len(tip_names)} sequences in FASTA")
    
    # Load behavior classifications from original simulation
    behavior_mapping = load_behavior_classifications(sim_dir)
    
    # Create mapping from nodename_samplingdate to tip info
    tip_mapping = {}
    for tip_name in tip_names:
        # Extract nodename_samplingdate from nodename_samplingdate_episodenumber
        parts = tip_name.split('_')
        if len(parts) >= 2:
            key = f"{parts[0]}_{parts[1]}"  # nodename_samplingdate
            tip_mapping[key] = tip_name
    
    # Process infections and extract behaviors
    tip_behaviors = []
    
    for _, infection in infections_df.iterrows():
        # Construct the matching key
        node_id = str(infection['infectee_node'])
        sampling_day = int(infection['day_of_sampling'])
        matching_key = f"{node_id}_{sampling_day}"
        transmission_id = infection['transmission_id']
        
        if matching_key in tip_mapping:
            tip_name = tip_mapping[matching_key]
            
            # Get dynamic behavior from classifications
            if transmission_id in behavior_mapping:
                dynamic_behavior = behavior_mapping[transmission_id]
            else:
                # Fallback to infectee behavior if no classification found
                dynamic_behavior = infection['behavior_infectee']
                print(f"Warning: No behavior classification for transmission {transmission_id}, using {dynamic_behavior}")
            
            tip_behaviors.append({
                'tip_name': tip_name,
                'transmission_id': transmission_id,
                'dynamic_behavior': dynamic_behavior,
                'node_id': node_id,
                'sampling_day': sampling_day
            })
    
    if not tip_behaviors:
        raise ValueError("No tip behaviors found! Check that infection data matches FASTA tips.")
    
    print(f"Extracted behaviors for {len(tip_behaviors)} tips")
    
    # Convert to DataFrame
    behavior_df = pd.DataFrame(tip_behaviors)
    
    # Add state mappings
    behavior_df['msmw_msm_state'] = behavior_df['dynamic_behavior'].apply(map_to_msmw_msm_states)
    behavior_df['msmw_msw_state'] = behavior_df['dynamic_behavior'].apply(map_to_msmw_msw_states)
    
    # Create output files
    create_behavior_files(behavior_df, output_prefix)
    
    return behavior_df

def map_to_msmw_msm_states(behavior):
    """MSMW+MSM analysis: MSMW and MSM vs WSM and MSW"""
    if behavior in ['MSMW', 'MSM']:
        return 'MSMW_MSM'
    elif behavior in ['WSM', 'MSW']:
        return 'WSM_MSW'
    else:
        return 'UNKNOWN'

def map_to_msmw_msw_states(behavior):
    """MSMW+MSW analysis: MSM vs everyone else"""
    if behavior == 'MSM':
        return 'MSM'
    elif behavior in ['MSMW', 'MSW', 'WSM']:
        return 'MSMW_MSW_WSM'
    else:
        return 'UNKNOWN'

def create_behavior_files(behavior_df, output_prefix):
    """Create the three TreeTime-compatible behavior files"""
    
    # Create output directory
    output_dir = os.path.dirname(output_prefix)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # 1. Dynamic behavior mapping (4-state)
    dynamic_file = f"{output_prefix}_dynamic_behaviors.csv"
    dynamic_output = behavior_df[['tip_name', 'dynamic_behavior']].copy()
    dynamic_output.columns = ['name', 'behavior']
    dynamic_output.to_csv(dynamic_file, index=False)
    print(f"✅ Dynamic behaviors saved to: {dynamic_file}")
    
    # Show distribution
    print("Dynamic behavior distribution:")
    for behavior, count in dynamic_output['behavior'].value_counts().items():
        print(f"  {behavior}: {count}")
    
    # 2. MSMW+MSM analysis (2-state)
    msmw_msm_file = f"{output_prefix}_msmw_msm_states.csv"
    msmw_msm_output = behavior_df[['tip_name', 'msmw_msm_state']].copy()
    msmw_msm_output.columns = ['name', 'state']
    msmw_msm_output.to_csv(msmw_msm_file, index=False)
    print(f"✅ MSMW+MSM states saved to: {msmw_msm_file}")
    
    # Show distribution
    print("MSMW+MSM state distribution:")
    for state, count in msmw_msm_output['state'].value_counts().items():
        print(f"  {state}: {count}")
    
    # 3. MSMW+MSW analysis (2-state)  
    msmw_msw_file = f"{output_prefix}_msmw_msw_states.csv"
    msmw_msw_output = behavior_df[['tip_name', 'msmw_msw_state']].copy()
    msmw_msw_output.columns = ['name', 'state']
    msmw_msw_output.to_csv(msmw_msw_file, index=False)
    print(f"✅ MSMW+MSW states saved to: {msmw_msw_file}")
    
    # Show distribution
    print("MSMW+MSW state distribution:")
    for state, count in msmw_msw_output['state'].value_counts().items():
        print(f"  {state}: {count}")

def main():
    parser = argparse.ArgumentParser(
        description='Extract behavior states for tree tips using dynamic behavior classifications',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python extract_tip_behaviors.py \\
      filtered_infections.csv \\
      sequences.fasta \\
      /path/to/original/sim_dir \\
      behaviors
  
  # Creates:
  # behaviors_dynamic_behaviors.csv      (4-state dynamic behaviors)
  # behaviors_msmw_msm_states.csv        (MSMW+MSM vs WSM+MSW)
  # behaviors_msmw_msw_states.csv        (MSM vs MSMW+MSW+WSM)
        """
    )
    
    parser.add_argument('infections_csv', help='Filtered infections CSV file')
    parser.add_argument('fasta_file', help='FASTA file with sequences (for tip names)')
    parser.add_argument('sim_dir', help='Original simulation directory (for behavior classifications)')
    parser.add_argument('output_prefix', help='Output file prefix (creates 3 files)')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.infections_csv):
        print(f"❌ Infections CSV not found: {args.infections_csv}")
        return 1
    
    if not os.path.exists(args.fasta_file):
        print(f"❌ FASTA file not found: {args.fasta_file}")
        return 1
    
    if not os.path.exists(args.sim_dir):
        print(f"❌ Simulation directory not found: {args.sim_dir}")
        return 1
    
    print("🧬 Starting tip behavior extraction...")
    
    try:
        # Extract behaviors and create files
        behavior_df = extract_tip_behaviors(
            args.infections_csv,
            args.fasta_file,
            args.sim_dir,
            args.output_prefix
        )
        
        print(f"\n🎉 Behavior extraction completed successfully!")
        print(f"Created 3 behavior files with prefix: {args.output_prefix}")
        print(f"Total tips processed: {len(behavior_df)}")
        
        return 0
        
    except Exception as e:
        print(f"❌ Behavior extraction failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())
