#!/usr/bin/env python3
"""
Create full_sequences downsamples for phylogenetic analysis
"""
import pandas as pd
import json
import os
import numpy as np
import random
import sys
from Bio import SeqIO

def main():
    if len(sys.argv) != 6:
        print("Usage: python create_full_sequences_downsamples.py <actual_sim_dir> <sim_id> <output_base> <pipeline_seed> <output_flag>")
        sys.exit(1)
    
    actual_sim_dir = sys.argv[1]
    sim_id = sys.argv[2]
    output_base = sys.argv[3]
    pipeline_seed = int(sys.argv[4])
    output_flag = sys.argv[5]
    
    # Set random seed
    random.seed(pipeline_seed)
    np.random.seed(pipeline_seed)
    
    print('Loading simulation data for full_sequences...')
    
    # Load post-burn-in, non-superseded, sampled infections
    transmission_df = pd.read_csv(os.path.join(actual_sim_dir, 'transmission_df.csv'))
    
    with open(os.path.join(actual_sim_dir, 'parameters_used.json'), 'r') as f:
        params = json.load(f)
    partnership_burnin = params['simulation']['partnership_burnin_days']
    transmission_burnin = params['simulation']['transmission_burnin_days']
    tracking_start_day = partnership_burnin + transmission_burnin
    
    base_infections = transmission_df[
        (transmission_df['day_of_transmission'] >= tracking_start_day) &
        (transmission_df['superseded_simultaneous'] == False) &
        (transmission_df['day_of_sampling'].notna())
    ].copy()
    
    print(f'Full sequences: {len(base_infections)} infections (should match FASTA)')
    
    # Load sequences for matching
    sequences_file = f'{output_base}/phylo/{sim_id}/sequences/sequences.fasta'
    sequences = list(SeqIO.parse(sequences_file, 'fasta'))
    seq_dict = {}
    for seq in sequences:
        seq_parts = seq.id.split('_')
        if len(seq_parts) >= 2:
            node_sampling_key = f'{seq_parts[0]}_{seq_parts[1]}'
            seq_dict[node_sampling_key] = seq
    
    print(f'Loaded {len(sequences)} sequences for matching')
    
    def save_full_sequences_dataset(infections_df, dataset_name):
        # Create directory structure
        if dataset_name == 'complete':
            dataset_dir = f'{output_base}/phylo/{sim_id}/filtering/full_sequences/complete'
        else:
            rate_part, rep_part = dataset_name.split('_')
            dataset_dir = f'{output_base}/phylo/{sim_id}/filtering/full_sequences/{rate_part}/{rep_part}'
        
        os.makedirs(dataset_dir, exist_ok=True)
        
        # Match sequences
        infections_df = infections_df.copy()
        infections_df['matching_key'] = (
            infections_df['infectee_node'].astype(str) + '_' +
            infections_df['day_of_sampling'].astype(int).astype(str)
        )
        
        matched_sequences = []
        matched_infection_ids = []
        
        for idx, row in infections_df.iterrows():
            matching_key = row['matching_key']
            if matching_key in seq_dict:
                matched_sequences.append(seq_dict[matching_key])
                matched_infection_ids.append(row['transmission_id'])
        
        # Save files
        infections_file = os.path.join(dataset_dir, 'infections.csv')
        fasta_file = os.path.join(dataset_dir, 'sequences.fasta')
        
        matched_infections_df = infections_df[infections_df['transmission_id'].isin(matched_infection_ids)]
        matched_infections_df.to_csv(infections_file, index=False)
        SeqIO.write(matched_sequences, fasta_file, 'fasta')
        
        print(f'  {dataset_name}: {len(matched_sequences)} sequences')
        
        return len(matched_sequences), infections_file, fasta_file
    
    # 1. Save complete full_sequences
    print('Creating complete full_sequences dataset...')
    n_seq, inf_file, fa_file = save_full_sequences_dataset(base_infections, 'complete')
    new_entries = []
    new_entries.append({
        'population': 'all',
        'dataset': 'full',
        'downsample_rate': '100pct',
        'replicate': 'full_sequences',
        'n_infections': len(base_infections),
        'n_sequences': n_seq,
        'infection_file': inf_file,
        'fasta_file': fa_file
    })
    
    # 2. Create temporal downsamples
    print('Creating full_sequences temporal downsamples...')
    
    # Calculate simulation years
    base_infections['sim_year'] = ((base_infections['day_of_transmission'] - tracking_start_day) // 365) + 1
    
    downsample_rates = [0.5, 0.25, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01]
    n_replicates = 10
    
    for rate in downsample_rates:
        print(f'  Downsampling full_sequences to {rate*100}%...')
        
        for replicate in range(n_replicates):
            # Sample by year
            yearly_samples = []
            
            for year in sorted(base_infections['sim_year'].unique()):
                year_infections = base_infections[base_infections['sim_year'] == year]
                n_sample = int(len(year_infections) * rate)
                
                if n_sample > 0:
                    sampled = year_infections.sample(
                        n=n_sample, 
                        random_state=42 + replicate * 1000 + int(rate * 1000)
                    )
                    yearly_samples.append(sampled)
            
            if yearly_samples:
                downsampled = pd.concat(yearly_samples, ignore_index=True)
                dataset_name = f'{rate*100:.0f}pct_rep{replicate+1:02d}'
                
                n_seq, inf_file, fa_file = save_full_sequences_dataset(downsampled, dataset_name)
                
                new_entries.append({
                    'population': 'all',
                    'dataset': 'full',
                    'downsample_rate': f'{rate*100:.0f}pct',
                    'replicate': f'rep{replicate+1:02d}',
                    'n_infections': len(downsampled),
                    'n_sequences': n_seq,
                    'infection_file': inf_file,
                    'fasta_file': fa_file
                })
    
    # 3. Add to filtering_summary.csv
    summary_file = f'{output_base}/phylo/{sim_id}/filtering/filtering_summary.csv'
    if os.path.exists(summary_file):
        summary_df = pd.read_csv(summary_file)
        
        # Add full_sequences entries at the TOP
        new_df = pd.DataFrame(new_entries)
        summary_df = pd.concat([new_df, summary_df], ignore_index=True)
        summary_df.to_csv(summary_file, index=False)
        
        print(f'Added {len(new_entries)} full_sequences datasets to filtering summary')
        print(f'Total datasets now: {len(summary_df)} (will all go through phylo pipeline!)')
    else:
        print('Warning: filtering_summary.csv not found, creating new one with just full_sequences')
        new_df = pd.DataFrame(new_entries)
        new_df.to_csv(summary_file, index=False)
        print(f'Created new filtering summary with {len(new_entries)} full_sequences datasets')
    
    # Create output flag
    with open(output_flag, 'w') as f:
        f.write("FULL_SEQUENCES_READY=True\n")
        f.write(f"ACTUAL_SIM_DIR={actual_sim_dir}\n")
    
    print('✅ Full sequences datasets created and added to pipeline!')

if __name__ == "__main__":
    main()
