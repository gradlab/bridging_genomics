#!/usr/bin/env python3
"""
Systematic Infection Filtering and Downsampling for Phylogenetic Analysis
Creates multiple filtered populations and temporal downsamples for bias studies
"""

import pandas as pd
import numpy as np
import argparse
import os
import json
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random
from collections import defaultdict

def load_simulation_data(sim_dir):
    """Load simulation outputs and metadata"""
    print(f"Loading simulation data from: {sim_dir}")
    
    # Load core data
    transmission_df = pd.read_csv(f"{sim_dir}/transmission_df.csv")
    nodes_df = pd.read_csv(f"{sim_dir}/nodes_df.csv", index_col=0)
    
    # Load parameters to get tracking start day
    with open(f"{sim_dir}/parameters_used.json", 'r') as f:
        params = json.load(f)
    
    partnership_burnin = params['simulation']['partnership_burnin_days']
    transmission_burnin = params['simulation']['transmission_burnin_days']  
    tracking_start_day = partnership_burnin + transmission_burnin
    
    print(f"Loaded data:")
    print(f"  Transmissions: {len(transmission_df)}")
    print(f"  Nodes: {len(nodes_df)}")
    print(f"  Tracking start day: {tracking_start_day}")
    
    return transmission_df, nodes_df, tracking_start_day

def create_base_populations(transmission_df, nodes_df, tracking_start_day, high_activity_subsample_fraction=0.25):
    """Create the 5 base filtered populations"""
    print("\n🔍 Creating base populations...")
    
    # Start with post-burn-in infections that actually occurred
    base_infections = transmission_df[
        (transmission_df['day_of_transmission'] >= tracking_start_day) &
        (transmission_df['superseded_simultaneous'] == False) &
        (transmission_df['day_of_sampling'].notna())
    ].copy()
    
    print(f"Base post-burn-in infections: {len(base_infections)}")
    
    populations = {}
    
    # 1. Detected only (remove natural clearance)
    detected_only = base_infections[
        base_infections['duration_source'] != 'natural_clearance'
    ].copy()
    populations['detected_only'] = detected_only
    print(f"  Detected only: {len(detected_only)}")
    
    # 2. Symptomatic only  
    symptomatic_only = detected_only[
        detected_only['infection_symptom_status'] == 'symptomatic'
    ].copy()
    populations['symptomatic_only'] = symptomatic_only
    print(f"  Symptomatic only: {len(symptomatic_only)}")
    
    # 3. Symptomatic men only (MSW, MSM, MSMW - no WSM)
    symptomatic_men = symptomatic_only[
        symptomatic_only['behavior_infectee'].isin(['MSW', 'MSM', 'MSMW'])
    ].copy()
    populations['symptomatic_men'] = symptomatic_men
    print(f"  Symptomatic men: {len(symptomatic_men)}")
    
    # 4. High activity symptomatic men
    high_activity_nodes = nodes_df[nodes_df['hi_risk'] == 1].index.tolist()
    
    high_activity_men = symptomatic_men[
        symptomatic_men['infectee_node'].isin(high_activity_nodes)
    ].copy()
    populations['high_activity_symptomatic_men'] = high_activity_men
    print(f"  High activity symptomatic men: {len(high_activity_men)}")
    
    # 5. Random subsample of high activity men (all their NON-NATURALLY-CLEARED infections)
    # Get all high activity men (MSW, MSM, MSMW)
    high_activity_men_nodes = nodes_df[
        (nodes_df['hi_risk'] == 1) &
        (nodes_df['behavior'].isin(['MSW', 'MSM', 'MSMW']))
    ].index.tolist()
    
    # Subsample nodes
    n_subsample = int(len(high_activity_men_nodes) * high_activity_subsample_fraction)
    if n_subsample == 0 and len(high_activity_men_nodes) > 0:
        n_subsample = 1  # Take at least 1 if any exist
    
    sampled_nodes = random.sample(high_activity_men_nodes, min(n_subsample, len(high_activity_men_nodes)))
    
    # Get ALL detected (non-naturally-cleared) infections from these sampled nodes
    random_subsample = detected_only[  # Changed from base_infections to detected_only
        detected_only['infectee_node'].isin(sampled_nodes)
    ].copy()
    populations['random_subsample_high_activity'] = random_subsample
    print(f"  Random subsample high activity ({high_activity_subsample_fraction*100}% of {len(high_activity_men_nodes)} nodes = {len(sampled_nodes)} nodes): {len(random_subsample)} infections")
    
    return populations

def apply_temporal_downsampling(infections_df, tracking_start_day, downsample_rates, n_replicates=10):
    """Apply yearly temporal downsampling to infection data"""
    print(f"\n⏰ Applying temporal downsampling...")
    
    # Calculate simulation years
    infections_df = infections_df.copy()
    infections_df['sim_year'] = ((infections_df['day_of_transmission'] - tracking_start_day) // 365) + 1
    
    downsampled_datasets = {}
    
    for rate in downsample_rates:
        print(f"  Downsampling to {rate*100}% per year...")
        
        for replicate in range(n_replicates):
            # Group by year and sample from each year
            yearly_samples = []
            
            for year in sorted(infections_df['sim_year'].unique()):
                year_infections = infections_df[infections_df['sim_year'] == year]
                n_sample = int(len(year_infections) * rate)
                
                if n_sample > 0:
                    # Use different random state for each replicate
                    sampled = year_infections.sample(n=n_sample, random_state=42 + replicate * 1000 + int(rate * 1000))
                    yearly_samples.append(sampled)
            
            if yearly_samples:
                downsampled = pd.concat(yearly_samples, ignore_index=True)
                dataset_name = f"{rate*100:.0f}pct_rep{replicate+1:02d}"
                downsampled_datasets[dataset_name] = downsampled
    
    print(f"  Created {len(downsampled_datasets)} downsampled datasets")
    return downsampled_datasets

def match_infections_to_sequences(infections_df, fasta_file, output_dir, dataset_name, population_name):
    """Match filtered infections to FASTA sequences with organized directory structure"""
    
    # Create nested directory structure
    if dataset_name == "complete":
        # population/complete/
        dataset_dir = os.path.join(output_dir, population_name, "complete")
    else:
        # Parse dataset_name like "50pct_rep01"
        rate_part, rep_part = dataset_name.split('_')
        # population/50pct/rep01/
        dataset_dir = os.path.join(output_dir, population_name, rate_part, rep_part)
    
    os.makedirs(dataset_dir, exist_ok=True)
    
    # Load original sequences
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    # Create a mapping from nodename_samplingdate to sequence
    seq_dict = {}
    for seq in sequences:
        # Extract nodename_samplingdate from nodename_samplingdate_episodenumber
        seq_parts = seq.id.split('_')
        if len(seq_parts) >= 2:
            node_sampling_key = f"{seq_parts[0]}_{seq_parts[1]}"  # nodename_samplingdate
            seq_dict[node_sampling_key] = seq
    
    # Construct matching keys from infection data: nodename_samplingdate
    infections_df = infections_df.copy()
    infections_df['matching_key'] = (
        infections_df['infectee_node'].astype(str) + '_' + 
        infections_df['day_of_sampling'].astype(int).astype(str)
    )
    
    # Filter sequences
    matched_sequences = []
    matched_infection_ids = []
    
    for idx, row in infections_df.iterrows():
        matching_key = row['matching_key']
        if matching_key in seq_dict:
            matched_sequences.append(seq_dict[matching_key])
            matched_infection_ids.append(row['transmission_id'])
        else:
            print(f"    Warning: No sequence found for {matching_key}")
    
    print(f"    Matched {len(matched_sequences)} sequences out of {len(infections_df)} infections")
    
    # Save files with simple names in organized directories
    infection_list_file = os.path.join(dataset_dir, "infections.csv")
    fasta_file_out = os.path.join(dataset_dir, "sequences.fasta")
    
    matched_infections_df = infections_df[infections_df['transmission_id'].isin(matched_infection_ids)]
    matched_infections_df.to_csv(infection_list_file, index=False)
    SeqIO.write(matched_sequences, fasta_file_out, "fasta")
    
    return len(matched_sequences), infection_list_file, fasta_file_out

def main():
    parser = argparse.ArgumentParser(
        description='Create filtered and downsampled infection datasets for phylogenetic analysis'
    )
    parser.add_argument('sim_dir', help='Simulation output directory')
    parser.add_argument('fasta_file', help='Input FASTA file from seq-gen')
    parser.add_argument('output_dir', help='Output directory for filtered datasets')
    parser.add_argument('--downsample-rates', nargs='+', type=float, 
                       default=[0.5, 0.25, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01],
                       help='Downsampling rates (default: 50%% to 1%%)')
    parser.add_argument('--n-replicates', type=int, default=10,
                       help='Number of replicates per downsampling rate (default: 10)')
    parser.add_argument('--high-activity-subsample', type=float, default=0.25,
                       help='Fraction of high activity nodes to subsample (default: 0.25)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility (default: 42)')
    
    args = parser.parse_args()
    
    # Set random seed
    random.seed(args.seed)
    np.random.seed(args.seed)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("🧬 Starting systematic infection filtering and downsampling...")
    print(f"Input simulation: {args.sim_dir}")
    print(f"Input FASTA: {args.fasta_file}")
    print(f"Output directory: {args.output_dir}")
    print(f"Downsample rates: {[f'{r*100}%' for r in args.downsample_rates]}")
    print(f"Replicates per rate: {args.n_replicates}")
    
    # Load simulation data
    transmission_df, nodes_df, tracking_start_day = load_simulation_data(args.sim_dir)
    
    # Create base populations
    populations = create_base_populations(
        transmission_df, nodes_df, tracking_start_day, args.high_activity_subsample
    )
    
    # Process each population
    summary_stats = []
    
    for pop_name, pop_infections in populations.items():
        print(f"\n📊 Processing population: {pop_name} ({len(pop_infections)} infections)")
        
        if len(pop_infections) == 0:
            print(f"  Skipping {pop_name} - no infections")
            continue
        # FIRST: Save the complete (100%) population 
        print(f"  Saving complete population...")
        n_sequences, infection_file, fasta_file = match_infections_to_sequences(
            pop_infections, args.fasta_file, args.output_dir, "complete", pop_name
        )
        
        summary_stats.append({
            'population': pop_name,
            'dataset': 'complete',
            'downsample_rate': '100pct',
            'replicate': 'all',
            'n_infections': len(pop_infections),
            'n_sequences': n_sequences,
            'infection_file': infection_file,
            'fasta_file': fasta_file
        })
        
        # Apply temporal downsampling
        downsampled_datasets = apply_temporal_downsampling(
            pop_infections, tracking_start_day, args.downsample_rates, args.n_replicates
        )
        
        # Process each downsampled dataset
        for dataset_name, dataset_infections in downsampled_datasets.items():
            n_sequences, infection_file, fasta_file = match_infections_to_sequences(
                dataset_infections, args.fasta_file, args.output_dir, dataset_name, pop_name
            )
            
            summary_stats.append({
                'population': pop_name,
                'dataset': dataset_name,
                'downsample_rate': dataset_name.split('_')[0],
                'replicate': dataset_name.split('_')[1],
                'n_infections': len(dataset_infections),
                'n_sequences': n_sequences,
                'infection_file': infection_file,
                'fasta_file': fasta_file
            })
    
    # Save summary
    summary_df = pd.DataFrame(summary_stats)
    summary_file = f"{args.output_dir}/filtering_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    
    print(f"\n🎉 Filtering complete!")
    print(f"Created {len(summary_stats)} filtered datasets")
    print(f"Summary saved to: {summary_file}")
    
    # Print summary by population
    for pop_name in populations.keys():
        pop_summary = summary_df[summary_df['population'] == pop_name]
        if len(pop_summary) > 0:
            print(f"\n{pop_name}: {len(pop_summary)} datasets")
            print(f"  Infections range: {pop_summary['n_infections'].min()}-{pop_summary['n_infections'].max()}")
            print(f"  Sequences range: {pop_summary['n_sequences'].min()}-{pop_summary['n_sequences'].max()}")

if __name__ == "__main__":
    main()
