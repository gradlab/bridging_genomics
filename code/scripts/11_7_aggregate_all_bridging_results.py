#!/usr/bin/env python3
"""
Aggregate bridging results across all 50 simulations
Creates 3 master files: phylogenetic, sampled, and ground truth bridging
"""

import pandas as pd
import os
import glob
from pathlib import Path

def extract_p_msmw_w(sim_id):
    """Extract p_msmw_w value from simulation ID (e.g., '0.05_016' -> 0.05)"""
    return float(sim_id.split('_')[0])

def aggregate_phylogenetic_bridging():
    """Aggregate phylogenetic bridging results"""
    print("🌳 Aggregating phylogenetic bridging results...")
    
    phylo_base = "../../output/baseline_phylo_bridging_complete_fresh_v2/phylo"
    all_phylo_data = []
    
    # Find all simulation directories
    sim_dirs = [d for d in os.listdir(phylo_base) if os.path.isdir(os.path.join(phylo_base, d))]
    sim_dirs.sort()
    
    for sim_id in sim_dirs:
        phylo_file = os.path.join(phylo_base, sim_id, "aggregated_phylogenetic_bridging_30_results.csv")
        
        if os.path.exists(phylo_file):
            try:
                df = pd.read_csv(phylo_file)
                df['sim_id'] = sim_id
                df['p_msmw_w'] = extract_p_msmw_w(sim_id)
                all_phylo_data.append(df)
                print(f"  ✅ Added {sim_id}: {len(df)} rows")
            except Exception as e:
                print(f"  ❌ Error reading {sim_id}: {e}")
        else:
            print(f"  ⚠️  Missing: {sim_id}")
    
    if all_phylo_data:
        master_phylo = pd.concat(all_phylo_data, ignore_index=True)
        
        # Reorder columns to put identifiers first
        cols = ['sim_id', 'p_msmw_w'] + [col for col in master_phylo.columns if col not in ['sim_id', 'p_msmw_w']]
        master_phylo = master_phylo[cols]
        
        master_phylo.to_csv("MASTER_phylogenetic_bridging_30_results.csv", index=False)
        print(f"🎉 Master phylogenetic bridging saved: {len(master_phylo)} total rows from {len(all_phylo_data)} simulations")
        return master_phylo
    else:
        print("❌ No phylogenetic bridging data found!")
        return None

def aggregate_sampled_bridging():
    """Aggregate sampled bridging results"""
    print("\n📊 Aggregating sampled bridging results...")
    
    phylo_base = "../../output/baseline_phylo_bridging_complete_fresh_v2/phylo"
    all_sampled_data = []
    
    # Find all simulation directories
    sim_dirs = [d for d in os.listdir(phylo_base) if os.path.isdir(os.path.join(phylo_base, d))]
    sim_dirs.sort()
    
    for sim_id in sim_dirs:
        sampled_file = os.path.join(phylo_base, sim_id, "aggregated_sampled_bridging_30_results.csv")
        
        if os.path.exists(sampled_file):
            try:
                df = pd.read_csv(sampled_file)
                df['sim_id'] = sim_id
                df['p_msmw_w'] = extract_p_msmw_w(sim_id)
                all_sampled_data.append(df)
                print(f"  ✅ Added {sim_id}: {len(df)} rows")
            except Exception as e:
                print(f"  ❌ Error reading {sim_id}: {e}")
        else:
            print(f"  ⚠️  Missing: {sim_id}")
    
    if all_sampled_data:
        master_sampled = pd.concat(all_sampled_data, ignore_index=True)
        
        # Reorder columns
        cols = ['sim_id', 'p_msmw_w'] + [col for col in master_sampled.columns if col not in ['sim_id', 'p_msmw_w']]
        master_sampled = master_sampled[cols]
        
        master_sampled.to_csv("MASTER_sampled_bridging_30_results.csv", index=False)
        print(f"🎉 Master sampled bridging saved: {len(master_sampled)} total rows from {len(all_sampled_data)} simulations")
        return master_sampled
    else:
        print("❌ No sampled bridging data found!")
        return None

def aggregate_ground_truth_bridging():
    """Aggregate ground truth bridging results"""
    print("\n🎯 Aggregating ground truth bridging results...")
    
    sim_base = "../../output/baseline_phylo_bridging_complete_fresh_v2/simulations"
    all_gt_data = []
    
    # Find all simulation directories
    sim_dirs = [d for d in os.listdir(sim_base) if os.path.isdir(os.path.join(sim_base, d))]
    sim_dirs.sort()
    
    for sim_id in sim_dirs:
        gt_file = os.path.join(sim_base, sim_id, "ground_truth_bridging.csv")
        
        if os.path.exists(gt_file):
            try:
                # Ground truth files are pivot tables (measure as index, methods as columns)
                df = pd.read_csv(gt_file, index_col=0)
                
                # Convert back to long format
                df_long = df.reset_index().melt(
                    id_vars=['measure'], 
                    var_name='method',
                    value_name='value'
                )
                
                df_long['sim_id'] = sim_id
                df_long['p_msmw_w'] = extract_p_msmw_w(sim_id)
                all_gt_data.append(df_long)
                print(f"  ✅ Added {sim_id}: {len(df_long)} rows")
            except Exception as e:
                print(f"  ❌ Error reading {sim_id}: {e}")
        else:
            print(f"  ⚠️  Missing: {sim_id}")
    
    if all_gt_data:
        master_gt = pd.concat(all_gt_data, ignore_index=True)
        
        # Reorder columns
        cols = ['sim_id', 'p_msmw_w', 'measure', 'method', 'value']
        master_gt = master_gt[cols]
        
        master_gt.to_csv("MASTER_ground_truth_bridging_results.csv", index=False)
        print(f"🎉 Master ground truth bridging saved: {len(master_gt)} total rows from {len(all_gt_data)} simulations")
        return master_gt
    else:
        print("❌ No ground truth bridging data found!")
        return None

def main():
    print("🚀 Starting aggregation of all bridging results across 50 simulations...")
    print("=" * 70)
    
    # Change to the right directory
    os.chdir("/n/netscratch/grad_lab/Lab/mkline/bridging_project/output/baseline_phylo_bridging_complete_fresh_v2")
    
    # Aggregate each type
    phylo_df = aggregate_phylogenetic_bridging()
    sampled_df = aggregate_sampled_bridging()
    gt_df = aggregate_ground_truth_bridging()
    
    print("\n" + "=" * 70)
    print("📈 AGGREGATION SUMMARY:")
    
    if phylo_df is not None:
        print(f"  🌳 Phylogenetic: {len(phylo_df)} rows from {phylo_df['sim_id'].nunique()} simulations")
        print(f"     Measures: {phylo_df['measure'].unique()}")
        print(f"     P_MSMW_W range: {phylo_df['p_msmw_w'].min():.2f} - {phylo_df['p_msmw_w'].max():.2f}")
    
    if sampled_df is not None:
        print(f"  📊 Sampled: {len(sampled_df)} rows from {sampled_df['sim_id'].nunique()} simulations")
        print(f"     Measures: {sampled_df['measure'].unique()}")
        print(f"     P_MSMW_W range: {sampled_df['p_msmw_w'].min():.2f} - {sampled_df['p_msmw_w'].max():.2f}")
    
    if gt_df is not None:
        print(f"  🎯 Ground Truth: {len(gt_df)} rows from {gt_df['sim_id'].nunique()} simulations")
        print(f"     Measures: {gt_df['measure'].unique()}")
        print(f"     P_MSMW_W range: {gt_df['p_msmw_w'].min():.2f} - {gt_df['p_msmw_w'].max():.2f}")
    
    print("\n🎉 Master aggregation completed!")
    print("📁 Output files:")
    print("   - MASTER_phylogenetic_bridging_30_results.csv")
    print("   - MASTER_sampled_bridging_30_results.csv") 
    print("   - MASTER_ground_truth_bridging_results.csv")

if __name__ == "__main__":
    main()
