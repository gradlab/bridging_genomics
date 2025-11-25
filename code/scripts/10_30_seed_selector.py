#!/usr/bin/env python3
"""
Simple Seed Selection from Successful Parameter Sweep
Reads aggregated results and copies parameter files for selected simulations
"""

import pandas as pd
import argparse
import os
import shutil
import json
from pathlib import Path

def load_successful_simulations(bridging_values_file):
    """Load successful simulations from aggregated bridging results"""
    
    print(f"📊 Loading successful simulations from: {bridging_values_file}")
    
    df = pd.read_csv(bridging_values_file)
    
    # Get unique combinations of p_msmw_w and replicate
    successful_sims = df[['p_msmw_w', 'replicate']].drop_duplicates()
    
    print(f"✅ Found {len(successful_sims)} successful simulations")
    
    # Group by p_msmw_w to see distribution
    counts = successful_sims.groupby('p_msmw_w').size()
    print(f"📈 Distribution by p_msmw_w:")
    for p_val, count in counts.items():
        print(f"  {p_val}: {count} successful simulations")
    
    return successful_sims

def select_simulations(successful_sims, n_per_param=10, random_seed=42):
    """Select specified number of simulations per parameter value"""
    
    print(f"\n🎯 Selecting {n_per_param} simulations per parameter value...")
    
    selected_sims = []
    
    for p_val in sorted(successful_sims['p_msmw_w'].unique()):
        # Get all successful sims for this parameter
        p_sims = successful_sims[successful_sims['p_msmw_w'] == p_val]
        
        if len(p_sims) >= n_per_param:
            # Randomly sample n_per_param
            selected = p_sims.sample(n=n_per_param, random_state=random_seed)
            print(f"  p_msmw_w = {p_val}: Selected {len(selected)} from {len(p_sims)} available")
        else:
            # Use all available
            selected = p_sims
            print(f"  p_msmw_w = {p_val}: Only {len(selected)} available (wanted {n_per_param})")
        
        selected_sims.append(selected)
    
    final_selection = pd.concat(selected_sims, ignore_index=True)
    
    print(f"\n📋 Final selection: {len(final_selection)} simulations")
    return final_selection

def copy_parameter_files(selected_sims, sweep_dir, output_dir):
    """Copy parameter files for selected simulations"""
    
    print(f"\n📁 Copying parameter files to: {output_dir}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    params_dir = Path(sweep_dir) / "params"
    copied_files = []
    
    for _, row in selected_sims.iterrows():
        p_val = row['p_msmw_w']
        rep = int(row['replicate'])
        
        # Format replicate as 3-digit string (matching original naming)
        rep_str = f"{rep:03d}"
        
        # Source parameter file
        source_file = params_dir / f"params_{p_val}_{rep_str}.json"
        
        # Destination file (simplified naming)
        dest_file = Path(output_dir) / f"params_{p_val}_{rep_str}.json"
        
        if source_file.exists():
            shutil.copy2(source_file, dest_file)
            copied_files.append({
                'p_msmw_w': p_val,
                'replicate': rep,
                'source': str(source_file),
                'destination': str(dest_file)
            })
            
            # Also extract and store the seed for reference
            try:
                with open(source_file, 'r') as f:
                    params = json.load(f)
                seed = params["simulation"]["seed"]
                copied_files[-1]['seed'] = seed
            except Exception as e:
                print(f"    Warning: Could not extract seed from {source_file}: {e}")
                copied_files[-1]['seed'] = None
        else:
            print(f"❌ Source file not found: {source_file}")
    
    print(f"✅ Copied {len(copied_files)} parameter files")
    
    return copied_files

def save_selection_summary(copied_files, output_file):
    """Save summary of selected simulations"""
    
    df = pd.DataFrame(copied_files)
    df.to_csv(output_file, index=False)
    
    print(f"💾 Saved selection summary to: {output_file}")
    
    # Print final summary
    print(f"\n📊 Final Selection Summary:")
    if not df.empty:
        counts = df.groupby('p_msmw_w').size()
        total_sims = len(df)
        for p_val, count in counts.items():
            print(f"  p_msmw_w = {p_val}: {count} simulations")
        print(f"  Total: {total_sims} simulations ready for phylogenetic analysis")
    
    return df

def main():
    parser = argparse.ArgumentParser(
        description='Select parameter files from successful parameter sweep simulations'
    )
    
    parser.add_argument('sweep_dir', help='Parameter sweep output directory')
    parser.add_argument('--n-per-param', '-n', type=int, default=10, 
                       help='Number of simulations to select per parameter value (default: 10)')
    parser.add_argument('--output-dir', '-o', default='selected_params',
                       help='Output directory for selected parameter files (default: selected_params)')
    parser.add_argument('--summary-file', '-s', default='selection_summary.csv',
                       help='Summary CSV file (default: selection_summary.csv)')
    parser.add_argument('--random-seed', type=int, default=42,
                       help='Random seed for selection (default: 42)')
    
    args = parser.parse_args()
    
    # Construct path to aggregated results
    bridging_file = Path(args.sweep_dir) / "aggregated" / "bridging_values.csv"
    
    if not bridging_file.exists():
        print(f"❌ Bridging values file not found: {bridging_file}")
        return 1
    
    print("🌱 Starting parameter file selection for phylogenetic analysis...")
    print(f"📂 Sweep directory: {args.sweep_dir}")
    print(f"🎯 Target per parameter: {args.n_per_param}")
    
    try:
        # 1. Load successful simulations from aggregated results
        successful_sims = load_successful_simulations(bridging_file)
        
        # 2. Select specified number per parameter
        selected_sims = select_simulations(
            successful_sims, 
            n_per_param=args.n_per_param,
            random_seed=args.random_seed
        )
        
        # 3. Copy parameter files
        copied_files = copy_parameter_files(selected_sims, args.sweep_dir, args.output_dir)
        
        # 4. Save summary
        df = save_selection_summary(copied_files, args.summary_file)
        
        print(f"\n🎉 Parameter selection complete!")
        print(f"📁 Parameter files ready in: {args.output_dir}")
        
        return 0
        
    except Exception as e:
        print(f"❌ Selection failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())
