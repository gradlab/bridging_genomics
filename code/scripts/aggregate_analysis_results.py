#!/usr/bin/env python3
"""
Aggregate interim ANALYSIS results from ongoing parameter sweep
Collects bridging and behavior classification results from completed analysis batches
"""

import pandas as pd
import numpy as np
import os
import glob
import json
from pathlib import Path
import argparse

def find_completed_analysis_batches(output_dir):
    """Find completed analysis batches"""
    print(f"Scanning for completed analysis results in: {output_dir}")
    
    analysis_flags = glob.glob(f"{output_dir}/batches/analysis_batch_*_complete.flag")
    
    print(f"Found {len(analysis_flags)} completed analysis batches")
    
    # Extract batch IDs
    completed_batches = []
    for flag in analysis_flags:
        batch_id = flag.split('analysis_batch_')[1].split('_complete')[0]
        completed_batches.append(int(batch_id))
    
    return sorted(completed_batches)

def load_batch_mapping(output_dir):
    """Load the batch mapping to understand parameter assignments"""
    mapping_file = f"{output_dir}/batch_mapping.txt"
    if not os.path.exists(mapping_file):
        print(f"Warning: No batch mapping found at {mapping_file}")
        return None
    
    mapping = []
    with open(mapping_file, 'r') as f:
        for line in f:
            batch_id, pmsmw_w, rep_id = line.strip().split(',')
            mapping.append({
                'batch_id': int(batch_id),
                'pmsmw_w': float(pmsmw_w),
                'rep_id': rep_id
            })
    
    return pd.DataFrame(mapping)

def aggregate_analysis_results(output_dir, completed_batches, batch_mapping, analysis_batch_size=10):
    """Aggregate results from completed analysis batches"""
    print("\n=== Aggregating Analysis Results ===")
    
    all_bridging_results = []
    all_bridging_detailed = []
    all_behavior_classifications = []
    all_enhanced_transmissions = []
    
    for batch_id in completed_batches:
        print(f"\nProcessing analysis batch {batch_id}...")
        
        # Calculate which combinations this batch processed
        start_idx = (batch_id - 1) * analysis_batch_size
        end_idx = start_idx + analysis_batch_size
        
        # Get the combinations for this batch
        batch_combos = batch_mapping.iloc[start_idx:end_idx]
        
        for _, combo in batch_combos.iterrows():
            pmsmw_w = combo['pmsmw_w']
            rep_id = combo['rep_id']
            
            # Directory structure: {pmsmw_w}/rep_{id}/
            sim_dir = f"{output_dir}/{pmsmw_w}/rep_{rep_id}"
            
            # 1. Load bridging results
            bridging_file = f"{sim_dir}/bridging_results.csv"
            if os.path.exists(bridging_file) and os.path.getsize(bridging_file) > 0:
                try:
                    df = pd.read_csv(bridging_file)
                    df['pmsmw_w'] = pmsmw_w
                    df['rep_id'] = rep_id
                    df['batch_id'] = batch_id
                    all_bridging_results.append(df)
                    print(f"  ✅ Loaded bridging results from {pmsmw_w}/rep_{rep_id}")
                except Exception as e:
                    print(f"  ❌ Error loading {bridging_file}: {e}")
            
            # 2. Load detailed bridging results
            bridging_detailed_file = f"{sim_dir}/bridging_detailed.csv"
            if os.path.exists(bridging_detailed_file) and os.path.getsize(bridging_detailed_file) > 0:
                try:
                    df = pd.read_csv(bridging_detailed_file)
                    df['pmsmw_w'] = pmsmw_w
                    df['rep_id'] = rep_id
                    df['batch_id'] = batch_id
                    all_bridging_detailed.append(df)
                    print(f"  ✅ Loaded detailed bridging from {pmsmw_w}/rep_{rep_id}")
                except Exception as e:
                    print(f"  ❌ Error loading {bridging_detailed_file}: {e}")
            
            # 3. Load behavior classifications
            behavior_file = f"{sim_dir}/infection_behavior_classifications.csv"
            if os.path.exists(behavior_file) and os.path.getsize(behavior_file) > 0:
                try:
                    df = pd.read_csv(behavior_file)
                    df['pmsmw_w'] = pmsmw_w
                    df['rep_id'] = rep_id
                    df['batch_id'] = batch_id
                    all_behavior_classifications.append(df)
                    print(f"  ✅ Loaded behavior classifications from {pmsmw_w}/rep_{rep_id}")
                except Exception as e:
                    print(f"  ❌ Error loading {behavior_file}: {e}")
            
            # 4. Load enhanced transmission data
            transmission_file = f"{sim_dir}/enhanced_transmission_df.csv"
            if os.path.exists(transmission_file) and os.path.getsize(transmission_file) > 0:
                try:
                    df = pd.read_csv(transmission_file)
                    df['pmsmw_w'] = pmsmw_w
                    df['rep_id'] = rep_id  
                    df['batch_id'] = batch_id
                    all_enhanced_transmissions.append(df)
                    print(f"  ✅ Loaded enhanced transmissions from {pmsmw_w}/rep_{rep_id}")
                except Exception as e:
                    print(f"  ❌ Error loading {transmission_file}: {e}")
    
    # Combine all results
    bridging_df = pd.concat(all_bridging_results, ignore_index=True) if all_bridging_results else pd.DataFrame()
    bridging_detailed_df = pd.concat(all_bridging_detailed, ignore_index=True) if all_bridging_detailed else pd.DataFrame()
    behavior_df = pd.concat(all_behavior_classifications, ignore_index=True) if all_behavior_classifications else pd.DataFrame()
    transmission_df = pd.concat(all_enhanced_transmissions, ignore_index=True) if all_enhanced_transmissions else pd.DataFrame()
    
    print(f"\n=== Aggregation Summary ===")
    print(f"Bridging results: {len(bridging_df)} rows from {bridging_df['rep_id'].nunique() if not bridging_df.empty else 0} replicates")
    print(f"Detailed bridging: {len(bridging_detailed_df)} rows from {bridging_detailed_df['rep_id'].nunique() if not bridging_detailed_df.empty else 0} replicates") 
    print(f"Behavior classifications: {len(behavior_df)} rows from {behavior_df['rep_id'].nunique() if not behavior_df.empty else 0} replicates")
    print(f"Enhanced transmissions: {len(transmission_df)} rows from {transmission_df['rep_id'].nunique() if not transmission_df.empty else 0} replicates")
    
    return bridging_df, bridging_detailed_df, behavior_df, transmission_df

def create_analysis_summary(bridging_df, behavior_df, transmission_df, output_dir):
    """Create summary of analysis results"""
    print("\n=== Creating Analysis Summary ===")
    
    summary = {
        'timestamp': pd.Timestamp.now().isoformat(),
        'analysis_results_summary': {}
    }
    
    if not bridging_df.empty:
        # Bridging results by parameter value
        bridging_summary = bridging_df.groupby('pmsmw_w').agg({
            'rep_id': 'nunique'
        }).rename(columns={'rep_id': 'completed_reps'})
        
        summary['analysis_results_summary']['bridging_data_by_pmsmw_w'] = bridging_summary.to_dict()
        
        # Overall bridging stats if there are numeric columns
        numeric_cols = bridging_df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 3:  # More than just pmsmw_w, rep_id, batch_id
            summary['analysis_results_summary']['bridging_columns_available'] = list(numeric_cols)
    
    if not behavior_df.empty:
        summary['analysis_results_summary']['behavior_analysis'] = {
            'total_infections_classified': len(behavior_df),
            'replicates_with_behavior_data': behavior_df['rep_id'].nunique(),
            'behavior_change_rate': float(behavior_df['behavior_changed'].mean()) if 'behavior_changed' in behavior_df.columns else None
        }
    
    if not transmission_df.empty:
        summary['analysis_results_summary']['transmission_analysis'] = {
            'total_transmissions': len(transmission_df),
            'replicates_with_transmission_data': transmission_df['rep_id'].nunique()
        }
        
        if 'dynamic_behavior' in transmission_df.columns:
            behavior_dist = transmission_df['dynamic_behavior'].value_counts().to_dict()
            summary['analysis_results_summary']['dynamic_behavior_distribution'] = behavior_dist
    
    # Save summary
    summary_file = f"{output_dir}/interim_analysis_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Analysis summary saved to: {summary_file}")
    return summary

def main():
    parser = argparse.ArgumentParser(description='Aggregate interim analysis results from ongoing parameter sweep')
    parser.add_argument('output_dir', help='Parameter sweep output directory')
    parser.add_argument('--analysis-batch-size', type=int, default=10,
                       help='Analysis batch size (default: 10)')
    parser.add_argument('--save-aggregated', action='store_true',
                       help='Save aggregated results to CSV files')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        print(f"❌ Output directory does not exist: {args.output_dir}")
        return 1
    
    print("📊 Aggregating interim ANALYSIS results from ongoing parameter sweep...")
    
    # Find completed analysis batches
    completed_batches = find_completed_analysis_batches(args.output_dir)
    
    if not completed_batches:
        print("❌ No completed analysis batches found yet. Check back later!")
        return 1
    
    print(f"Found completed analysis batches: {completed_batches}")
    
    # Load batch mapping
    batch_mapping = load_batch_mapping(args.output_dir)
    if batch_mapping is None:
        print("❌ Could not load batch mapping")
        return 1
    
    # Aggregate analysis results
    bridging_df, bridging_detailed_df, behavior_df, transmission_df = aggregate_analysis_results(
        args.output_dir, completed_batches, batch_mapping, args.analysis_batch_size
    )
    
    # Create summary
    summary = create_analysis_summary(bridging_df, behavior_df, transmission_df, args.output_dir)
    
    # Save aggregated data if requested
    if args.save_aggregated:
        if not bridging_df.empty:
            bridging_file = f"{args.output_dir}/interim_bridging_results.csv"
            bridging_df.to_csv(bridging_file, index=False)
            print(f"✅ Interim bridging results saved to: {bridging_file}")
        
        if not bridging_detailed_df.empty:
            detailed_file = f"{args.output_dir}/interim_bridging_detailed.csv"
            bridging_detailed_df.to_csv(detailed_file, index=False)
            print(f"✅ Interim detailed bridging saved to: {detailed_file}")
        
        if not behavior_df.empty:
            behavior_file = f"{args.output_dir}/interim_behavior_classifications.csv"
            behavior_df.to_csv(behavior_file, index=False)
            print(f"✅ Interim behavior classifications saved to: {behavior_file}")
        
        if not transmission_df.empty:
            transmission_file = f"{args.output_dir}/interim_enhanced_transmissions.csv"
            transmission_df.to_csv(transmission_file, index=False)
            print(f"✅ Interim enhanced transmissions saved to: {transmission_file}")
    
    print(f"\n🎉 Interim analysis aggregation complete!")
    print(f"Check {args.output_dir}/interim_analysis_summary.json for details")
    
    return 0

if __name__ == "__main__":
    exit(main())
