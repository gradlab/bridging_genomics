#!/usr/bin/env python3
"""
Infection Behavior Classification Analysis
Analyzes simulation outputs to classify behavior based on partnerships in year prior to infection
"""

import pandas as pd
import numpy as np
import argparse
import os
import json
from pathlib import Path

def classify_infection_behaviors(transmission_df, edge_df, nodes_df, tracking_start_day):
    """
    For each infection, classify behavior based on partnerships in past 365 days
    
    Parameters:
    -----------
    transmission_df : DataFrame with columns [transmission_id, day_of_transmission, infectee_node, superseded_simultaneous, ...]
    edge_df : DataFrame with columns [day, node1, node2, behavior1, behavior2, duration, ...]
    nodes_df : DataFrame with node behavior assignments
    tracking_start_day : int, day when post-burnin period started
    
    Returns:
    --------
    DataFrame with infection behavior classifications
    """
    
    print("Classifying infection behaviors based on past-year partnerships...")
    
    results = []
    
    # Only look at post-burnin infections that actually occurred (winners)
    sim_infections = transmission_df[
        (transmission_df['day_of_transmission'] >= tracking_start_day) &
        (transmission_df['superseded_simultaneous'] == False) &
        (transmission_df['day_of_sampling'].notna())  # ADD THIS - only successful infections
    ].copy()
    
    print(f"Found {len(sim_infections)} post-burn-in infections to classify")
    
    for idx, infection in sim_infections.iterrows():
        infectee = infection['infectee_node']
        infection_day = infection['day_of_transmission']
        lookback_start = max(tracking_start_day, infection_day - 365)
        
        # Find partnerships involving this infectee in the lookback window
        infectee_partnerships = edge_df[
            ((edge_df['node1'] == infectee) | (edge_df['node2'] == infectee)) &
            (edge_df['day'] < infection_day) &  # Partnership started before infection
            ((edge_df['day'] + edge_df['duration']) > lookback_start)  # Partnership ended after lookback start
        ].copy()
        
        # Count partnerships by partner behavior type
        partner_counts = {'MSM': 0, 'MSMW': 0, 'WSM': 0, 'MSW': 0}
        
        for _, partnership in infectee_partnerships.iterrows():
            # Determine partner's behavior
            if partnership['node1'] == infectee:
                partner_behavior = partnership['behavior2']
            else:
                partner_behavior = partnership['behavior1']
            
            partner_counts[partner_behavior] += 1
        
        # Get original behavior
        if infectee in nodes_df.index:
            original_behavior = nodes_df.loc[infectee, 'behavior']
        else:
            original_behavior = 'UNKNOWN'
        
        # Apply classification rules (only reclassify MSMW)
        if original_behavior == 'MSMW':
            has_male_partners = partner_counts['MSM'] + partner_counts['MSMW'] > 0
            has_female_partners = partner_counts['WSM'] > 0
            
            if has_male_partners and not has_female_partners:
                classified_behavior = 'MSM'
            elif has_female_partners and not has_male_partners:
                classified_behavior = 'MSW'
            else:
                classified_behavior = 'MSMW'  # Had both, neither, or no partnerships
        else:
            classified_behavior = original_behavior  # Non-MSMW don't get reclassified
        
        # Calculate simulation year for analysis
        sim_year = (infection_day - tracking_start_day) // 365 + 1
        
        results.append({
            'transmission_id': infection['transmission_id'],
            'infectee_node': infectee,
            'infection_day': infection_day,
            'lookback_start_day': lookback_start,
            'lookback_window_days': infection_day - lookback_start,
            'original_behavior': original_behavior,
            'classified_behavior': classified_behavior,
            'partnerships_past_year': len(infectee_partnerships),
            'MSM_partners_past_year': partner_counts['MSM'],
            'MSMW_partners_past_year': partner_counts['MSMW'],
            'WSM_partners_past_year': partner_counts['WSM'],
            'MSW_partners_past_year': partner_counts['MSW'],
            'male_partners_past_year': partner_counts['MSM'] + partner_counts['MSMW'],
            'female_partners_past_year': partner_counts['WSM'] + partner_counts['MSW'],
            'sim_year': sim_year,
            'behavior_changed': classified_behavior != original_behavior
        })
    
    results_df = pd.DataFrame(results)
    
    # Print summary statistics
    if not results_df.empty:
        print(f"\nClassification Summary:")
        print(f"Total post-burn-in infections classified: {len(results_df)}")
        print(f"MSMW infections: {(results_df['original_behavior'] == 'MSMW').sum()}")
        print(f"MSMW -> MSM reclassifications: {((results_df['original_behavior'] == 'MSMW') & (results_df['classified_behavior'] == 'MSM')).sum()}")
        print(f"MSMW -> MSW reclassifications: {((results_df['original_behavior'] == 'MSMW') & (results_df['classified_behavior'] == 'MSW')).sum()}")
        print(f"MSMW -> MSMW (no change): {((results_df['original_behavior'] == 'MSMW') & (results_df['classified_behavior'] == 'MSMW')).sum()}")
        
        print(f"\nInfections by simulation year:")
        year_counts = results_df['sim_year'].value_counts().sort_index()
        for year, count in year_counts.items():
            print(f"  Year {year}: {count} infections")
            
        behavior_changes = results_df['behavior_changed'].sum()
        print(f"\nTotal behavior changes: {behavior_changes} ({100*behavior_changes/len(results_df):.1f}%)")
    
    return results_df

def create_enhanced_transmission_log(transmission_df, classification_df):
    """
    Add behavior classification columns to the original transmission log
    
    Parameters:
    -----------
    transmission_df : Original transmission DataFrame
    classification_df : Behavior classification results
    
    Returns:
    --------
    Enhanced transmission DataFrame with classification columns
    """
    print("Creating enhanced transmission log with behavior classifications...")
    
    # Create a mapping from transmission_id to classification info
    classification_map = classification_df.set_index('transmission_id')[
        ['classified_behavior', 'behavior_changed', 'partnerships_past_year', 
         'male_partners_past_year', 'female_partners_past_year']
    ].to_dict('index')
    
    # Add classification columns to transmission_df
    enhanced_df = transmission_df.copy()
    
    # Initialize new columns
    enhanced_df['dynamic_behavior'] = enhanced_df['behavior_infectee']  # Default to original
    enhanced_df['behavior_changed'] = False
    enhanced_df['partnerships_past_year'] = np.nan
    enhanced_df['male_partners_past_year'] = np.nan
    enhanced_df['female_partners_past_year'] = np.nan
    
    # Update with classifications for post-burn-in infections
    for tx_id, classification in classification_map.items():
        if tx_id in enhanced_df['transmission_id'].values:
            mask = enhanced_df['transmission_id'] == tx_id
            enhanced_df.loc[mask, 'dynamic_behavior'] = classification['classified_behavior']
            enhanced_df.loc[mask, 'behavior_changed'] = classification['behavior_changed']
            enhanced_df.loc[mask, 'partnerships_past_year'] = classification['partnerships_past_year']
            enhanced_df.loc[mask, 'male_partners_past_year'] = classification['male_partners_past_year']
            enhanced_df.loc[mask, 'female_partners_past_year'] = classification['female_partners_past_year']
    
    print(f"Enhanced transmission log: {len(enhanced_df)} total transmissions")
    print(f"  {(enhanced_df['dynamic_behavior'] != enhanced_df['behavior_infectee']).sum()} with behavior changes")
    
    return enhanced_df

def load_simulation_data(output_dir):
    """Load core simulation outputs with robust parameter detection"""
    print(f"Loading simulation data from: {output_dir}")
    
    # Load core files
    transmission_df = pd.read_csv(os.path.join(output_dir, "transmission_df.csv"))
    edge_df = pd.read_csv(os.path.join(output_dir, "edge_df.csv"))
    nodes_df = pd.read_csv(os.path.join(output_dir, "nodes_df.csv"), index_col=0)
    
    # Load parameters to get tracking_start_day
    tracking_start_day = None
    
    # Try parameters_used.json first (most reliable)
    try:
        with open(os.path.join(output_dir, "parameters_used.json"), 'r') as f:
            params = json.load(f)
        partnership_burnin = params['simulation']['partnership_burnin_days']
        transmission_burnin = params['simulation']['transmission_burnin_days']
        tracking_start_day = partnership_burnin + transmission_burnin
        print(f"Loaded parameters: partnership_burnin={partnership_burnin}, transmission_burnin={transmission_burnin}")
    
    except Exception as e:
        print(f"Warning: Could not load parameters_used.json: {e}")
        
        # Fallback: try to infer from transmission data
        try:
            # Look for a gap in transmission timing that might indicate burn-in end
            tx_days = transmission_df['day_of_transmission'].value_counts().sort_index()
            # This is a rough heuristic - adjust based on your simulation structure
            tracking_start_day = 12000  # Common default in your simulations
            print(f"Warning: Using default tracking_start_day = {tracking_start_day}")
        except:
            tracking_start_day = 12000
            print(f"Warning: Using fallback tracking_start_day = {tracking_start_day}")
    
    print(f"Loaded data:")
    print(f"  Nodes: {len(nodes_df)}")
    print(f"  Partnerships: {len(edge_df)}")
    print(f"  Transmissions: {len(transmission_df)}")
    print(f"  Tracking start day: {tracking_start_day}")
    
    return transmission_df, edge_df, nodes_df, tracking_start_day

def main():
    parser = argparse.ArgumentParser(description='Analyze infection behavior classifications from simulation outputs')
    parser.add_argument('output_dir', help='Directory containing simulation outputs')
    parser.add_argument('--output-classifications', '-c', 
                       help='Output file for behavior classifications (default: infection_behavior_classifications.csv)')
    parser.add_argument('--output-enhanced-transmission', '-e', 
                       help='Output file for enhanced transmission log with classifications (default: none)')
    parser.add_argument('--both', action='store_true', 
                       help='Output both classifications and enhanced transmission log')
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.exists(args.output_dir):
        print(f"❌ Output directory does not exist: {args.output_dir}")
        return 1
    
    # Set default output files
    if args.output_classifications is None:
        args.output_classifications = os.path.join(args.output_dir, "infection_behavior_classifications.csv")
    
    if args.both and args.output_enhanced_transmission is None:
        args.output_enhanced_transmission = os.path.join(args.output_dir, "enhanced_transmission_df.csv")
    
    try:
        print("🧬 Starting infection behavior classification analysis...")
        
        # Load data
        transmission_df, edge_df, nodes_df, tracking_start_day = load_simulation_data(args.output_dir)
        
        # Perform classification analysis
        classification_df = classify_infection_behaviors(transmission_df, edge_df, nodes_df, tracking_start_day)
        
        if classification_df.empty:
            print("❌ No post-burn-in infections found to classify")
            return 1
        
        # Save behavior classifications
        classification_df.to_csv(args.output_classifications, index=False)
        print(f"\n✅ Behavior classifications saved to: {args.output_classifications}")
        
        # Optionally create enhanced transmission log
        if args.output_enhanced_transmission or args.both:
            enhanced_transmission_df = create_enhanced_transmission_log(transmission_df, classification_df)
            
            output_file = args.output_enhanced_transmission or os.path.join(args.output_dir, "enhanced_transmission_df.csv")
            enhanced_transmission_df.to_csv(output_file, index=False)
            print(f"✅ Enhanced transmission log saved to: {output_file}")
        
        return 0
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())