#!/usr/bin/env python3
"""
Calculate sampled bridging metrics from filtered infection data
Matches 10_30 ground truth bridging analysis format and measures
"""

import pandas as pd
import argparse
import os

def get_dynamic_behavior_map(classifications_df, transmission_df):
    """Create mapping from transmission_id to dynamic behavior"""
    behavior_map = {}
    
    # Map from classifications
    for _, row in classifications_df.iterrows():
        behavior_map[row['transmission_id']] = row['classified_behavior']
    
    # For transmissions not in classifications, use original behavior
    for _, row in transmission_df.iterrows():
        if row['transmission_id'] not in behavior_map:
            behavior_map[row['transmission_id']] = row['behavior_infectee']
    
    return behavior_map

def get_infector_transmission_id(infector_node, transmission_day, post_burnin):
    """Find the transmission_id for when the infector was infected (and still active)"""
    
    # Find all infections of the infector
    infector_infections = post_burnin[post_burnin['infectee_node'] == infector_node]
    
    # Find the infection that was active when they infected someone else
    active_infections = infector_infections[
        (infector_infections['day_of_transmission'] <= transmission_day) &
        (infector_infections['day_of_sampling'] >= transmission_day)
    ]
    
    if len(active_infections) > 0:
        # Return the most recent active infection
        return active_infections.sort_values('day_of_transmission').iloc[-1]['transmission_id']
    else:
        return None

def get_network_groups(behavior, grouping):
    """Map behaviors to consistent WSM_group vs MSM_group (matching 10_30 format)"""
    
    if grouping == "MSM+MSMW":
        # MSM_group = MSM + MSMW, WSM_group = WSM + MSW
        if behavior in ['MSM', 'MSMW']:
            return 'MSM_group'
        elif behavior in ['WSM', 'MSW']:
            return 'WSM_group'
    
    elif grouping == "MSMW+MSW":
        # MSM_group = MSM only, WSM_group = MSMW + MSW + WSM  
        if behavior == 'MSM':
            return 'MSM_group'
        elif behavior in ['MSMW', 'MSW', 'WSM']:
            return 'WSM_group'
    
    return 'UNKNOWN'

def calculate_between_network_proportions(behavior_df, grouping):
    """Calculate between network proportions (matching 10_30 format)"""
    
    # Map behaviors to consistent network groups
    behavior_df['infector_network'] = behavior_df['infector_dynamic'].apply(
        lambda x: get_network_groups(x, grouping)
    )
    behavior_df['infectee_network'] = behavior_df['infectee_dynamic'].apply(
        lambda x: get_network_groups(x, grouping)
    )
    
    # Count transitions
    between_network_count = 0
    WSM_to_MSM_count = 0  
    MSM_to_WSM_count = 0  
    
    for _, row in behavior_df.iterrows():
        infector_network = row['infector_network']
        infectee_network = row['infectee_network']
        
        if infector_network != infectee_network:
            between_network_count += 1
            
            # Track consistent directionality (matching 10_30)
            if infector_network == "WSM_group" and infectee_network == "MSM_group":
                WSM_to_MSM_count += 1  # WSM_group → MSM_group
            elif infector_network == "MSM_group" and infectee_network == "WSM_group":
                MSM_to_WSM_count += 1  # MSM_group → WSM_group
    
    # Calculate proportions
    total = len(behavior_df)
    overall_proportion = between_network_count / total if total > 0 else 0
    WSM_to_MSM_proportion = WSM_to_MSM_count / total if total > 0 else 0
    MSM_to_WSM_proportion = MSM_to_WSM_count / total if total > 0 else 0
    
    return {
        'overall_proportion': overall_proportion,
        'overall_count': between_network_count,
        'WSM_to_MSM_proportion': WSM_to_MSM_proportion,
        'WSM_to_MSM_count': WSM_to_MSM_count,
        'MSM_to_WSM_proportion': MSM_to_WSM_proportion,
        'MSM_to_WSM_count': MSM_to_WSM_count,
        'total': total
    }

def calculate_sampled_bridging(infections_csv, sim_dir):
    """Calculate sampled bridging metrics (matching 10_30 format)"""
    
    print(f"Loading filtered infections from: {infections_csv}")
    infections_df = pd.read_csv(infections_csv)
    
    # Load behavior classifications 
    print(f"Loading behavior classifications from: {sim_dir}")
    classifications_file = os.path.join(sim_dir, "infection_behavior_classifications.csv")
    classifications_df = pd.read_csv(classifications_file)
    
    # Create behavior mapping
    behavior_mapping = get_dynamic_behavior_map(classifications_df, infections_df)
    print(f"Created behavior mapping for {len(behavior_mapping)} transmissions")
    
    # Get dynamic behaviors for infectors and infectees
    dynamic_behaviors = []
    infector_lookup_stats = {'found': 0, 'not_found': 0}
    
    for _, row in infections_df.iterrows():
        transmission_id = row['transmission_id']
        infector_node = row['infector_node']
        infectee_node = row['infectee_node']
        transmission_day = row['day_of_transmission']
        
        # Get dynamic behavior for infectee
        infectee_dynamic = behavior_mapping.get(transmission_id, row['behavior_infectee'])
        
        # Get dynamic behavior for infector
        infector_transmission_id = get_infector_transmission_id(
            infector_node, transmission_day, infections_df
        )
        
        if infector_transmission_id is not None:
            infector_lookup_stats['found'] += 1
            infector_dynamic = behavior_mapping.get(infector_transmission_id, row['behavior_infector'])
        else:
            infector_lookup_stats['not_found'] += 1
            infector_dynamic = row['behavior_infector']  # Fallback to static
        
        dynamic_behaviors.append({
            'transmission_id': transmission_id,
            'infector_node': infector_node,
            'infectee_node': infectee_node,
            'infector_dynamic': infector_dynamic,
            'infectee_dynamic': infectee_dynamic,
            'infector_transmission_id': infector_transmission_id,
            'infector_static': row['behavior_infector'],
            'infectee_static': row['behavior_infectee']
        })
    
    print(f"Infector lookup stats: {infector_lookup_stats['found']} found, {infector_lookup_stats['not_found']} not found in filtered data")
    
    behavior_df = pd.DataFrame(dynamic_behaviors)
    
    # Initialize results (matching 10_30 format)
    results = []
    
    # Calculate for both groupings (matching 10_30)
    for grouping in ["MSM+MSMW", "MSMW+MSW"]:
        print(f"\nCalculating sampled bridging for {grouping} grouping...")
        
        # Calculate between network proportions
        between_results = calculate_between_network_proportions(behavior_df, grouping)
        
        # Add results (matching 10_30 measure names exactly)
        results.extend([
            {"measure": "total_between_network_proportion", "grouping": grouping, 
             "value": between_results['overall_proportion'], 
             "count": between_results['overall_count'], 
             "total": between_results['total']},
            
            {"measure": "WSM_group_to_MSM_group_proportion", "grouping": grouping,
             "value": between_results['WSM_to_MSM_proportion'],
             "count": between_results['WSM_to_MSM_count'],
             "total": between_results['total']},
            
            {"measure": "MSM_group_to_WSM_group_proportion", "grouping": grouping,
             "value": between_results['MSM_to_WSM_proportion'],
             "count": between_results['MSM_to_WSM_count'],
             "total": between_results['total']}
        ])
        
        # Print summary for this grouping
        print(f"  Between network - Total: {between_results['overall_proportion']:.4f} ({between_results['overall_count']}/{between_results['total']})")
        print(f"  Between network - WSM→MSM: {between_results['WSM_to_MSM_proportion']:.4f} ({between_results['WSM_to_MSM_count']}/{between_results['total']})")
        print(f"  Between network - MSM→WSM: {between_results['MSM_to_WSM_proportion']:.4f} ({between_results['MSM_to_WSM_count']}/{between_results['total']})")
    
    return results, behavior_df

def main():
    parser = argparse.ArgumentParser(
        description='Calculate sampled bridging metrics matching 10_30 ground truth format'
    )
    
    parser.add_argument('infections_csv', help='Filtered infections CSV file')
    parser.add_argument('sim_dir', help='Original simulation directory')
    parser.add_argument('--output', '-o', default='sampled_bridging_30.csv', 
                       help='Output file for results table')
    parser.add_argument('--detailed-output', '-d', default='sampled_bridging_30_detailed.csv',
                       help='Output file for detailed behaviors')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.infections_csv):
        print(f"❌ Infections file not found: {args.infections_csv}")
        return 1
    
    if not os.path.exists(args.sim_dir):
        print(f"❌ Simulation directory not found: {args.sim_dir}")
        return 1
    
    print("📊 Starting sampled bridging analysis (10_30 format)...")
    
    try:
        # Calculate sampled bridging metrics
        results, behavior_df = calculate_sampled_bridging(args.infections_csv, args.sim_dir)
        
        # Create output tables (matching 10_30 format)
        results_df = pd.DataFrame(results)
        
        # Create separate pivot tables for values and counts
        value_pivot = results_df.pivot(index='measure', columns='grouping', values='value')
        count_pivot = results_df.pivot(index='measure', columns='grouping', values='count')
        total_pivot = results_df.pivot(index='measure', columns='grouping', values='total')
        
        # Save all tables
        value_pivot.to_csv(args.output)
        count_pivot.to_csv(args.output.replace('.csv', '_counts.csv'))
        total_pivot.to_csv(args.output.replace('.csv', '_totals.csv'))
        
        # Save detailed behaviors
        behavior_df.to_csv(args.detailed_output, index=False)
        
        print(f"\n✅ Results saved to: {args.output}")
        print(f"✅ Counts saved to: {args.output.replace('.csv', '_counts.csv')}")
        print(f"✅ Totals saved to: {args.output.replace('.csv', '_totals.csv')}")  
        print(f"✅ Detailed behaviors saved to: {args.detailed_output}")
        
        # Print summary tables
        print(f"\n📊 Final Results Summary (3 measures × 2 groupings = 6 total):")
        print(value_pivot)
        
        return 0
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())
