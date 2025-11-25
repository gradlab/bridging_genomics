#!/usr/bin/env python3
"""
Ground Truth Bridging Analysis - Combined Best Methods
Calculates bridging in simulation data using 9 measures × 2 groupings = 18 total measurements
"""

import pandas as pd
import numpy as np
import argparse
import os
from pathlib import Path
from collections import defaultdict
import json

def load_data(output_dir, classifications_file):
    """Load all required data files"""
    print("Loading data files...")
    
    # Load transmission data
    transmission_df = pd.read_csv(os.path.join(output_dir, "transmission_df.csv"))
    
    # Load behavior classifications
    classifications_df = pd.read_csv(classifications_file)
    
    # Load parameters for tracking start day
    with open(os.path.join(output_dir, "parameters_used.json"), 'r') as f:
        params = json.load(f)
    tracking_start_day = params['simulation']['partnership_burnin_days'] + params['simulation']['transmission_burnin_days']
    
    print(f"Loaded: {len(transmission_df)} transmissions, {len(classifications_df)} classifications")
    return transmission_df, classifications_df, tracking_start_day

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
    """Map behaviors to consistent WSM_group vs MSM_group"""
    
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

def calculate_between_network_proportion(transmission_df, behavior_map, tracking_start_day, grouping):
    """Calculate between network transmission proportions with consistent WSM/MSM directionality"""
    
    # Filter to post-burn-in successful transmissions
    post_burnin = transmission_df[
        (transmission_df['day_of_transmission'] >= tracking_start_day) &
        (transmission_df['superseded_simultaneous'] == False) &
        (transmission_df['day_of_sampling'].notna())
    ]
    
    between_network_count = 0
    WSM_to_MSM_count = 0  
    MSM_to_WSM_count = 0  
    
    for _, row in post_burnin.iterrows():
        # Get dynamic behavior for infector
        infector_tx_id = get_infector_transmission_id(row['infector_node'], row['day_of_transmission'], post_burnin)
        
        if infector_tx_id is not None and infector_tx_id in behavior_map:
            infector_behavior = behavior_map[infector_tx_id]  # Dynamic behavior
        else:
            infector_behavior = row['behavior_infector']  # Fallback to behavior at transmission
        
        # Get dynamic behavior for infectee
        if row['transmission_id'] in behavior_map:
            infectee_behavior = behavior_map[row['transmission_id']]  # Dynamic behavior
        else:
            infectee_behavior = row['behavior_infectee']  # Fallback
        
        # Get consistent network groups
        infector_group = get_network_groups(infector_behavior, grouping)
        infectee_group = get_network_groups(infectee_behavior, grouping)
        
        if infector_group != infectee_group:
            between_network_count += 1
            
            # Track consistent directionality
            if infector_group == "WSM_group" and infectee_group == "MSM_group":
                WSM_to_MSM_count += 1  # WSM_group → MSM_group
            elif infector_group == "MSM_group" and infectee_group == "WSM_group":
                MSM_to_WSM_count += 1  # MSM_group → WSM_group
    
    # Calculate proportions
    total = len(post_burnin)
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

def calculate_triplet_proportion(transmission_df, behavior_map, tracking_start_day, grouping):
    """Calculate triplet (flanked between-network) proportions with consistent WSM/MSM directionality"""
    
    # Filter to post-burn-in successful transmissions
    post_burnin = transmission_df[
        (transmission_df['day_of_transmission'] >= tracking_start_day) &
        (transmission_df['superseded_simultaneous'] == False) &
        (transmission_df['day_of_sampling'].notna())
    ]
    
    # Build transmission network maps
    infector_map = {}  # infectee -> infector info
    infectee_map = defaultdict(list)  # infector -> list of infectees
    
    for _, row in post_burnin.iterrows():
        # Get dynamic behavior for infector
        infector_tx_id = get_infector_transmission_id(
            row['infector_node'], 
            row['day_of_transmission'], 
            post_burnin
        )
        if infector_tx_id is not None:
            infector_behavior = behavior_map.get(infector_tx_id, row['behavior_infector'])
        else:
            infector_behavior = row['behavior_infector']  # Fallback
        
        # Get dynamic behavior for infectee
        infectee_behavior = behavior_map.get(row['transmission_id'], row['behavior_infectee'])
        
        infector_node = row['infector_node']
        infectee_node = row['infectee_node']
        
        infector_map[infectee_node] = {
            'behavior': infector_behavior,
            'tx_id': row['transmission_id']
        }
        
        infectee_map[infector_node].append({
            'behavior': infectee_behavior,
            'tx_id': row['transmission_id']
        })
    
    # Find between-network transmissions that are flanked by within-group transmissions
    flanked_between_network_count = 0
    flanked_WSM_to_MSM_count = 0  
    flanked_MSM_to_WSM_count = 0  
    
    for _, row in post_burnin.iterrows():
        # Get dynamic behavior for infector
        infector_tx_id = get_infector_transmission_id(
            row['infector_node'], 
            row['day_of_transmission'], 
            post_burnin
        )
        if infector_tx_id is not None:
            infector_behavior = behavior_map.get(infector_tx_id, row['behavior_infector'])
        else:
            infector_behavior = row['behavior_infector']  # Fallback
        
        # Get dynamic behavior for infectee
        infectee_behavior = behavior_map.get(row['transmission_id'], row['behavior_infectee'])
        
        infector_group = get_network_groups(infector_behavior, grouping)
        infectee_group = get_network_groups(infectee_behavior, grouping)
        
        # Check if this is a between-network transmission
        if infector_group != infectee_group:
            infector_node = row['infector_node']
            infectee_node = row['infectee_node']
            
            # Check for within-group transmission BEFORE
            infector_source = infector_map.get(infector_node)
            within_group_before = (infector_source and 
                                 get_network_groups(infector_source['behavior'], grouping) == infector_group)
            
            # Check for within-group transmission AFTER
            within_group_after = any(
                get_network_groups(target['behavior'], grouping) == infectee_group
                for target in infectee_map.get(infectee_node, [])
            )
            
            # Only count if flanked by within-group transmissions
            if within_group_before and within_group_after:
                flanked_between_network_count += 1
                
                # Track consistent directionality
                if infector_group == "WSM_group" and infectee_group == "MSM_group":
                    flanked_WSM_to_MSM_count += 1  # WSM_group → MSM_group
                elif infector_group == "MSM_group" and infectee_group == "WSM_group":
                    flanked_MSM_to_WSM_count += 1  # MSM_group → WSM_group
    
    # Calculate proportions
    total = len(post_burnin)
    overall_proportion = flanked_between_network_count / total if total > 0 else 0
    WSM_to_MSM_proportion = flanked_WSM_to_MSM_count / total if total > 0 else 0
    MSM_to_WSM_proportion = flanked_MSM_to_WSM_count / total if total > 0 else 0
    
    return {
        'overall_proportion': overall_proportion,
        'overall_count': flanked_between_network_count,
        'WSM_to_MSM_proportion': WSM_to_MSM_proportion,
        'WSM_to_MSM_count': flanked_WSM_to_MSM_count,
        'MSM_to_WSM_proportion': MSM_to_WSM_proportion,
        'MSM_to_WSM_count': flanked_MSM_to_WSM_count,
        'total': total
    }

def calculate_average_time_in_network(transmission_df, behavior_map, tracking_start_day, grouping):
    """Calculate average time spent in networks before transitions using transmission chains"""
    
    print(f"Analyzing network residence times for {grouping}...")
    
    # Get post-burn-in successful transmissions
    post_burnin = transmission_df[
        (transmission_df['day_of_transmission'] >= tracking_start_day) &
        (transmission_df['superseded_simultaneous'] == False) &
        (transmission_df['day_of_sampling'].notna())
    ].copy()
    
    print(f"  Working with {len(post_burnin)} post-burn-in transmissions")
    
    # Create lookup for faster searching
    infectee_lookup = {}
    for _, row in post_burnin.iterrows():
        infectee = row['infectee_node']
        if infectee not in infectee_lookup:
            infectee_lookup[infectee] = []
        infectee_lookup[infectee].append(row)
    
    # Sort each infectee's infections by day
    for infectee in infectee_lookup:
        infectee_lookup[infectee].sort(key=lambda x: x['day_of_transmission'])
    
    network_durations = []
    chains_processed = 0
    
    # For each infection, build transmission chain going backwards
    for _, infection in post_burnin.iterrows():
        chain = []
        current_infection = infection
        
        # Build chain by following infector relationships backwards
        max_chain_length = 100  # Prevent infinite loops
        for step in range(max_chain_length):
            
            # Get behavior for current infection
            if current_infection['transmission_id'] not in behavior_map:
                # Skip if no dynamic behavior available
                break
            current_behavior = behavior_map[current_infection['transmission_id']]  # Dynamic only
            
            # Get consistent network group
            current_network = get_network_groups(current_behavior, grouping)
            
            if current_network == 'UNKNOWN':
                break
            
            # Add to chain
            chain.append({
                'node': current_infection['infectee_node'],
                'day': current_infection['day_of_transmission'],
                'sampling_day': current_infection['day_of_sampling'],
                'behavior': current_behavior,
                'network': current_network,
                'transmission_id': current_infection['transmission_id']
            })
            
            # Find who infected the current infector
            infector_node = current_infection['infector_node']
            
            # Look for infection of the infector
            infector_infections = infectee_lookup.get(infector_node, [])
            
            # Find the most recent infection of the infector before current transmission
            next_infection = None
            current_tx_day = current_infection['day_of_transmission']
            
            for inf_infection in reversed(infector_infections):  # Most recent first
                if inf_infection['day_of_transmission'] < current_tx_day:
                    # Check if this infection was still active when current transmission occurred
                    if (pd.notna(inf_infection['day_of_sampling']) and 
                        inf_infection['day_of_sampling'] >= current_tx_day):
                        next_infection = inf_infection
                        break
            
            if next_infection is None:
                break  # End of chain
            
            current_infection = next_infection
        
        # Analyze network transitions in this chain
        if len(chain) > 1:
            # Reverse chain so it goes chronologically (oldest to newest)
            chain.reverse()
            
            # Find network transitions
            current_network = None
            network_start_day = None
            
            for i, link in enumerate(chain):
                if current_network is None:
                    current_network = link['network']
                    network_start_day = link['day']
                elif link['network'] != current_network:
                    # Network transition occurred
                    duration = link['day'] - network_start_day
                    if duration > 0:
                        network_durations.append({
                            'from_network': current_network,
                            'to_network': link['network'],
                            'duration': duration,
                            'start_day': network_start_day,
                            'end_day': link['day'],
                            'chain_id': chains_processed,
                            'chain_length': len(chain),
                            'grouping': grouping,
                            'transition_position': i  # Position in chain where transition occurred
                        })
                    
                    current_network = link['network']
                    network_start_day = link['day']
        
        chains_processed += 1
        
        if chains_processed % 1000 == 0:
            print(f"    Processed {chains_processed} chains...")
    
    # Calculate average durations
    all_durations = [d['duration'] for d in network_durations]
    WSM_network_durations = [d['duration'] for d in network_durations if d['from_network'] == 'WSM_group']
    MSM_network_durations = [d['duration'] for d in network_durations if d['from_network'] == 'MSM_group']
    
    avg_overall = np.mean(all_durations) if all_durations else 0
    avg_WSM_network = np.mean(WSM_network_durations) if WSM_network_durations else 0
    avg_MSM_network = np.mean(MSM_network_durations) if MSM_network_durations else 0
    
    print(f"  Processed {chains_processed} transmission chains")
    print(f"  Found {len(network_durations)} network transitions")
    print(f"  Average time in network overall: {avg_overall:.1f} days")
    print(f"  Average time in WSM_group: {avg_WSM_network:.1f} days ({len(WSM_network_durations)} transitions)")
    print(f"  Average time in MSM_group: {avg_MSM_network:.1f} days ({len(MSM_network_durations)} transitions)")
    
    return {
        'overall_avg': avg_overall,
        'overall_count': len(all_durations),
        'WSM_network_avg': avg_WSM_network,
        'WSM_network_count': len(WSM_network_durations),
        'MSM_network_avg': avg_MSM_network,
        'MSM_network_count': len(MSM_network_durations),
        'network_transitions': network_durations
    }

def main():
    parser = argparse.ArgumentParser(description='Calculate comprehensive ground truth bridging measures')
    parser.add_argument('output_dir', help='Simulation output directory')
    parser.add_argument('classifications_file', help='Infection behavior classifications file (.csv)')
    parser.add_argument('--output', '-o', default='bridging_analysis.csv', 
                       help='Output file for results table')
    parser.add_argument('--detailed-output', '-d', default='bridging_detailed.csv',
                       help='Output file for detailed network transitions')
    args = parser.parse_args()

    try:
        print("🌉 Starting comprehensive ground truth bridging analysis...")
        
        # Load data
        transmission_df, classifications_df, tracking_start_day = load_data(
            args.output_dir, args.classifications_file
        )
        
        # Create behavior mapping
        behavior_map = get_dynamic_behavior_map(classifications_df, transmission_df)
        
        # Initialize results table
        results = []
        detailed_results = []
        
        # Calculate all 9 measures for both groupings
        for grouping in ["MSM+MSMW", "MSMW+MSW"]:
            print(f"\nCalculating 9 measures for {grouping} grouping...")
            
            # Measures 1-3: Between network proportions
            between_results = calculate_between_network_proportion(
                transmission_df, behavior_map, tracking_start_day, grouping
            )
            
            # Measures 4-6: Triplet proportions
            triplet_results = calculate_triplet_proportion(
                transmission_df, behavior_map, tracking_start_day, grouping
            )
            
            # Measures 7-9: Average time in network
            time_results = calculate_average_time_in_network(
                transmission_df, behavior_map, tracking_start_day, grouping
            )
            
            # Add all 9 results
            results.extend([
                # Between network measures (1-3)
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
                 "total": between_results['total']},
                
                # Triplet measures (4-6)
                {"measure": "total_triplet_proportion", "grouping": grouping,
                 "value": triplet_results['overall_proportion'],
                 "count": triplet_results['overall_count'],
                 "total": triplet_results['total']},
                 
                {"measure": "WSM_group_to_MSM_group_triplet_proportion", "grouping": grouping,
                 "value": triplet_results['WSM_to_MSM_proportion'],
                 "count": triplet_results['WSM_to_MSM_count'],
                 "total": triplet_results['total']},
                 
                {"measure": "MSM_group_to_WSM_group_triplet_proportion", "grouping": grouping,
                 "value": triplet_results['MSM_to_WSM_proportion'],
                 "count": triplet_results['MSM_to_WSM_count'],
                 "total": triplet_results['total']},
                
                # Time in network measures (7-9)
                {"measure": "average_time_in_network_overall", "grouping": grouping,
                 "value": time_results['overall_avg'],
                 "count": time_results['overall_count'],
                 "total": time_results['overall_count']},
                 
                {"measure": "average_time_in_WSM_network", "grouping": grouping,
                 "value": time_results['WSM_network_avg'],
                 "count": time_results['WSM_network_count'],
                 "total": time_results['WSM_network_count']},
                 
                {"measure": "average_time_in_MSM_network", "grouping": grouping,
                 "value": time_results['MSM_network_avg'],
                 "count": time_results['MSM_network_count'],
                 "total": time_results['MSM_network_count']}
            ])
            
            # Extend detailed results with network transition data
            detailed_results.extend(time_results['network_transitions'])
            
            # Print summary for this grouping
            print(f"  Between network - Total: {between_results['overall_proportion']:.4f} ({between_results['overall_count']}/{between_results['total']})")
            print(f"  Between network - WSM→MSM: {between_results['WSM_to_MSM_proportion']:.4f} ({between_results['WSM_to_MSM_count']}/{between_results['total']})")
            print(f"  Between network - MSM→WSM: {between_results['MSM_to_WSM_proportion']:.4f} ({between_results['MSM_to_WSM_count']}/{between_results['total']})")
            
            print(f"  Triplet - Total: {triplet_results['overall_proportion']:.4f} ({triplet_results['overall_count']}/{triplet_results['total']})")
            print(f"  Triplet - WSM→MSM: {triplet_results['WSM_to_MSM_proportion']:.4f} ({triplet_results['WSM_to_MSM_count']}/{triplet_results['total']})")
            print(f"  Triplet - MSM→WSM: {triplet_results['MSM_to_WSM_proportion']:.4f} ({triplet_results['MSM_to_WSM_count']}/{triplet_results['total']})")
            
            print(f"  Time in network - Overall: {time_results['overall_avg']:.1f} days ({time_results['overall_count']} transitions)")
            print(f"  Time in network - WSM: {time_results['WSM_network_avg']:.1f} days ({time_results['WSM_network_count']} transitions)")
            print(f"  Time in network - MSM: {time_results['MSM_network_avg']:.1f} days ({time_results['MSM_network_count']} transitions)")
        
        # Create output tables
        results_df = pd.DataFrame(results)
        
        # Create separate pivot tables for values and counts
        value_pivot = results_df.pivot(index='measure', columns='grouping', values='value')
        count_pivot = results_df.pivot(index='measure', columns='grouping', values='count')
        total_pivot = results_df.pivot(index='measure', columns='grouping', values='total')
        
        # Save all tables
        value_pivot.to_csv(args.output)
        count_pivot.to_csv(args.output.replace('.csv', '_counts.csv'))
        total_pivot.to_csv(args.output.replace('.csv', '_totals.csv'))
        
        print(f"\n✅ Results saved to: {args.output}")
        print(f"✅ Counts saved to: {args.output.replace('.csv', '_counts.csv')}")
        print(f"✅ Totals saved to: {args.output.replace('.csv', '_totals.csv')}")
        
        # Save detailed network transitions
        if detailed_results:
            detailed_df = pd.DataFrame(detailed_results)
            detailed_df.to_csv(args.detailed_output, index=False)
            print(f"✅ Detailed network transitions saved to: {args.detailed_output}")
        
        # Print final summary tables
        print(f"\n📊 Final Results Summary (9 measures × 2 groupings = 18 total):")
        print(value_pivot)
        
        return 0
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())
