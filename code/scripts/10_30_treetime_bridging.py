#!/usr/bin/env python3
"""
TreeTime Bridging Analysis - Matching Ground Truth Methods
Calculates bridging from TreeTime output using 6 measures × 2 groupings = 12 total measurements
"""

import pandas as pd
import numpy as np
import argparse
import os
import re

def extract_lsd2_rate(lsd2_result_file):
    """Extract substitution rate from LSD2 output file"""
    
    if not os.path.exists(lsd2_result_file):
        raise FileNotFoundError(f"LSD2 result file not found: {lsd2_result_file}")
    
    with open(lsd2_result_file, 'r') as f:
        content = f.read()
    
    # Look for rate in format: "rate 2.11241e-06"
    rate_match = re.search(r'rate\s+([0-9\.e\-\+]+)', content)
    
    if not rate_match:
        raise ValueError("Could not find substitution rate in LSD2 output")
    
    rate = float(rate_match.group(1))
    print(f"Extracted LSD2 substitution rate: {rate:.2e} substitutions per day")
    
    return rate

def parse_treetime_nexus(tree_file):
    """Parse TreeTime Nexus output with state annotations - IMPROVED"""
    print(f"Parsing TreeTime Nexus file: {tree_file}")
    
    with open(tree_file, 'r') as f:
        content = f.read()
    
    # Look for the tree line more carefully
    lines = content.strip().split('\n')
    tree_string = None
    
    for line in lines:
        line = line.strip()
        # Skip comments and headers
        if line.startswith('#') or line.startswith('Begin') or line.startswith('End') or not line:
            continue
        # Look for the actual tree (contains parentheses and ends with semicolon)
        if '(' in line and line.endswith(';'):
            tree_string = line
            break
    
    if not tree_string:
        # Fallback: try to find tree in original regex approach
        tree_match = re.search(r'([^;]*\([^;]*\);)', content)
        if tree_match:
            tree_string = tree_match.group(1)
    
    print(f"Extracted tree string: {len(tree_string) if tree_string else 0} characters")
    return tree_string

def extract_node_info_from_nexus(tree_string, substitution_rate):
    """Extract node states and convert branch lengths to calendar time"""
    
    # Find all nodes with states: nodename:branchlength[&state="STATE"]
    pattern = r'([^,\(\):]+):([0-9\.]+)\[&state="([^"]+)"\]'
    matches = re.findall(pattern, tree_string)
    
    node_states = {}
    node_times_days = {}
    
    for node_name, branch_length_evo, state in matches:
        # Convert evolutionary time to calendar time
        branch_length_days = float(branch_length_evo) / substitution_rate
        
        node_states[node_name] = state
        node_times_days[node_name] = branch_length_days
    
    print(f"Extracted states for {len(node_states)} nodes")
    print(f"Converted branch lengths from evolutionary units to days")
    
    return node_states, node_times_days

def get_network_groups(state, grouping):
    """Map TreeTime combined states to consistent WSM_group vs MSM_group"""
    
    if grouping == "MSM+MSMW":
        # For MSMW+MSM vs WSM+MSW analysis
        if state == 'MSMW_MSM':  # This represents the MSM+MSMW group
            return 'MSM_group'
        elif state == 'WSM_MSW':  # This represents the WSM+MSW group  
            return 'WSM_group'
    
    elif grouping == "MSMW+MSW":
        # For MSM vs MSMW+MSW+WSM analysis  
        if state == 'MSM':  # Pure MSM group
            return 'MSM_group'
        elif state == 'MSMW_MSW_WSM':  # The MSMW+MSW+WSM group
            return 'WSM_group'
    
    return 'UNKNOWN'

def calculate_between_network_proportions(transitions_df):
    """Calculate between network proportions (matching ground truth format)"""
    
    if transitions_df.empty:
        return {
            'overall_proportion': 0,
            'overall_count': 0,
            'WSM_to_MSM_proportion': 0,
            'WSM_to_MSM_count': 0,
            'MSM_to_WSM_proportion': 0,
            'MSM_to_WSM_count': 0,
            'total': 0
        }
    
    # Count transitions
    total_branches = len(transitions_df)
    between_network_transitions = transitions_df['transition'].sum()
    
    # Count directional transitions
    WSM_to_MSM = len(transitions_df[
        (transitions_df['parent_network'] == 'WSM_group') & 
        (transitions_df['child_network'] == 'MSM_group')
    ])
    
    MSM_to_WSM = len(transitions_df[
        (transitions_df['parent_network'] == 'MSM_group') & 
        (transitions_df['child_network'] == 'WSM_group')
    ])
    
    # Calculate proportions
    overall_proportion = between_network_transitions / total_branches if total_branches > 0 else 0
    WSM_to_MSM_proportion = WSM_to_MSM / total_branches if total_branches > 0 else 0
    MSM_to_WSM_proportion = MSM_to_WSM / total_branches if total_branches > 0 else 0
    
    return {
        'overall_proportion': overall_proportion,
        'overall_count': between_network_transitions,
        'WSM_to_MSM_proportion': WSM_to_MSM_proportion,
        'WSM_to_MSM_count': WSM_to_MSM,
        'MSM_to_WSM_proportion': MSM_to_WSM_proportion,
        'MSM_to_WSM_count': MSM_to_WSM,
        'total': total_branches
    }

def calculate_average_time_in_network(network_segments):
    """Calculate average time spent in networks before transitions (matching ground truth)"""
    
    if not network_segments:
        return {
            'overall_avg': 0,
            'overall_count': 0,
            'WSM_network_avg': 0,
            'WSM_network_count': 0,
            'MSM_network_avg': 0,
            'MSM_network_count': 0
        }
    
    # Extract time durations by network
    all_durations = [seg['time_days'] for seg in network_segments if seg['time_days'] > 0]
    WSM_durations = [seg['time_days'] for seg in network_segments 
                     if seg['network'] == 'WSM_group' and seg['time_days'] > 0]
    MSM_durations = [seg['time_days'] for seg in network_segments 
                     if seg['network'] == 'MSM_group' and seg['time_days'] > 0]
    
    # Calculate means (matching ground truth approach)
    avg_overall = np.mean(all_durations) if all_durations else 0
    avg_WSM_network = np.mean(WSM_durations) if WSM_durations else 0
    avg_MSM_network = np.mean(MSM_durations) if MSM_durations else 0
    
    return {
        'overall_avg': avg_overall,
        'overall_count': len(all_durations),
        'WSM_network_avg': avg_WSM_network,
        'WSM_network_count': len(WSM_durations),
        'MSM_network_avg': avg_MSM_network,
        'MSM_network_count': len(MSM_durations)
    }

def parse_transitions_from_raw_tree(tree_string, node_states, node_times_days, grouping):
    """Parse transitions directly from the raw tree string using regex"""
    
    print(f"Parsing transitions directly from tree string for {grouping}...")
    
    transitions = []
    network_segments = []
    
    # Extract all node information with states from the tree string
    # Pattern: nodename:branchlength[&state="STATE"]
    pattern = r'([^,\(\):]+):([0-9\.]+)\[&state="([^"]+)"\]'
    matches = re.findall(pattern, tree_string)
    
    print(f"  Found {len(matches)} nodes with states in tree")
    
    # Create network segments
    for node_name, branch_length, state in matches:
        network = get_network_groups(state, grouping)
        time_days = node_times_days.get(node_name, 0)
        
        if network != 'UNKNOWN':
            # Determine if it's a terminal or internal node
            node_type = 'internal' if node_name.startswith('NODE_') else 'terminal'
            
            network_segments.append({
                'node': node_name,
                'state': state,
                'network': network,
                'time_days': time_days,
                'node_type': node_type
            })
    
    # Find parent-child relationships by analyzing the tree structure
    # Look for patterns like: (child1,child2)parent or child)parent
    
    # Split on major delimiters and find parent-child patterns
    parent_child_pattern = r'\(([^)]+)\)([^,\(\):]+:[0-9\.]+\[&state="[^"]+"\])'
    parent_matches = re.findall(parent_child_pattern, tree_string)
    
    print(f"  Found {len(parent_matches)} potential parent-child groups")
    
    for children_str, parent_str in parent_matches:
        # Extract parent info
        parent_match = re.search(r'([^:]+):[0-9\.]+\[&state="([^"]+)"\]', parent_str)
        if not parent_match:
            continue
            
        parent_name = parent_match.group(1)
        parent_state = parent_match.group(2)
        parent_network = get_network_groups(parent_state, grouping)
        
        if parent_network == 'UNKNOWN':
            continue
        
        # Extract children info
        child_matches = re.findall(r'([^,\(\):]+):[0-9\.]+\[&state="([^"]+)"\]', children_str)
        
        for child_name, child_state in child_matches:
            child_network = get_network_groups(child_state, grouping)
            child_time = node_times_days.get(child_name, 0)
            
            if child_network != 'UNKNOWN':
                transition_occurred = (parent_network != child_network)
                transitions.append({
                    'parent_node': parent_name,
                    'child_node': child_name,
                    'parent_state': parent_state,
                    'child_state': child_state,
                    'parent_network': parent_network,
                    'child_network': child_network,
                    'transition': transition_occurred,
                    'time_days': child_time,
                    'grouping': grouping
                })
    
    print(f"  Found {len(transitions)} parent-child relationships")
    if transitions:
        transition_count = sum(1 for t in transitions if t['transition'])
        print(f"  Actual transitions: {transition_count}/{len(transitions)}")
    
    # Debug: show some examples
    if transitions:
        print(f"  Sample transitions:")
        for i, t in enumerate(transitions[:3]):
            print(f"    {t['parent_node']}({t['parent_network']}) -> {t['child_node']}({t['child_network']}) = {t['transition']}")
    
    return pd.DataFrame(transitions), network_segments

def analyze_single_treetime_file(nexus_file, substitution_rate, analysis_name, grouping):
    """Analyze a single TreeTime file - USING DIRECT PARSING"""
    
    print(f"\n=== Analyzing {analysis_name} for {grouping} grouping ===")
    
    # Parse TreeTime output
    tree_string = parse_treetime_nexus(nexus_file)
    node_states, node_times_days = extract_node_info_from_nexus(tree_string, substitution_rate)
    
    if not node_states:
        print(f"No states found in {nexus_file}")
        return []
    
    # Use direct parsing instead of dendropy
    transitions_df, network_segments = parse_transitions_from_raw_tree(
        tree_string, node_states, node_times_days, grouping
    )
    
    # Rest stays the same...
    between_results = calculate_between_network_proportions(transitions_df)
    time_results = calculate_average_time_in_network(network_segments)
    
    # [Rest of function unchanged - the results formatting part]
    results = [
        {"measure": "total_between_network_proportion", "grouping": grouping, 
         "value": between_results['overall_proportion'], 
         "count": between_results['overall_count'], 
         "total": between_results['total'],
         "analysis": analysis_name},
        
        {"measure": "WSM_group_to_MSM_group_proportion", "grouping": grouping,
         "value": between_results['WSM_to_MSM_proportion'],
         "count": between_results['WSM_to_MSM_count'],
         "total": between_results['total'],
         "analysis": analysis_name},
        
        {"measure": "MSM_group_to_WSM_group_proportion", "grouping": grouping,
         "value": between_results['MSM_to_WSM_proportion'],
         "count": between_results['MSM_to_WSM_count'],
         "total": between_results['total'],
         "analysis": analysis_name},
        
        {"measure": "average_time_in_network_overall", "grouping": grouping,
         "value": time_results['overall_avg'],
         "count": time_results['overall_count'],
         "total": time_results['overall_count'],
         "analysis": analysis_name},
         
        {"measure": "average_time_in_WSM_network", "grouping": grouping,
         "value": time_results['WSM_network_avg'],
         "count": time_results['WSM_network_count'],
         "total": time_results['WSM_network_count'],
         "analysis": analysis_name},
         
        {"measure": "average_time_in_MSM_network", "grouping": grouping,
         "value": time_results['MSM_network_avg'],
         "count": time_results['MSM_network_count'],
         "total": time_results['MSM_network_count'],
         "analysis": analysis_name}
    ]
    
    # Print summary
    print(f"  Between network - Total: {between_results['overall_proportion']:.4f} ({between_results['overall_count']}/{between_results['total']})")
    print(f"  Between network - WSM→MSM: {between_results['WSM_to_MSM_proportion']:.4f} ({between_results['WSM_to_MSM_count']}/{between_results['total']})")
    print(f"  Between network - MSM→WSM: {between_results['MSM_to_WSM_proportion']:.4f} ({between_results['MSM_to_WSM_count']}/{between_results['total']})")
    
    print(f"  Time in network - Overall: {time_results['overall_avg']:.1f} days ({time_results['overall_count']} segments)")
    print(f"  Time in network - WSM: {time_results['WSM_network_avg']:.1f} days ({time_results['WSM_network_count']} segments)")
    print(f"  Time in network - MSM: {time_results['MSM_network_avg']:.1f} days ({time_results['MSM_network_count']} segments)")
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description='Calculate bridging metrics from TreeTime matching ground truth methods'
    )
    
    parser.add_argument('msmw_msm_nexus', help='TreeTime Nexus file for MSMW+MSM analysis')
    parser.add_argument('msmw_msw_nexus', help='TreeTime Nexus file for MSMW+MSW analysis')
    parser.add_argument('lsd2_result_file', help='LSD2 result file (for substitution rate)')
    parser.add_argument('output_prefix', help='Output file prefix')
    
    args = parser.parse_args()
    
    # Validate inputs
    for file_path in [args.msmw_msm_nexus, args.msmw_msw_nexus, args.lsd2_result_file]:
        if not os.path.exists(file_path):
            print(f"❌ File not found: {file_path}")
            return 1
    
    print("🌳 Starting TreeTime bridging analysis matching ground truth methods...")
    
    try:
        # Extract substitution rate from LSD2
        substitution_rate = extract_lsd2_rate(args.lsd2_result_file)
        
        all_results = []
        
        # FIXED: Each TreeTime file corresponds to one specific grouping
        treetime_analyses = [
            (args.msmw_msm_nexus, "MSMW+MSM_vs_WSM+MSW", "MSM+MSMW"),
            (args.msmw_msw_nexus, "MSM_vs_MSMW+MSW+WSM", "MSMW+MSW")
        ]
        
        for nexus_file, analysis_name, grouping in treetime_analyses:
            results = analyze_single_treetime_file(nexus_file, substitution_rate, analysis_name, grouping)
            all_results.extend(results)
        
        if not all_results:
            raise ValueError("No valid results from either analysis")
        
        # Create results DataFrame
        results_df = pd.DataFrame(all_results)
        
        # Create pivot tables (matching ground truth format)
        value_pivot = results_df.pivot_table(
            index='measure', 
            columns='grouping', 
            values='value'
        )
        count_pivot = results_df.pivot_table(
            index='measure', 
            columns='grouping', 
            values='count'
        )
        total_pivot = results_df.pivot_table(
            index='measure', 
            columns='grouping', 
            values='total'
        )
        
        # Save results
        value_file = f"{args.output_prefix}_treetime_bridging.csv"
        count_file = f"{args.output_prefix}_treetime_bridging_counts.csv"
        total_file = f"{args.output_prefix}_treetime_bridging_totals.csv"
        
        value_pivot.to_csv(value_file)
        count_pivot.to_csv(count_file)
        total_pivot.to_csv(total_file)
        
        print(f"\n🎉 TreeTime bridging analysis completed!")
        print(f"✅ Results saved to: {value_file}")
        print(f"✅ Counts saved to: {count_file}")
        print(f"✅ Totals saved to: {total_file}")
        
        print(f"\n📊 Final Results Summary (6 measures × 2 groupings = 12 total):")
        print(value_pivot)
        
        return 0
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())
