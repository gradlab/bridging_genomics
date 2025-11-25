#!/usr/bin/env python3
"""
Parameter Interpretation Script
Reads a simulation parameters JSON file and creates an interpretable summary
"""

import json
import pandas as pd
import numpy as np
from scipy.stats import gamma
import argparse
import sys

def shifted_gamma_median(shape, scale, min_rate):
    """Calculate median of shifted gamma distribution"""
    return gamma.ppf(0.5, a=shape, scale=scale) + min_rate

def regular_gamma_median(shape, scale):
    """Calculate median of regular gamma distribution"""
    return gamma.ppf(0.5, a=shape, scale=scale)

def interpret_parameters(json_file):
    """Read JSON parameters and create interpretable summary"""
    
    # Load parameters
    with open(json_file, 'r') as f:
        params = json.load(f)
    
    print(f"📊 Interpreting parameters from: {json_file}")
    
    # Extract and calculate interpretable values
    interpretable_data = []
    
    # =================================================================
    # BEHAVIOR COMPOSITION
    # =================================================================
    behaviors = params['behaviors']
    interpretable_data.extend([
        ('WSM proportion', behaviors['WSM']),
        ('MSW proportion', behaviors['MSW']),
        ('MSM proportion', behaviors['MSM']),
        ('MSMW proportion', behaviors['MSMW'])
    ])
    
    # =================================================================
    # RISK ACTIVITY LEVELS & PARTNER ACQUISITION RATES
    # =================================================================
    risk_dist = params['risk_distribution']
    
    for behavior in ['WSM', 'MSW', 'MSM', 'MSMW']:
        behavior_risk = risk_dist[behavior]
        
        # High activity percentage
        interpretable_data.append((f'{behavior} % high activity', behavior_risk['prop_hi'] * 100))
        
        # Low activity median (shifted gamma)
        lo_median = shifted_gamma_median(
            behavior_risk['lo_par']['shape'],
            behavior_risk['lo_par']['scale'], 
            behavior_risk['lo_par']['min_rate']
        )
        interpretable_data.append((f'{behavior} low activity median', lo_median))
        
        # High activity median (shifted gamma)
        hi_median = shifted_gamma_median(
            behavior_risk['hi_par']['shape'],
            behavior_risk['hi_par']['scale'],
            behavior_risk['hi_par']['min_rate']
        )
        interpretable_data.append((f'{behavior} high activity median', hi_median))
    
    # =================================================================
    # PARTNERSHIP PARAMETERS
    # =================================================================
    partnerships = params['partnerships']
    
    # Basic partnership parameters
    interpretable_data.extend([
        ('Casual rate', partnerships['casual_rate']),
        ('p_MSMW_W', partnerships['p_MSMW_W'])
    ])
    
    # Steady partnership duration (regular gamma)
    steady_duration_median = regular_gamma_median(
        partnerships['steady_shape'],
        partnerships['steady_scale']
    )
    interpretable_data.append(('Steady duration median (days)', steady_duration_median))
    
    # Capacity mix
    capacity_mix = partnerships['steady_capacity_mix']
    interpretable_data.extend([
        ('Hi activity: 0 steady partners', capacity_mix['hi']['0']),
        ('Hi activity: 1 steady partner', capacity_mix['hi']['1']),
        ('Hi activity: 2 steady partners', capacity_mix['hi']['2']),
        ('Lo activity: 0 steady partners', capacity_mix['lo']['0']),
        ('Lo activity: 1 steady partner', capacity_mix['lo']['1']),
        ('Lo activity: 2 steady partners', capacity_mix['lo']['2'])
    ])
    
    # =================================================================
    # TRANSMISSION PROBABILITIES
    # =================================================================
    transmission_probs = params['infections']['transmission_probs_per_partner_per_day']
    
    # Organize transmission probabilities
    steady_probs = [(k, v) for k, v in transmission_probs.items() if 'steady' in k]
    casual_probs = [(k, v) for k, v in transmission_probs.items() if 'casual' in k]
    
    for prob_name, prob_value in steady_probs:
        interpretable_data.append((f'{prob_name} transmission prob', prob_value))
    
    for prob_name, prob_value in casual_probs:
        interpretable_data.append((f'{prob_name} transmission prob', prob_value))
    
    # =================================================================
    # SYMPTOMATIC PROBABILITIES
    # =================================================================
    p_sx = params['infections']['p_sx']
    for behavior, prob in p_sx.items():
        interpretable_data.append((f'p_sx {behavior}', prob))
    
    # =================================================================
    # INFECTION DURATIONS
    # =================================================================
    duration_params = params['infections']['infxn_duration_params_dict']
    
    for behavior in ['WSM', 'MSW', 'MSM', 'MSMW']:
        for activity in ['lo', 'hi']:
            for status in ['symptomatic', 'asymptomatic', 'natural_clearance']:
                duration_data = duration_params[behavior][activity][status]
                
                # Regular gamma (no minimum)
                duration_median = regular_gamma_median(
                    duration_data['shape'],
                    duration_data['scale']
                )
                
                interpretable_data.append((
                    f'{behavior} {activity} {status} duration (days)',
                    duration_median
                ))
    
    # =================================================================
    # SIMULATION PARAMETERS
    # =================================================================
    sim_params = params['simulation']
    interpretable_data.extend([
        ('Population size (N)', sim_params['N']),
        ('Partnership burnin days', sim_params['partnership_burnin_days']),
        ('Transmission burnin days', sim_params['transmission_burnin_days']),
        ('Simulation days', sim_params['sim_days']),
        ('Random seed', sim_params['seed'])
    ])
    
    return interpretable_data

def create_summary_table(interpretable_data, output_file=None):
    """Create and display organized parameter summary"""
    
    # Create DataFrame
    df = pd.DataFrame(interpretable_data, columns=['Parameter', 'Value'])
    df['Value'] = pd.to_numeric(df['Value'], errors='coerce').round(4)
    
    # Save to CSV if requested
    if output_file:
        df.to_csv(output_file, index=False)
        print(f"💾 Saved complete summary to: {output_file}")
    
    # Print organized sections
    print("\n" + "="*70)
    print("INTERPRETABLE PARAMETER SUMMARY")
    print("="*70)
    
    # Define section patterns
    sections = [
        ("BEHAVIOR COMPOSITION", ["WSM proportion", "MSW proportion", "MSM proportion", "MSMW proportion"]),
        ("RISK ACTIVITY LEVELS", [p for p in df['Parameter'] if '% high activity' in p]),
        ("MEDIAN PARTNER ACQUISITION RATES - LOW ACTIVITY", [p for p in df['Parameter'] if 'low activity median' in p]),
        ("MEDIAN PARTNER ACQUISITION RATES - HIGH ACTIVITY", [p for p in df['Parameter'] if 'high activity median' in p]),
        ("PARTNERSHIP PARAMETERS", ["Casual rate", "p_MSMW_W", "Steady duration median (days)"]),
        ("STEADY PARTNERSHIP CAPACITY MIX", [p for p in df['Parameter'] if 'steady partners' in p]),
        ("TRANSMISSION PROBABILITIES - STEADY", [p for p in df['Parameter'] if 'steady transmission prob' in p]),
        ("TRANSMISSION PROBABILITIES - CASUAL", [p for p in df['Parameter'] if 'casual transmission prob' in p]),
        ("SYMPTOMATIC PROBABILITIES", [p for p in df['Parameter'] if 'p_sx' in p]),
        ("INFECTION DURATIONS - WSM", [p for p in df['Parameter'] if 'WSM' in p and 'duration (days)' in p]),
        ("INFECTION DURATIONS - MSW", [p for p in df['Parameter'] if 'MSW' in p and 'duration (days)' in p]),
        ("INFECTION DURATIONS - MSM", [p for p in df['Parameter'] if 'MSM ' in p and 'duration (days)' in p]),
        ("INFECTION DURATIONS - MSMW", [p for p in df['Parameter'] if 'MSMW' in p and 'duration (days)' in p]),
        ("SIMULATION PARAMETERS", [p for p in df['Parameter'] if any(x in p for x in ['Population', 'burnin', 'Simulation', 'seed'])])
    ]
    
    for section_name, param_patterns in sections:
        section_df = df[df['Parameter'].isin(param_patterns)]
        if len(section_df) > 0:
            print(f"\n{section_name}")
            print("-" * len(section_name))
            for _, row in section_df.iterrows():
                print(f"  {row['Parameter']}: {row['Value']}")
    
    return df

def main():
    parser = argparse.ArgumentParser(description='Interpret simulation parameters from JSON file')
    parser.add_argument('json_file', help='Path to JSON parameters file')
    parser.add_argument('--output', '-o', help='Output CSV file (optional)')
    parser.add_argument('--quiet', '-q', action='store_true', help='Only save to file, don\'t print')
    
    args = parser.parse_args()
    
    try:
        # Interpret parameters
        interpretable_data = interpret_parameters(args.json_file)
        
        # Create summary
        if not args.quiet:
            summary_df = create_summary_table(interpretable_data, args.output)
        else:
            summary_df = pd.DataFrame(interpretable_data, columns=['Parameter', 'Value'])
            summary_df['Value'] = pd.to_numeric(summary_df['Value'], errors='coerce').round(4)
            if args.output:
                summary_df.to_csv(args.output, index=False)
                print(f"💾 Saved summary to: {args.output}")
        
        print(f"\n✅ Successfully interpreted {len(interpretable_data)} parameters!")
        
    except Exception as e:
        print(f"❌ Error interpreting parameters: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
