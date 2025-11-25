#!/usr/bin/env python3
"""
Extract tip dates from infection CSV for LSD2 input
Creates properly formatted dates file from infection data
"""

import pandas as pd
import argparse
import os
import json
from Bio import SeqIO

def load_burnin_period(sim_dir):
    """Load burn-in period from simulation parameters"""
    params_file = os.path.join(sim_dir, "parameters_used.json")
    
    if os.path.exists(params_file):
        with open(params_file, 'r') as f:
            params = json.load(f)
        
        partnership_burnin = params['simulation']['partnership_burnin_days']
        transmission_burnin = params['simulation']['transmission_burnin_days']
        total_burnin = partnership_burnin + transmission_burnin
        
        print(f"Loaded burn-in from parameters: {total_burnin} days")
        return total_burnin
    else:
        # Fallback to burnin_period.txt if it exists
        burnin_file = os.path.join(os.path.dirname(sim_dir), "burnin_period.txt")
        if os.path.exists(burnin_file):
            with open(burnin_file, 'r') as f:
                burnin = int(f.read().strip())
            print(f"Loaded burn-in from burnin_period.txt: {burnin} days")
            return burnin
        else:
            raise FileNotFoundError(f"Could not find burn-in period in {params_file} or {burnin_file}")

def extract_tip_dates_from_csv(infections_csv, fasta_file, sim_dir, output_file):
    """Extract tip dates from infections CSV and create LSD2-formatted dates file"""
    
    print(f"Reading infections from: {infections_csv}")
    infections_df = pd.read_csv(infections_csv)
    
    print(f"Reading sequences from: {fasta_file}")
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    # Get actual tip names from FASTA file
    tip_names = [seq.id for seq in sequences]
    print(f"Found {len(tip_names)} sequences in FASTA")
    
    # Load burn-in period
    burnin_period = load_burnin_period(sim_dir)
    
    # Create mapping from nodename_samplingdate to full tip name
    tip_mapping = {}
    for tip_name in tip_names:
        # Extract nodename_samplingdate from nodename_samplingdate_episodenumber
        parts = tip_name.split('_')
        if len(parts) >= 2:
            key = f"{parts[0]}_{parts[1]}"  # nodename_samplingdate
            tip_mapping[key] = tip_name
    
    # Process infections and create dates
    tip_dates = []
    
    for _, infection in infections_df.iterrows():
        # Construct the matching key
        node_id = str(infection['infectee_node'])
        sampling_day = int(infection['day_of_sampling'])
        matching_key = f"{node_id}_{sampling_day}"
        
        if matching_key in tip_mapping:
            tip_name = tip_mapping[matching_key]
            
            # Adjust sampling day by subtracting burn-in
            adjusted_date = sampling_day - burnin_period
            
            tip_dates.append((tip_name, adjusted_date))
            
        else:
            print(f"Warning: No tip found for infection {matching_key}")
    
    # Sort by tip name for consistency
    tip_dates.sort(key=lambda x: x[0])
    
    print(f"Extracted dates for {len(tip_dates)} tips")
    
    # Write LSD2-formatted output
    with open(output_file, 'w') as f:
        # First line: number of tips
        f.write(f"{len(tip_dates)}\n")
        
        # Following lines: tip_name\tdate
        for tip_name, date in tip_dates:
            f.write(f"{tip_name}\t{date}\n")
    
    print(f"✅ Dates file created: {output_file}")
    print(f"Format: {len(tip_dates)} entries with dates adjusted by burn-in period ({burnin_period} days)")
    
    return len(tip_dates)

def main():
    parser = argparse.ArgumentParser(
        description='Extract tip dates from infection CSV for LSD2',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python extract_tip_dates.py infections.csv sequences.fasta /path/to/sim_dir dates.txt
  
  # If files are in organized directory structure:
  python extract_tip_dates.py detected_only/50pct/rep01/infections.csv \\
                               detected_only/50pct/rep01/sequences.fasta \\
                               /path/to/original/sim_dir \\
                               detected_only/50pct/rep01/dates.txt
        """
    )
    
    parser.add_argument('infections_csv', help='CSV file with infection data')
    parser.add_argument('fasta_file', help='FASTA file with sequences (for tip names)')
    parser.add_argument('sim_dir', help='Original simulation directory (for burn-in period)')
    parser.add_argument('output_file', help='Output dates file for LSD2')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.infections_csv):
        print(f"❌ Infections CSV not found: {args.infections_csv}")
        return 1
    
    if not os.path.exists(args.fasta_file):
        print(f"❌ FASTA file not found: {args.fasta_file}")
        return 1
    
    if not os.path.exists(args.sim_dir):
        print(f"❌ Simulation directory not found: {args.sim_dir}")
        return 1
    
    print("🕐 Starting tip date extraction...")
    
    try:
        # Create output directory if needed
        os.makedirs(os.path.dirname(args.output_file), exist_ok=True)
        
        # Extract dates
        n_dates = extract_tip_dates_from_csv(
            args.infections_csv,
            args.fasta_file, 
            args.sim_dir,
            args.output_file
        )
        
        print(f"\n🎉 Date extraction completed successfully!")
        print(f"Created {args.output_file} with {n_dates} tip dates")
        
        return 0
        
    except Exception as e:
        print(f"❌ Date extraction failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())
