#!/usr/bin/env python3
"""
Comprehensive Simulation Plotting Script
Creates plots and analyses from core simulation outputs with organized category structure
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import json
from collections import defaultdict, Counter

# ========== HELPER FUNCTIONS ==========

def _prep_daily_df(list_of_dicts):
    """Convert list of dicts to DataFrame indexed by day"""
    if not list_of_dicts:
        return pd.DataFrame()
    df = pd.DataFrame(list_of_dicts).fillna(0)
    if 'day' not in df.columns:
        df['day'] = np.arange(len(df))
    if 'phase' not in df.columns:
        df['phase'] = 'post_tx'
    df = df.sort_values('day').set_index('day')
    return df

def _shade_burnin(ax, partnership_burnin_days: int, tracking_start_day: int):
    """Shade the transmission burn-in, draw a vertical line at the end"""
    if partnership_burnin_days is not None and tracking_start_day is not None:
        if tracking_start_day > partnership_burnin_days:
            ax.axvspan(partnership_burnin_days, tracking_start_day - 1, alpha=0.10, label="Tx burn-in")
        ax.axvline(tracking_start_day, linestyle='--', linewidth=1.5, label="Burn-in over")

def _line_plot_multi(df, columns, title, ylabel, path, partnership_burnin_days, tracking_start_day, legend_cols=1):
    """Create multi-line time series plot"""
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 6))
    for col in columns:
        if col in df.columns:
            ax.plot(df.index, df[col], label=str(col), alpha=0.7, linewidth=1.8)
    ax.set_title(title)
    ax.set_xlabel("Day")
    ax.set_ylabel(ylabel)
    _shade_burnin(ax, partnership_burnin_days, tracking_start_day)
    ax.legend(ncol=legend_cols, fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)

# ========== DATA RECONSTRUCTION FROM CORE OUTPUTS ==========

def reconstruct_daily_epidemiology(transmission_df, nodes_dict, partnership_burnin_days, transmission_burnin_days, sim_days):
    """Reconstruct daily incidence and prevalence from transmission events"""
    
    total_days = partnership_burnin_days + transmission_burnin_days + sim_days
    
    # Create behavior and risk mappings
    node_behavior = {nid: data['behavior'] for nid, data in nodes_dict.items()}
    node_risk = {nid: 'hi' if data['hi_risk'] else 'lo' for nid, data in nodes_dict.items()}
    node_group = {nid: f"{data['behavior']}_{'hi' if data['hi_risk'] else 'lo'}" 
                  for nid, data in nodes_dict.items()}
    
    # Add mapped columns to transmission_df
    tx_df = transmission_df.copy()
    tx_df['behavior_infectee'] = tx_df['infectee_node'].map(node_behavior)
    tx_df['risk_infectee'] = tx_df['infectee_node'].map(node_risk)
    tx_df['group_infectee'] = tx_df['infectee_node'].map(node_group)
    tx_df['symptom_infectee'] = tx_df.get('infection_symptom_status', 'unknown')
    
    # Daily incidence by various groupings
    daily_incidence_behavior = tx_df.groupby(['day_of_transmission', 'behavior_infectee']).size().unstack(fill_value=0)
    daily_incidence_risk = tx_df.groupby(['day_of_transmission', 'risk_infectee']).size().unstack(fill_value=0)
    daily_incidence_group = tx_df.groupby(['day_of_transmission', 'group_infectee']).size().unstack(fill_value=0)
    daily_incidence_symptom = tx_df.groupby(['day_of_transmission', 'symptom_infectee']).size().unstack(fill_value=0)
    
    # Calculate prevalence (active infections each day)
    print("Computing daily prevalence...")
    
    # Filter to successful infections
    successful_tx = tx_df[tx_df['day_of_sampling'].notna()].copy()
    
    prevalence_data = []
    for day in range(total_days):
        if day % 500 == 0:
            print(f"  Day {day}/{total_days}")
            
        # Active infections on this day
        active = successful_tx[
            (successful_tx['day_of_transmission'] <= day) & 
            (successful_tx['day_of_sampling'] > day)
        ]
        
        # Count by groupings
        prev_behavior = active['behavior_infectee'].value_counts()
        prev_risk = active['risk_infectee'].value_counts()
        prev_group = active['group_infectee'].value_counts()
        prev_symptom = active['symptom_infectee'].value_counts()
        
        # Combine into row
        row = {'day': day}
        for behavior in ['WSM', 'MSW', 'MSM', 'MSMW']:
            row[behavior] = prev_behavior.get(behavior, 0)
        for risk in ['hi', 'lo']:
            row[f'risk_{risk}'] = prev_risk.get(risk, 0)
        for group in [f"{b}_{r}" for b in ['WSM', 'MSW', 'MSM', 'MSMW'] for r in ['hi', 'lo']]:
            row[f'group_{group}'] = prev_group.get(group, 0)
        for symptom in ['symptomatic', 'asymptomatic']:
            row[f'symptom_{symptom}'] = prev_symptom.get(symptom, 0)
            
        prevalence_data.append(row)
    
    prevalence_df = pd.DataFrame(prevalence_data).set_index('day')
    
    return {
        'daily_incidence_behavior': daily_incidence_behavior,
        'daily_incidence_risk': daily_incidence_risk,
        'daily_incidence_group': daily_incidence_group,
        'daily_incidence_symptom': daily_incidence_symptom,
        'daily_prevalence': prevalence_df
    }

def compute_concurrency_metrics(edge_df, nodes_dict):
    """Compute concurrency metrics from edge data"""
    
    if edge_df.empty:
        return {}
    
    # Same-day concurrency
    max_conc_by_node = {}
    conc_hist_all = Counter()
    conc_hist_gt1 = Counter()
    conc_type_counts = {'steady_steady': 0, 'casual_casual': 0, 'mixed': 0}
    same_day_ever_by_node = {}
    
    start_day = int(edge_df['day'].min())
    end_day = int((edge_df['day'] + edge_df['duration'] - 1).max())
    day_node_counts = defaultdict(lambda: defaultdict(lambda: [0, 0]))  # [steady, casual]
    
    for _, r in edge_df.iterrows():
        s = int(r['day'])
        e = int(r['day'] + r['duration'] - 1)
        cflag = int(r.get('casual', 0))
        u, v = int(r['node1']), int(r['node2'])
        for d in range(s, e + 1):
            if cflag == 1:
                day_node_counts[d][u][1] += 1
                day_node_counts[d][v][1] += 1
            else:
                day_node_counts[d][u][0] += 1
                day_node_counts[d][v][0] += 1
    
    for d in range(start_day, end_day + 1):
        nodes_today = day_node_counts.get(d, {})
        for n, (steady_ct, casual_ct) in nodes_today.items():
            total_ct = steady_ct + casual_ct
            conc_hist_all[total_ct] += 1
            if total_ct >= 2:
                conc_hist_gt1[total_ct] += 1
                if casual_ct == 0:
                    conc_type_counts['steady_steady'] += 1
                elif steady_ct == 0:
                    conc_type_counts['casual_casual'] += 1
                else:
                    conc_type_counts['mixed'] += 1
                same_day_ever_by_node[n] = True
            prev = max_conc_by_node.get(n, 0)
            if total_ct > prev:
                max_conc_by_node[n] = total_ct
    
    # Fill missing nodes
    for n in nodes_dict.keys():
        same_day_ever_by_node.setdefault(int(n), False)
    
    return {
        'max_conc_by_node': max_conc_by_node,
        'conc_hist_all': conc_hist_all,
        'conc_hist_gt1': conc_hist_gt1,
        'conc_type_counts': conc_type_counts,
        'same_day_ever_by_node': same_day_ever_by_node
    }

# ========== PLOTTING CATEGORIES ==========

def plot_basic_epidemiology(epi_data, nodes_dict, transmission_df, partnership_burnin_days, tracking_start_day, plots_dir, additional_outputs_dir):
    """Category 1: Basic epidemiology plots - always essential"""
    
    print("Creating basic epidemiology plots...")
    
    # Population denominators
    pop_by_behavior = defaultdict(int)
    pop_by_risk = defaultdict(int)
    pop_by_group = defaultdict(int)
    total_pop = len(nodes_dict)
    
    for _, d in nodes_dict.items():
        pop_by_behavior[d['behavior']] += 1
        risk_key = 'hi' if d['hi_risk'] else 'lo'
        pop_by_risk[risk_key] += 1
        group_key = f"{d['behavior']}_{risk_key}"
        pop_by_group[group_key] += 1
    
    # Save incidence CSVs
    epi_data['daily_incidence_behavior'].to_csv(
        os.path.join(additional_outputs_dir, "daily_incidence_by_behavior.csv"))
    epi_data['daily_incidence_group'].to_csv(
        os.path.join(additional_outputs_dir, "daily_incidence_by_group.csv"))
    epi_data['daily_prevalence'].to_csv(
        os.path.join(additional_outputs_dir, "current_infected_by_group_counts.csv"))
    
    # 1. NEW: True incidence over time (total)
    behavior_cols = [col for col in epi_data['daily_incidence_behavior'].columns]
    if behavior_cols:
        total_incidence = epi_data['daily_incidence_behavior'][behavior_cols].sum(axis=1) / total_pop
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(total_incidence.index, total_incidence.values, label="Total incidence", 
                linewidth=2)
        ax.set_title("True Incidence Over Time (Total)")
        ax.set_xlabel("Day")
        ax.set_ylabel("Incidence")
        ax.grid(True, alpha=0.3)
        _shade_burnin(ax, partnership_burnin_days, tracking_start_day)
        fig.tight_layout()
        fig.savefig(os.path.join(plots_dir, "true_incidence_over_time.png"), dpi=300)
        plt.close(fig)
    
    # 2. Incidence plots by behavior (both counts and rates)
    if not epi_data['daily_incidence_behavior'].empty:
        # Counts
        _line_plot_multi(
            epi_data['daily_incidence_behavior'], 
            epi_data['daily_incidence_behavior'].columns,
            "Daily New Infections by Behavior (Counts)", "New infections",
            os.path.join(plots_dir, "infections_by_behavior_counts.png"),
            partnership_burnin_days, tracking_start_day, legend_cols=2
        )
        
        # NEW: Incidence rates (per person)
        inc_behavior = epi_data['daily_incidence_behavior'].copy()
        for col in inc_behavior.columns:
            denom = pop_by_behavior.get(col, 0)
            inc_behavior[col] = inc_behavior[col] / denom if denom > 0 else 0.0
        _line_plot_multi(
            inc_behavior, inc_behavior.columns,
            "Daily Incidence by Behavior (per person)", "Incidence",
            os.path.join(plots_dir, "incidence_by_behavior.png"),
            partnership_burnin_days, tracking_start_day, legend_cols=2
        )
    
    # 3. NEW: Incidence by risk (both counts and rates)
    if not epi_data['daily_incidence_risk'].empty:
        _line_plot_multi(
            epi_data['daily_incidence_risk'], 
            epi_data['daily_incidence_risk'].columns,
            "Daily New Infections by Risk (Counts)", "New infections",
            os.path.join(plots_dir, "infections_by_risk_counts.png"),
            partnership_burnin_days, tracking_start_day
        )
        
        # NEW: Risk incidence rates
        inc_risk = epi_data['daily_incidence_risk'].copy()
        for col in inc_risk.columns:
            denom = pop_by_risk.get(col, 0)
            inc_risk[col] = inc_risk[col] / denom if denom > 0 else 0.0
        _line_plot_multi(
            inc_risk, inc_risk.columns,
            "Daily Incidence by Risk (per person)", "Incidence",
            os.path.join(plots_dir, "incidence_by_risk.png"),
            partnership_burnin_days, tracking_start_day
        )
    
    # 4. NEW: Incidence by symptom (both counts and rates)
    if not epi_data['daily_incidence_symptom'].empty:
        _line_plot_multi(
            epi_data['daily_incidence_symptom'],
            epi_data['daily_incidence_symptom'].columns,
            "Daily New Infections by Symptom Status (Counts)", "New infections",
            os.path.join(plots_dir, "infections_by_symptom_counts.png"),
            partnership_burnin_days, tracking_start_day
        )
        
        # NEW: Symptom incidence rates  
        inc_symptom = epi_data['daily_incidence_symptom'].copy()
        for col in inc_symptom.columns:
            inc_symptom[col] = inc_symptom[col] / total_pop
        _line_plot_multi(
            inc_symptom, inc_symptom.columns,
            "Daily Incidence by Symptom Status (per person)", "Incidence",
            os.path.join(plots_dir, "incidence_by_symptom.png"),
            partnership_burnin_days, tracking_start_day
        )
    
    # 5. Prevalence plots (rest of function unchanged...)
    prev_df = epi_data['daily_prevalence']
    
    # Overall prevalence
    behavior_cols = [col for col in prev_df.columns if col in ['WSM', 'MSW', 'MSM', 'MSMW']]
    if behavior_cols:
        overall_counts = prev_df[behavior_cols].sum(axis=1)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(overall_counts.index, overall_counts.values / total_pop)
        ax.set_title("Overall Prevalence Over Time")
        ax.set_xlabel("Day")
        ax.set_ylabel("Prevalence")
        ax.grid(True, alpha=0.3)
        _shade_burnin(ax, partnership_burnin_days, tracking_start_day)
        fig.tight_layout()
        fig.savefig(os.path.join(plots_dir, "prevalence_overall.png"), dpi=300)
        plt.close(fig)
        
        # Prevalence by behavior
        prev_beh = prev_df[behavior_cols].copy()
        for col in prev_beh.columns:
            denom = pop_by_behavior.get(col, 0)
            prev_beh[col] = prev_beh[col] / denom if denom > 0 else 0.0
        _line_plot_multi(
            prev_beh, prev_beh.columns,
            "Prevalence by Behavior Over Time", "Prevalence",
            os.path.join(plots_dir, "prevalence_by_behavior.png"),
            partnership_burnin_days, tracking_start_day, legend_cols=2
        )
    
    # 6. Prevalence by risk
    risk_cols = [col for col in prev_df.columns if col.startswith('risk_')]
    if risk_cols:
        prev_risk = prev_df[risk_cols].copy()
        for col in prev_risk.columns:
            risk_key = col.replace('risk_', '')
            denom = pop_by_risk.get(risk_key, 0)
            prev_risk[col] = prev_risk[col] / denom if denom > 0 else 0.0
        
        # Clean column names for legend
        prev_risk.columns = [col.replace('risk_', '') for col in prev_risk.columns]
        _line_plot_multi(
            prev_risk, prev_risk.columns,
            "Prevalence by Risk Over Time", "Prevalence",
            os.path.join(plots_dir, "prevalence_by_risk.png"),
            partnership_burnin_days, tracking_start_day
        )
    
    # 7. NEW: Prevalence by behavior × risk (your special plot)
    group_cols = [col for col in prev_df.columns if col.startswith('group_')]
    if group_cols:
        prev_group = prev_df[group_cols].copy()
        for col in prev_group.columns:
            group_key = col.replace('group_', '')
            denom = pop_by_group.get(group_key, 0)
            prev_group[col] = prev_group[col] / denom if denom > 0 else 0.0
        
        # Create the behavior × risk plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        behavior_colors = {'WSM': '#E74C3C', 'MSW': '#3498DB', 'MSM': '#2ECC71', 'MSMW': '#F39C12'}
        risk_styles = {'hi': '-', 'lo': '--'}
        
        for col in prev_group.columns:
            group_key = col.replace('group_', '')
            parts = group_key.split('_')
            if len(parts) == 2:
                behavior, risk = parts
                color = behavior_colors.get(behavior, 'gray')
                linestyle = risk_styles.get(risk, '-')
                
                ax.plot(prev_group.index, prev_group[col],
                       color=color, linestyle=linestyle, linewidth=2,
                       label=f"{behavior} ({risk} risk)")
        
        ax.set_title("Prevalence Over Time by Behavior Group × Risk Level")
        ax.set_xlabel("Day")
        ax.set_ylabel("Prevalence")
        ax.legend(ncol=2, fontsize=10)
        ax.grid(True, alpha=0.3)
        _shade_burnin(ax, partnership_burnin_days, tracking_start_day)
        
        fig.tight_layout()
        fig.savefig(os.path.join(plots_dir, "prevalence_by_behavior_risk.png"), dpi=300)
        plt.close(fig)

    # 8. NEW: Additional missing plots - ADD THIS SECTION:
    
    # Infections per person annualized (needs transmission_df)
    if not transmission_df.empty and 'infectee_node' in transmission_df.columns:
        min_day = int(transmission_df['day_of_transmission'].min())
        max_day = int(transmission_df['day_of_transmission'].max())
        total_days = max_day - min_day + 1
        years_obs = max(total_days / 365.0, 1e-9)
        
        trans_counts = transmission_df['infectee_node'].value_counts()
        trans_counts_per_year = trans_counts / years_obs
        
        plt.figure(figsize=(10, 6))
        plt.hist(trans_counts_per_year, bins=30, edgecolor='black', alpha=0.7, color='steelblue')
        plt.title("Infections Per Person (Annualized)")
        plt.xlabel("Infections per Year")
        plt.ylabel("Count")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "infections_per_person_annualized.png"), dpi=300)
        plt.close()
        
        # Infection count distribution (two-panel plot)
        repeat_threshold = 10
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Left panel: Overall distribution with log scale
        ax1.hist(trans_counts.values, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
        ax1.axvline(repeat_threshold, color='red', linestyle='--', 
                    label=f'Repeat threshold ({repeat_threshold})')
        ax1.set_xlabel('Number of Infections')
        ax1.set_ylabel('Number of People')
        ax1.set_title('Distribution of Infection Counts per Person')
        ax1.set_yscale('log')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Right panel: Detail view of repeat infectors
        repeat_infectors = trans_counts[trans_counts >= repeat_threshold]
        if len(repeat_infectors) > 0:
            ax2.hist(repeat_infectors.values, bins=min(20, len(repeat_infectors)), 
                     edgecolor='black', alpha=0.7, color='orange')
            ax2.set_xlabel('Number of Infections')
            ax2.set_ylabel('Number of People')
            ax2.set_title(f'Repeat Infectors ({repeat_threshold}+ infections)')
            ax2.grid(True, alpha=0.3)
        else:
            ax2.text(0.5, 0.5, f'No individuals with {repeat_threshold}+ infections', 
                     ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title(f'Repeat Infectors ({repeat_threshold}+ infections)')

        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "infection_count_distribution.png"), dpi=300)
        plt.close()

    # Infections by behavior × symptom
    if not transmission_df.empty and 'infection_symptom_status' in transmission_df.columns:
        tx_df = transmission_df.copy()
        tx_df['behavior_infectee'] = tx_df['infectee_node'].map({nid: data['behavior'] for nid, data in nodes_dict.items()})
        tx_df['behavior_symptom'] = tx_df['behavior_infectee'] + '_' + tx_df['infection_symptom_status']
        
        daily_beh_sym = tx_df.groupby(['day_of_transmission', 'behavior_symptom']).size().unstack(fill_value=0)
        
        if not daily_beh_sym.empty:
            _line_plot_multi(
                daily_beh_sym, daily_beh_sym.columns,
                "Daily New Infections by Behavior × Symptom", "New infections",
                os.path.join(plots_dir, "infections_by_behavior_symptom_counts.png"),
                partnership_burnin_days, tracking_start_day, legend_cols=2
            )
    
    # Infections by behavior × symptom  
    if 'infection_symptom_status' in transmission_df.columns:
            tx_df = transmission_df.copy()
            tx_df['behavior_infectee'] = tx_df['infectee_node'].map({nid: data['behavior'] for nid, data in nodes_dict.items()})
            tx_df['behavior_symptom'] = tx_df['behavior_infectee'] + '_' + tx_df['infection_symptom_status']
            
            daily_beh_sym = tx_df.groupby(['day_of_transmission', 'behavior_symptom']).size().unstack(fill_value=0)
            
            if not daily_beh_sym.empty:
                _line_plot_multi(
                    daily_beh_sym, daily_beh_sym.columns,
                    "Daily New Infections by Behavior × Symptom", "New infections",
                    os.path.join(plots_dir, "infections_by_behavior_symptom_counts.png"),
                    partnership_burnin_days, tracking_start_day, legend_cols=2
                )
    
    # 9. Prevalence by symptom (this one we can do now)
    symptom_cols = [col for col in prev_df.columns if col.startswith('symptom_')]
    if symptom_cols:
        prev_symptom = prev_df[symptom_cols].copy()
        for col in prev_symptom.columns:
            prev_symptom[col] = prev_symptom[col] / total_pop
        
        # Clean column names
        prev_symptom.columns = [col.replace('symptom_', '') for col in prev_symptom.columns]
        _line_plot_multi(
            prev_symptom, prev_symptom.columns,
            "Prevalence by Symptom Status Over Time", "Prevalence",
            os.path.join(plots_dir, "prevalence_by_symptom.png"),
            partnership_burnin_days, tracking_start_day
        )
   

def plot_partnerships(edge_df, nodes_dict, partnership_burnin_days, tracking_start_day, plots_dir, additional_outputs_dir):
    """Category 2: Partnership formation & dynamics"""
    
    print("Creating partnership plots...")
    
    if edge_df.empty:
        print("No edge data available for partnership plots")
        return
    
    # Active partnerships timeline
    start_day = int(edge_df['day'].min())
    end_day = int(edge_df['day'].max())

    daily_counts = defaultdict(lambda: defaultdict(int))

    for _, row in edge_df.iterrows():
        partnership_start = int(row['day'])
        duration = int(row['duration'])
        partnership_end = partnership_start + duration
        
        partnership_type = 'casual' if row.get('casual', 0) == 1 else 'steady'
        
        for day in range(partnership_start, partnership_end):
            daily_counts[day][partnership_type] += 1

    # Convert to DataFrame (KEEP this - CSV generation)
    all_days = range(start_day, end_day + 1)
    timeline_data = []
    for day in all_days:
        timeline_data.append({
            'day': day,
            'casual': daily_counts[day]['casual'],
            'steady': daily_counts[day]['steady'],
            'total': daily_counts[day]['casual'] + daily_counts[day]['steady']
        })

    timeline_df = pd.DataFrame(timeline_data)
    timeline_df.to_csv(os.path.join(additional_outputs_dir, "active_partnerships_timeline.csv"), index=False)

    # REPLACE the plotting section with this:

    # Plot 1: 5-panel faceted active partnerships by partnership type
    if 'edge_type' in edge_df.columns:
        # Create daily counts for each partnership type
        daily_partnership_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        
        for _, row in edge_df.iterrows():
            partnership_start = int(row['day'])
            duration = int(row['duration'])
            partnership_end = partnership_start + duration
            edge_type = row['edge_type']
            partnership_category = 'casual' if row.get('casual', 0) == 1 else 'steady'
            
            for day in range(partnership_start, partnership_end):
                daily_partnership_counts[day][edge_type][partnership_category] += 1
        
        # Get top 5 partnership types by frequency
        type_counts = edge_df['edge_type'].value_counts()
        top_5_types = type_counts.head(5).index.tolist()
        
        # Create 5-panel plot
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        axes = axes.flatten()
        
        for i, p_type in enumerate(top_5_types):
            if i >= 5:  # Only plot first 5
                break
                
            ax = axes[i]
            
            # Collect data for this partnership type
            days = []
            casual_counts = []
            steady_counts = []
            
            for day in range(start_day, end_day + 1):
                days.append(day)
                casual_counts.append(daily_partnership_counts[day][p_type]['casual'])
                steady_counts.append(daily_partnership_counts[day][p_type]['steady'])
            
            # Plot
            ax.plot(days, casual_counts, label='casual', color='#E74C3C', alpha=0.8)
            ax.plot(days, steady_counts, label='steady', color='#3498DB', alpha=0.8)
            ax.set_title(f'{p_type} Partnerships')
            ax.set_ylabel('Active Partnerships')
            ax.legend()
            ax.grid(True, alpha=0.3)
            _shade_burnin(ax, partnership_burnin_days, tracking_start_day)
        
        # Hide the 6th subplot
        axes[5].set_visible(False)
        
        # Common x-label
        fig.text(0.5, 0.02, 'Day', ha='center', fontsize=12)
        fig.suptitle('Active Partnerships Over Time by Type', fontsize=14)
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.96])
        plt.savefig(os.path.join(plots_dir, "active_partnerships_timeline.png"), dpi=300)
        plt.close()

    # Plot 2: Total partnerships timeline (separate plot)
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(timeline_df['day'], timeline_df['casual'], label='Casual', color='#E74C3C', linewidth=2)
    ax.plot(timeline_df['day'], timeline_df['steady'], label='Steady', color='#3498DB', linewidth=2)
    ax.plot(timeline_df['day'], timeline_df['total'], label='Total', color='black', linewidth=2, linestyle='--')

    ax.set_title('Total Active Partnerships Over Time')
    ax.set_xlabel('Day')
    ax.set_ylabel('Active Partnerships')
    ax.legend()
    ax.grid(True, alpha=0.3)
    _shade_burnin(ax, partnership_burnin_days, tracking_start_day)

    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "total_partnerships_timeline.png"), dpi=300)
    plt.close()
    
    # Partnership durations by type and risk
    if {'duration', 'casual', 'risk_level'}.issubset(edge_df.columns):
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        
        # Casual partnerships
        sub = edge_df[edge_df['casual'] == 1]
        axes[0].set_title("Casual Partnerships")
        for risk_label in ['low', 'high']:
            ds = sub[sub['risk_level'] == risk_label]['duration']
            if not ds.empty:
                axes[0].hist(ds, bins=30, edgecolor='black', alpha=0.7, label=f"{risk_label} risk")
        axes[0].set_xlabel("Duration (days)")
        axes[0].set_ylabel("Count")
        axes[0].legend()
        
        # Steady partnerships
        sub = edge_df[edge_df['casual'] == 0]
        axes[1].set_title("Steady Partnerships")
        for risk_label in ['low', 'high']:
            ds = sub[sub['risk_level'] == risk_label]['duration']
            if not ds.empty:
                axes[1].hist(ds, bins=30, edgecolor='black', alpha=0.7, label=f"{risk_label} risk")
        axes[1].set_xlabel("Duration (days)")
        axes[1].set_ylabel("Count")
        axes[1].legend()
        
        fig.suptitle("Partnership Durations by Type and Risk")
        fig.tight_layout()
        plt.savefig(os.path.join(plots_dir, "partnership_duration_facets.png"), dpi=300)
        plt.close()
        # Calculate years for annualization (post-partnership-burnin period)
    post_partnership_days = (end_day - partnership_burnin_days + 1) if end_day > partnership_burnin_days else 1
    years = max(post_partnership_days / 365.0, 1e-9)
    
    # 1. Edge type frequencies (annualized)
    if 'edge_type' in edge_df.columns:
        # Casual edge types
        casual_edges = edge_df[edge_df['casual'] == 1]
        if not casual_edges.empty:
            casual_counts = casual_edges['edge_type'].value_counts().sort_values(ascending=False)
            annualized_casual = casual_counts / years
            
            plt.figure(figsize=(12, 8))
            annualized_casual.plot(kind='bar', edgecolor='black', alpha=0.7)
            plt.title("Edge Type Frequencies (Casual) — Annualized")
            plt.xlabel("Edge Type")
            plt.ylabel("Partnerships per Year")
            plt.xticks(rotation=45, ha='right')
            plt.grid(True, alpha=0.3, axis='y')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "edge_type_freq_casual_annualized.png"), dpi=300)
            plt.close()
        
        # Steady edge types
        steady_edges = edge_df[edge_df['casual'] == 0]
        if not steady_edges.empty:
            steady_counts = steady_edges['edge_type'].value_counts().sort_values(ascending=False)
            annualized_steady = steady_counts / years
            
            plt.figure(figsize=(12, 8))
            annualized_steady.plot(kind='bar', edgecolor='black', alpha=0.7)
            plt.title("Edge Type Frequencies (Steady) — Annualized")
            plt.xlabel("Edge Type")
            plt.ylabel("Partnerships per Year")
            plt.xticks(rotation=45, ha='right')
            plt.grid(True, alpha=0.3, axis='y')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "edge_type_freq_steady_annualized.png"), dpi=300)
            plt.close()
    
    # 2. Partner count distributions by behavior and risk (4-panel plots)
    if {'node1', 'node2'}.issubset(edge_df.columns):
        # Count partnerships per node
        node_partner_counts = defaultdict(int)
        for _, row in edge_df.iterrows():
            node_partner_counts[row['node1']] += 1
            node_partner_counts[row['node2']] += 1
        
        # Annualize partner counts
        partner_counts_annualized = {node: count / years for node, count in node_partner_counts.items()}
        
        # Organize by behavior and risk
        partner_data = {'WSM': {'hi': [], 'lo': []}, 'MSW': {'hi': [], 'lo': []}, 
                    'MSM': {'hi': [], 'lo': []}, 'MSMW': {'hi': [], 'lo': []}}
        
        for node_id, count in partner_counts_annualized.items():
            if node_id in nodes_dict:
                behavior = nodes_dict[node_id]['behavior']
                risk = 'hi' if nodes_dict[node_id]['hi_risk'] else 'lo'
                partner_data[behavior][risk].append(count)
        
        colors = {'lo': 'skyblue', 'hi': 'orange'}
        behaviors = ['WSM', 'MSW', 'MSM', 'MSMW']
        
        # Plot all (4-panel by behavior)
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        axes = axes.flatten()
        
        for i, behavior in enumerate(behaviors):
            ax = axes[i]
            for risk in ['lo', 'hi']:
                data = partner_data[behavior][risk]
                if data:
                    ax.hist(data, bins=30, alpha=0.7, label=f'hi_risk={1 if risk=="hi" else 0}',
                        color=colors[risk], edgecolor='black')
            ax.set_title(behavior)
            ax.set_xlabel('Partners per Year')
            ax.set_ylabel('Count')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        fig.suptitle('Partner Counts by Behavior and Risk (all)')
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "network_partner_counts_all.png"), dpi=300)
        plt.close()
        
        # Plot high risk only (4-panel by behavior)
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        axes = axes.flatten()
        
        for i, behavior in enumerate(behaviors):
            ax = axes[i]
            data_hi = partner_data[behavior]['hi']
            if data_hi:
                ax.hist(data_hi, bins=30, alpha=0.7, color=colors['hi'], edgecolor='black')
            ax.set_title(f'{behavior} (High Risk)')
            ax.set_xlabel('Partners per Year')
            ax.set_ylabel('Count')
            ax.grid(True, alpha=0.3)
        
        fig.suptitle('Partner Counts by Behavior (High Risk Only)')
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "network_partner_counts_high.png"), dpi=300)
        plt.close()
        
        # Plot low risk only (4-panel by behavior)
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        axes = axes.flatten()
        
        for i, behavior in enumerate(behaviors):
            ax = axes[i]
            data_lo = partner_data[behavior]['lo']
            if data_lo:
                ax.hist(data_lo, bins=30, alpha=0.7, color=colors['lo'], edgecolor='black')
            ax.set_title(f'{behavior} (Low Risk)')
            ax.set_xlabel('Partners per Year')
            ax.set_ylabel('Count')
            ax.grid(True, alpha=0.3)
        
        fig.suptitle('Partner Counts by Behavior (Low Risk Only)')
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "network_partner_counts_low.png"), dpi=300)
        plt.close()
    
    # 3. Partnerships per node by type (casual vs steady)
    if 'edge_type' in edge_df.columns:
        # Calculate population by behavior for normalization
        pop_by_behavior = defaultdict(int)
        for _, d in nodes_dict.items():
            pop_by_behavior[d['behavior']] += 1
        
        # Casual partnerships per node
        casual_counts = edge_df[edge_df['casual'] == 1]['edge_type'].value_counts()
        casual_per_node = {}
        
        for edge_type, count in casual_counts.items():
            behaviors = edge_type.split('-')
            if len(behaviors) == 2:
                b1, b2 = behaviors
                total_nodes = pop_by_behavior.get(b1, 0) + pop_by_behavior.get(b2, 0)
                if total_nodes > 0:
                    casual_per_node[edge_type] = (count / years) / total_nodes
        
        if casual_per_node:
            plt.figure(figsize=(12, 7))
            edge_types = list(casual_per_node.keys())
            values = list(casual_per_node.values())
            plt.bar(edge_types, values, edgecolor='black', alpha=0.7)
            plt.title("Annualized Casual Partnerships per Node by Edge Type")
            plt.xlabel("Edge Type")
            plt.ylabel("Partnerships per Person per Year")
            plt.xticks(rotation=45, ha='right')
            plt.grid(True, alpha=0.3, axis='y')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "partnerships_per_node_casual.png"), dpi=300)
            plt.close()
        
        # Steady partnerships per node
        steady_counts = edge_df[edge_df['casual'] == 0]['edge_type'].value_counts()
        steady_per_node = {}
        
        for edge_type, count in steady_counts.items():
            behaviors = edge_type.split('-')
            if len(behaviors) == 2:
                b1, b2 = behaviors
                total_nodes = pop_by_behavior.get(b1, 0) + pop_by_behavior.get(b2, 0)
                if total_nodes > 0:
                    steady_per_node[edge_type] = (count / years) / total_nodes
        
        if steady_per_node:
            plt.figure(figsize=(12, 7))
            edge_types = list(steady_per_node.keys())
            values = list(steady_per_node.values())
            plt.bar(edge_types, values, edgecolor='black', alpha=0.7)
            plt.title("Annualized Steady Partnerships per Node by Edge Type")
            plt.xlabel("Edge Type")
            plt.ylabel("Partnerships per Person per Year")
            plt.xticks(rotation=45, ha='right')
            plt.grid(True, alpha=0.3, axis='y')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "partnerships_per_node_steady.png"), dpi=300)
            plt.close()
    
    # 4. Unpartnered duration analysis
    if {'node1', 'node2'}.issubset(edge_df.columns):
        # Calculate partnered days for each node
        partnered_days_by_node = defaultdict(set)
        
        for _, row in edge_df.iterrows():
            start = int(row['day'])
            duration = int(row['duration'])
            end = start + duration - 1
            
            node1 = int(row['node1'])
            node2 = int(row['node2'])
            
            for day in range(max(start, start_day), min(end + 1, end_day + 1)):
                partnered_days_by_node[node1].add(day)
                partnered_days_by_node[node2].add(day)
        
        # Calculate unpartnered durations
        unpartnered_durations = []
        gaps_by_group = defaultdict(list)
        period_days = end_day - start_day + 1
        
        for node_id in nodes_dict.keys():
            node_id = int(node_id)
            partnered_days = partnered_days_by_node.get(node_id, set())
            node_info = nodes_dict[node_id]
            group = f"{node_info['behavior']}_{'hi' if node_info['hi_risk'] else 'lo'}"
            
            if not partnered_days:
                # Never partnered
                gap = period_days
                unpartnered_durations.append(gap)
                gaps_by_group[group].append(gap)
                continue
            
            # Calculate gaps between partnerships
            all_possible_days = set(range(start_day, end_day + 1))
            unpartnered_days = all_possible_days - partnered_days
            
            if unpartnered_days:
                # Calculate consecutive gaps
                unpartnered_sorted = sorted(unpartnered_days)
                gap_length = 1
                
                for i in range(1, len(unpartnered_sorted)):
                    if unpartnered_sorted[i] == unpartnered_sorted[i-1] + 1:
                        gap_length += 1
                    else:
                        unpartnered_durations.append(gap_length)
                        gaps_by_group[group].append(gap_length)
                        gap_length = 1
                
                # Don't forget last gap
                unpartnered_durations.append(gap_length)
                gaps_by_group[group].append(gap_length)
        
        # Plot overall unpartnered durations
        if unpartnered_durations:
            plt.figure(figsize=(10, 6))
            plt.hist(unpartnered_durations, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
            plt.title("Histogram of Unpartnered Durations (days)")
            plt.xlabel("Days Unpartnered")
            plt.ylabel("Count of Gaps")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "network_unpartnered_duration_hist.png"), dpi=300)
            plt.close()
            
            # Stacked histogram by group
            bins = [(0, 0), (1, 7), (8, 30), (31, 90), (91, 365), (366, float('inf'))]
            bin_labels = ["0", "1-7", "8-30", "31-90", "91-365", "366+"]
            
            all_groups = set(gaps_by_group.keys())
            if not all_groups:
                behaviors = ['WSM', 'MSW', 'MSM', 'MSMW']
                risks = ['hi', 'lo']
                all_groups = {f"{b}_{r}" for b in behaviors for r in risks}
            
            groups = sorted(all_groups)
            data = {g: [0] * len(bins) for g in groups}
            
            for group, gaps in gaps_by_group.items():
                for gap in gaps:
                    for i, (lo, hi) in enumerate(bins):
                        if lo <= gap <= hi:
                            data[group][i] += 1
                            break
            
            gap_stack_df = pd.DataFrame(data, index=bin_labels)
            
            ax = gap_stack_df.plot(kind='bar', stacked=True, figsize=(12, 6), edgecolor='black')
            ax.set_title("Unpartnered Duration Distribution by Behavior×Risk")
            ax.set_xlabel("Gap duration bin (days)")
            ax.set_ylabel("Count of Gaps")
            ax.set_xticklabels(bin_labels, rotation=0)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "network_unpartnered_duration_stacked_by_group.png"), dpi=300)
            plt.close()
            
            # Save CSV
            gap_stack_df.to_csv(os.path.join(additional_outputs_dir, "network_gap_durations_by_group_bins.csv"))

def plot_concurrency(edge_df, nodes_dict, plots_dir, additional_outputs_dir):
    """Category 3: Concurrency analysis"""
    
    print("Creating concurrency plots...")
    
    concurrency_metrics = compute_concurrency_metrics(edge_df, nodes_dict)
    
    if not concurrency_metrics:
        print("No concurrency data available")
        return
    
    max_conc_by_node = concurrency_metrics['max_conc_by_node']
    conc_type_counts = concurrency_metrics['conc_type_counts']
    same_day_ever_by_node = concurrency_metrics['same_day_ever_by_node']
    
    # Save concurrency summary
    never = sum(1 for v in max_conc_by_node.values() if v == 0)
    mono_only = sum(1 for v in max_conc_by_node.values() if v == 1)
    concurrent = sum(1 for v in max_conc_by_node.values() if v >= 2)
    
    conc_summary = pd.DataFrame({
        'category': ['never_partnered', 'monogamous_only', 'concurrent_ever'],
        'count': [never, mono_only, concurrent]
    })
    conc_summary.to_csv(os.path.join(additional_outputs_dir, "concurrency_summary.csv"), index=False)
    
    # Plot max concurrency per node
    if len(max_conc_by_node) > 0:
        plt.figure(figsize=(10, 6))
        plt.hist(list(max_conc_by_node.values()),
                bins=range(0, max(max_conc_by_node.values()) + 2), 
                align='left', edgecolor='black', alpha=0.7)
        plt.title("Distribution of Max Concurrent Partnerships per Node")
        plt.xlabel("Max concurrent partners")
        plt.ylabel("Count of nodes")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "concurrency_max_per_node_hist.png"), dpi=300)
        plt.close()
        
        # Concurrency > 1 only
        vals_gt1 = [v for v in max_conc_by_node.values() if v > 1]
        if len(vals_gt1) > 0:
            plt.figure(figsize=(10, 6))
            plt.hist(vals_gt1, bins=range(2, max(vals_gt1) + 2), 
                    align='left', edgecolor='black', alpha=0.7)
            plt.title("Distribution of Max Concurrency (Nodes with >1 at any time)")
            plt.xlabel("Max concurrent partners")
            plt.ylabel("Count of nodes")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "concurrency_max_gt1_hist.png"), dpi=300)
            plt.close()
    
    # Concurrency type composition
    if conc_type_counts:
        plt.figure(figsize=(10, 6))
        labels = list(conc_type_counts.keys())
        values = [conc_type_counts[k] for k in labels]
        colors = ['#3498DB', '#E74C3C', '#F39C12']
        
        plt.bar(labels, values, edgecolor='black', alpha=0.7, color=colors[:len(labels)])
        plt.title("Composition of Concurrency by Partnership Types\n(node-days with ≥2 partners)")
        plt.ylabel("Count of node-days")
        plt.xticks(rotation=45)
        plt.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "concurrency_type_composition.png"), dpi=300)
        plt.close()

    # Rolling window concurrency analysis (3, 4, 5 day windows)
    print("Computing rolling window concurrency...")
    
    if not edge_df.empty:
        # Get simulation time parameters
        start_day = int(edge_df['day'].min())
        end_day = int(edge_df['day'].max())
        total_days = end_day - start_day + 1
        
        node_ids = sorted(nodes_dict.keys())
        id_to_idx = {nid: i for i, nid in enumerate(node_ids)}
        n = len(node_ids)
        
        # Build daily partner sets for each node
        daily_partners = [defaultdict(set) for _ in range(total_days)]
        
        for _, row in edge_df.iterrows():
            partnership_start = int(row['day']) - start_day  # Convert to 0-based indexing
            duration = int(row['duration'])
            partnership_end = partnership_start + duration - 1
            
            node1 = int(row['node1'])
            node2 = int(row['node2'])
            
            for day_idx in range(max(0, partnership_start), min(total_days, partnership_end + 1)):
                daily_partners[day_idx][node1].add(node2)
                daily_partners[day_idx][node2].add(node1)
        
        # Compute rolling concurrency for different windows
        for W in [3, 4, 5]:
            print(f"  Processing {W}-day rolling window...")
            
            # Track rolling partner counts
            node_rolling_counts = np.zeros((total_days, n), dtype=int)
            
            for day_idx in range(total_days):
                # For each node, count distinct partners in past W days
                for node_id in node_ids:
                    partners_in_window = set()
                    
                    # Look back W days (including current day)
                    for lookback in range(max(0, day_idx - W + 1), day_idx + 1):
                        partners_in_window.update(daily_partners[lookback].get(node_id, set()))
                    
                    node_idx = id_to_idx[node_id]
                    node_rolling_counts[day_idx, node_idx] = len(partners_in_window)
            
            # Calculate metrics
            # Average rolling concurrency per node over all days
            avg_rolling_conc = node_rolling_counts.mean(axis=0)
            
            # Identify nodes that ever had rolling concurrency >= 2
            ever_concurrent_mask = (node_rolling_counts >= 2).any(axis=0)
            
            # Average rolling concurrency on days when >= 2 (for nodes that ever had >= 2)
            avg_on_concurrent_days = np.full(n, 0.0)
            for i in range(n):
                concurrent_days = node_rolling_counts[:, i] >= 2
                if concurrent_days.any():
                    avg_on_concurrent_days[i] = node_rolling_counts[concurrent_days, i].mean()
            
            # Save CSVs
            rolling_metrics = pd.DataFrame({
                'node_id': node_ids,
                f'avg_rolling_concurrency_{W}': avg_rolling_conc,
                f'avg_on_concurrent_days_{W}': avg_on_concurrent_days,
                f'ever_rolling_concurrent_{W}': ever_concurrent_mask
            })
            rolling_metrics.to_csv(
                os.path.join(additional_outputs_dir, f"rolling_concurrency_node_metrics_{W}.csv"), index=False)
            
            # Plot 1: Average rolling concurrency (all nodes)
            plt.figure(figsize=(10, 6))
            plt.hist(avg_rolling_conc, bins=30, edgecolor='black', alpha=0.7)
            plt.title(f"Rolling Window Concurrency (W={W}) — Average per Node")
            plt.xlabel(f"Average # distinct partners in past {W} days")
            plt.ylabel("Count of nodes")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"rolling_concurrency_{W}_all.png"), dpi=300)
            plt.close()
            
            # Plot 2: Average on concurrent days (nodes that ever had >= 2)
            concurrent_nodes_avg = avg_on_concurrent_days[ever_concurrent_mask]
            
            plt.figure(figsize=(10, 6))
            if len(concurrent_nodes_avg) > 0:
                plt.hist(concurrent_nodes_avg, bins=30, edgecolor='black', alpha=0.7)
                plt.title(f"Rolling Window Concurrency (W={W}) — Degree on Concurrent Days")
                plt.xlabel(f"Average window-degree on concurrent days")
                plt.ylabel("Count of nodes")
            else:
                plt.text(0.5, 0.5, f"No nodes reached ≥2 partners in {W}-day window",
                         ha="center", va="center", transform=plt.gca().transAxes)
                plt.xlim(0, 1)
                plt.title(f"Rolling Window Concurrency (W={W}) — Degree on Concurrent Days")
                plt.xlabel(f"Average window-degree on concurrent days")
                plt.ylabel("Count of nodes")
            
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"rolling_concurrency_{W}_if_concurrent.png"), dpi=300)
            plt.close()
            
            # Save additional CSVs for detailed analysis
            # Daily overall means
            daily_means = node_rolling_counts.mean(axis=1)
            daily_overall_df = pd.DataFrame({
                'day_since_start': np.arange(total_days),
                f'mean_rolling_concurrency_{W}': daily_means
            })
            daily_overall_df.to_csv(
                os.path.join(additional_outputs_dir, f"rolling_concurrency_daily_overall_{W}.csv"), index=False)
            
            # Events (days when nodes had >= 2 partners in window)
            events_data = []
            for day_idx in range(total_days):
                for node_idx, node_id in enumerate(node_ids):
                    degree = node_rolling_counts[day_idx, node_idx]
                    if degree >= 2:
                        events_data.append({
                            'day_since_start': day_idx,
                            'node_id': node_id,
                            'window': W,
                            'window_degree': degree
                        })
            
            if events_data:
                events_df = pd.DataFrame(events_data)
                events_df.to_csv(
                    os.path.join(additional_outputs_dir, f"rolling_concurrency_events_{W}.csv"), index=False)
            
            # Daily histogram
            hist_data = []
            for day_idx in range(total_days):
                degree_counts = Counter(node_rolling_counts[day_idx, :])
                for degree, count in degree_counts.items():
                    hist_data.append({
                        'day_since_start': day_idx,
                        'degree': degree,
                        'count': count,
                        'window': W
                    })
            
            if hist_data:
                hist_df = pd.DataFrame(hist_data)
                hist_df.to_csv(
                    os.path.join(additional_outputs_dir, f"rolling_concurrency_daily_hist_{W}.csv"), index=False)
    
    print("Rolling concurrency analysis complete!")

def plot_transmission(transmission_df, edge_df, partnership_burnin_days, tracking_start_day, plots_dir, additional_outputs_dir):
    """Category 4: Transmission patterns"""
    
    print("Creating transmission analysis...")
    
    if transmission_df.empty:
        return
    
    # Bridging transmission analysis
    if {'behavior_infector', 'behavior_infectee'}.issubset(transmission_df.columns):
        total_events = len(transmission_df)
        
        def prop(count): 
            return (count / total_events) if total_events > 0 else 0.0
        
        cats = {
            'WSM_to_MSMW': ((transmission_df['behavior_infector'] == 'WSM') & 
                          (transmission_df['behavior_infectee'] == 'MSMW')).sum(),
            'MSMW_to_WSM': ((transmission_df['behavior_infector'] == 'MSMW') & 
                          (transmission_df['behavior_infectee'] == 'WSM')).sum(),
            'MSMW_to_MSMW': ((transmission_df['behavior_infector'] == 'MSMW') & 
                           (transmission_df['behavior_infectee'] == 'MSMW')).sum(),
            'MSMW_to_MSM': ((transmission_df['behavior_infector'] == 'MSMW') & 
                          (transmission_df['behavior_infectee'] == 'MSM')).sum(),
            'MSM_to_MSMW': ((transmission_df['behavior_infector'] == 'MSM') & 
                          (transmission_df['behavior_infectee'] == 'MSMW')).sum(),
        }
        
        other_mask = ((transmission_df['behavior_infector'] == 'MSMW') | 
                     (transmission_df['behavior_infectee'] == 'MSMW'))
        known_mask = (
            ((transmission_df['behavior_infector'] == 'WSM') & (transmission_df['behavior_infectee'] == 'MSMW')) |
            ((transmission_df['behavior_infector'] == 'MSMW') & (transmission_df['behavior_infectee'] == 'WSM')) |
            ((transmission_df['behavior_infector'] == 'MSMW') & (transmission_df['behavior_infectee'] == 'MSMW')) |
            ((transmission_df['behavior_infector'] == 'MSMW') & (transmission_df['behavior_infectee'] == 'MSM')) |
            ((transmission_df['behavior_infector'] == 'MSM') & (transmission_df['behavior_infectee'] == 'MSMW'))
        )
        cats['other_MSMW_involving'] = (other_mask & ~known_mask).sum()
        cats['any_MSMW_involved'] = other_mask.sum()
        
        bridging_data = [{'category': k, 'count': int(v), 'proportion_of_all_transmissions': prop(v)} 
                        for k, v in cats.items()]
        bridging_df = pd.DataFrame(bridging_data)
        bridging_df.to_csv(os.path.join(additional_outputs_dir, "bridging_transmission_proportions.csv"), index=False)
    
    # Transmission proportions by partnership type
    if not edge_df.empty and 'day_of_transmission' in transmission_df.columns and 'partnership_category' in transmission_df.columns:
        # Calculate active partnerships by day and type
        active_counts = defaultdict(lambda: defaultdict(int))
        
        for _, row in edge_df.iterrows():
            start = int(row['day'])
            duration = int(row['duration'])
            p_type = 'casual' if row.get('casual', 0) == 1 else 'steady'
            
            for day in range(start, start + duration):
                active_counts[day][p_type] += 1
        
        # Calculate transmissions by day and type
        daily_trans = transmission_df.groupby(['day_of_transmission', 'partnership_category']).size().unstack(fill_value=0)
        
        # Calculate rate
        rate_data = []
        for day in daily_trans.index:
            row_data = {'day': day}
            for p_type in ['casual', 'steady']:
                trans_count = daily_trans.loc[day, p_type] if p_type in daily_trans.columns else 0
                active_count = active_counts[day][p_type]
                rate = trans_count / active_count if active_count > 0 else 0
                row_data[p_type] = rate
            rate_data.append(row_data)
        
        rate_df = pd.DataFrame(rate_data).set_index('day')
        
        _line_plot_multi(
            rate_df, rate_df.columns,
            "Transmissions per Active Partnership (Daily)", "Events per active partnership",
            os.path.join(plots_dir, "transmissions_per_active_partnership_by_type_daily.png"),
            partnership_burnin_days, tracking_start_day
        )
        
        # Daily transmission proportions
        if 'day_of_transmission' in transmission_df.columns:
            daily_type = transmission_df.groupby(['day_of_transmission', 'partnership_category']).size().unstack(fill_value=0)
            denom = daily_type.sum(axis=1).replace(0, np.nan)
            prop_df = daily_type.div(denom, axis=0).fillna(0)
            
            fig, ax = plt.subplots(figsize=(10, 6))
            for col in ['casual', 'steady']:
                if col in prop_df.columns:
                    ax.plot(prop_df.index, prop_df[col], label=col, linewidth=2)
            ax.set_title("Proportion of Transmissions by Partnership Type (Daily)")
            ax.set_xlabel("Day")
            ax.set_ylabel("Proportion")
            ax.legend()
            ax.grid(True, alpha=0.3)
            _shade_burnin(ax, partnership_burnin_days, tracking_start_day)
            fig.tight_layout()
            fig.savefig(os.path.join(plots_dir, "proportion_transmissions_by_type_daily.png"), dpi=300)
            plt.close()
            # 1. Proportion transmissions by detailed pair type (using partnership_type)
    if 'partnership_type' in transmission_df.columns:
        trans_counts = transmission_df['partnership_type'].value_counts()
        trans_prop = (trans_counts / trans_counts.sum()).sort_values(ascending=False)
        
        plt.figure(figsize=(14, 8))
        trans_prop.plot(kind='bar', color='steelblue', edgecolor='black', alpha=0.7)
        plt.title("Proportion of Transmission Events by Partnership Type")
        plt.xlabel("Partnership Type")
        plt.ylabel("Proportion")
        plt.xticks(rotation=45, ha='right')
        plt.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "proportion_transmissions_by_pairtype_overall.png"), dpi=300)
        plt.close()
    
    # 2. Normalized by partnership formation rates
    if not edge_df.empty and 'partnership_type' in transmission_df.columns:
        # Create partnership_type column in edge_df to match transmission_df
        edge_df_copy = edge_df.copy()
        
        # Map edge_df to partnership_type format
        edge_df_copy['partnership_type'] = edge_df_copy.apply(lambda row: 
            f"{row.get('edge_type', 'Unknown')}, {'casual' if row.get('casual', 0) == 1 else 'steady'}", axis=1)
        
        # Calculate partnership formation proportions
        edge_counts = edge_df_copy['partnership_type'].value_counts(normalize=True)
        
        # Calculate transmission proportions  
        trans_counts = transmission_df['partnership_type'].value_counts(normalize=True)
        
        # Find common partnership types
        common_types = set(edge_counts.index).intersection(set(trans_counts.index))
        
        if common_types:
            ratios = []
            labels = []
            
            for p_type in sorted(common_types):
                edge_prop = edge_counts[p_type]
                trans_prop = trans_counts[p_type]
                ratio = trans_prop / edge_prop if edge_prop > 0 else 0
                ratios.append(ratio)
                labels.append(p_type)
            
            plt.figure(figsize=(14, 8))
            plt.bar(labels, ratios, color='orange', edgecolor='black', alpha=0.7)
            plt.title("Transmissions / Partnership Share by Detailed Pair Type")
            plt.xlabel("Partnership Type")
            plt.ylabel("Ratio (1 = proportional to partnerships)")
            plt.xticks(rotation=45, ha='right')
            plt.grid(True, alpha=0.3, axis='y')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "proportion_transmissions_normalized_by_pairtype_overall.png"), dpi=300)
            plt.close()

def plot_infection_recovery(transmission_df, plots_dir, additional_outputs_dir):
    """Category 5: Infection recovery & duration"""
    
    print("Creating infection recovery analysis...")
    
    if transmission_df.empty:
        return
    
    # Calculate durations
    df = transmission_df.copy()
    df = df[df['day_of_sampling'].notna()].copy()
    df['duration'] = df['day_of_sampling'] - df['day_of_transmission']
    df = df[df['duration'] > 0]
    
    # Duration by symptom status
    if 'infection_symptom_status' in df.columns:
        df_sym = df[df['infection_symptom_status'].isin(['symptomatic', 'asymptomatic'])]
        
        if not df_sym.empty:
            # Summary stats
            summary_sym = df_sym.groupby('infection_symptom_status')['duration'].describe()
            summary_sym.to_csv(os.path.join(additional_outputs_dir, "infection_duration_by_symptom_summary.csv"))
            
            # Plot
            plt.figure(figsize=(10, 6))
            bins = np.linspace(0, df_sym['duration'].max(), 50)
            
            for status, color in [('symptomatic', 'red'), ('asymptomatic', 'blue')]:
                sub = df_sym[df_sym['infection_symptom_status'] == status]['duration']
                if not sub.empty:
                    plt.hist(sub, bins=bins, alpha=0.6, label=status, 
                           edgecolor='black', color=color)
            
            plt.title("Distribution of Infection Durations by Symptom Status")
            plt.xlabel("Duration (days)")
            plt.ylabel("Count")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "infection_duration_by_symptom.png"), dpi=300)
            plt.close()
    
    # Duration by recovery method
    if 'duration_source' in df.columns:
        df_recovery = df[df['duration_source'].isin(['symptomatic', 'asymptomatic', 'natural_clearance'])]
        
        if not df_recovery.empty:
            # Summary stats
            summary_recovery = df_recovery.groupby('duration_source')['duration'].describe()
            summary_recovery.to_csv(os.path.join(additional_outputs_dir, "infection_duration_by_recovery_method_summary.csv"))
            
            # Plot
            plt.figure(figsize=(10, 6))
            bins = np.linspace(0, df_recovery['duration'].max(), 50)
            
            colors = {'symptomatic': 'red', 'asymptomatic': 'blue', 'natural_clearance': 'green'}
            
            for method in ['symptomatic', 'asymptomatic', 'natural_clearance']:
                sub = df_recovery[df_recovery['duration_source'] == method]['duration']
                if not sub.empty:
                    plt.hist(sub, bins=bins, alpha=0.6, label=method.replace('_', ' ').title(),
                           edgecolor='black', color=colors[method])
            
            plt.title("Distribution of Infection Durations by Recovery Method")
            plt.xlabel("Duration (days)")
            plt.ylabel("Count")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, "infection_duration_by_recovery_method.png"), dpi=300)
            plt.close()
        
        # Asymptomatic recovery analysis
        if 'behavior_infectee' in df.columns:
            asymptomatic_df = df[df['infection_symptom_status'] == 'asymptomatic'].copy()
            asymptomatic_df = asymptomatic_df[asymptomatic_df['duration_source'].isin(['asymptomatic', 'natural_clearance'])]
            
            if not asymptomatic_df.empty:
                # Summary by behavior
                behavior_summary = []
                for behavior in asymptomatic_df['behavior_infectee'].unique():
                    behavior_data = asymptomatic_df[asymptomatic_df['behavior_infectee'] == behavior]
                    
                    total = len(behavior_data)
                    natural_count = len(behavior_data[behavior_data['duration_source'] == 'natural_clearance'])
                    asymptomatic_count = len(behavior_data[behavior_data['duration_source'] == 'asymptomatic'])
                    
                    behavior_summary.append({
                        'behavior': behavior,
                        'total_asymptomatic_infections': total,
                        'natural_clearance_count': natural_count,
                        'asymptomatic_recovery_count': asymptomatic_count,
                        'proportion_natural_clearance': natural_count / total if total > 0 else 0,
                        'proportion_asymptomatic_recovery': asymptomatic_count / total if total > 0 else 0
                    })
                
                summary_df = pd.DataFrame(behavior_summary)
                summary_df.to_csv(os.path.join(additional_outputs_dir, "asymptomatic_recovery_proportions_by_behavior.csv"), index=False)
        # Missing Plot 1: Asymptomatic recovery by duration (stacked histogram)
        if 'duration_source' in df.columns and 'infection_symptom_status' in df.columns:
            asymptomatic_only = df[df['infection_symptom_status'] == 'asymptomatic'].copy()
            asymptomatic_only = asymptomatic_only[asymptomatic_only['duration_source'].isin(['asymptomatic', 'natural_clearance'])]
            
            if not asymptomatic_only.empty:
                plt.figure(figsize=(12, 7))
                
                # Create duration bins
                max_duration = asymptomatic_only['duration'].max()
                bins = np.linspace(0, max_duration, 40)
                
                # Get data for each recovery method
                natural_clearance_durations = asymptomatic_only[asymptomatic_only['duration_source'] == 'natural_clearance']['duration']
                asymptomatic_recovery_durations = asymptomatic_only[asymptomatic_only['duration_source'] == 'asymptomatic']['duration']
                
                # Create stacked histogram
                plt.hist([natural_clearance_durations, asymptomatic_recovery_durations], 
                        bins=bins, 
                        label=['Natural Clearance', 'Asymptomatic Recovery'],
                        color=['green', 'blue'],
                        alpha=0.7,
                        edgecolor='black',
                        stacked=True)
                
                plt.title("Recovery Method for Asymptomatic Infections by Duration")
                plt.xlabel("Duration (days)")
                plt.ylabel("Count of Infections")
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                plt.savefig(os.path.join(plots_dir, "asymptomatic_recovery_by_duration.png"), dpi=300)
                plt.close()
        
        # Missing Plot 2: Asymptomatic recovery proportions by behavior (bar chart)
        if 'behavior_infectee' in df.columns and 'duration_source' in df.columns and 'infection_symptom_status' in df.columns:
            asymptomatic_df = df[df['infection_symptom_status'] == 'asymptomatic'].copy()
            asymptomatic_df = asymptomatic_df[asymptomatic_df['duration_source'].isin(['asymptomatic', 'natural_clearance'])]
            
            if not asymptomatic_df.empty:
                # Calculate proportions by behavior (reuse the logic you already have)
                behavior_summary = []
                for behavior in asymptomatic_df['behavior_infectee'].unique():
                    behavior_data = asymptomatic_df[asymptomatic_df['behavior_infectee'] == behavior]
                    
                    total = len(behavior_data)
                    natural_count = len(behavior_data[behavior_data['duration_source'] == 'natural_clearance'])
                    asymptomatic_count = len(behavior_data[behavior_data['duration_source'] == 'asymptomatic'])
                    
                    if total > 0:
                        behavior_summary.append({
                            'behavior': behavior,
                            'proportion_natural_clearance': natural_count / total,
                            'proportion_asymptomatic_recovery': asymptomatic_count / total
                        })
                
                if behavior_summary:
                    summary_df = pd.DataFrame(behavior_summary).sort_values('behavior')
                    
                    # Create bar chart
                    plt.figure(figsize=(10, 6))
                    
                    behaviors = summary_df['behavior']
                    natural_props = summary_df['proportion_natural_clearance']
                    asymptomatic_props = summary_df['proportion_asymptomatic_recovery']
                    
                    x = np.arange(len(behaviors))
                    width = 0.35
                    
                    plt.bar(x - width/2, natural_props, width, label='Natural Clearance', 
                        color='green', alpha=0.7, edgecolor='black')
                    plt.bar(x + width/2, asymptomatic_props, width, label='Asymptomatic Recovery', 
                        color='blue', alpha=0.7, edgecolor='black')
                    
                    plt.xlabel('Behavior Group')
                    plt.ylabel('Proportion of Asymptomatic Infections')
                    plt.title('Recovery Method Proportions for Asymptomatic Infections by Behavior')
                    plt.xticks(x, behaviors)
                    plt.legend()
                    plt.grid(True, alpha=0.3, axis='y')
                    plt.tight_layout()
                    plt.savefig(os.path.join(plots_dir, "asymptomatic_recovery_proportions_by_behavior.png"), dpi=300)
                    plt.close()
        
        print(f"Infection recovery analysis complete!")

def plot_repeat_infections(transmission_df, nodes_dict, plots_dir, additional_outputs_dir, repeat_threshold=10):
    """Category 6: Repeat infections"""
    
    print("Creating repeat infection analysis...")
    
    if transmission_df.empty:
        return
    
    # Filter successful transmissions
    df = transmission_df[transmission_df['superseded_simultaneous'] == False].copy()
    df = df[df['day_of_sampling'].notna()].copy()
    
    # Count infections per person
    infection_counts = df['infectee_node'].value_counts()
    repeat_infectors = infection_counts[infection_counts >= repeat_threshold]
    
    print(f"Found {len(repeat_infectors)} people with {repeat_threshold}+ infections")
    
    if len(repeat_infectors) == 0:
        return
    
    # Infection count distribution
    plt.figure(figsize=(12, 8))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Overall distribution
    ax1.hist(infection_counts.values, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.axvline(repeat_threshold, color='red', linestyle='--', 
                label=f'Repeat threshold ({repeat_threshold})')
    ax1.set_xlabel('Number of Infections')
    ax1.set_ylabel('Number of People')
    ax1.set_title('Distribution of Infection Counts per Person')
    ax1.set_yscale('log')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Detail view of repeat infectors
    if len(repeat_infectors) > 0:
        ax2.hist(repeat_infectors.values, bins=min(20, len(repeat_infectors)), 
                 edgecolor='black', alpha=0.7, color='orange')
        ax2.set_xlabel('Number of Infections')
        ax2.set_ylabel('Number of People')
        ax2.set_title(f'Repeat Infectors ({repeat_threshold}+ infections)')
        ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "infection_count_distribution.png"), dpi=300)
    plt.close()

    # 1. Partner diversity analysis
    print("Analyzing transmission partner diversity...")
    partner_analysis_data = []
    
    for node_id in repeat_infectors.index:
        node_infections = df[df['infectee_node'] == node_id]
        
        total_infections = len(node_infections)
        unique_partners = node_infections['infector_node'].nunique()
        most_common_partner_count = node_infections['infector_node'].value_counts().iloc[0]
        most_common_partner_id = node_infections['infector_node'].value_counts().index[0]
        
        partner_analysis_data.append({
            'node_id': node_id,
            'total_infections': total_infections,
            'unique_partners': unique_partners,
            'partner_diversity': unique_partners / total_infections,
            'most_common_partner_count': most_common_partner_count,
            'most_common_partner_id': most_common_partner_id,
            'most_common_partner_prop': most_common_partner_count / total_infections
        })
    
    partner_df = pd.DataFrame(partner_analysis_data)
    partner_df.to_csv(os.path.join(additional_outputs_dir, "repeat_infections_partner_diversity.csv"), index=False)
    
    # Partner diversity plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Partner diversity distribution
    ax1.hist(partner_df['partner_diversity'], bins=20, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.set_xlabel('Partner Diversity (unique partners / total infections)')
    ax1.set_ylabel('Number of People')
    ax1.set_title('Transmission Partner Diversity')
    ax1.grid(True, alpha=0.3)
    
    # Unique partners vs total infections
    ax2.scatter(partner_df['total_infections'], partner_df['unique_partners'], alpha=0.6, color='orange')
    max_inf = partner_df['total_infections'].max()
    ax2.plot([0, max_inf], [0, max_inf], 'r--', label='Perfect diversity line')
    ax2.set_xlabel('Total Infections')
    ax2.set_ylabel('Unique Partners')
    ax2.set_title('Partner Count vs Infection Count')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Most common partner proportion
    ax3.hist(partner_df['most_common_partner_prop'], bins=20, edgecolor='black', alpha=0.7, color='green')
    ax3.set_xlabel('Proportion from Most Common Partner')
    ax3.set_ylabel('Number of People')
    ax3.set_title('Dependence on Single Partner')
    ax3.grid(True, alpha=0.3)
    
    # Summary statistics
    summary_stats = [
        f"Avg partner diversity: {partner_df['partner_diversity'].mean():.3f}",
        f"People with single partner: {(partner_df['unique_partners'] == 1).sum()}",
        f"People with >50% from one partner: {(partner_df['most_common_partner_prop'] > 0.5).sum()}"
    ]
    ax4.text(0.1, 0.5, '\n'.join(summary_stats), transform=ax4.transAxes, 
             fontsize=12, verticalalignment='center')
    ax4.set_title('Summary Statistics')
    ax4.axis('off')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "repeat_infections_partner_analysis.png"), dpi=300)
    plt.close()
    
    # 2. Spacing impact analysis
    print("Analyzing infection spacing impact...")
    
    def calculate_distinct_infections(node_infections, spacing_days):
        """Calculate distinct infections with minimum spacing requirement"""
        if len(node_infections) <= 1:
            return len(node_infections)
        
        node_infections_sorted = node_infections.sort_values('day_of_transmission')
        distinct_infections = []
        
        for _, infection in node_infections_sorted.iterrows():
            tx_day = infection['day_of_transmission']
            sample_day = infection['day_of_sampling']
            
            # Check if this infection is distinct from previous ones
            is_distinct = True
            
            for prev_infection in distinct_infections:
                prev_sample = prev_infection['day_of_sampling']
                
                # If current infection starts within spacing_days of previous recovery
                if tx_day < prev_sample + spacing_days:
                    is_distinct = False
                    break
            
            if is_distinct:
                distinct_infections.append({
                    'day_of_transmission': tx_day,
                    'day_of_sampling': sample_day
                })
        
        return len(distinct_infections)
    
    spacing_analysis_data = []
    
    for node_id in repeat_infectors.index:
        node_infections = df[df['infectee_node'] == node_id].copy()
        total_infections = len(node_infections)
        
        # Calculate distinct infections with different spacing requirements
        distinct_30_day = calculate_distinct_infections(node_infections, spacing_days=30)
        distinct_60_day = calculate_distinct_infections(node_infections, spacing_days=60)
        
        spacing_analysis_data.append({
            'node_id': node_id,
            'total_infections': total_infections,
            'distinct_infections_30_day': distinct_30_day,
            'distinct_infections_60_day': distinct_60_day,
            'reduction_30_day': (total_infections - distinct_30_day) / total_infections,
            'reduction_60_day': (total_infections - distinct_60_day) / total_infections
        })
    
    spacing_df = pd.DataFrame(spacing_analysis_data)
    spacing_df.to_csv(os.path.join(additional_outputs_dir, "repeat_infections_spacing_analysis.csv"), index=False)
    
    # Spacing impact plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Comparison: total vs distinct (30-day)
    ax1.scatter(spacing_df['total_infections'], spacing_df['distinct_infections_30_day'], alpha=0.6)
    max_val = max(spacing_df['total_infections'].max(), spacing_df['distinct_infections_30_day'].max())
    ax1.plot([0, max_val], [0, max_val], 'r--', label='No reduction line')
    ax1.set_xlabel('Total Infections')
    ax1.set_ylabel('Distinct Infections (30-day spacing)')
    ax1.set_title('Impact of 30-day Spacing Requirement')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Comparison: total vs distinct (60-day)
    ax2.scatter(spacing_df['total_infections'], spacing_df['distinct_infections_60_day'], 
                alpha=0.6, color='orange')
    ax2.plot([0, max_val], [0, max_val], 'r--', label='No reduction line')
    ax2.set_xlabel('Total Infections')
    ax2.set_ylabel('Distinct Infections (60-day spacing)')
    ax2.set_title('Impact of 60-day Spacing Requirement')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Reduction proportions (30-day)
    ax3.hist(spacing_df['reduction_30_day'], bins=20, edgecolor='black', alpha=0.7, color='steelblue')
    ax3.set_xlabel('Proportion of Infections Eliminated')
    ax3.set_ylabel('Number of People')
    ax3.set_title('Infection Reduction with 30-day Spacing')
    ax3.grid(True, alpha=0.3)
    
    # Reduction proportions (60-day)
    ax4.hist(spacing_df['reduction_60_day'], bins=20, edgecolor='black', alpha=0.7, color='orange')
    ax4.set_xlabel('Proportion of Infections Eliminated')
    ax4.set_ylabel('Number of People')
    ax4.set_title('Infection Reduction with 60-day Spacing')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "repeat_infections_spacing_impact.png"), dpi=300)
    plt.close()
    
    # 3. Symptom and recovery patterns analysis
    print("Analyzing symptom and recovery patterns...")
    
    symptom_analysis_data = []
    recovery_analysis_data = []
    
    for node_id in repeat_infectors.index:
        node_infections = df[df['infectee_node'] == node_id]
        
        total_infections = len(node_infections)
        
        # Symptom analysis
        if 'infection_symptom_status' in df.columns:
            symptomatic_count = (node_infections['infection_symptom_status'] == 'symptomatic').sum()
            asymptomatic_count = (node_infections['infection_symptom_status'] == 'asymptomatic').sum()
            
            symptom_analysis_data.append({
                'node_id': node_id,
                'total_infections': total_infections,
                'symptomatic_count': symptomatic_count,
                'asymptomatic_count': asymptomatic_count,
                'prop_symptomatic': symptomatic_count / total_infections,
                'prop_asymptomatic': asymptomatic_count / total_infections
            })
        
        # Recovery analysis
        if 'duration_source' in df.columns:
            recovery_counts = node_infections['duration_source'].value_counts()
            natural_clearance = recovery_counts.get('natural_clearance', 0)
            symptomatic_recovery = recovery_counts.get('symptomatic', 0)
            asymptomatic_recovery = recovery_counts.get('asymptomatic', 0)
            
            recovery_analysis_data.append({
                'node_id': node_id,
                'total_infections': total_infections,
                'natural_clearance': natural_clearance,
                'symptomatic_recovery': symptomatic_recovery,
                'asymptomatic_recovery': asymptomatic_recovery,
                'prop_natural_clearance': natural_clearance / total_infections,
                'prop_symptomatic_recovery': symptomatic_recovery / total_infections,
                'prop_asymptomatic_recovery': asymptomatic_recovery / total_infections
            })
    
    # Save CSVs
    if symptom_analysis_data:
        symptom_df = pd.DataFrame(symptom_analysis_data)
        symptom_df.to_csv(os.path.join(additional_outputs_dir, "repeat_infections_symptom_patterns.csv"), index=False)
    
    if recovery_analysis_data:
        recovery_df = pd.DataFrame(recovery_analysis_data)
        recovery_df.to_csv(os.path.join(additional_outputs_dir, "repeat_infections_recovery_patterns.csv"), index=False)
    
    # Symptom and recovery plots
    if symptom_analysis_data and recovery_analysis_data:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Proportion symptomatic distribution
        ax1.hist(symptom_df['prop_symptomatic'], bins=20, edgecolor='black', alpha=0.7, color='red')
        ax1.set_xlabel('Proportion of Infections that were Symptomatic')
        ax1.set_ylabel('Number of People')
        ax1.set_title('Symptom Pattern Distribution (Repeat Infectors)')
        ax1.grid(True, alpha=0.3)
        
        # Recovery method proportions
        ax2.hist(recovery_df['prop_natural_clearance'], bins=20, edgecolor='black', 
                 alpha=0.7, color='green', label='Natural Clearance')
        ax2.set_xlabel('Proportion of Infections with Natural Clearance')
        ax2.set_ylabel('Number of People')
        ax2.set_title('Natural Clearance Pattern Distribution')
        ax2.grid(True, alpha=0.3)
        
        # Scatter: total infections vs proportion symptomatic
        ax3.scatter(symptom_df['total_infections'], symptom_df['prop_symptomatic'], alpha=0.6)
        ax3.set_xlabel('Total Infections')
        ax3.set_ylabel('Proportion Symptomatic')
        ax3.set_title('Infections vs Symptom Rate')
        ax3.grid(True, alpha=0.3)
        
        # Recovery method stacked bar (summary)
        recovery_summary = recovery_df[['prop_natural_clearance', 'prop_symptomatic_recovery', 
                                       'prop_asymptomatic_recovery']].mean()
        ax4.bar(['Natural\nClearance', 'Symptomatic\nRecovery', 'Asymptomatic\nRecovery'], 
                recovery_summary.values, color=['green', 'red', 'blue'], alpha=0.7, edgecolor='black')
        ax4.set_ylabel('Average Proportion')
        ax4.set_title('Average Recovery Method Distribution')
        ax4.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "repeat_infections_symptom_recovery_analysis.png"), dpi=300)
        plt.close()
    
    # Create comprehensive summary
    print("Creating comprehensive summary...")
    summary_data = []
    
    for node_id in repeat_infectors.index:
        node_infections = df[df['infectee_node'] == node_id]
        node_info = nodes_dict.get(node_id, {})
        
        # Basic info
        total_infections = len(node_infections)
        behavior = node_info.get('behavior', 'Unknown')
        hi_risk = node_info.get('hi_risk', 'Unknown')
        
        # Partner diversity
        unique_partners = node_infections['infector_node'].nunique()
        
        # Time span
        first_infection = node_infections['day_of_transmission'].min()
        last_infection = node_infections['day_of_transmission'].max()
        time_span = last_infection - first_infection
        
        # Distinct infections with spacing
        distinct_30 = calculate_distinct_infections(node_infections.sort_values('day_of_transmission'), 30)
        distinct_60 = calculate_distinct_infections(node_infections.sort_values('day_of_transmission'), 60)
        
        summary_data.append({
            'node_id': node_id,
            'behavior': behavior,
            'hi_risk': hi_risk,
            'total_infections': total_infections,
            'unique_partners': unique_partners,
            'partner_diversity': unique_partners / total_infections,
            'first_infection_day': first_infection,
            'last_infection_day': last_infection,
            'infection_time_span_days': time_span,
            'distinct_infections_30_day': distinct_30,
            'distinct_infections_60_day': distinct_60,
            'reduction_30_day': (total_infections - distinct_30) / total_infections,
            'reduction_60_day': (total_infections - distinct_60) / total_infections
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('total_infections', ascending=False)
    summary_df.to_csv(os.path.join(additional_outputs_dir, "repeat_infections_comprehensive_summary.csv"), index=False)
    
    print(f"Repeat infection analysis complete. Summary saved for {len(summary_df)} individuals.")

def check_optional_files(output_dir):
    """Check which optional files exist"""
    
    optional_files = {
        'unmatched': os.path.join(output_dir, "unmatched_df.csv"),
        'seekers': os.path.join(output_dir, "seeker_df.csv")  # CHANGED: from seeker_time_series.csv
    }
    
    available = {key: os.path.exists(path) for key, path in optional_files.items()}
    return available

def plot_unmatched_dynamics(output_dir, nodes_dict, partnership_burnin_days, tracking_start_day, plots_dir):
    """Category 8: Unmatched/seeking dynamics (requires optional files)"""
    
    print("Creating unmatched dynamics plots...")
    
    try:
        unmatched_df = pd.read_csv(os.path.join(output_dir, "unmatched_df.csv"))
        seeker_df = pd.read_csv(os.path.join(output_dir, "seeker_df.csv"))
    except FileNotFoundError as e:
        print(f"Optional file missing: {e}")
        return
    
    # Calculate population by behavior
    pop_by_behavior = defaultdict(int)
    for _, d in nodes_dict.items():
        pop_by_behavior[d['behavior']] += 1
    
    print(f"Population by behavior: {dict(pop_by_behavior)}")
    
    # Filter to post-burnin period
    post_unmatched = unmatched_df[unmatched_df['day'] >= partnership_burnin_days].copy()
    
    # 1. Absolute counts (unchanged)
    fig, ax = plt.subplots(figsize=(10, 6))
    for behavior in ['WSM', 'MSW', 'MSM', 'MSMW']:
        if behavior in post_unmatched.columns:
            ax.plot(post_unmatched['day'], post_unmatched[behavior], 
                   label=behavior, alpha=0.7, linewidth=2)
    ax.set_title("Unmatched Individuals Over Time (post partnership burn-in)")
    ax.set_xlabel("Day")
    ax.set_ylabel("Unmatched count")
    ax.legend()
    ax.grid(True, alpha=0.3)
    _shade_burnin(ax, partnership_burnin_days, tracking_start_day)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "unmatched_counts_postburnin.png"), dpi=300)
    plt.close()
    
    # 2. FIXED: Actually normalized by population
    fig, ax = plt.subplots(figsize=(10, 6))
    for behavior in ['WSM', 'MSW', 'MSM', 'MSMW']:
        if behavior in post_unmatched.columns:
            denom = pop_by_behavior.get(behavior, 0)
            if denom > 0:
                # Actually normalize by dividing by population
                normalized_counts = post_unmatched[behavior] / denom
                ax.plot(post_unmatched['day'], normalized_counts, 
                       label=behavior, alpha=0.7, linewidth=2)
    ax.set_title("Normalized Unmatched Individuals Over Time (post partnership burn-in)")
    ax.set_xlabel("Day")
    ax.set_ylabel("Proportion unmatched")  # Now actually proportions!
    ax.legend()
    ax.grid(True, alpha=0.3)
    _shade_burnin(ax, partnership_burnin_days, tracking_start_day)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "unmatched_normalized_postburnin.png"), dpi=300)
    plt.close()
    
    # 3. Proportion of seekers who were unmatched - by behavior (lines)
    post_seekers = seeker_df[seeker_df['day'] >= partnership_burnin_days].copy()
    merged = post_unmatched.merge(post_seekers, on='day', how='inner')
    
    fig, ax = plt.subplots(figsize=(12, 6))
    behaviors = ['WSM', 'MSW', 'MSM', 'MSMW']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for i, behavior in enumerate(behaviors):
        unmatched_col = behavior
        seeker_col = f'{behavior}_seekers'
        if unmatched_col in merged.columns and seeker_col in merged.columns:
            # Calculate proportion: unmatched / seekers
            proportion = merged[unmatched_col] / merged[seeker_col].replace(0, np.nan)
            proportion = proportion.fillna(0)
            
            ax.plot(merged['day'], proportion, 
                    label=behavior, color=colors[i], alpha=0.7, linewidth=1.5)
    
    ax.set_title("Proportion of Seekers Who Were Unmatched by Behavior")
    ax.set_xlabel("Day")
    ax.set_ylabel("Proportion Unmatched Among Seekers")
    ax.set_ylim(0, 1)
    ax.legend()
    ax.grid(True, alpha=0.3)
    _shade_burnin(ax, partnership_burnin_days, tracking_start_day)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "unmatched_proportion_by_behavior_lines.png"), dpi=300)
    plt.close()
    
    # 4. Faceted plot (each behavior in own subplot)
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True, sharey=True)
    axes = axes.flatten()
    
    for i, behavior in enumerate(behaviors):
        ax = axes[i]
        unmatched_col = behavior
        seeker_col = f'{behavior}_seekers'
        
        if unmatched_col in merged.columns and seeker_col in merged.columns:
            proportion = merged[unmatched_col] / merged[seeker_col].replace(0, np.nan)
            proportion = proportion.fillna(0)
            
            ax.plot(merged['day'], proportion, color=colors[i], alpha=0.8, linewidth=1.5)
            
            # Add smoothing for trend visibility
            if len(proportion) > 50:
                smoothed = proportion.rolling(window=7, center=True).mean()
                ax.plot(merged['day'], smoothed, color='black', alpha=0.5, 
                       linewidth=2, linestyle='--')
        
        ax.set_title(f'{behavior} - Proportion Unmatched')
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)
        _shade_burnin(ax, partnership_burnin_days, tracking_start_day)
    
    fig.text(0.5, 0.04, 'Day', ha='center', fontsize=12)
    fig.text(0.04, 0.5, 'Proportion of Seekers Unmatched', va='center', rotation='vertical', fontsize=12)
    fig.suptitle('Unmatched Proportions by Behavior Group (Faceted)', fontsize=14)
    
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    fig.savefig(os.path.join(plots_dir, "unmatched_proportion_by_behavior_faceted.png"), dpi=300)
    plt.close()
    
    # 5. Weekly aggregated proportions  
    merged['week'] = (merged['day'] - merged['day'].min()) // 7
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True, sharey=True)
    axes = axes.flatten()
    
    for i, behavior in enumerate(behaviors):
        ax = axes[i]
        unmatched_col = behavior
        seeker_col = f'{behavior}_seekers'
        
        if unmatched_col in merged.columns and seeker_col in merged.columns:
            # Weekly aggregation
            weekly = merged.groupby('week').agg({
                unmatched_col: 'sum',
                seeker_col: 'sum',
                'day': 'mean'  # Average day for x-axis
            }).reset_index()
            
            # Calculate weekly proportion
            weekly['proportion'] = weekly[unmatched_col] / weekly[seeker_col].replace(0, np.nan)
            weekly['proportion'] = weekly['proportion'].fillna(0)
            
            # Only plot weeks with sufficient seekers
            sufficient_data = weekly[weekly[seeker_col] >= 5]
            
            if len(sufficient_data) > 0:
                ax.plot(sufficient_data['day'], sufficient_data['proportion'], 
                       color=colors[i], alpha=0.8, linewidth=2, marker='o', markersize=3)
        
        ax.set_title(f'{behavior} - Weekly Proportion Unmatched')
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)
    
    fig.text(0.5, 0.04, 'Day', ha='center', fontsize=12)
    fig.text(0.04, 0.5, 'Weekly Proportion of Seekers Unmatched', va='center', rotation='vertical', fontsize=12)
    fig.suptitle('Weekly Unmatched Proportions by Behavior Group', fontsize=14)
    
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    fig.savefig(os.path.join(plots_dir, "unmatched_proportion_weekly.png"), dpi=300)
    plt.close()
    
    print("Unmatched dynamics analysis complete!")

# ========== MAIN FUNCTIONS ==========

def load_simulation_data(output_dir):
    """Load simulation data from output directory"""
    
    print(f"Loading simulation data from: {output_dir}")
    
    # Core files
    transmission_df = pd.read_csv(os.path.join(output_dir, "transmission_df.csv"))
    nodes_df = pd.read_csv(os.path.join(output_dir, "nodes_df.csv"), index_col=0)
    nodes_dict = nodes_df.to_dict('index')
    
    # Optional files
    try:
        edge_df = pd.read_csv(os.path.join(output_dir, "edge_df.csv"))
    except FileNotFoundError:
        edge_df = pd.DataFrame()
    
    # Parameters
    try:
        with open(os.path.join(output_dir, "parameters_used.json"), 'r') as f:
            params = json.load(f)
        partnership_burnin_days = params['simulation']['partnership_burnin_days']
        transmission_burnin_days = params['simulation']['transmission_burnin_days']  # ADD THIS
        tracking_start_day = partnership_burnin_days + transmission_burnin_days
        sim_days = params['simulation']['sim_days']
    except FileNotFoundError:
        raise FileNotFoundError(f"parameters_used.json not found in {output_dir}")
    
    print(f"Loaded:")
    print(f"  Nodes: {len(nodes_dict)}")
    print(f"  Transmissions: {len(transmission_df)}")
    print(f"  Partnerships: {len(edge_df)}")
    print(f"  Partnership burnin: {partnership_burnin_days}")
    print(f"  Transmission burnin: {transmission_burnin_days}")  # ADD THIS
    print(f"  Tracking start: {tracking_start_day}")
    print(f"  Sim days: {sim_days}")
    
    return {
        'nodes_dict': nodes_dict,
        'transmission_df': transmission_df,
        'edge_df': edge_df,
        'partnership_burnin_days': partnership_burnin_days,
        'transmission_burnin_days': transmission_burnin_days,  # ADD THIS
        'tracking_start_day': tracking_start_day,
        'sim_days': sim_days
    }

def create_all_plots(data, plot_categories, output_dir):
    """Create all requested plot categories"""
    
    # Create output directories
    plots_dir = os.path.join(output_dir, "plots")
    additional_outputs_dir = os.path.join(output_dir, "additional_outputs")
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(additional_outputs_dir, exist_ok=True)
    
    # Always include basic epidemiology (essential)
    print("\n" + "="*60)
    print("RECONSTRUCTING EPIDEMIOLOGICAL DATA FROM CORE OUTPUTS")
    print("="*60)
    
    epi_data = reconstruct_daily_epidemiology(
    data['transmission_df'], data['nodes_dict'],
    data['partnership_burnin_days'], data['transmission_burnin_days'], data['sim_days']  # ADD transmission_burnin_days
)
    
    if 'basic_epi' in plot_categories or 'all' in plot_categories:
        plot_basic_epidemiology(
    epi_data, data['nodes_dict'], data['transmission_df'],  # Added transmission_df
    data['partnership_burnin_days'], data['tracking_start_day'],
    plots_dir, additional_outputs_dir
)
    
    if 'partnerships' in plot_categories or 'all' in plot_categories:
        plot_partnerships(
            data['edge_df'], data['nodes_dict'],
            data['partnership_burnin_days'], data['tracking_start_day'],
            plots_dir, additional_outputs_dir
        )
    
    if 'concurrency' in plot_categories or 'all' in plot_categories:
        plot_concurrency(
            data['edge_df'], data['nodes_dict'],
            plots_dir, additional_outputs_dir
        )
    
    if 'transmission' in plot_categories or 'all' in plot_categories:
        plot_transmission(
            data['transmission_df'], data['edge_df'],
            data['partnership_burnin_days'], data['tracking_start_day'],
            plots_dir, additional_outputs_dir
        )
    
    if 'infection_recovery' in plot_categories or 'all' in plot_categories:
        plot_infection_recovery(
            data['transmission_df'], plots_dir, additional_outputs_dir
        )
    
    if 'repeat_infections' in plot_categories or 'all' in plot_categories:
        plot_repeat_infections(
            data['transmission_df'], data['nodes_dict'],
            plots_dir, additional_outputs_dir
        )
    
    # Optional: unmatched dynamics
    if 'unmatched' in plot_categories or 'all' in plot_categories:
        available_optional = check_optional_files(output_dir)
        if available_optional['unmatched'] and available_optional['seekers']:
            plot_unmatched_dynamics(
                output_dir, data['nodes_dict'],  # ADD nodes_dict here
                data['partnership_burnin_days'], data['tracking_start_day'],
                plots_dir
            )
        else:
            print("Skipping unmatched dynamics - optional files not available")

def main():
    parser = argparse.ArgumentParser(description='Create comprehensive plots from simulation outputs')
    parser.add_argument('output_dir', help='Directory containing simulation outputs')
    parser.add_argument('--categories', nargs='+', 
                       choices=['basic_epi', 'partnerships', 'concurrency', 'transmission', 
                               'infection_recovery', 'repeat_infections', 'unmatched', 'all'],
                       default=['all'], help='Plot categories to create')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        print(f"❌ Directory does not exist: {args.output_dir}")
        return 1
    
    # Expand 'all' category
    if 'all' in args.categories:
        plot_categories = ['basic_epi', 'partnerships', 'concurrency', 'transmission', 
                          'infection_recovery', 'repeat_infections', 'unmatched']
    else:
        plot_categories = args.categories
    
    try:
        print(f"🎨 Creating plots for categories: {', '.join(plot_categories)}")
        
        # Load data
        data = load_simulation_data(args.output_dir)
        
        # Create plots
        create_all_plots(data, plot_categories, args.output_dir)
        
        print(f"\n✅ Plots saved to: {os.path.join(args.output_dir, 'plots')}")
        print(f"✅ Additional outputs saved to: {os.path.join(args.output_dir, 'additional_outputs')}")
        return 0
        
    except Exception as e:
        print(f"❌ Plotting failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())