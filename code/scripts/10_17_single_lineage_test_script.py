#!/usr/bin/env python3
"""
Lineage Assessment Script
Analyzes simulation outputs to determine if single lineage dominates post-burnin period
"""

import pandas as pd
import numpy as np
import argparse
import os
import json
import itertools
import matplotlib.pyplot as plt
from collections import defaultdict

# ========== HELPER PLOTTING FUNCTIONS ==========

def _prep_daily_df(list_of_dicts):
    """
    list_of_dicts -> tidy DataFrame indexed by day, with 'phase' column preserved.
    Non-count columns: 'day', 'phase'. Everything else treated as count columns.
    """
    if not list_of_dicts:
        return pd.DataFrame()
    df = pd.DataFrame(list_of_dicts).fillna(0)
    if 'day' not in df.columns:
        # fall back to implicit day index if needed
        df['day'] = np.arange(len(df))
    if 'phase' not in df.columns:
        df['phase'] = 'post_tx'
    df = df.sort_values('day').set_index('day')
    return df

def _shade_burnin(ax, partnership_burnin_days: int, tracking_start_day: int):
    """Shade the transmission burn-in, draw a vertical line at the end (tracking_start_day)."""
    if partnership_burnin_days is not None and tracking_start_day is not None:
        if tracking_start_day > partnership_burnin_days:
            ax.axvspan(partnership_burnin_days, tracking_start_day - 1, alpha=0.10, label="Tx burn-in")
        ax.axvline(tracking_start_day, linestyle='--', linewidth=1.5, label="Burn-in over")

def _line_plot_multi(df, columns, title, ylabel, path, partnership_burnin_days, tracking_start_day, legend_cols=1):
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
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)

# ========== CORE LINEAGE FUNCTIONS ==========

def build_lineage_episodes(initial_infectors, initial_infectors_day, transmission_df):
    """
    initial_infectors: list[(node_id, remaining_days)] from the sim result
    initial_infectors_day: int (the seeding day; e.g. partnership_burnin_days)
    transmission_df: must include day_of_transmission, day_of_sampling, infector_node, infectee_node
                     and losers have day_of_sampling == NaN (ignored)
    Returns: list of dicts with keys: lineage, node, start, end
    """
    # 1) Seeded episodes create their own lineages
    episodes = []
    active_by_node = defaultdict(list)  # node -> list of (start,end,lineage)
    for node_id, remaining in initial_infectors:
        start = int(initial_infectors_day)
        end   = start + int(remaining)
        lineage = int(node_id)          # lineage id = root node id
        ep = {"lineage": lineage, "node": int(node_id), "start": start, "end": end}
        episodes.append(ep)
        active_by_node[int(node_id)].append((start, end, lineage))

    # 2) Transmissions in time order → inherit infector lineage
    tx = transmission_df.copy()
    tx = tx[pd.notna(tx["day_of_sampling"])].copy()  # ignore losers (NA end time)
    tx["day_of_transmission"] = tx["day_of_transmission"].astype(int)
    tx["day_of_sampling"]     = tx["day_of_sampling"].astype(int)
    tx = tx.sort_values("day_of_transmission")

    def lineage_of_infector_at(infector, day):
        # find an active episode that covers 'day' (start <= day < end)
        for s,e,L in active_by_node.get(int(infector), []):
            if s <= day < e:
                return L
        return None  # should not happen if logs are consistent

    for _, r in tx.iterrows():
        t   = int(r["day_of_transmission"])
        end = int(r["day_of_sampling"])
        inf = int(r["infector_node"])
        sus = int(r["infectee_node"])

        L = lineage_of_infector_at(inf, t)
        if L is None:
            # Fallback: if something's off, treat infectee as starting a micro-lineage
            L = inf

        ep = {"lineage": L, "node": sus, "start": t, "end": end}
        episodes.append(ep)
        active_by_node[sus].append((t, end, L))

    return episodes

def lineage_prevalence_timeseries(episodes, total_days=None):
    """
    Returns a DataFrame indexed by day, columns=lineage_id, values=#currently infected in that lineage.
    """
    if not episodes:
        return pd.DataFrame()

    if total_days is None:
        total_days = max(e["end"] for e in episodes)

    # per-lineage delta array: +1 at start, -1 at end (end is exclusive)
    deltas = defaultdict(lambda: np.zeros(total_days+1, dtype=np.int32))
    for e in episodes:
        L = int(e["lineage"])
        s = int(e["start"])
        eend = int(e["end"])
        s = max(0, s)
        if s <= total_days:
            deltas[L][s] += 1
        if eend <= total_days:
            deltas[L][eend] -= 1

    # cumulative sum → prevalence per lineage
    data = {}
    for L, arr in deltas.items():
        data[L] = np.cumsum(arr)[:total_days]  # days [0..total_days-1]

    df = pd.DataFrame(data)
    df.index.name = "day"
    df = df.fillna(0).astype(int)
    return df  

def plot_lineage_prevalence(
    df_lineage_prev,
    initial_infectors_day,
    tracking_start_day,
    top_k=10,
    highlight=None,           # optional iterable of lineage IDs to color & label
    smooth=None,              # e.g., 7 for rolling mean
    save_path=None,
    other_color="gray",
    other_alpha=0.35,
    other_lw=1.0,
    top_alpha=0.9,
    top_lw=2.0,
):
    if df_lineage_prev.empty:
        return

    # Make sure we start with a fresh figure
    fig, ax = plt.subplots(figsize=(12, 7))

    # optional smoothing
    plot_df = df_lineage_prev.copy()
    if smooth and smooth > 1:
        plot_df = plot_df.rolling(window=int(smooth), min_periods=1, center=True).mean()

    # choose which lineages to highlight
    if highlight is None:
        # Use burn-in period for highlighting (where lineages compete before extinction)
        burnin_period = plot_df.loc[(plot_df.index >= initial_infectors_day) & 
                                    (plot_df.index < tracking_start_day)]
        if not burnin_period.empty:
            burnin_max = burnin_period.max()
            highlight = list(burnin_max.sort_values(ascending=False).head(top_k).index)
        else:
            # Fallback to overall max if no burn-in data
            maxvals = plot_df.max(axis=0)
            highlight = list(maxvals.sort_values(ascending=False).head(top_k).index)
    else:
        # keep only those that exist as columns
        highlight = [h for h in highlight if h in plot_df.columns]

    # cast to a common type for robust membership checks
    hset = set(map(int, highlight))
    cols_int = [int(c) for c in plot_df.columns]
    others = [c for c in cols_int if c not in hset]

    # palette for highlighted lineages
    colors = plt.rcParams['axes.prop_cycle'].by_key().get('color', ['C0','C1','C2','C3','C4','C5','C6','C7'])
    color_cycle = itertools.cycle(colors)

    # plot the "other" lineages first (thin gray, suppressed legend)
    for col in others:
        y = plot_df[int(col)].values
        ax.plot(plot_df.index, y, color=other_color, alpha=other_alpha, linewidth=other_lw,
                label=f"_{col}")  # leading underscore: ignore in legend

    # plot highlighted lineages (colored + legend)
    for col in highlight:
        y = plot_df[int(col)].values
        ax.plot(plot_df.index, y, color=next(color_cycle), alpha=top_alpha, linewidth=top_lw,
                label=f"Lineage {col}")

    # burn-in shading & markers (single artists)
    if tracking_start_day > initial_infectors_day:
        ax.axvspan(initial_infectors_day, tracking_start_day - 1,
                   alpha=0.12, color="lightsteelblue", label="Tx burn-in")
    ax.axvline(initial_infectors_day, linestyle=":",  linewidth=1.5, color="tab:blue",  alpha=0.7, label="Seeding")
    ax.axvline(tracking_start_day,    linestyle="--", linewidth=1.5, color="tab:blue",  alpha=0.7, label="Burn-in over")

    ax.set_title("Lineage prevalence over time (active infections per lineage)")
    ax.set_xlabel("Day")
    ax.set_ylabel("# currently infected in lineage")

    # de-duplicate legend entries while preserving order
    handles, labels = ax.get_legend_handles_labels()
    seen = set()
    uniq_handles, uniq_labels = [], []
    for h, l in zip(handles, labels):
        if l.startswith('_'):  # matplotlib's "ignore" convention
            continue
        if l not in seen:
            uniq_handles.append(h); uniq_labels.append(l); seen.add(l)
    ax.legend(uniq_handles, uniq_labels, ncol=2, fontsize=8)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300)
        plt.close(fig)

# ========== DATA LOADING ==========

def load_simulation_data(output_dir):
    """Load simulation data from output directory"""
    print(f"Loading simulation data from: {output_dir}")
    
    # Load core files
    transmission_df = pd.read_csv(os.path.join(output_dir, "transmission_df.csv"))
    
    # Load initial infectors
    try:
        initial_infectors_df = pd.read_csv(os.path.join(output_dir, "initial_infectors.csv"))
        initial_infectors = list(zip(initial_infectors_df['node_id'], 
                                   initial_infectors_df['remaining_days_infectious']))
    except FileNotFoundError:
        print("Warning: initial_infectors.csv not found")
        initial_infectors = []
    
    # Get parameters
    try:
        with open(os.path.join(output_dir, "parameters_used.json"), 'r') as f:
            params = json.load(f)
        
        initial_infectors_day = params['simulation']['partnership_burnin_days']
        tracking_start_day = (params['simulation']['partnership_burnin_days'] + 
                            params['simulation']['transmission_burnin_days'])
        total_days = (params['simulation']['partnership_burnin_days'] + 
                     params['simulation']['transmission_burnin_days'] + 
                     params['simulation']['sim_days'])
    except FileNotFoundError:
        raise FileNotFoundError(f"parameters_used.json not found in {output_dir}. Cannot determine simulation timing parameters.")
    
    print(f"Loaded data:")
    print(f"  Transmissions: {len(transmission_df)}")
    print(f"  Initial infectors: {len(initial_infectors)}")
    print(f"  Initial infectors day: {initial_infectors_day}")
    print(f"  Tracking start day: {tracking_start_day}")
    print(f"  Total days: {total_days}")
    
    return transmission_df, initial_infectors, initial_infectors_day, tracking_start_day, total_days

def analyze_lineage_dominance(lineage_df, tracking_start_day, total_days, lineage_threshold=0.8, min_prevalence=10):
    """Analyze lineage dominance (from original Steps 2-3)"""
    
    print("\n🔍 Analyzing lineage dominance...")
    
    if lineage_df.empty:
        return {
            'success': False,
            'reason': 'No lineages found - simulation may have died out'
        }
    
    # Look at prevalence right after burn-in ends
    post_burnin_start = tracking_start_day
    post_burnin_end = min(tracking_start_day + 365, total_days - 1)  # First year post-burn-in
    
    if post_burnin_start >= len(lineage_df):
        return {
            'success': False, 
            'reason': f'Burn-in period ({tracking_start_day}) exceeds simulation length ({len(lineage_df)})'
        }
    
    # Average prevalence in first post-burn-in period
    post_burnin_prevalence = lineage_df.iloc[post_burnin_start:post_burnin_end].mean()
    post_burnin_prevalence = post_burnin_prevalence[post_burnin_prevalence > 0].sort_values(ascending=False)
    
    if len(post_burnin_prevalence) == 0:
        return {
            'success': False,
            'reason': 'No active lineages found in post-burn-in period - epidemic died out'
        }
    
    total_post_burnin_infections = post_burnin_prevalence.sum()
    dominant_lineage = post_burnin_prevalence.index[0]
    dominant_prevalence = post_burnin_prevalence.iloc[0]
    dominant_fraction = dominant_prevalence / total_post_burnin_infections
    
    print(f"📊 Post-burn-in analysis:")
    print(f"   Total active lineages: {len(post_burnin_prevalence)}")
    print(f"   Dominant lineage: {dominant_lineage}")
    print(f"   Dominant prevalence: {dominant_prevalence:.1f}")
    print(f"   Dominant fraction: {dominant_fraction:.3f}")
    print(f"   Threshold for single lineage: {lineage_threshold}")
    
    # Check if single lineage dominates
    single_lineage = (dominant_fraction >= lineage_threshold and dominant_prevalence >= min_prevalence)
    
    if single_lineage:
        print(f"✅ Single dominant lineage identified: {dominant_lineage} ({dominant_fraction:.3f} of infections)")
    else:
        print(f"❌ Multiple lineages still competing (dominant = {dominant_fraction:.3f}, threshold = {lineage_threshold})")
    
    return {
        'success': bool(True),
        'single_lineage': bool(single_lineage),
        'dominant_lineage': int(dominant_lineage),
        'dominant_fraction': float(dominant_fraction),
        'total_lineages': int(len(post_burnin_prevalence)),
        'post_burnin_prevalence': {str(k): float(v) for k, v in post_burnin_prevalence.to_dict().items()},
        'reason': str('Single lineage dominates' if single_lineage else f'Multiple lineages competing - dominant lineage {dominant_lineage} has only {dominant_fraction:.3f} of infections (threshold: {lineage_threshold})')
    }

# ========== MAIN SCRIPT ==========

def main():
    parser = argparse.ArgumentParser(description='Analyze lineage dominance in simulation outputs')
    parser.add_argument('output_dir', help='Directory containing simulation outputs')
    parser.add_argument('--threshold', type=float, default=1.0, 
                       help='Fraction threshold for single lineage dominance (default: 1.0)')
    parser.add_argument('--min-prevalence', type=float, default=10,
                       help='Minimum prevalence for viable lineage (default: 10)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        print(f"❌ Directory does not exist: {args.output_dir}")
        return 1
    
    try:
        print(f"🧬 Starting lineage assessment for: {args.output_dir}")
        
        # Step 1: Load simulation data
        transmission_df, initial_infectors, initial_infectors_day, tracking_start_day, total_days = \
            load_simulation_data(args.output_dir)
        
        # Step 2: Build lineage episodes
        print("\n🧬 Building lineage episodes...")
        episodes = build_lineage_episodes(
            initial_infectors=initial_infectors,
            initial_infectors_day=initial_infectors_day,
            transmission_df=transmission_df
        )
        print(f"Built {len(episodes)} lineage episodes")
        
        # Step 3: Create lineage prevalence timeseries
        print("\n📊 Creating lineage prevalence timeseries...")
        lineage_df = lineage_prevalence_timeseries(episodes, total_days=total_days)
        print(f"Created timeseries with {len(lineage_df.columns)} lineages over {len(lineage_df)} days")
        
        # Step 4: Create lineage plot
        print("\n📈 Creating lineage prevalence plot...")
        lineage_plot_path = os.path.join(args.output_dir, "lineage_analysis.png")
        plot_lineage_prevalence(
            lineage_df,
            initial_infectors_day=initial_infectors_day,
            tracking_start_day=tracking_start_day,
            top_k=10,
            smooth=7,
            save_path=lineage_plot_path
        )
        print(f"Lineage plot saved: {lineage_plot_path}")
        
        # Step 5: Analyze dominance
        dominance_results = analyze_lineage_dominance(
            lineage_df, tracking_start_day, total_days, 
            lineage_threshold=args.threshold,
            min_prevalence=args.min_prevalence
        )
        
        # Step 6: Save results
        results_path = os.path.join(args.output_dir, "lineage_test_results.json")
        with open(results_path, 'w') as f:
            json.dump(dominance_results, f, indent=2)
        
        # Step 7: Save episodes for tree building (if needed)
        episodes_path = os.path.join(args.output_dir, "lineage_episodes.json")
        with open(episodes_path, 'w') as f:
            json.dump(episodes, f, indent=2)
        
        print(f"\n📁 Results saved:")
        print(f"   - Lineage test: {results_path}")
        print(f"   - Episodes data: {episodes_path}")
        print(f"   - Lineage plot: {lineage_plot_path}")
        
        # Final result
        if dominance_results['success'] and dominance_results['single_lineage']:
            print(f"\n✅ SINGLE LINEAGE DETECTED: Lineage {dominance_results['dominant_lineage']}")
            print(f"   Dominance: {dominance_results['dominant_fraction']:.3f} (threshold: {args.threshold})")
        else:
            print(f"\n❌ MULTIPLE LINEAGES OR FAILURE")
            print(f"   Reason: {dominance_results['reason']}")
        
        return 0
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())