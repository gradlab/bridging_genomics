#!/usr/bin/env python3
"""
Core Simulation Runner - JSON Parameter Input
Runs simulation with parameters from JSON file
"""

import json
import argparse
import os
import sys
from datetime import datetime
from collections import defaultdict
from collections import Counter
import numpy as np
import pandas as pd

### Core simulation functions ###

# Generate distribution sampler (for exponential casual duration distribution)
def make_exponential_duration_sampler(rate=0.1, seed=None):
    """
    Returns a sampler for integer durations based on an Exponential(rate) distribution.
    - mean = 1/rate
    - durations are rounded up (ceil) to ensure >=1 day.
    
    Parameters
    ----------
    rate : float
        Rate λ > 0 (mean = 1/λ).
    seed : int or None
        Optional RNG seed.
    """
    if rate <= 0:
        raise ValueError("Rate λ must be positive.")
    
    rng = np.random.default_rng(seed)

    def sampler():
        x = rng.exponential(1.0 / rate)   # mean = 1/rate
        return max(1, int(np.ceil(x)))    # ceil to integer days, at least 1
    
    sampler.rate = rate
    sampler.mean_continuous = 1.0 / rate  # expected value of continuous distribution
    return sampler


# Node initialization and risk characteristic assignments
def initialize_nodes(N, behavior_proportions, seed):
    if abs(sum(behavior_proportions.values()) - 1.0) > 1e-6:
        raise ValueError("Behavior proportions must sum to 1")
    nodes = {}
    node_id = 1
    rng = np.random.default_rng(seed)
    group_counts = {k: int(N * v) for k, v in behavior_proportions.items()}
    assigned_total = sum(group_counts.values())
    if assigned_total != N:
        raise ValueError(f"Assigned {assigned_total} nodes but expected {N}. Adjust N or your behavior proportions.")
    for behavior, count in group_counts.items():
        for _ in range(count):
            nodes[node_id] = {"behavior": behavior}
            node_id += 1
    return nodes

def assign_risk_chars(
    nodes_dict,
    risk_dist_dict,
    seed,
    steady_capacity_mix,
):
    """
    Assign per-node risk characteristics with shifted gamma distributions.
    
    Now supports minimum partnership seeking rates via shifted gamma:
    value = gamma(shape, scale) + min_rate
    
    risk_dist_dict structure:
    {
        behavior: {
            "prop_hi": float,
            "lo_par": {"shape": float, "scale": float, "min_rate": float},
            "hi_par": {"shape": float, "scale": float, "min_rate": float}
        }
    }
    """

    def _strict_probs(mix_for_risk, risk_label):
        # [same validation logic as before]
        expected_keys = {0, 1, 2}
        missing = expected_keys - set(mix_for_risk.keys())
        extra   = set(mix_for_risk.keys()) - expected_keys
        if missing:
            raise ValueError(f"steady_capacity_mix['{risk_label}'] missing keys: {sorted(missing)}")
        if extra:
            raise ValueError(f"steady_capacity_mix['{risk_label}'] has unexpected keys: {sorted(extra)}")

        probs = np.array([float(mix_for_risk[k]) for k in (0, 1, 2)], dtype=float)

        if np.any((probs < 0) | (probs > 1)):
            raise ValueError(f"steady_capacity_mix['{risk_label}'] entries must be in [0,1]; got {probs.tolist()}")

        s = probs.sum()
        if abs(s - 1.0) > 1e-9:
            raise ValueError(
                f"steady_capacity_mix['{risk_label}'] probabilities must sum to 1; got sum={s}"
            )
        return probs

    # Validation
    for key in ("hi", "lo"):
        if key not in steady_capacity_mix:
            raise ValueError(f"steady_capacity_mix is missing key '{key}'")

    probs_hi = _strict_probs(steady_capacity_mix["hi"], "hi")
    probs_lo = _strict_probs(steady_capacity_mix["lo"], "lo")

    rng = np.random.default_rng(seed)
    states = np.array([0, 1, 2], dtype=int)

    for node_id, data in nodes_dict.items():
        behavior = data['behavior']
        risk_info = risk_dist_dict[behavior]

        # hi vs lo risk
        is_hi = rng.random() < risk_info['prop_hi']
        data['hi_risk'] = int(is_hi)

        # partner-change rate from SHIFTED gamma
        dist_par = risk_info['hi_par'] if is_hi else risk_info['lo_par']
        
        # Sample from gamma and shift by minimum rate
        raw_value = rng.gamma(dist_par['shape'], dist_par['scale'])
        value = raw_value + dist_par['min_rate']
        
        data['yearly_partner_change'] = value
        data['daily_partner_change'] = value / 365.0

        # steady capacity (unchanged)
        probs = probs_hi if is_hi else probs_lo
        data['max_steady_partners'] = int(rng.choice(states, p=probs))
        data['current_steady'] = 0

    return nodes_dict


# Function to get starting prevalence table from inputs 
def calculate_prevalences_and_create_table(initial_prevalences, risk_dist_dict, nodes_dict):
    """
    Calculate implied low-risk prevalences and create summary table.
    
    Parameters:
    - initial_prevalences: dict with format {"behavior": {"all": float, "hi": float}}
    - risk_dist_dict: dict with risk distribution info including prop_hi
    - nodes_dict: dict of nodes to get actual population counts
    
    Returns:
    - calculated_prevalences: dict in the old format {"behavior_risk": prevalence}
    - prevalence_table: pandas DataFrame with the breakdown
    """
    
    calculated_prevalences = {}
    table_rows = []
    
    # Get population counts by behavior and risk
    pop_counts = defaultdict(lambda: {"hi": 0, "lo": 0, "all": 0})
    for node_id, data in nodes_dict.items():
        behavior = data['behavior']
        risk_label = "hi" if data['hi_risk'] == 1 else "lo"
        pop_counts[behavior][risk_label] += 1
        pop_counts[behavior]["all"] += 1
    
    for behavior in initial_prevalences:
        P_all = initial_prevalences[behavior]["all"]
        P_hi = initial_prevalences[behavior]["hi"]
        
        # Get the proportion in hi-risk group for this behavior
        prop_hi = risk_dist_dict[behavior]["prop_hi"]
        prop_lo = 1 - prop_hi
        
        # Calculate implied low-risk prevalence
        # P_all = P_hi * prop_hi + P_lo * prop_lo
        # So: P_lo = (P_all - P_hi * prop_hi) / prop_lo
        
        if prop_lo > 0:
            P_lo = (P_all - P_hi * prop_hi) / prop_lo
        else:
            P_lo = 0  # Edge case: if everyone is high-risk
            
        # Ensure P_lo is non-negative
        P_lo = max(0, P_lo)
        
        # Store in old format for compatibility
        calculated_prevalences[f"{behavior}_hi"] = P_hi
        calculated_prevalences[f"{behavior}_lo"] = P_lo
        
        # Add to table
        table_rows.append({
            "behavior": behavior,
            "group": "all",
            "target_prevalence": P_all,
            "calculated_prevalence": P_all,  # This is what we specified
            "population": pop_counts[behavior]["all"],
            "expected_infected": P_all * pop_counts[behavior]["all"]
        })
        
        table_rows.append({
            "behavior": behavior,
            "group": "hi",
            "target_prevalence": P_hi,
            "calculated_prevalence": P_hi,  # This is what we specified
            "population": pop_counts[behavior]["hi"],
            "expected_infected": P_hi * pop_counts[behavior]["hi"]
        })
        
        table_rows.append({
            "behavior": behavior,
            "group": "lo", 
            "target_prevalence": None,  # We didn't specify this
            "calculated_prevalence": P_lo,  # This is what we calculated
            "population": pop_counts[behavior]["lo"],
            "expected_infected": P_lo * pop_counts[behavior]["lo"]
        })
    
    prevalence_table = pd.DataFrame(table_rows)
    
    return calculated_prevalences, prevalence_table

# Main partnership formation function
def form_daily_partnerships_with_unmatched_two_tier(
    nodes_dict,
    allowed_partners,
    p_MSMW_W,
    rng,
    casual_duration_sampler,
    steady_duration_sampler
):
    """
    Two-tier priority matching with consistent MSMW target selection.
    MSMW individuals roll once at the beginning to determine male vs female seeking.
    """
    df = pd.DataFrame.from_dict(nodes_dict, orient='index')

    # Who seeks today?
    seeking_mask = rng.random(len(df)) < df['daily_partner_change'].values
    seekers_df = df[seeking_mask].copy()
    if seekers_df.empty:
        return {}, {}, {}

    # PRE-DETERMINE MSMW targets for the entire day 
    msmw_targets = {}
    for node_id in seekers_df.index:
        if nodes_dict[node_id]['behavior'] == 'MSMW':
            seeks_female = rng.random() < p_MSMW_W
            msmw_targets[node_id] = 'female' if seeks_female else 'male'

    new_edges = {}
    all_matched = set()
    unmatched = set(seekers_df.index.tolist())

    def decide_casual(u, v):
        u_at_cap = nodes_dict[u]['current_steady'] >= nodes_dict[u]['max_steady_partners']
        v_at_cap = nodes_dict[v]['current_steady'] >= nodes_dict[v]['max_steady_partners']
        return 1 if (u_at_cap or v_at_cap) else 0

    def get_target_behaviors_consistent(node_id):
        """Get target behaviors using pre-determined MSMW choices"""
        behavior = nodes_dict[node_id]['behavior']
        if behavior == 'MSMW':
            # Use the pre-determined target instead of rolling again
            if msmw_targets[node_id] == 'female':
                return {'WSM'}
            else:
                return {'MSM', 'MSMW'}
        else:
            return allowed_partners[behavior]

    def try_match_group(seeker_ids):
        """Try to match a group of seekers"""
        seeker_list = list(seeker_ids)
        rng.shuffle(seeker_list)
        
        for i, u in enumerate(seeker_list):
            if u not in unmatched:
                continue
            
            u_behavior = nodes_dict[u]['behavior']
            acceptable_for_u = get_target_behaviors_consistent(u)
            
            candidates = []
            for v in unmatched:
                if v == u:
                    continue
                if nodes_dict[v]['behavior'] not in acceptable_for_u:
                    continue
                acceptable_for_v = get_target_behaviors_consistent(v)
                if u_behavior in acceptable_for_v:
                    candidates.append(v)
            
            if not candidates:
                # No one to match with - continue to next seeker
                continue
            else:
                # SELECT PARTNER AND CREATE EDGE 
                v = int(rng.choice(candidates))
                
                casual_flag = decide_casual(u, v)
                if casual_flag == 0:
                    duration = int(round(steady_duration_sampler()))
                else:
                    duration = int(round(casual_duration_sampler()))

                # Record and remove from pool
                unmatched.discard(u)
                unmatched.discard(v)
                all_matched.update([u, v])

                new_edges[(u, v)] = {
                    'casual': casual_flag,
                    'duration': max(1, duration)
                }

    # Separate seekers into priority groups
    high_priority = []  # MSM + MSMW seeking male partners (most constrained)
    low_priority = []   # Everyone else

    for u in unmatched:
        behavior = nodes_dict[u]['behavior']
        
        # HIGH PRIORITY: MSM (all) + MSMW seeking male partners
        if behavior == 'MSM':
            high_priority.append(u)
        elif behavior == 'MSMW' and msmw_targets.get(u) == 'male':
            high_priority.append(u)
        else:
            # LOW PRIORITY: WSM, MSW, MSMW seeking female
            low_priority.append(u)

    # Match high priority first (those seeking scarce MSM/MSMW)
    try_match_group(high_priority)
    
    # Then match everyone else
    try_match_group(low_priority)

    # Create unmatched diagnostics
    all_seekers = set(seekers_df.index.tolist())
    unmatched_ids = list(all_seekers - all_matched)
    unmatched_attrs = seekers_df.loc[unmatched_ids].to_dict(orient='index')
    seeker_attrs = seekers_df.to_dict(orient='index')

    return new_edges, unmatched_attrs, seeker_attrs


# Simulation Run Function

def run_partnership_and_infection_simulation(
    N,
    behavior_props,
    risk_dist_dict,
    allowed_partners,
    p_MSMW_W,
    rng,
    casual_duration_sampler,
    steady_duration_sampler,
    initial_prevalences,
    transmission_probs_per_partner_per_day,
    infxn_duration_params_dict,
    prop_init_sx,
    p_sx,
    partnership_burnin_days,
    transmission_burnin_days,
    sim_days,
    seed,
    steady_capacity_mix,
    optional_outputs=None,  # NEW: Set/list of optional outputs to track
):
    if steady_capacity_mix is None:
        raise ValueError("steady_capacity_mix must be provided.")
    
    # Handle optional outputs
    if optional_outputs is None:
        optional_outputs = set()
    optional_outputs = set(optional_outputs)  # Convert to set for fast lookup
    
    # Define what each optional output includes
    AVAILABLE_OPTIONAL_OUTPUTS = {
        'daily_infection_series',      # All daily infection tracking
        'partnership_metrics',         # yearly_partners, avg_concurrency  
        'unmatched_seeker_data',      # unmatched_df, seeker_df
        'infection_duration_records',  # infection_duration_records
        'edge_type_counts',           # edge_type_counts_df
        'prevalence_calculations',    # prevalence_table, calculated_prevalences
        'concurrency_tracking',       # concurrency_tracker
        'prevalence_time_series',     #track prevalences after burn-ins    
    }
    
    # Validate optional outputs
    invalid_outputs = optional_outputs - AVAILABLE_OPTIONAL_OUTPUTS
    if invalid_outputs:
        raise ValueError(f"Invalid optional outputs: {invalid_outputs}. Available: {AVAILABLE_OPTIONAL_OUTPUTS}")

    def _gamma_round(shape, scale):
        return max(1, int(round(rng.gamma(shape, scale))))

    def _draw_duration_from_params(params_dict, behavior, risk_label, status):
        par = params_dict[behavior][risk_label][status]
        return _gamma_round(par["shape"], par["scale"])

    SX, ASX, NAT = "symptomatic", "asymptomatic", "natural_clearance"

    # ---------- init nodes ----------
    nodes_dict = assign_risk_chars(
        nodes_dict=initialize_nodes(N, behavior_props, seed),
        risk_dist_dict=risk_dist_dict,
        seed=seed,
        steady_capacity_mix=steady_capacity_mix
    )
    for node_id in nodes_dict:
        nodes_dict[node_id].update({
            'currently_infected': 0,
            'infected_time_remaining': 0,
            'infection_day': -1,
            'infection_status': None,
            'duration_source': None
        })

    total_days = partnership_burnin_days + transmission_burnin_days + sim_days
    tracking_start = partnership_burnin_days + transmission_burnin_days

    # ---------- CORE OUTPUTS (always tracked) ----------
    active_edges = {}
    transmission_log = []
    all_edge_data = []
    initial_infectors = []

    # ---------- OPTIONAL OUTPUTS (only if requested) ----------
    
    # Daily infection series
    if 'daily_infection_series' in optional_outputs:
        new_infections_by_day = []
        current_infected_by_day = []
        new_infections_by_day_symptom = []
        new_infections_by_day_behavior = []
        new_infections_by_day_risk = []
        new_infections_by_day_behavior_symptom = []
        new_infections_by_day_behavior_risk = []
        new_infections_by_day_behavior_risk_symptom = []
        current_infected_by_day_symptom = []
        current_infected_by_day_behavior = []
        current_infected_by_day_risk = []
        current_infected_by_day_behavior_symptom = []
        current_infected_by_day_behavior_risk = []
        current_infected_by_day_behavior_risk_symptom = []

    # Partnership metrics (includes concurrency)
    if 'partnership_metrics' in optional_outputs:
        node_partner_counts = {node: 0 for node in nodes_dict}
        concurrency_tracker = {node: [0] * total_days for node in nodes_dict}  # Always needed for avg_concurrency

    # Concurrency tracking (extra detail beyond just averages)  
    if 'concurrency_tracking' in optional_outputs or 'partnership_metrics' in optional_outputs:
        if 'partnership_metrics' not in optional_outputs:
            concurrency_tracker = {node: [0] * total_days for node in nodes_dict}

    # Unmatched/seeker data
    if 'unmatched_seeker_data' in optional_outputs:
        unmatched_time_series = []
        seeker_time_series = []

    # Edge type counts
    if 'edge_type_counts' in optional_outputs:
        edge_type_counter = Counter()

    # Infection duration records
    if 'infection_duration_records' in optional_outputs:
        infection_duration_records = []

    # Prevalence calculations (will be set during seeding if requested)
    if 'prevalence_calculations' in optional_outputs:
        calculated_prevalences = None
        prevalence_table = None
    
    # Prevalence time series (every 100 steps post burn-in)
    if 'prevalence_time_series' in optional_outputs:
        prevalence_time_series = []

    # ---------- MAIN LOOP ----------
    for t in range(total_days):
        # === Seed infections at START of transmission burn-in ===
        if t == partnership_burnin_days:
            # Calculate prevalences
            calculated_prevalences_temp, prevalence_table_temp = calculate_prevalences_and_create_table(
                initial_prevalences, risk_dist_dict, nodes_dict
            )
            
            # Store if requested
            if 'prevalence_calculations' in optional_outputs:
                calculated_prevalences = calculated_prevalences_temp
                prevalence_table = prevalence_table_temp
            
            for node_id, data in nodes_dict.items():
                behavior = data['behavior']
                risk_label = "hi" if data['hi_risk'] == 1 else "lo"
                prevalence = float(calculated_prevalences_temp[f"{behavior}_{risk_label}"])
                
                if rng.random() < prevalence:
                    is_symptomatic = (rng.random() < float(prop_init_sx[behavior]))
                    status_choice = SX if is_symptomatic else ASX

                    dur_status = _draw_duration_from_params(infxn_duration_params_dict, behavior, risk_label, status_choice)
                    dur_nat = _draw_duration_from_params(infxn_duration_params_dict, behavior, risk_label, NAT)

                    if dur_nat < dur_status:
                        total_duration, duration_source = dur_nat, NAT
                    else:
                        total_duration, duration_source = dur_status, status_choice

                    elapsed_days = rng.integers(0, total_duration)
                    remaining_duration = max(1, total_duration - elapsed_days)

                    data.update({
                        'currently_infected': 1,
                        'infected_time_remaining': remaining_duration,
                        'infection_day': t - elapsed_days,
                        'infection_status': status_choice,
                        'duration_source': duration_source
                    })

                    # Track infection duration records if requested
                    if 'infection_duration_records' in optional_outputs:
                        infection_duration_records.append({
                            "day": t, "behavior": behavior, "risk": risk_label,
                            "symptom_status": status_choice,
                            "total_duration": total_duration,
                            "elapsed_duration": elapsed_days,
                            "remaining_duration": remaining_duration,
                            "duration_source": duration_source,
                            "seeded": True
                        })

            # SNAPSHOT initial infectors (CORE - always tracked)
            initial_infectors = [
                (node_id, d['infected_time_remaining'])
                for node_id, d in nodes_dict.items()
                if d['currently_infected'] == 1 and d['infected_time_remaining'] > 0
            ]

        # === Recovery step ===
        if t >= partnership_burnin_days:
            expired_nodes = []
            for node_id, d in nodes_dict.items():
                if d['currently_infected'] == 1:
                    d['infected_time_remaining'] -= 1
                    if d['infected_time_remaining'] <= 0:
                        expired_nodes.append(node_id)
            for node_id in expired_nodes:
                nodes_dict[node_id].update({
                    'currently_infected': 0,
                    'infected_time_remaining': 0,
                    'infection_day': -1,
                    'infection_status': None,
                    'duration_source': None
                })

        # === Dissolve expired partnerships ===
        expired_edges = []
        for edge, attr in active_edges.items():
            attr['duration'] -= 1
            if attr['duration'] <= 0:
                expired_edges.append(edge)
        for edge in expired_edges:
            n1, n2 = edge
            if active_edges[edge]['casual'] == 0:
                nodes_dict[n1]['current_steady'] = max(0, nodes_dict[n1]['current_steady'] - 1)
                nodes_dict[n2]['current_steady'] = max(0, nodes_dict[n2]['current_steady'] - 1)
            del active_edges[edge]

        # === Form new partnerships ===
        new_edges, unmatched_attrs, seeker_attrs = form_daily_partnerships_with_unmatched_two_tier(
            nodes_dict, allowed_partners, p_MSMW_W, rng,
            casual_duration_sampler, steady_duration_sampler
        )

        # Track unmatched/seeker data if requested
        if 'unmatched_seeker_data' in optional_outputs:
            unmatched_summary = {'day': t}
            for b in ['WSM', 'MSW', 'MSM', 'MSMW']:
                unmatched_summary[b] = sum(1 for u in unmatched_attrs.values() if u['behavior'] == b)
            unmatched_time_series.append(unmatched_summary)

            seeker_summary = {'day': t}
            for b in ['WSM', 'MSW', 'MSM', 'MSMW']:
                seeker_summary[f'{b}_seekers'] = sum(1 for s in seeker_attrs.values() if s['behavior'] == b)
                for r in ['hi', 'lo']:
                    risk_key = f'{b}_{r}_seekers'
                    seeker_summary[risk_key] = sum(1 for s in seeker_attrs.values()
                                                if s['behavior'] == b and
                                                   (s['hi_risk'] == 1 if r == 'hi' else s['hi_risk'] == 0))
            seeker_time_series.append(seeker_summary)

        # Add new edges and record partnership metrics
        for (n1, n2), attr in new_edges.items():
            active_edges[(n1, n2)] = attr
            if attr['casual'] == 0:
                nodes_dict[n1]['current_steady'] += 1
                nodes_dict[n2]['current_steady'] += 1

            if t >= partnership_burnin_days:
                # Partnership metrics (optional)
                if 'partnership_metrics' in optional_outputs:
                    node_partner_counts[n1] += 1
                    node_partner_counts[n2] += 1

                # Edge data (CORE - always tracked)
                b1 = nodes_dict[n1]['behavior']; b2 = nodes_dict[n2]['behavior']
                r1 = 'hi' if nodes_dict[n1]['hi_risk'] else 'lo'
                r2 = 'hi' if nodes_dict[n2]['hi_risk'] else 'lo'
                edge_type = '-'.join(sorted([b1, b2]))
                br_pairs = sorted([(b1, r1), (b2, r2)], key=lambda x: x[0])
                edge_type_byrisk = f"{br_pairs[0][0]}_{br_pairs[0][1]}-{br_pairs[1][0]}_{br_pairs[1][1]}"
                
                # Edge type counter (optional)
                if 'edge_type_counts' in optional_outputs:
                    edge_type_counter[edge_type] += 1

                hi1 = nodes_dict[n1]['hi_risk']; hi2 = nodes_dict[n2]['hi_risk']
                risk_level_any_hi = 'high' if (hi1 or hi2) else 'low'

                all_edge_data.append({
                    'day': t, 'node1': n1, 'node2': n2,
                    'behavior1': b1, 'behavior2': b2,
                    'risk1': r1, 'risk2': r2,
                    'edge_type': edge_type,
                    'edge_type_byrisk': edge_type_byrisk,
                    'casual': attr['casual'],
                    'duration': attr['duration'],
                    'risk_level': risk_level_any_hi
                })

                # Concurrency tracking (optional)
                if 'concurrency_tracking' in optional_outputs:
                    end_day = min(t + attr['duration'], total_days)
                    for day in range(t, end_day):
                        concurrency_tracker[n1][day] += 1
                        concurrency_tracker[n2][day] += 1

        # === Transmission step ===
        if t >= partnership_burnin_days:
            phase = "burnin_tx" if t < tracking_start else "post_tx"
            candidates_by_infectee = defaultdict(list)

            for (n1, n2), attr in active_edges.items():
                d1 = nodes_dict[n1]; d2 = nodes_dict[n2]
                inf1 = d1['currently_infected'] == 1 and d1['infection_day'] < t
                inf2 = d2['currently_infected'] == 1 and d2['infection_day'] < t
                if inf1 == inf2:
                    continue

                infector, infectee = (n1, n2) if inf1 else (n2, n1)

                if nodes_dict[infectee]['currently_infected'] == 1:
                    continue

                b_inf = nodes_dict[infector]['behavior']
                b_sus = nodes_dict[infectee]['behavior']
                partnership_category = 'casual' if attr['casual'] == 1 else 'steady'
                pair_type = f"{b_inf}-{b_sus}, {partnership_category}"
                prob = float(transmission_probs_per_partner_per_day[pair_type])

                rv = rng.random()
                did_tx = (rv < prob)
                if did_tx:
                    is_symptomatic = (rng.random() < float(p_sx[b_sus]))
                    status_choice = SX if is_symptomatic else ASX
                    risk_label = "hi" if nodes_dict[infectee]['hi_risk'] == 1 else "lo"

                    dur_status = _draw_duration_from_params(infxn_duration_params_dict, b_sus, risk_label, status_choice)
                    dur_nat = _draw_duration_from_params(infxn_duration_params_dict, b_sus, risk_label, NAT)
                    duration, duration_source = (dur_nat, NAT) if dur_nat < dur_status else (dur_status, status_choice)

                    candidates_by_infectee[infectee].append({
                        "infector": infector,
                        "duration": duration,
                        "duration_source": duration_source,
                        "status_choice": status_choice,
                        "pair_type": pair_type,
                        "partnership_category": partnership_category,
                        "behavior_infector": b_inf,
                        "behavior_infectee": b_sus,
                        "phase": phase
                    })

            # Process transmission candidates
            winners_today, losers_today = [], []
            for infectee, cand_list in candidates_by_infectee.items():
                win_idx = rng.integers(len(cand_list))
                for j, cand in enumerate(cand_list):
                    if j == win_idx:
                        nodes_dict[infectee].update({
                            'currently_infected': 1,
                            'infected_time_remaining': cand["duration"],
                            'infection_day': t,
                            'infection_status': cand["status_choice"],
                            'duration_source': cand["duration_source"]
                        })
                        winners_today.append((infectee, cand))

                        # Track infection duration records if requested
                        if 'infection_duration_records' in optional_outputs:
                            risk_label = "hi" if nodes_dict[infectee]['hi_risk'] == 1 else "lo"
                            infection_duration_records.append({
                                "day": t, "behavior": nodes_dict[infectee]['behavior'], "risk": risk_label,
                                "symptom_status": cand["status_choice"], "duration": cand["duration"],
                                "duration_source": cand["duration_source"], "seeded": False
                            })
                    else:
                        losers_today.append((infectee, cand))

            # Log transmission events (CORE - always tracked)
            for infectee, cand in winners_today:
                b_inf = cand["behavior_infector"]
                b_sus = cand["behavior_infectee"]
                part_cat = cand["partnership_category"]
                status_choice = cand["status_choice"]
                duration_source = cand["duration_source"]
                duration = cand["duration"]
                ph = cand["phase"]

                transmission_log.append({
                    'transmission_id': len(transmission_log) + 1,
                    'day_of_transmission': t,
                    'phase': ph,
                    'infector_node': cand["infector"],
                    'infectee_node': infectee,
                    'behavior_infector': b_inf,
                    'behavior_infectee': b_sus,
                    'partnership_type': cand["pair_type"],
                    'day_of_sampling': t + duration,
                    'partnership_category': part_cat,
                    'behavior_pair': cand["pair_type"].replace(",", "_"),
                    'superseded_simultaneous': False,
                    'infection_symptom_status': status_choice,
                    'duration_source': duration_source
                })

            for infectee, cand in losers_today:
                transmission_log.append({
                    'transmission_id': len(transmission_log) + 1,
                    'day_of_transmission': t,
                    'phase': cand["phase"],
                    'infector_node': cand["infector"],
                    'infectee_node': infectee,
                    'behavior_infector': cand["behavior_infector"],
                    'behavior_infectee': cand["behavior_infectee"],
                    'partnership_type': cand["pair_type"],
                    'day_of_sampling': None,
                    'partnership_category': cand["partnership_category"],
                    'behavior_pair': cand["pair_type"].replace(",", "_"),
                    'superseded_simultaneous': True,
                    'infection_symptom_status': None,
                    'duration_source': None
                })
            # === Prevalence tracking (every 100 steps post burn-in) ===
            if ('prevalence_time_series' in optional_outputs and 
                t >= tracking_start and 
                (t - tracking_start) % 100 == 0):
                
                prevalence_snapshot = {'day': t, 'day_post_burnin': t - tracking_start}
                
                # Calculate prevalence for each behavior group
                for behavior in ['WSM', 'MSW', 'MSM', 'MSMW']:
                    # Overall prevalence for this behavior
                    total_behavior = sum(1 for node_data in nodes_dict.values() 
                                    if node_data['behavior'] == behavior)
                    infected_behavior = sum(1 for node_data in nodes_dict.values() 
                                        if node_data['behavior'] == behavior and 
                                            node_data['currently_infected'] == 1)
                    prevalence_snapshot[f'{behavior}_prevalence'] = (
                        infected_behavior / total_behavior if total_behavior > 0 else 0.0
                    )
                    
                    # High activity prevalence for this behavior
                    total_behavior_hi = sum(1 for node_data in nodes_dict.values() 
                                        if node_data['behavior'] == behavior and 
                                            node_data['hi_risk'] == 1)
                    infected_behavior_hi = sum(1 for node_data in nodes_dict.values() 
                                            if node_data['behavior'] == behavior and 
                                                node_data['hi_risk'] == 1 and 
                                                node_data['currently_infected'] == 1)
                    prevalence_snapshot[f'{behavior}_hi_prevalence'] = (
                        infected_behavior_hi / total_behavior_hi if total_behavior_hi > 0 else 0.0
                    )
                
                prevalence_time_series.append(prevalence_snapshot)


            # Daily infection series tracking (optional)
            if 'daily_infection_series' in optional_outputs:
                # Build daily tallies
                new_group_counts = defaultdict(int)
                new_symptom_counts = defaultdict(int)
                new_behavior_counts = defaultdict(int)
                new_risk_counts = defaultdict(int)
                new_behavior_symptom = defaultdict(int)
                new_behavior_risk = defaultdict(int)
                new_behavior_risk_symptom = defaultdict(int)

                for infectee, cand in winners_today:
                    b_sus = cand["behavior_infectee"]
                    risk_label = "hi" if nodes_dict[infectee]['hi_risk'] == 1 else "lo"
                    status_choice = cand["status_choice"]

                    new_group_counts[f"{b_sus}_{risk_label}"] += 1
                    new_symptom_counts[status_choice] += 1
                    new_behavior_counts[b_sus] += 1
                    new_risk_counts[risk_label] += 1
                    new_behavior_symptom[f"{b_sus}_{status_choice}"] += 1
                    new_behavior_risk[f"{b_sus}_{risk_label}"] += 1
                    new_behavior_risk_symptom[f"{b_sus}_{risk_label}_{status_choice}"] += 1

                # Current infected tallies
                cur_group = defaultdict(int)
                cur_symptom = defaultdict(int)
                cur_behavior = defaultdict(int)
                cur_risk = defaultdict(int)
                cur_behavior_symptom = defaultdict(int)
                cur_behavior_risk = defaultdict(int)
                cur_behavior_risk_symptom = defaultdict(int)

                for node_id, d in nodes_dict.items():
                    if d['currently_infected'] == 1:
                        behavior = d['behavior']
                        risk_label = "hi" if d['hi_risk'] == 1 else "lo"
                        status_choice = d['infection_status'] if d['infection_status'] in (SX, ASX) else ASX
                        cur_group[f"{behavior}_{risk_label}"] += 1
                        cur_symptom[status_choice] += 1
                        cur_behavior[behavior] += 1
                        cur_risk[risk_label] += 1
                        cur_behavior_symptom[f"{behavior}_{status_choice}"] += 1
                        cur_behavior_risk[f"{behavior}_{risk_label}"] += 1
                        cur_behavior_risk_symptom[f"{behavior}_{risk_label}_{status_choice}"] += 1

                def _with_meta(dct):
                    dct = dict(dct)
                    dct["day"] = t
                    dct["phase"] = phase
                    return dct

                current_infected_by_day.append(_with_meta(cur_group))
                new_infections_by_day.append(_with_meta(new_group_counts))
                new_infections_by_day_symptom.append(_with_meta(new_symptom_counts))
                new_infections_by_day_behavior.append(_with_meta(new_behavior_counts))
                new_infections_by_day_risk.append(_with_meta(new_risk_counts))
                new_infections_by_day_behavior_symptom.append(_with_meta(new_behavior_symptom))
                new_infections_by_day_behavior_risk.append(_with_meta(new_behavior_risk))
                new_infections_by_day_behavior_risk_symptom.append(_with_meta(new_behavior_risk_symptom))
                current_infected_by_day_symptom.append(_with_meta(cur_symptom))
                current_infected_by_day_behavior.append(_with_meta(cur_behavior))
                current_infected_by_day_risk.append(_with_meta(cur_risk))
                current_infected_by_day_behavior_symptom.append(_with_meta(cur_behavior_symptom))
                current_infected_by_day_behavior_risk.append(_with_meta(cur_behavior_risk))
                current_infected_by_day_behavior_risk_symptom.append(_with_meta(cur_behavior_risk_symptom))

    # ---------- BUILD RETURN DICTIONARY ----------
    results = {
        # CORE OUTPUTS (always included)
        "nodes_dict": nodes_dict,
        "edge_df": pd.DataFrame(all_edge_data),
        "transmission_df": pd.DataFrame(transmission_log),
        "initial_infectors": initial_infectors,
        "initial_infectors_day": partnership_burnin_days,
        
        # Simulation parameters (always useful)
        "partnership_burnin_days": partnership_burnin_days,
        "transmission_burnin_days": transmission_burnin_days,
        "tracking_start_day": tracking_start,
        "sim_days": sim_days,
    }

    # Add optional outputs
    if 'daily_infection_series' in optional_outputs:
        results.update({
            "new_infections_by_day": new_infections_by_day,
            "current_infected_by_day": current_infected_by_day,
            "new_infections_by_day_symptom": new_infections_by_day_symptom,
            "new_infections_by_day_behavior": new_infections_by_day_behavior,
            "new_infections_by_day_risk": new_infections_by_day_risk,
            "new_infections_by_day_behavior_symptom": new_infections_by_day_behavior_symptom,
            "new_infections_by_day_behavior_risk": new_infections_by_day_behavior_risk,
            "new_infections_by_day_behavior_risk_symptom": new_infections_by_day_behavior_risk_symptom,
            "current_infected_by_day_symptom": current_infected_by_day_symptom,
            "current_infected_by_day_behavior": current_infected_by_day_behavior,
            "current_infected_by_day_risk": current_infected_by_day_risk,
            "current_infected_by_day_behavior_symptom": current_infected_by_day_behavior_symptom,
            "current_infected_by_day_behavior_risk": current_infected_by_day_behavior_risk,
            "current_infected_by_day_behavior_risk_symptom": current_infected_by_day_behavior_risk_symptom,
        })

    if 'partnership_metrics' in optional_outputs:
        results["yearly_partners"] = pd.Series(node_partner_counts)
        avg_concurrency = {
            node: np.mean(concurrency_tracker[node][partnership_burnin_days:])
            for node in nodes_dict
        }
        results["avg_concurrency"] = pd.Series(avg_concurrency)

    if 'unmatched_seeker_data' in optional_outputs:
        results.update({
            "unmatched_df": pd.DataFrame(unmatched_time_series),
            "seeker_df": pd.DataFrame(seeker_time_series),
        })

    if 'infection_duration_records' in optional_outputs:
        results["infection_duration_records"] = pd.DataFrame(infection_duration_records)

    if 'edge_type_counts' in optional_outputs:
        if len(edge_type_counter) > 0:
            edge_type_counts_df = (
                pd.DataFrame([{"edge_type": et, "count": c} for et, c in edge_type_counter.items()])
                .sort_values("count", ascending=False)
                .reset_index(drop=True)
            )
        else:
            edge_type_counts_df = pd.DataFrame(columns=["edge_type", "count"])
        results["edge_type_counts_df"] = edge_type_counts_df

    if 'prevalence_calculations' in optional_outputs:
        results.update({
            "prevalence_table": prevalence_table,
            "calculated_prevalences": calculated_prevalences,
        })

    if 'concurrency_tracking' in optional_outputs:
        results["concurrency_tracker"] = concurrency_tracker
    
    if 'prevalence_time_series' in optional_outputs:
        prevalence_df = pd.DataFrame(prevalence_time_series)
        
        # Calculate average prevalences AND standard errors in DataFrame format
        avg_prevalences_data = []
        for behavior in ['WSM', 'MSW', 'MSM', 'MSMW']:
            # Overall behavior group
            if f'{behavior}_prevalence' in prevalence_df.columns:
                values = prevalence_df[f'{behavior}_prevalence']
                avg_prevalences_data.append({
                    'group': f'{behavior}_prevalence',
                    'mean': values.mean(),
                    'se': values.std() / np.sqrt(len(values)) if len(values) > 1 else 0.0
                })
            
            # High activity group  
            if f'{behavior}_hi_prevalence' in prevalence_df.columns:
                values_hi = prevalence_df[f'{behavior}_hi_prevalence']
                avg_prevalences_data.append({
                    'group': f'{behavior}_hi_prevalence', 
                    'mean': values_hi.mean(),
                    'se': values_hi.std() / np.sqrt(len(values_hi)) if len(values_hi) > 1 else 0.0
                })
        
        results["prevalence_time_series"] = prevalence_df
        results["average_prevalences"] = pd.DataFrame(avg_prevalences_data).set_index('group')

    return results


### Load in Parameters and Run Simulation ###


def load_parameters(json_path):
    """Load simulation parameters from JSON file"""
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"Parameter file not found: {json_path}")
    
    with open(json_path, 'r') as f:
        params = json.load(f)
    
    return params

def create_output_directory(params):
    """Create output directory"""
    base_dir = params['output']['output_dir']
    run_name = params['output']['run_name']
    
    output_dir = os.path.join(base_dir, run_name)
    os.makedirs(output_dir, exist_ok=True)
    
    return output_dir

def run_simulation_from_json(json_path):
    """Main function to run simulation from JSON parameters"""
    
    print(f"Loading parameters from: {json_path}")
    params = load_parameters(json_path)
    
    # Create output directory
    output_dir = create_output_directory(params)
    print(f"Output directory: {output_dir}")
    
    # Save parameters to output directory for reproducibility
    with open(os.path.join(output_dir, "parameters_used.json"), 'w') as f:
        json.dump(params, f, indent=2)
    
    # Set up RNG
    seed = params['simulation']['seed']
    rng = np.random.default_rng(seed)
    
    # Create duration samplers
    casual_duration_sampler = make_exponential_duration_sampler(
        rate=params['partnerships']['casual_rate'], 
        seed=seed
    )
    steady_duration_sampler = lambda: rng.gamma(
        shape=params['partnerships']['steady_shape'],
        scale=params['partnerships']['steady_scale']
    )
    
    # Convert JSON parameters to expected format
    behavior_props = params['behaviors']
    risk_dist_dict = params['risk_distribution']
    allowed_partners = {k: set(v) for k, v in params['partnerships']['allowed_partners'].items()}
    
    # Extract infection parameters
    initial_prevalences = params['infections']['initial_prevalences']
    transmission_probs_per_partner_per_day = params['infections']['transmission_probs_per_partner_per_day']
    infxn_duration_params_dict = params['infections']['infxn_duration_params_dict']
    prop_init_sx = params['infections']['prop_init_sx']
    p_sx = params['infections']['p_sx']
    
    # Extract optional outputs (default to empty list if not specified)
    optional_outputs = params.get('optional_outputs', [])
    
    # Convert steady capacity mix
    steady_capacity_mix = {}
    for risk_level, mix in params['partnerships']['steady_capacity_mix'].items():
        steady_capacity_mix[risk_level] = {int(k): v for k, v in mix.items()}
    
    print("Starting simulation...")
    print(f"Optional outputs requested: {optional_outputs}")
    
    # Run simulation with optional outputs
    sim_results = run_partnership_and_infection_simulation(
        N=params['simulation']['N'],
        behavior_props=behavior_props,
        risk_dist_dict=risk_dist_dict,
        allowed_partners=allowed_partners,
        p_MSMW_W=params['partnerships']['p_MSMW_W'],
        rng=rng,
        casual_duration_sampler=casual_duration_sampler,
        steady_duration_sampler=steady_duration_sampler,
        initial_prevalences=initial_prevalences,
        transmission_probs_per_partner_per_day=transmission_probs_per_partner_per_day,
        infxn_duration_params_dict=infxn_duration_params_dict,
        prop_init_sx=prop_init_sx,
        p_sx=p_sx,
        partnership_burnin_days=params['simulation']['partnership_burnin_days'],
        transmission_burnin_days=params['simulation']['transmission_burnin_days'],
        sim_days=params['simulation']['sim_days'],
        seed=seed,
        steady_capacity_mix=steady_capacity_mix,
        optional_outputs=optional_outputs  # PASS THE OPTIONAL OUTPUTS
    )
    
    print("Simulation complete! Saving outputs...")
    
    # SAVE CORE OUTPUTS (always available)
    nodes_df = pd.DataFrame.from_dict(sim_results["nodes_dict"], orient='index')
    nodes_df.to_csv(os.path.join(output_dir, "nodes_df.csv"), index=True)
    sim_results["edge_df"].to_csv(os.path.join(output_dir, "edge_df.csv"), index=False)
    sim_results["transmission_df"].to_csv(os.path.join(output_dir, "transmission_df.csv"), index=False)
    initial_infectors_df = pd.DataFrame(sim_results["initial_infectors"], columns=["node_id", "remaining_days_infectious"])
    initial_infectors_df.to_csv(os.path.join(output_dir, "initial_infectors.csv"), index=False)
    
    # SAVE OPTIONAL OUTPUTS (only if they exist)
    
    # Unmatched/seeker data
    if "unmatched_df" in sim_results:
        sim_results["unmatched_df"].to_csv(os.path.join(output_dir, "unmatched_df.csv"), index=False)
    if "seeker_df" in sim_results:
        sim_results["seeker_df"].to_csv(os.path.join(output_dir, "seeker_df.csv"), index=False)
    
    # Partnership metrics  
    if "yearly_partners" in sim_results:
        sim_results["yearly_partners"].to_csv(os.path.join(output_dir, "yearly_partners.csv"))
    if "avg_concurrency" in sim_results:
        sim_results["avg_concurrency"].to_csv(os.path.join(output_dir, "avg_concurrency.csv"))
    
    # Daily infection series
    if "new_infections_by_day" in sim_results:
        pd.DataFrame(sim_results["new_infections_by_day"]).to_csv(os.path.join(output_dir, "new_infections_by_day.csv"), index=False)
    if "current_infected_by_day" in sim_results:
        pd.DataFrame(sim_results["current_infected_by_day"]).to_csv(os.path.join(output_dir, "current_infected_by_day.csv"), index=False)
    if "new_infections_by_day_behavior" in sim_results:
        pd.DataFrame(sim_results["new_infections_by_day_behavior"]).to_csv(os.path.join(output_dir, "new_infections_by_day_behavior.csv"), index=False)
    if "current_infected_by_day_behavior" in sim_results:
        pd.DataFrame(sim_results["current_infected_by_day_behavior"]).to_csv(os.path.join(output_dir, "current_infected_by_day_behavior.csv"), index=False)
    
    # Infection duration records
    if "infection_duration_records" in sim_results:
        sim_results["infection_duration_records"].to_csv(os.path.join(output_dir, "infection_duration_records.csv"), index=False)
    
    # Edge type counts
    if "edge_type_counts_df" in sim_results:
        sim_results["edge_type_counts_df"].to_csv(os.path.join(output_dir, "edge_type_counts.csv"), index=False)
    
    # Prevalence calculations
    if "prevalence_table" in sim_results:
        sim_results["prevalence_table"].to_csv(os.path.join(output_dir, "prevalence_table.csv"), index=False)
    if "calculated_prevalences" in sim_results:
        pd.Series(sim_results["calculated_prevalences"]).to_csv(os.path.join(output_dir, "calculated_prevalences.csv"))
    
    # Concurrency tracker (large file - save as pickle)
    if "concurrency_tracker" in sim_results:
        import pickle
        with open(os.path.join(output_dir, "concurrency_tracker.pkl"), 'wb') as f:
            pickle.dump(sim_results["concurrency_tracker"], f)
    
    # Prevalence time series
    if "prevalence_time_series" in sim_results:
        sim_results["prevalence_time_series"].to_csv(os.path.join(output_dir, "prevalence_time_series.csv"), index=False)
    if "average_prevalences" in sim_results:
        sim_results["average_prevalences"].to_csv(os.path.join(output_dir, "average_prevalences.csv"))
    
    # Save metadata
    metadata = {
        'simulation_complete': True,
        'total_nodes': len(nodes_df),
        'total_partnerships': len(sim_results["edge_df"]),
        'total_transmissions': len(sim_results["transmission_df"]),
        'optional_outputs_used': optional_outputs,
        'output_directory': output_dir,
        'parameter_file': json_path
    }
    
    with open(os.path.join(output_dir, "metadata.json"), 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"✅ Simulation complete! Results saved to: {output_dir}")
    return output_dir


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(description='Run epidemic simulation from JSON parameters')
    parser.add_argument('json_file', help='Path to JSON parameter file')
    parser.add_argument('--validate-only', action='store_true', 
                       help='Only validate parameters without running simulation')
    
    args = parser.parse_args()
    
    if args.validate_only:
        try:
            params = load_parameters(args.json_file)
            print("✅ Parameter file is valid!")
            return
        except Exception as e:
            print(f"❌ Parameter file validation failed: {e}")
            sys.exit(1)
    
    try:
        output_dir = run_simulation_from_json(args.json_file)
        print(f"Success! Results in: {output_dir}")
    except Exception as e:
        print(f"❌ Simulation failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
