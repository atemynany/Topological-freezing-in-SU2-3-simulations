#!/usr/bin/env python3
"""
==============================================================================
analyze_topcharge.py
==============================================================================
Statistical analysis of topological charge measurements.
Computes mean, standard deviation, and jackknife errors with bias estimation.

Usage:
    python3 analysis/analyze_topcharge.py

Author: Alexander de Barros Noll
Date: January 2026
==============================================================================
"""

import os
import numpy as np

from Q_estimator import QEstimator


# ==============================================================================
# Paths
# ==============================================================================

base_dir = "/home/alex/Desktop/workspace/master_thesis/lattice_qcd_topolgical_charge/Topological_charge_my_implementation"
output_dir = os.path.join(base_dir, "output")
topcharge_file = os.path.join(output_dir, "topcharge.dat")


# ==============================================================================
# Data Loading
# ==============================================================================

def load_topcharge_data(filename: str) -> dict:
    """Load topological charge data from output file."""
    data = {
        'smear_steps': [],
        'config': [],
        'Q': [],
        'plaquette': []
    }
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) >= 4:
                try:
                    data['smear_steps'].append(int(parts[0]))
                    data['config'].append(int(parts[1]))
                    data['Q'].append(float(parts[2]))
                    data['plaquette'].append(float(parts[3]))
                except ValueError:
                    continue
    
    for key in data:
        data[key] = np.array(data[key])
    
    return data


# ==============================================================================
# Jackknife Analysis
# ==============================================================================

def jackknife_analysis(values: np.ndarray, func=np.mean) -> dict:
    """
    Perform jackknife resampling to estimate error and bias.
    
    The jackknife method:
    1. Compute the estimator on the full sample: theta_full = func(values)
    2. For each i, compute theta_i = func(values with i-th element removed)
    3. Jackknife estimate: theta_jack = n * theta_full - (n-1) * mean(theta_i)
    4. Jackknife error: sigma_jack = sqrt((n-1)/n * sum((theta_i - mean(theta_i))^2))
    5. Bias estimate: bias = (n-1) * (mean(theta_i) - theta_full)
    
    Parameters:
        values: array of measurements
        func: function to apply (default: np.mean)
    
    Returns:
        dict with: estimate, error, bias, jackknife_samples
    """
    n = len(values)
    
    # Full sample estimate
    theta_full = func(values)
    
    # Jackknife resamples (leave-one-out)
    theta_jack = np.zeros(n)
    for i in range(n):
        # Create sample with i-th element removed
        sample_i = np.delete(values, i)
        theta_jack[i] = func(sample_i)
    
    # Mean of jackknife estimates
    theta_jack_mean = np.mean(theta_jack)
    
    # Jackknife error estimate
    error = np.sqrt((n - 1) / n * np.sum((theta_jack - theta_jack_mean)**2))
    
    # Bias estimate
    bias = (n - 1) * (theta_jack_mean - theta_full)
    
    # Bias-corrected estimate
    estimate_corrected = theta_full - bias
    
    return {
        'estimate': theta_full,
        'estimate_corrected': estimate_corrected,
        'error': error,
        'bias': bias,
        'jackknife_samples': theta_jack
    }


def jackknife_Q2(values: np.ndarray) -> dict:
    """Jackknife analysis for <Q^2> (topological susceptibility numerator)."""
    return jackknife_analysis(values, func=lambda x: np.mean(x**2))


# ==============================================================================
# Main Analysis
# ==============================================================================

def main():
    print("="*60)
    print("Topological Charge Analysis")
    print("="*60)
    
    # Load data
    print(f"\nLoading data from: {topcharge_file}")
    data = load_topcharge_data(topcharge_file)
    
    if len(data['Q']) == 0:
        print("Error: No data found")
        return
    
    # Get unique smearing levels
    smear_levels = np.unique(data['smear_steps'])
    print(f"Found smearing levels: {smear_levels}")
    
    # Analyze each smearing level
    print("\n" + "-"*60)
    print(f"{'Smear':>6} {'N':>5} {'<Q>':>10} {'std':>10} {'JK error':>10} {'JK bias':>10}")
    print("-"*60)
    
    for smear in smear_levels:
        mask = data['smear_steps'] == smear
        Q_values = data['Q'][mask]
        n = len(Q_values)
        
        if n < 2:
            continue
        
        # Basic statistics
        Q_mean = np.mean(Q_values)
        Q_std = np.std(Q_values, ddof=1)
        
        # Jackknife analysis for mean
        jk_mean = jackknife_analysis(Q_values, func=np.mean)
        
        print(f"{smear:>6} {n:>5} {Q_mean:>10.4f} {Q_std:>10.4f} "
              f"{jk_mean['error']:>10.4f} {jk_mean['bias']:>10.6f}")
    
    print("-"*60)
    
    # Detailed analysis for highest smearing level
    max_smear = smear_levels[-1]
    mask = data['smear_steps'] == max_smear
    Q_values = data['Q'][mask]
    
    print(f"\n\nDetailed analysis at smearing = {max_smear}:")
    print("="*60)
    
    # Mean
    jk_mean = jackknife_analysis(Q_values, func=np.mean)
    print(f"\n<Q>:")
    print(f"  Estimate:          {jk_mean['estimate']:.6f}")
    print(f"  Jackknife error:   {jk_mean['error']:.6f}")
    print(f"  Jackknife bias:    {jk_mean['bias']:.8f}")
    print(f"  Bias-corrected:    {jk_mean['estimate_corrected']:.6f}")
    
    # <Q^2> for susceptibility
    jk_Q2 = jackknife_Q2(Q_values)
    print(f"\n<Q²>:")
    print(f"  Estimate:          {jk_Q2['estimate']:.6f}")
    print(f"  Jackknife error:   {jk_Q2['error']:.6f}")
    print(f"  Jackknife bias:    {jk_Q2['bias']:.8f}")
    print(f"  Bias-corrected:    {jk_Q2['estimate_corrected']:.6f}")
    
    # Standard deviation via jackknife
    jk_std = jackknife_analysis(Q_values, func=lambda x: np.std(x, ddof=1))
    print(f"\nstd(Q):")
    print(f"  Estimate:          {jk_std['estimate']:.6f}")
    print(f"  Jackknife error:   {jk_std['error']:.6f}")
    print(f"  Jackknife bias:    {jk_std['bias']:.8f}")
    
    # ==================================================================
    # Q Estimator Analysis: Rescaled charge with optimal rescaling
    # ==================================================================
    print("\n" + "="*60)
    print("Q Estimator: Rescaled Charge with Optimal Rescaling")
    print("="*60)
    print("\nMethod: Q_rescaled = round(α * Q̂)")
    print("α minimizes ⟨(α * Q̂ - Q_rescaled)²⟩")
    
    # Compute Q estimator using the module
    est = QEstimator(Q_values)
    
    print(f"\nOptimal rescaling factor:")
    print(f"  α = {est.alpha:.6f}")
    
    print(f"\nDeviation from integer:")
    print(f"  ⟨(αQ̂ - Q_rescaled)²⟩ = {est.result.mean_deviation_squared:.6f}")
    print(f"  RMS deviation = {est.result.rms_deviation:.6f}")
    
    print(f"\nRescaled charge statistics:")
    print(f"  <Q_rescaled>  = {est.result.Q_mean:.4f}")
    print(f"  std           = {est.result.Q_std:.4f}")
    print(f"  <Q_rescaled²> = {est.result.Q2_mean:.4f}")
    
    # Distribution of rescaled charges
    unique_Q, counts = est.result.histogram()
    print(f"\nRescaled Q distribution:")
    for q, c in zip(unique_Q, counts):
        pct = 100.0 * c / len(est.Q_rescaled)
        print(f"  Q = {q:+2d}: {c:4d} ({pct:5.1f}%)")
    
    # Jackknife on rescaled charges
    Q_int = est.Q_rescaled
    if len(Q_int) >= 2:
        jk_Qi = jackknife_analysis(Q_int.astype(float), func=np.mean)
        jk_Qi2 = jackknife_analysis(Q_int.astype(float), func=lambda x: np.mean(x**2))
        print(f"\nJackknife analysis on integer charges:")
        print(f"  <Qi>  = {jk_Qi['estimate']:.4f} ± {jk_Qi['error']:.4f}")
        print(f"  <Qi²> = {jk_Qi2['estimate']:.4f} ± {jk_Qi2['error']:.4f}")
    
    print("\n" + "="*60)
    print("Done!")


if __name__ == '__main__':
    main()
