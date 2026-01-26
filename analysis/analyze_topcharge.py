#!/usr/bin/env python3
"""
==============================================================================
analyze_topcharge.py
==============================================================================
Statistical analysis of topological charge measurements.
Computes topological susceptibility, autocorrelation times, and performs
binning analysis for error estimation.

Usage:
    python3 analysis/analyze_topcharge.py --input <data_file> --output <output_file>

Author: Alexander de Barros Noll
Date: January 2026
==============================================================================
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np


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
# Statistical Analysis
# ==============================================================================

def compute_mean_error(values: np.ndarray) -> tuple:
    """Compute mean and standard error of the mean."""
    n = len(values)
    mean = np.mean(values)
    std = np.std(values, ddof=1)
    error = std / np.sqrt(n)
    return mean, error


def compute_autocorrelation(values: np.ndarray, max_lag: int = None) -> np.ndarray:
    """Compute normalized autocorrelation function."""
    n = len(values)
    if max_lag is None:
        max_lag = n // 2
    
    values_centered = values - np.mean(values)
    var = np.var(values, ddof=0)
    
    if var < 1e-15:
        return np.zeros(max_lag)
    
    autocorr = np.correlate(values_centered, values_centered, mode='full')[n-1:]
    autocorr = autocorr[:max_lag] / (var * n)
    
    return autocorr


def compute_integrated_autocorr_time(autocorr: np.ndarray, c: float = 5.0) -> tuple:
    """
    Compute integrated autocorrelation time with automatic window selection.
    
    Uses the Madras-Sokal automatic windowing procedure:
    tau_int = 0.5 + sum_{t=1}^{W} rho(t)
    where W is chosen as the first point where W >= c * tau_int(W)
    """
    n = len(autocorr)
    tau_int = 0.5
    
    for w in range(1, n):
        tau_int += autocorr[w]
        
        if w >= c * tau_int:
            break
    
    # Error estimate
    tau_int_error = tau_int * np.sqrt(2 * (2*w + 1) / n)
    
    return tau_int, tau_int_error, w


def binning_analysis(values: np.ndarray, max_binsize: int = None) -> dict:
    """
    Perform binning analysis to estimate true errors accounting for autocorrelation.
    
    Returns dictionary with:
        - bin_sizes: array of bin sizes tested
        - errors: error estimates for each bin size
        - effective_samples: effective number of independent samples
    """
    n = len(values)
    if max_binsize is None:
        max_binsize = n // 10
    
    bin_sizes = []
    errors = []
    
    for binsize in range(1, max_binsize + 1):
        n_bins = n // binsize
        if n_bins < 2:
            break
        
        # Compute binned values
        binned = np.array([np.mean(values[i*binsize:(i+1)*binsize]) 
                          for i in range(n_bins)])
        
        # Error of binned mean
        mean = np.mean(binned)
        std = np.std(binned, ddof=1)
        error = std / np.sqrt(n_bins)
        
        bin_sizes.append(binsize)
        errors.append(error)
    
    bin_sizes = np.array(bin_sizes)
    errors = np.array(errors)
    
    # Estimate effective number of samples
    # When errors plateau, that gives the true error
    if len(errors) > 0:
        naive_error = errors[0]  # Error without binning
        plateau_error = np.max(errors)  # Approximate plateau error
        # Effective samples: n_eff = n * (naive_error / plateau_error)^2
        if plateau_error > 0:
            n_eff = n * (naive_error / plateau_error) ** 2
        else:
            n_eff = n
    else:
        n_eff = n
    
    return {
        'bin_sizes': bin_sizes,
        'errors': errors,
        'effective_samples': n_eff
    }


def compute_topological_susceptibility(Q_values: np.ndarray, volume: int) -> tuple:
    """
    Compute topological susceptibility chi_t = <Q^2> / V
    
    Returns: (chi_t, chi_t_error)
    """
    Q2_values = Q_values ** 2
    Q2_mean, Q2_error = compute_mean_error(Q2_values)
    
    chi_t = Q2_mean / volume
    chi_t_error = Q2_error / volume
    
    return chi_t, chi_t_error


# ==============================================================================
# Analysis
# ==============================================================================

def analyze_smearing_level(Q_values: np.ndarray, smear_level: int, volume: int) -> dict:
    """Perform complete analysis for one smearing level."""
    results = {
        'smear_level': smear_level,
        'n_configs': len(Q_values)
    }
    
    # Basic statistics
    results['Q_mean'], results['Q_mean_error'] = compute_mean_error(Q_values)
    results['Q_std'] = np.std(Q_values, ddof=1)
    results['Q_min'] = np.min(Q_values)
    results['Q_max'] = np.max(Q_values)
    
    # Topological susceptibility
    results['chi_t'], results['chi_t_error'] = compute_topological_susceptibility(Q_values, volume)
    
    # Autocorrelation
    autocorr = compute_autocorrelation(Q_values)
    tau_int, tau_int_error, window = compute_integrated_autocorr_time(autocorr)
    results['tau_int'] = tau_int
    results['tau_int_error'] = tau_int_error
    results['tau_window'] = window
    
    # Binning analysis
    binning = binning_analysis(Q_values)
    results['effective_samples'] = binning['effective_samples']
    results['binning_error_plateau'] = np.max(binning['errors']) if len(binning['errors']) > 0 else 0
    
    # Corrected error (accounting for autocorrelation)
    results['Q_mean_error_corrected'] = results['Q_mean_error'] * np.sqrt(2 * tau_int)
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Statistical analysis of topological charge data'
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Input data file from meas_topcharge')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file for analysis results')
    parser.add_argument('--T', type=int, default=16,
                        help='Temporal lattice extent')
    parser.add_argument('--L', type=int, default=16,
                        help='Spatial lattice extent')
    
    args = parser.parse_args()
    
    # Load data
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    print(f"Loading data from: {args.input}")
    data = load_topcharge_data(args.input)
    
    if len(data['Q']) == 0:
        print("Error: No data found")
        sys.exit(1)
    
    volume = args.T * args.L ** 3
    print(f"Lattice volume: {args.T} x {args.L}^3 = {volume}")
    
    # Analyze each smearing level
    smear_levels = np.unique(data['smear_steps'])
    all_results = []
    
    for smear in smear_levels:
        mask = data['smear_steps'] == smear
        Q_values = data['Q'][mask]
        
        if len(Q_values) < 2:
            continue
        
        results = analyze_smearing_level(Q_values, smear, volume)
        all_results.append(results)
        
        print(f"\nSmearing level {smear}:")
        print(f"  <Q> = {results['Q_mean']:.4f} +/- {results['Q_mean_error_corrected']:.4f}")
        print(f"  chi_t = {results['chi_t']:.6f} +/- {results['chi_t_error']:.6f}")
        print(f"  tau_int = {results['tau_int']:.2f} +/- {results['tau_int_error']:.2f}")
    
    # Write results
    with open(args.output, 'w') as f:
        f.write("# Topological Charge Analysis Results\n")
        f.write(f"# Lattice: {args.T} x {args.L}^3\n")
        f.write(f"# Volume: {volume}\n")
        f.write("#\n")
        f.write("# Columns: smear_level  n_configs  Q_mean  Q_mean_err  Q_std  chi_t  chi_t_err  tau_int  tau_int_err\n")
        
        for r in all_results:
            f.write(f"{r['smear_level']:5d}  {r['n_configs']:5d}  "
                    f"{r['Q_mean']:10.6f}  {r['Q_mean_error_corrected']:10.6f}  {r['Q_std']:10.6f}  "
                    f"{r['chi_t']:12.8f}  {r['chi_t_error']:12.8f}  "
                    f"{r['tau_int']:8.3f}  {r['tau_int_error']:8.3f}\n")
    
    print(f"\nResults written to: {args.output}")


if __name__ == '__main__':
    main()
