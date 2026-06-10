#!/usr/bin/env python3
"""
==============================================================================
plot_topcharge.py
==============================================================================
Plotting script for SU(2) topological charge analysis.
Creates publication-quality plots of Q vs smearing steps and Monte Carlo time.

Usage:
    python3 analysis/plot_topcharge.py --input <data_file> --output <output_dir>

Author: Alexander de Barros Noll
Date: January 2026
==============================================================================
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from calculations import find_optimal_alpha

# ==============================================================================
# Plot Configuration
# ==============================================================================

# Publication-quality plot settings
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 12
rcParams['axes.labelsize'] = 14
rcParams['axes.titlesize'] = 14
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12
rcParams['legend.fontsize'] = 11
rcParams['figure.figsize'] = (8, 6)
rcParams['figure.dpi'] = 150
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'


# ==============================================================================
# Data Loading
# ==============================================================================

def load_topcharge_data(filename: str) -> dict:
    """
    Load topological charge data from output file.
    
    Expected format:
        # Header comments
        smear_steps  config_number  Q  plaquette
        
    Returns:
        Dictionary with arrays: smear_steps, configs, Q, plaquette
    """
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


def parse_input_file(input_file: str) -> dict:
    """Parse simple key-value input files from run directories."""
    info = {}
    if not input_file or not os.path.isfile(input_file):
        return info
    with open(input_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                info[parts[0]] = parts[1]
    return info


def infer_run_dir(input_path: str) -> str:
    """Infer run directory from .../output/topcharge.dat paths."""
    path = os.path.abspath(input_path)
    parent = os.path.dirname(path)
    if os.path.basename(parent) in {"output", "topcharge"}:
        return os.path.dirname(parent)
    return parent


def timeslice_file_for_input(input_path: str, explicit_path: str = None) -> str:
    if explicit_path:
        return explicit_path
    run_dir = infer_run_dir(input_path)
    for candidate in [
        os.path.join(run_dir, "output", "topcharge_timeslice.dat"),
        os.path.join(run_dir, "topcharge", "topcharge_timeslice.dat"),
    ]:
        if os.path.isfile(candidate):
            return candidate
    return ""


def load_temporal_subvolume_data(timeslice_file: str, reference_data: dict,
                                 T: int, exclude_boundary_slices: int,
                                 alpha_rescale: bool = True) -> dict:
    """Load continuous temporal-subvolume Q from topcharge_timeslice.dat."""
    plaquette_by_key = {
        (int(s), int(c)): float(p)
        for s, c, p in zip(reference_data['smear_steps'],
                           reference_data['config'],
                           reference_data['plaquette'])
    }
    t_min = max(0, exclude_boundary_slices)
    t_max = T - max(0, exclude_boundary_slices)
    by_key = {}

    with open(timeslice_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                smear = int(parts[0])
                conf = int(parts[1])
                t = int(parts[2])
                q_t = float(parts[3])
            except ValueError:
                continue
            if t_min <= t < t_max:
                by_key[(smear, conf)] = by_key.get((smear, conf), 0.0) + q_t

    if not by_key:
        return reference_data

    rows = []
    alpha_by_smear = {}
    for smear in sorted({key[0] for key in by_key}):
        configs = sorted(conf for s, conf in by_key if s == smear)
        q_raw = np.array([by_key[(smear, conf)] for conf in configs], dtype=float)
        alpha = find_optimal_alpha(q_raw) if alpha_rescale and len(q_raw) else 1.0
        alpha_by_smear[smear] = alpha
        for conf, q_val in zip(configs, alpha * q_raw):
            rows.append((smear, conf, q_val, plaquette_by_key.get((smear, conf), np.nan)))

    data = {
        'smear_steps': np.array([row[0] for row in rows]),
        'config': np.array([row[1] for row in rows]),
        'Q': np.array([row[2] for row in rows], dtype=float),
        'plaquette': np.array([row[3] for row in rows], dtype=float),
        'continuous_q': True,
        'alpha_by_smear': alpha_by_smear,
        'source': 'temporal_subvolume',
    }
    return data


def extract_metadata(filename: str) -> dict:
    """Extract metadata from header comments."""
    metadata = {}
    
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                break
            line = line[1:].strip()
            if ':' in line:
                key, value = line.split(':', 1)
                metadata[key.strip()] = value.strip()
    
    return metadata


# ==============================================================================
# Plotting Functions
# ==============================================================================

def plot_Q_vs_smearing(data: dict, output_dir: str, title: str = None):
    """
    Plot topological charge Q as a function of smearing steps.
    Shows Q at different smearing levels for each configuration.
    """
    fig, ax = plt.subplots()
    
    # Get unique configurations and smearing steps
    configs = np.unique(data['config'])
    smear_steps = np.unique(data['smear_steps'])
    
    # Plot Q vs smearing for each configuration
    cmap = plt.cm.viridis
    colors = cmap(np.linspace(0, 1, len(configs)))
    
    for i, conf in enumerate(configs):
        mask = data['config'] == conf
        ax.plot(data['smear_steps'][mask], data['Q'][mask], 
                'o-', color=colors[i], alpha=0.7, markersize=4,
                label=f'conf {conf}' if len(configs) <= 10 else None)
    
    # Add horizontal lines at integer values
    Q_range = (np.floor(data['Q'].min()) - 1, np.ceil(data['Q'].max()) + 1)
    for q_int in range(int(Q_range[0]), int(Q_range[1]) + 1):
        ax.axhline(y=q_int, color='gray', linestyle='--', alpha=0.3, linewidth=0.5)
    
    ax.set_xlabel('APE Smearing Steps')
    ylabel = r'Topological Charge $\alpha Q$' if data.get('continuous_q') else 'Topological Charge Q'
    ax.set_ylabel(ylabel)
    ax.set_title(title or 'Topological Charge vs Smearing Steps')
    
    if len(configs) <= 10:
        ax.legend(loc='best', ncol=2)
    
    ax.grid(True, alpha=0.3)
    
    # Save
    output_path = os.path.join(output_dir, 'Q_vs_smearing.pdf')
    plt.savefig(output_path)
    plt.savefig(output_path.replace('.pdf', '.png'))
    plt.close()
    
    print(f"Saved: {output_path}")


def plot_Q_vs_mctime(data: dict, output_dir: str, smear_level: int = None, title: str = None):
    """
    Plot topological charge Q as a function of Monte Carlo time.
    Shows Q after a specific number of smearing steps.
    """
    fig, ax = plt.subplots()
    
    smear_steps = np.unique(data['smear_steps'])
    
    # If smear_level not specified, use the largest available
    if smear_level is None:
        smear_level = smear_steps[-1] if len(smear_steps) > 0 else 0
    
    # Filter data for this smearing level
    mask = data['smear_steps'] == smear_level
    configs = data['config'][mask]
    Q_values = data['Q'][mask]
    
    ax.plot(configs, Q_values, 'o', color='blue', markersize=6)
    
    # Add horizontal lines at integer values
    Q_range = (np.floor(Q_values.min()) - 1, np.ceil(Q_values.max()) + 1)
    for q_int in range(int(Q_range[0]), int(Q_range[1]) + 1):
        ax.axhline(y=q_int, color='gray', linestyle='--', alpha=0.3, linewidth=0.5)
    
    ax.set_xlabel('Monte Carlo Time (Configuration Number)')
    ylabel = r'Topological Charge $\alpha Q$' if data.get('continuous_q') else 'Topological Charge Q'
    ax.set_ylabel(ylabel)
    ax.set_title(title or f'Topological Charge vs MC Time (smear={smear_level})')
    ax.grid(True, alpha=0.3)
    
    # Add statistics
    Q_mean = np.mean(Q_values)
    Q_std = np.std(Q_values)
    stats_text = f'$\\langle Q \\rangle = {Q_mean:.3f} \\pm {Q_std:.3f}$'
    ax.text(0.95, 0.95, stats_text, transform=ax.transAxes,
            ha='right', va='top', fontsize=11,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Save
    output_path = os.path.join(output_dir, f'Q_vs_mctime_smear{smear_level}.pdf')
    plt.savefig(output_path)
    plt.savefig(output_path.replace('.pdf', '.png'))
    plt.close()
    
    print(f"Saved: {output_path}")


def plot_Q_histogram(data: dict, output_dir: str, smear_level: int = None, title: str = None):
    """
    Plot histogram of topological charge values.
    """
    fig, ax = plt.subplots()
    
    smear_steps = np.unique(data['smear_steps'])
    
    if smear_level is None:
        smear_level = smear_steps[-1] if len(smear_steps) > 0 else 0
    
    mask = data['smear_steps'] == smear_level
    Q_values = data['Q'][mask]
    
    if data.get('continuous_q'):
        Q_min = np.floor(Q_values.min()) - 0.5
        Q_max = np.ceil(Q_values.max()) + 0.5
        bins = np.linspace(Q_min, Q_max, 40)
    else:
        Q_min = np.floor(Q_values.min()) - 0.5
        Q_max = np.ceil(Q_values.max()) + 0.5
        bins = np.arange(Q_min, Q_max + 1, 1)
    
    ax.hist(Q_values, bins=bins, color='steelblue', edgecolor='black', alpha=0.7)
    
    xlabel = r'Topological Charge $\alpha Q$' if data.get('continuous_q') else 'Topological Charge Q'
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Count')
    ax.set_title(title or f'Topological Charge Distribution (smear={smear_level})')
    
    # Add statistics
    Q_mean = np.mean(Q_values)
    Q_std = np.std(Q_values)
    chi_susc = Q_std**2  # Topological susceptibility (unnormalized)
    
    stats_text = f'$\\langle Q \\rangle = {Q_mean:.3f}$\n$\\sigma_Q = {Q_std:.3f}$\n$\\chi_t = {chi_susc:.3f}$'
    ax.text(0.95, 0.95, stats_text, transform=ax.transAxes,
            ha='right', va='top', fontsize=11,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.grid(True, alpha=0.3, axis='y')
    
    # Save
    output_path = os.path.join(output_dir, f'Q_histogram_smear{smear_level}.pdf')
    plt.savefig(output_path)
    plt.savefig(output_path.replace('.pdf', '.png'))
    plt.close()
    
    print(f"Saved: {output_path}")


def plot_plaquette_vs_smearing(data: dict, output_dir: str, title: str = None):
    """
    Plot average plaquette as a function of smearing steps.
    """
    fig, ax = plt.subplots()
    
    smear_steps = np.unique(data['smear_steps'])
    
    # Compute mean and std for each smearing level
    plaq_mean = []
    plaq_std = []
    
    for smear in smear_steps:
        mask = data['smear_steps'] == smear
        plaq_mean.append(np.mean(data['plaquette'][mask]))
        plaq_std.append(np.std(data['plaquette'][mask]))
    
    plaq_mean = np.array(plaq_mean)
    plaq_std = np.array(plaq_std)
    
    ax.errorbar(smear_steps, plaq_mean, yerr=plaq_std, 
                fmt='o-', color='green', markersize=6, capsize=3)
    
    ax.set_xlabel('APE Smearing Steps')
    ax.set_ylabel('Average Plaquette $\\langle P \\rangle$')
    ax.set_title(title or 'Plaquette vs Smearing Steps')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0, top=1.05)
    
    # Save
    output_path = os.path.join(output_dir, 'plaquette_vs_smearing.pdf')
    plt.savefig(output_path)
    plt.savefig(output_path.replace('.pdf', '.png'))
    plt.close()
    
    print(f"Saved: {output_path}")


def plot_Q_smearing_heatmap(data: dict, output_dir: str, title: str = None):
    """
    Create a heatmap of Q values: x-axis = MC time, y-axis = smearing steps.
    """
    configs = np.unique(data['config'])
    smear_steps = np.unique(data['smear_steps'])
    
    # Create 2D array
    Q_matrix = np.zeros((len(smear_steps), len(configs)))
    
    for i, smear in enumerate(smear_steps):
        for j, conf in enumerate(configs):
            mask = (data['smear_steps'] == smear) & (data['config'] == conf)
            if np.any(mask):
                Q_matrix[i, j] = data['Q'][mask][0]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x_min, x_max = configs[0], configs[-1]
    y_min, y_max = smear_steps[0], smear_steps[-1]
    if x_min == x_max:
        x_min, x_max = x_min - 0.5, x_max + 0.5
    if y_min == y_max:
        y_min, y_max = y_min - 0.5, y_max + 0.5

    im = ax.imshow(Q_matrix, aspect='auto', origin='lower',
                   extent=[x_min, x_max, y_min, y_max],
                   cmap='RdBu_r')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar_label = r'Topological Charge $\alpha Q$' if data.get('continuous_q') else 'Topological Charge Q'
    cbar.set_label(cbar_label)
    
    ax.set_xlabel('Monte Carlo Time (Configuration Number)')
    ax.set_ylabel('APE Smearing Steps')
    ax.set_title(title or 'Topological Charge Heatmap')
    
    # Save
    output_path = os.path.join(output_dir, 'Q_heatmap.pdf')
    plt.savefig(output_path)
    plt.savefig(output_path.replace('.pdf', '.png'))
    plt.close()
    
    print(f"Saved: {output_path}")


def plot_Q_autocorrelation(data: dict, output_dir: str, smear_level: int = None, title: str = None):
    """
    Plot autocorrelation function of topological charge.
    """
    smear_steps = np.unique(data['smear_steps'])
    
    if smear_level is None:
        smear_level = smear_steps[-1] if len(smear_steps) > 0 else 0
    
    mask = data['smear_steps'] == smear_level
    Q_values = data['Q'][mask]
    
    # Compute autocorrelation
    n = len(Q_values)
    Q_centered = Q_values - np.mean(Q_values)
    autocorr = np.correlate(Q_centered, Q_centered, mode='full')[n-1:]
    autocorr = autocorr / autocorr[0]
    
    fig, ax = plt.subplots()
    
    max_lag = min(len(autocorr), 50)
    lags = np.arange(max_lag)
    
    ax.bar(lags, autocorr[:max_lag], color='steelblue', edgecolor='black', alpha=0.7)
    ax.axhline(y=0, color='black', linewidth=0.5)
    
    # 95% confidence bounds for white noise
    conf_bound = 1.96 / np.sqrt(n)
    ax.axhline(y=conf_bound, color='red', linestyle='--', alpha=0.5)
    ax.axhline(y=-conf_bound, color='red', linestyle='--', alpha=0.5)
    
    ax.set_xlabel('Lag (Configuration Separation)')
    ax.set_ylabel('Autocorrelation')
    ax.set_title(title or f'Q Autocorrelation (smear={smear_level})')
    ax.set_xlim(-0.5, max_lag - 0.5)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Save
    output_path = os.path.join(output_dir, f'Q_autocorrelation_smear{smear_level}.pdf')
    plt.savefig(output_path)
    plt.savefig(output_path.replace('.pdf', '.png'))
    plt.close()
    
    print(f"Saved: {output_path}")


# ==============================================================================
# Main
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Plot topological charge results from SU(2) lattice simulation'
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Input data file from meas_topcharge')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory for plots')
    parser.add_argument('--title', '-t', default=None,
                        help='Base title for plots')
    parser.add_argument('--smear-level', type=int, default=None,
                        help='Specific smearing level for MC time plots')
    parser.add_argument('--timeslice-input', default=None,
                        help='Optional topcharge_timeslice.dat for temporal-subvolume Q reconstruction')
    parser.add_argument('--T', type=int, default=None,
                        help='Temporal lattice extent; inferred from input.txt when possible')
    parser.add_argument('--exclude-boundary-slices', type=int, default=None,
                        help='Use continuous temporal-subvolume Q by excluding this many slices at each edge')
    parser.add_argument('--no-alpha-rescale', action='store_true',
                        help='For temporal subvolume plots, show raw cropped Q instead of alpha*Q')
    
    args = parser.parse_args()
    
    # Validate input
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Load data
    print(f"Loading data from: {args.input}")
    data = load_topcharge_data(args.input)
    data['continuous_q'] = False
    
    if len(data['Q']) == 0:
        print("Error: No data found in input file")
        sys.exit(1)
    
    print(f"Loaded {len(data['Q'])} data points")
    print(f"Configurations: {len(np.unique(data['config']))}")
    print(f"Smearing levels: {len(np.unique(data['smear_steps']))}")

    run_dir = infer_run_dir(args.input)
    input_info = parse_input_file(os.path.join(run_dir, "input.txt"))
    T = args.T if args.T is not None else int(input_info.get("T", 0) or 0)
    exclude_boundary_slices = (
        args.exclude_boundary_slices
        if args.exclude_boundary_slices is not None
        else int(input_info.get("exclude_boundary_slices", 0) or 0)
    )
    timeslice_file = timeslice_file_for_input(args.input, args.timeslice_input)
    if exclude_boundary_slices > 0:
        if T <= 0:
            print("Error: --T is required for temporal-subvolume plots when input.txt is unavailable")
            sys.exit(1)
        if not timeslice_file or not os.path.isfile(timeslice_file):
            print("Error: topcharge_timeslice.dat is required for temporal-subvolume plots")
            sys.exit(1)
        print(f"Using temporal subvolume from: {timeslice_file}")
        print(f"Excluding t < {exclude_boundary_slices} and t >= {T - exclude_boundary_slices}")
        data = load_temporal_subvolume_data(
            timeslice_file,
            data,
            T=T,
            exclude_boundary_slices=exclude_boundary_slices,
            alpha_rescale=not args.no_alpha_rescale,
        )
        print("Plotting continuous alpha*Q values (not integer-rounded)")
        for smear, alpha in data.get('alpha_by_smear', {}).items():
            print(f"  smear {smear}: alpha={alpha:.6f}")
    
    # Generate plots
    print("\nGenerating plots...")
    
    plot_Q_vs_smearing(data, args.output, args.title)
    plot_Q_vs_mctime(data, args.output, args.smear_level, args.title)
    plot_Q_histogram(data, args.output, args.smear_level, args.title)
    plot_plaquette_vs_smearing(data, args.output, args.title)
    plot_Q_smearing_heatmap(data, args.output, args.title)
    plot_Q_autocorrelation(data, args.output, args.smear_level, args.title)
    
    print(f"\nAll plots saved to: {args.output}")


if __name__ == '__main__':
    main()
