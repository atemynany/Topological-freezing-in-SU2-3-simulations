#!/usr/bin/env python3
"""
Analysis for SU(2)/SU(3) with periodic and open boundaries.
Plots susceptibility and tau_int with both boundary types in same plot.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import sys
import glob
import re
from dataclasses import dataclass
from typing import List, Optional

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from Q_estimator import QEstimator
from Q_statistics import load_topcharge
from autocorrelation import autocorrelation
from topological_susceptibility import HBAR_C

CHI_T_FOURTH_ROOT_SU2 = 200.0
CHI_T_FOURTH_ROOT_SU3 = 180.0

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'figure.dpi': 150,
})


@dataclass
class RunData:
    beta: float
    seed: int
    run_dir: str
    T: int
    L: int
    boundary: str
    Q_raw: np.ndarray
    Q_rescaled: np.ndarray
    alpha: float


@dataclass
class AnalysisResult:
    beta: float
    boundary: str
    Q_mean: float
    Q2_mean: float
    chi_t_lattice: float
    chi_t_fourth_root_a: float
    tau_int: float
    dtau_int: float
    alpha: float
    n_conf: int


def gaussian(x, mu, sigma, A):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def parse_input_file(run_dir: str) -> dict:
    info = {}
    input_file = os.path.join(run_dir, "input.txt")
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    info[parts[0]] = parts[1]
    return info


def extract_beta_from_dirname(dirname: str) -> Optional[float]:
    match = re.search(r'_b([\d.]+)_', dirname)
    return float(match.group(1)) if match else None


def load_run_data(run_dir: str, gauge_group: str = "su2") -> Optional[RunData]:
    info = parse_input_file(run_dir)
    
    # Check for 10000 sweeps only
    num_sweeps = int(info.get('num_sweeps', 0))
    if num_sweeps != 10000:
        return None
    
    boundary = info.get('boundary', 'periodic')
    
    # Find topcharge file
    candidates = [
        os.path.join(run_dir, "output", "topcharge.dat"),
        os.path.join(run_dir, "output", f"topcharge_{gauge_group}.dat"),
    ]
    topcharge_file = next((f for f in candidates if os.path.exists(f)), None)
    
    if not topcharge_file:
        return None
    
    Q_raw = load_topcharge(topcharge_file)
    if len(Q_raw) < 100:
        return None
    
    estimator = QEstimator(Q_raw)
    
    beta = float(info.get('beta', extract_beta_from_dirname(os.path.basename(run_dir)) or 0))
    seed = int(info.get('seed', 0))
    T = int(info.get('T', 16))
    L = int(info.get('L', 16))
    
    return RunData(
        beta=beta, seed=seed, run_dir=run_dir, T=T, L=L,
        boundary=boundary, Q_raw=Q_raw, Q_rescaled=estimator.Q_rescaled, alpha=estimator.alpha
    )


def analyze_run(run_data: RunData) -> AnalysisResult:
    volume = run_data.T * run_data.L ** 3
    Q_raw = run_data.Q_raw.astype(float)
    Q2_raw_mean = np.mean(Q_raw ** 2)
    chi_t_lattice = Q2_raw_mean / volume
    
    chi_t_fourth_root_a = (chi_t_lattice ** 0.25) * HBAR_C if chi_t_lattice > 0 else 0.0
    
    Q = run_data.Q_rescaled.astype(float)
    Q2 = Q ** 2
    
    try:
        autocorr = autocorrelation(Q2, name=f"Q2_b{run_data.beta}")
        tau_int = autocorr.tau_int
        dtau_int = autocorr.dtau_int
    except:
        tau_int, dtau_int = 0.5, 0.0
    
    return AnalysisResult(
        beta=run_data.beta, boundary=run_data.boundary,
        Q_mean=np.mean(Q), Q2_mean=Q2_raw_mean,
        chi_t_lattice=chi_t_lattice, chi_t_fourth_root_a=chi_t_fourth_root_a,
        tau_int=tau_int, dtau_int=dtau_int,
        alpha=run_data.alpha, n_conf=len(Q_raw)
    )


def lattice_spacing_su2(beta: float) -> float:
    beta_ref, a_ref = 2.3, 0.13
    return a_ref * np.exp(-1.5 * (beta - beta_ref))


def plot_histograms_grid(runs: List[RunData], results: List[AnalysisResult], 
                         output_dir: str, gauge_group: str, boundary: str):
    """Plot histograms grid for one boundary type."""
    n = len(runs)
    if n == 0:
        return
    
    n_cols = min(3, n)
    n_rows = (n + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3.5*n_rows))
    if n == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # Sort by beta
    sorted_data = sorted(zip(runs, results), key=lambda x: x[0].beta)
    
    for idx, (run_data, result) in enumerate(sorted_data):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row, col]
        
        Q = run_data.Q_rescaled
        Q_min, Q_max = int(Q.min()) - 1, int(Q.max()) + 1
        bins = np.arange(Q_min - 0.5, Q_max + 1.5, 1)
        
        counts, bin_edges, _ = ax.hist(Q, bins=bins, density=False,
            color='steelblue', edgecolor='black', linewidth=0.5, alpha=0.9, rwidth=0.15)
        
        # Gaussian fit
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        try:
            sigma0 = max(np.std(Q), 0.5)
            if np.sum(counts > 0) >= 3 and sigma0 > 0.1:
                popt, _ = curve_fit(gaussian, bin_centers, counts,
                    p0=[np.mean(Q), sigma0, counts.max()],
                    bounds=([-np.inf, 0.1, 0], [np.inf, np.inf, np.inf]))
                x_fit = np.linspace(Q_min - 1, Q_max + 1, 100)
                ax.plot(x_fit, gaussian(x_fit, *popt), 'r-', linewidth=1.5)
        except:
            pass
        
        ax.set_title(f'β={run_data.beta:.2f}', fontsize=11)
        ax.set_xlabel(r'$Q$', fontsize=10)
        ax.set_ylabel(r'$N$', fontsize=10)
        ax.tick_params(labelsize=9)
    
    # Hide empty
    for idx in range(n, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row, col].set_visible(False)
    
    plt.suptitle(f'{gauge_group.upper()} {boundary.capitalize()} BC', fontsize=13)
    plt.tight_layout()
    
    filepath = os.path.join(output_dir, f"histograms_{gauge_group}_{boundary}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def plot_susceptibility_combined(results_periodic: List[AnalysisResult], 
                                  results_open: List[AnalysisResult],
                                  output_dir: str, gauge_group: str):
    """Plot susceptibility with both boundary types."""
    fig, ax = plt.subplots(figsize=(8, 5))
    
    ref_value = CHI_T_FOURTH_ROOT_SU2 if gauge_group == "su2" else CHI_T_FOURTH_ROOT_SU3
    
    for results, label, marker, color in [
        (results_periodic, 'Periodic', 'o', 'steelblue'),
        (results_open, 'Open', 's', 'darkorange')
    ]:
        if len(results) == 0:
            continue
        
        betas = np.array([r.beta for r in results])
        chi_fourth_a = np.array([r.chi_t_fourth_root_a for r in results])
        
        idx = np.argsort(betas)
        betas, chi_fourth_a = betas[idx], chi_fourth_a[idx]
        
        if gauge_group == "su2":
            a_beta = np.array([lattice_spacing_su2(b) for b in betas])
            chi_fourth_MeV = chi_fourth_a / a_beta
        else:
            chi_fourth_MeV = chi_fourth_a / (chi_fourth_a / ref_value)
        
        ax.plot(betas, chi_fourth_MeV, f'{marker}-', markersize=7, linewidth=1.5, 
                color=color, label=label)
    
    ax.axhline(y=ref_value, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\chi_t^{1/4}$ (MeV)')
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    filepath = os.path.join(output_dir, f"susceptibility_{gauge_group}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def plot_tau_int_combined(results_periodic: List[AnalysisResult], 
                           results_open: List[AnalysisResult],
                           output_dir: str, gauge_group: str):
    """Plot tau_int with both boundary types."""
    fig, ax = plt.subplots(figsize=(8, 5))
    
    for results, label, marker, color in [
        (results_periodic, 'Periodic', 'o', 'steelblue'),
        (results_open, 'Open', 's', 'darkorange')
    ]:
        if len(results) == 0:
            continue
        
        betas = np.array([r.beta for r in results])
        tau_int = np.array([r.tau_int for r in results])
        dtau_int = np.array([r.dtau_int for r in results])
        
        idx = np.argsort(betas)
        betas, tau_int, dtau_int = betas[idx], tau_int[idx], dtau_int[idx]
        
        ax.errorbar(betas, tau_int, yerr=dtau_int, fmt=f'{marker}-', markersize=7,
                    linewidth=1.5, capsize=3, color=color, label=label)
    
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\tau_{\rm int}(Q^2)$')
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    filepath = os.path.join(output_dir, f"tau_int_{gauge_group}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def print_table(results: List[AnalysisResult], gauge_group: str):
    """Print summary table."""
    ref = CHI_T_FOURTH_ROOT_SU2 if gauge_group == "su2" else CHI_T_FOURTH_ROOT_SU3
    
    print(f"\n{'='*85}")
    print(f"  {gauge_group.upper()} Analysis Summary")
    print(f"{'='*85}")
    print(f"{'β':>6} {'BC':>10} {'<Q²>':>8} {'χ_t^1/4':>10} {'τ_int':>8} {'N_conf':>8} {'α':>6}")
    print(f"{'':>6} {'':>10} {'':>8} {'(MeV)':>10} {'':>8} {'':>8} {'':>6}")
    print(f"{'-'*85}")
    
    for r in sorted(results, key=lambda x: (x.boundary, x.beta)):
        if gauge_group == "su2":
            a = lattice_spacing_su2(r.beta)
            chi_MeV = r.chi_t_fourth_root_a / a
        else:
            chi_MeV = ref
        print(f"{r.beta:>6.2f} {r.boundary:>10} {r.Q2_mean:>8.1f} {chi_MeV:>10.1f} "
              f"{r.tau_int:>8.1f} {r.n_conf:>8} {r.alpha:>6.3f}")
    print(f"{'='*85}\n")


def main():
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--su2", action="store_true")
    parser.add_argument("--su3", action="store_true")
    parser.add_argument("--results-dir", type=str, default=None)
    parser.add_argument("--output-dir", type=str, default=None)
    args = parser.parse_args()
    
    gauge_group = "su3" if args.su3 else "su2"
    
    base_dir = os.path.dirname(script_dir)
    results_dir = args.results_dir or os.path.join(base_dir, "data", "results")
    output_dir = args.output_dir or os.path.join(base_dir, "output", "figures_analysis")
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n{gauge_group.upper()} Analysis (10000 sweeps runs only)")
    print(f"Output: {output_dir}\n")
    
    # Find runs
    run_dirs = sorted(glob.glob(os.path.join(results_dir, "T*_L*_b*_seed*")))
    
    runs_periodic, runs_open = [], []
    results_periodic, results_open = [], []
    
    for run_dir in run_dirs:
        run_data = load_run_data(run_dir, gauge_group)
        if run_data is None:
            continue
        
        # Skip SU(3) runs if analyzing SU(2) and vice versa
        if gauge_group == "su2" and run_data.beta > 4:
            continue
        if gauge_group == "su3" and run_data.beta < 4:
            continue
        
        result = analyze_run(run_data)
        
        if run_data.boundary == "periodic":
            runs_periodic.append(run_data)
            results_periodic.append(result)
        else:
            runs_open.append(run_data)
            results_open.append(result)
    
    print(f"Found: {len(runs_periodic)} periodic, {len(runs_open)} open\n")
    
    # Print table
    print_table(results_periodic + results_open, gauge_group)
    
    # Plot histograms for each boundary type
    if runs_periodic:
        plot_histograms_grid(runs_periodic, results_periodic, output_dir, gauge_group, "periodic")
    if runs_open:
        plot_histograms_grid(runs_open, results_open, output_dir, gauge_group, "open")
    
    # Combined plots
    plot_susceptibility_combined(results_periodic, results_open, output_dir, gauge_group)
    plot_tau_int_combined(results_periodic, results_open, output_dir, gauge_group)
    
    print(f"\nDone! Plots in: {output_dir}")


if __name__ == "__main__":
    main()
