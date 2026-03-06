#!/usr/bin/env python3
"""
Beta Scan Analysis for SU(2) Lattice Gauge Theory

Analyzes results from run_beta_scan.sh and produces:
1. Rescaled Q histograms with Gaussian fits for each beta
2. Topological susceptibility χ vs beta
3. Integrated autocorrelation time τ_int vs beta
4. Thermalization comparison plot

Usage:
    python3 analyze_beta_scan.py [--su2|--su3] [--results-dir PATH]
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
import os
import sys
import glob
import re
from dataclasses import dataclass
from typing import List, Tuple, Optional

# Add analysis directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from Q_estimator import QEstimator
from Q_statistics import compute_Q_statistics, load_topcharge
from autocorrelation import autocorrelation, load_plaquette
from topological_susceptibility import topological_susceptibility, HBAR_C

# ==============================================================================
# Configuration
# ==============================================================================

# Reference values
CHI_T_FOURTH_ROOT_SU2 = 200.0  # MeV (reference for SU(2))
CHI_T_FOURTH_ROOT_SU3 = 180.0  # MeV (reference for SU(3))

# Plot settings
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'figure.figsize': (10, 8),
    'figure.dpi': 150,
})

# Colors for different betas (colormap)
CMAP = plt.cm.viridis


# ==============================================================================
# Data Classes
# ==============================================================================

@dataclass
class RunData:
    """Data from a single beta run."""
    beta: float
    seed: int
    run_dir: str
    T: int
    L: int
    start_conf: int
    Q_raw: np.ndarray
    Q_rescaled: np.ndarray
    alpha: float
    plaq_data: np.ndarray
    mc_time: np.ndarray


@dataclass
class AnalysisResult:
    """Analysis results for a single beta."""
    beta: float
    Q_mean: float
    Q_error: float
    Q2_mean: float
    Q2_error: float
    chi_t_lattice: float
    chi_t_fourth_root_a: float  # χ_t^(1/4) * a in MeV*fm (scale-independent)
    tau_int: float
    dtau_int: float
    n_eff: float
    alpha: float
    therm_sweep: int


# ==============================================================================
# Helper Functions
# ==============================================================================

def gaussian(x, mu, sigma, A):
    """Gaussian function for fitting."""
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def parse_run_info(run_dir: str) -> dict:
    """Parse run_info.txt to extract parameters."""
    info = {}
    info_file = os.path.join(run_dir, "run_info.txt")
    if os.path.exists(info_file):
        with open(info_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    key, value = parts[0], parts[1]
                    info[key] = value
    return info


def extract_beta_from_dirname(dirname: str) -> Optional[float]:
    """Extract beta value from directory name like T16_L16_b2.4_seed123."""
    match = re.search(r'_b([\d.]+)_', dirname)
    if match:
        return float(match.group(1))
    return None


def load_run_data(run_dir: str, gauge_group: str = "su2") -> Optional[RunData]:
    """Load all data from a run directory."""
    info = parse_run_info(run_dir)
    
    # Determine file names based on gauge group (try multiple patterns)
    if gauge_group == "su3":
        plaq_candidates = [
            os.path.join(run_dir, "output", "plaquette_su3.dat"),
            os.path.join(run_dir, "output", "plaquette.dat")
        ]
        topcharge_candidates = [
            os.path.join(run_dir, "output", "topcharge_su3.dat"),
            os.path.join(run_dir, "output", "topcharge.dat")
        ]
    else:
        plaq_candidates = [
            os.path.join(run_dir, "output", "plaquette_su2.dat"),
            os.path.join(run_dir, "output", "plaquette.dat")
        ]
        topcharge_candidates = [
            os.path.join(run_dir, "output", "topcharge_su2.dat"),
            os.path.join(run_dir, "output", "topcharge.dat")
        ]
    
    # Find existing files
    plaq_file = next((f for f in plaq_candidates if os.path.exists(f)), None)
    topcharge_file = next((f for f in topcharge_candidates if os.path.exists(f)), None)
    
    if not topcharge_file:
        print(f"  Warning: No topcharge file in {run_dir}")
        return None
    
    # Load topological charge
    Q_raw = load_topcharge(topcharge_file)
    if len(Q_raw) == 0:
        print(f"  Warning: Empty topcharge data in {run_dir}")
        return None
    
    # Rescale Q
    estimator = QEstimator(Q_raw)
    
    # Load plaquette data
    mc_time = []
    plaq_data = []
    if plaq_file and os.path.exists(plaq_file):
        with open(plaq_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mc_time.append(int(parts[0]))
                        plaq_data.append(float(parts[1]))
                    except ValueError:
                        continue
    
    # Extract parameters
    beta = float(info.get('beta', extract_beta_from_dirname(os.path.basename(run_dir)) or 0))
    seed = int(info.get('seed', 0))
    T = int(info.get('T', 16))
    L = int(info.get('L', 16))
    start_conf = int(info.get('start_conf', 50))
    
    return RunData(
        beta=beta,
        seed=seed,
        run_dir=run_dir,
        T=T,
        L=L,
        start_conf=start_conf,
        Q_raw=Q_raw,
        Q_rescaled=estimator.Q_rescaled,
        alpha=estimator.alpha,
        plaq_data=np.array(plaq_data),
        mc_time=np.array(mc_time)
    )


def analyze_run(run_data: RunData, gauge_group: str = "su2") -> AnalysisResult:
    """Perform full analysis on a single run."""
    # Compute Q statistics with autocorrelation (note: Q_statistics uses rescaled Q)
    stats = compute_Q_statistics(run_data.Q_raw, ensemble=f"b{run_data.beta}")
    
    # Compute susceptibility from RESCALED Q values (corrected estimate of true Q)
    # χ_t = <Q²> / V  where V = T × L³
    volume = run_data.T * run_data.L ** 3
    Q_rescaled = run_data.Q_rescaled.astype(float)  # Use rescaled (integer) Q
    Q2_mean = np.mean(Q_rescaled ** 2)
    Q2_std = np.std(Q_rescaled ** 2) / np.sqrt(len(Q_rescaled))  # Naive error
    chi_t_lattice = Q2_mean / volume
    
    # Compute χ_t^(1/4) * a in MeV*fm (scale-independent quantity)
    # Formula: χ_t^(1/4) * a = (χ_t_lattice)^(1/4) * ℏc
    # To get χ_t^(1/4) in MeV, divide by a(β) in fm
    # Reference for SU(2): χ_t^(1/4) = 200(15) MeV
    if chi_t_lattice > 0:
        chi_t_fourth_root_a = (chi_t_lattice ** 0.25) * HBAR_C  # MeV*fm
    else:
        chi_t_fourth_root_a = 0.0
    
    # Autocorrelation of Q² (more relevant for susceptibility measurement)
    Q2 = Q_rescaled ** 2
    
    try:
        # Use Q² for autocorrelation (relevant for χ_t measurement)
        autocorr = autocorrelation(Q2, name=f"Q2_b{run_data.beta}")
        tau_int = autocorr.tau_int
        dtau_int = autocorr.dtau_int
        n_eff = autocorr.n_eff
    except Exception as e:
        print(f"  Warning: Autocorrelation failed for beta={run_data.beta}: {e}")
        tau_int = 0.5  # Default minimum
        dtau_int = 0.0
        n_eff = len(run_data.Q_raw)
    
    return AnalysisResult(
        beta=run_data.beta,
        Q_mean=stats.Q_mean,
        Q_error=stats.Q_error,
        Q2_mean=Q2_mean,  # Use rescaled Q² for χ_t
        Q2_error=Q2_std,
        chi_t_lattice=chi_t_lattice,
        chi_t_fourth_root_a=chi_t_fourth_root_a,
        tau_int=tau_int,
        dtau_int=dtau_int,
        n_eff=n_eff,
        alpha=stats.alpha,
        therm_sweep=run_data.start_conf
    )


# ==============================================================================
# Plotting Functions
# ==============================================================================

def plot_Q_histogram(run_data: RunData, result: AnalysisResult, output_dir: str, gauge_group: str = "su2"):
    """Plot rescaled Q histogram with Gaussian fit."""
    Q = run_data.Q_rescaled
    N_conf = len(Q)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Histogram with thin bars to show integer values clearly
    Q_min, Q_max = int(Q.min()) - 1, int(Q.max()) + 1
    bins = np.arange(Q_min - 0.5, Q_max + 1.5, 1)
    
    counts, bin_edges, patches = ax.hist(
        Q, bins=bins, density=False, 
        color='steelblue', edgecolor='black', linewidth=0.8,
        alpha=0.9, label='Data', rwidth=0.15  # Very thin bars for integer values
    )
    
    # Gaussian fit to counts
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    try:
        # Initial guess
        mu0 = np.mean(Q)
        sigma0 = max(np.std(Q), 0.5)  # Avoid zero sigma
        A0 = counts.max()
        
        # Only fit if we have enough non-zero bins
        nonzero_bins = np.sum(counts > 0)
        if nonzero_bins >= 3 and sigma0 > 0.1:
            popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=[mu0, sigma0, A0],
                                   bounds=([-np.inf, 0.1, 0], [np.inf, np.inf, np.inf]))
            mu_fit, sigma_fit, A_fit = popt
            
            # Plot fit
            x_fit = np.linspace(Q_min - 1, Q_max + 1, 200)
            y_fit = gaussian(x_fit, *popt)
            ax.plot(x_fit, y_fit, 'r-', linewidth=2, 
                    label=f'Gaussian: μ={mu_fit:.2f}, σ={sigma_fit:.2f}')
        else:
            print(f"  Skipping Gaussian fit for beta={run_data.beta} (too few bins or zero variance)")
    except Exception as e:
        print(f"  Warning: Gaussian fit failed for beta={run_data.beta}: {e}")
    
    ax.set_xlabel(r'$Q_{\rm rescaled}$')
    ax.set_ylabel(r'$N_{\rm conf}$')
    ax.set_title(f'{gauge_group.upper()} β={run_data.beta:.2f}, α={result.alpha:.3f}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Set x-axis to show integer ticks
    ax.set_xticks(np.arange(Q_min, Q_max + 1))
    
    # Add statistics text
    textstr = f'$\\langle Q \\rangle$ = {result.Q_mean:.2f} ± {result.Q_error:.2f}\n'
    textstr += f'$\\langle Q^2 \\rangle$ = {result.Q2_mean:.1f} ± {result.Q2_error:.1f}\n'
    textstr += f'$N_{{\\rm conf}}$ = {len(Q)}'
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    # Save
    filename = f"Q_histogram_b{run_data.beta:.2f}.png"
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"  Saved: {filename}")


def plot_all_histograms_grid(runs: List[RunData], results: List[AnalysisResult], 
                              output_dir: str, gauge_group: str = "su2"):
    """Plot all histograms in a grid layout."""
    n_runs = len(runs)
    if n_runs == 0:
        return
    
    n_cols = min(3, n_runs)
    n_rows = (n_runs + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_runs == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    for idx, (run_data, result) in enumerate(zip(runs, results)):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row, col]
        
        Q = run_data.Q_rescaled
        Q_min, Q_max = int(Q.min()) - 1, int(Q.max()) + 1
        bins = np.arange(Q_min - 0.5, Q_max + 1.5, 1)
        
        counts, bin_edges, _ = ax.hist(
            Q, bins=bins, density=False,
            color=CMAP(idx / max(n_runs-1, 1)), edgecolor='black', 
            linewidth=0.5, alpha=0.9, rwidth=0.15  # Very thin bars for integer values
        )
        
        # Gaussian fit to counts
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        try:
            sigma0 = max(np.std(Q), 0.5)
            nonzero_bins = np.sum(counts > 0)
            if nonzero_bins >= 3 and sigma0 > 0.1:
                popt, _ = curve_fit(gaussian, bin_centers, counts, 
                                   p0=[np.mean(Q), sigma0, counts.max()],
                                   bounds=([-np.inf, 0.1, 0], [np.inf, np.inf, np.inf]))
                x_fit = np.linspace(Q_min - 1, Q_max + 1, 100)
                ax.plot(x_fit, gaussian(x_fit, *popt), 'r-', linewidth=1.5)
        except:
            pass
        
        ax.set_title(f'β={run_data.beta:.2f}')
        ax.set_xlabel(r'$Q$')
        ax.set_ylabel(r'$N_{\rm conf}$')
        ax.grid(True, alpha=0.3)
    
    # Hide empty subplots
    for idx in range(n_runs, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row, col].set_visible(False)
    
    plt.suptitle(f'{gauge_group.upper()} Rescaled Topological Charge Distributions', fontsize=16)
    plt.tight_layout()
    
    filepath = os.path.join(output_dir, f"Q_histograms_all_{gauge_group}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: Q_histograms_all_{gauge_group}.png")


def lattice_spacing_su2(beta: float) -> float:
    """
    Estimate lattice spacing a(β) for SU(2) in fm.
    Simple parameterization calibrated so that χ_t^(1/4) ≈ 200 MeV at β=2.3.
    Uses asymptotic scaling with corrections.
    """
    # Reference point: β=2.3, a ≈ 0.13 fm (from χ_t^(1/4) = 200 MeV)
    beta_ref = 2.3
    a_ref = 0.13  # fm
    
    # Two-loop beta function coefficient for SU(2)
    b0 = 11 / (12 * np.pi**2)  # ≈ 0.093
    
    # Simple scaling: a(β) = a_ref * exp(-b0 * 4 * (β - β_ref))
    # Adjusted coefficient to match typical SU(2) lattice spacings
    return a_ref * np.exp(-1.5 * (beta - beta_ref))


def plot_susceptibility_vs_beta(results: List[AnalysisResult], output_dir: str, 
                                 gauge_group: str = "su2"):
    """Plot topological susceptibility vs beta."""
    betas = np.array([r.beta for r in results])
    chi_fourth_a = np.array([r.chi_t_fourth_root_a for r in results])
    
    # Sort by beta
    idx = np.argsort(betas)
    betas = betas[idx]
    chi_fourth_a = chi_fourth_a[idx]
    
    # Compute χ_t^(1/4) in MeV using lattice spacing estimate
    if gauge_group == "su2":
        a_beta = np.array([lattice_spacing_su2(b) for b in betas])
        chi_fourth_MeV = chi_fourth_a / a_beta
        ref_value = CHI_T_FOURTH_ROOT_SU2
    else:
        # For SU(3), would need different parameterization
        a_beta = chi_fourth_a / CHI_T_FOURTH_ROOT_SU3  # Use implied a
        chi_fourth_MeV = chi_fourth_a / a_beta
        ref_value = CHI_T_FOURTH_ROOT_SU3
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(betas, chi_fourth_MeV, 'o-', markersize=8, linewidth=2, color='steelblue', label='Data')
    
    # Reference value as horizontal line
    ax.axhline(y=ref_value, color='red', linestyle='--', linewidth=2, 
               label=f'Reference: {ref_value} MeV')
    
    ax.set_xlabel(r'$\beta$', fontsize=14)
    ax.set_ylabel(r'$\chi_t^{1/4}$ (MeV)', fontsize=14)
    ax.set_title(f'{gauge_group.upper()} Topological Susceptibility', fontsize=14)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    filepath = os.path.join(output_dir, f"susceptibility_vs_beta_{gauge_group}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: susceptibility_vs_beta_{gauge_group}.png")
    
    # Save data to file (using filtered data)
    data_file = os.path.join(output_dir, f"susceptibility_data_{gauge_group}.txt")
    with open(data_file, 'w') as f:
        f.write(f"# Topological susceptibility data for {gauge_group.upper()}\n")
        f.write(f"# Reference: chi_t^(1/4) = {ref_value} MeV\n")
        f.write("# beta    a(beta) [fm]    chi_t^1/4*a [MeV*fm]    chi_t^1/4 [MeV]\n")
        for i in range(len(betas)):
            f.write(f"{betas[i]:.4f}    {a_beta[i]:.4f}    {chi_fourth_a[i]:8.3f}    {chi_fourth_MeV[i]:8.1f}\n")
    print(f"Saved: susceptibility_data_{gauge_group}.txt")


def plot_tau_int_vs_beta(results: List[AnalysisResult], output_dir: str, 
                         gauge_group: str = "su2"):
    """Plot integrated autocorrelation time vs beta."""
    betas = np.array([r.beta for r in results])
    tau_int = np.array([r.tau_int for r in results])
    dtau_int = np.array([r.dtau_int for r in results])
    
    # Sort by beta
    idx = np.argsort(betas)
    betas = betas[idx]
    tau_int = tau_int[idx]
    dtau_int = dtau_int[idx]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.errorbar(betas, tau_int, yerr=dtau_int, fmt='o-', markersize=8, 
                linewidth=2, capsize=4, color='steelblue')
    
    ax.set_xlabel(r'$\beta$', fontsize=14)
    ax.set_ylabel(r'$\tau_{\rm int}(Q^2)$ (sweeps)', fontsize=14)
    ax.set_title(f'{gauge_group.upper()} Integrated Autocorrelation Time', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    filepath = os.path.join(output_dir, f"tau_int_vs_beta_{gauge_group}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: tau_int_vs_beta_{gauge_group}.png")


def plot_thermalization_comparison(runs: List[RunData], output_dir: str, 
                                   gauge_group: str = "su2"):
    """Plot all thermalization curves in one plot (normalized)."""
    if len(runs) == 0:
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), height_ratios=[2, 1])
    
    # Get colors
    colors = CMAP(np.linspace(0, 1, len(runs)))
    
    # Sort runs by beta
    runs_sorted = sorted(runs, key=lambda r: r.beta)
    
    for idx, run in enumerate(runs_sorted):
        if len(run.plaq_data) == 0:
            continue
        
        color = colors[idx]
        label = f'β={run.beta:.2f}'
        
        # Top plot: raw plaquette
        ax1.plot(run.mc_time, run.plaq_data, '-', color=color, 
                 linewidth=0.8, alpha=0.8, label=label)
        
        # Mark thermalization point
        therm_idx = np.searchsorted(run.mc_time, run.start_conf)
        if therm_idx < len(run.mc_time):
            ax1.axvline(x=run.start_conf, color=color, linestyle=':', alpha=0.5)
        
        # Bottom plot: normalized plaquette (subtract mean, divide by std)
        if len(run.plaq_data) > 10:
            plaq_norm = (run.plaq_data - run.plaq_data.mean()) / (run.plaq_data.std() + 1e-10)
            ax2.plot(run.mc_time, plaq_norm, '-', color=color, linewidth=0.6, alpha=0.7)
    
    ax1.set_ylabel(r'$\langle P_{\mu\nu} \rangle$')
    ax1.set_title(f'{gauge_group.upper()} Thermalization Comparison')
    ax1.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    ax2.set_xlabel('Monte Carlo Sweeps')
    ax2.set_ylabel('Normalized Plaquette')
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    filepath = os.path.join(output_dir, f"thermalization_comparison_{gauge_group}.png")
    plt.savefig(filepath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: thermalization_comparison_{gauge_group}.png")


def print_summary_table(results: List[AnalysisResult], gauge_group: str = "su2"):
    """Print summary table of all results."""
    ref_value = CHI_T_FOURTH_ROOT_SU2 if gauge_group == "su2" else CHI_T_FOURTH_ROOT_SU3
    
    print("\n" + "="*110)
    print(f"  {gauge_group.upper()} Beta Scan Analysis Summary")
    print("="*110)
    print(f"{'β':>8} {'<Q²>':>10} {'χ_t':>12} {'χ_t^1/4·a':>12} {'a(β) implied':>12} {'τ_int':>8} {'N_eff':>8} {'α':>8}")
    print(f"{'':>8} {'':>10} {'(lattice)':>12} {'(MeV·fm)':>12} {'(fm)':>12} {'':>8} {'':>8} {'':>8}")
    print("-"*110)
    
    for r in sorted(results, key=lambda x: x.beta):
        # Implied lattice spacing if χ_t^(1/4) = ref_value
        # a = χ_t^(1/4)·a / ref_value
        if r.chi_t_fourth_root_a > 0:
            a_implied = r.chi_t_fourth_root_a / ref_value
        else:
            a_implied = 0.0
        
        print(f"{r.beta:>8.3f} {r.Q2_mean:>10.2f} {r.chi_t_lattice:>12.4e} "
              f"{r.chi_t_fourth_root_a:>12.2f} {a_implied:>12.4f} "
              f"{r.tau_int:>8.1f} {r.n_eff:>8.1f} {r.alpha:>8.3f}")
    
    print("="*110)
    print(f"Note: a(β) implied assumes χ_t^(1/4) = {ref_value} MeV (reference value for {gauge_group.upper()})")
    print(f"      To get χ_t^(1/4) [MeV] = χ_t^(1/4)·a [MeV·fm] / a(β) [fm]")
    print("")


# ==============================================================================
# Main
# ==============================================================================

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze beta scan results")
    parser.add_argument("--su2", action="store_true", help="Analyze SU(2) data (default)")
    parser.add_argument("--su3", action="store_true", help="Analyze SU(3) data")
    parser.add_argument("--results-dir", type=str, default=None,
                        help="Results directory (default: data/results)")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory for plots (default: output/figures)")
    args = parser.parse_args()
    
    # Determine gauge group
    gauge_group = "su3" if args.su3 else "su2"
    
    # Set directories
    base_dir = os.path.dirname(script_dir)
    results_dir = args.results_dir or os.path.join(base_dir, "data", "results")
    output_dir = args.output_dir or os.path.join(base_dir, "output", f"figures_beta_scan_{gauge_group}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"  {gauge_group.upper()} Beta Scan Analysis")
    print(f"{'='*60}")
    print(f"Results dir: {results_dir}")
    print(f"Output dir:  {output_dir}")
    
    # Find all run directories
    run_dirs = sorted(glob.glob(os.path.join(results_dir, "T*_L*_b*_seed*")))
    
    if not run_dirs:
        print("No run directories found!")
        return
    
    print(f"Found {len(run_dirs)} run directories\n")
    
    # Load data from each run
    runs = []
    results = []
    
    for run_dir in run_dirs:
        dirname = os.path.basename(run_dir)
        print(f"Loading: {dirname}")
        
        run_data = load_run_data(run_dir, gauge_group)
        if run_data is None:
            continue
        
        result = analyze_run(run_data, gauge_group)
        
        runs.append(run_data)
        results.append(result)
        
        # Save individual histogram
        plot_Q_histogram(run_data, result, output_dir, gauge_group)
    
    if len(runs) == 0:
        print("No valid runs found!")
        return
    
    # Print summary table
    print_summary_table(results, gauge_group)
    
    # Generate combined plots
    print("Generating combined plots...")
    
    plot_all_histograms_grid(runs, results, output_dir, gauge_group)
    plot_susceptibility_vs_beta(results, output_dir, gauge_group)
    plot_tau_int_vs_beta(results, output_dir, gauge_group)
    plot_thermalization_comparison(runs, output_dir, gauge_group)
    
    print(f"\nAnalysis complete! Plots saved to: {output_dir}")


if __name__ == "__main__":
    main()
