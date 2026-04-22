#!/usr/bin/env python3
"""
Plotting functions for SU(2)/SU(3) topological charge analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
from typing import List

from calculations import (
    RunData, AnalysisResult, gaussian,
    lattice_spacing_su2, lattice_spacing_su3,
    CHI_T_FOURTH_ROOT_SU2, CHI_T_FOURTH_ROOT_SU3,
    find_optimal_alpha
)

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'figure.dpi': 150,
})


def _smooth(y: np.ndarray, window: int = None) -> np.ndarray:
    """Rolling mean for a guiding line. Window defaults to ~10% of length, min 3."""
    w = max(3, len(y) // 10) if window is None else window
    return np.convolve(y, np.ones(w) / w, mode='same')


def plot_Q_vs_mctime_grid(runs: List[RunData], output_dir: str, gauge_group: str, boundary: str):
    """Plot Q vs MC time grid for one boundary type."""
    n = len(runs)
    if n == 0:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    n_cols = 1
    n_rows = n

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7, 3 * n_rows))
    axes = np.array(axes).reshape(-1, 1)

    sorted_runs = sorted(runs, key=lambda x: x.beta)

    # Global y-range across all runs so panels are comparable.
    # Use floor/ceil (not int, which truncates toward zero and clips negative open-BC values).
    raw_min = min(r.Q_rescaled.min() for r in sorted_runs)
    raw_max = max(r.Q_rescaled.max() for r in sorted_runs)
    global_min = int(np.floor(raw_min))
    global_max = int(np.ceil(raw_max))
    pad = 0.5

    for idx, run_data in enumerate(sorted_runs):
        ax = axes[idx, 0]

        Q = run_data.Q_rescaled
        mc_time = np.arange(len(Q))

        ax.plot(mc_time, Q.astype(float), color='steelblue', lw=0.8, alpha=0.25, zorder=2)
        ax.scatter(mc_time, Q, s=3, color='steelblue', alpha=0.8, zorder=3)

        for q_int in range(global_min, global_max + 1):
            ax.axhline(y=q_int, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

        ax.set_ylim(global_min - pad, global_max + pad)

        a = lattice_spacing(run_data.beta)
        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel('MC time', fontsize=9)
        ax.set_ylabel(r'$Q_{\rm re}$', fontsize=9)
        ax.tick_params(labelsize=8)
        ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                fontsize=11, va='top', ha='left')
    
    plt.tight_layout()
    filepath = os.path.join(output_dir, f"Q_vs_mctime_{gauge_group}_{boundary}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def plot_histograms_grid(runs: List[RunData], output_dir: str, gauge_group: str, boundary: str):
    """Plot histograms grid for one boundary type."""
    n = len(runs)
    if n == 0:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    n_cols = min(3, n)
    n_rows = (n + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3.5*n_rows))
    if n == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    sorted_runs = sorted(runs, key=lambda x: x.beta)

    # Global x-range across all runs so panels are comparable.
    # For periodic, Q_rescaled is already rounded; for open, Q_rescaled is continuous.
    # Include the unrounded alpha*Q_raw in the range so rug lines are never clipped.
    q_mins = [min(r.Q_rescaled.min(), (r.alpha * r.Q_raw).min()) for r in sorted_runs]
    q_maxs = [max(r.Q_rescaled.max(), (r.alpha * r.Q_raw).max()) for r in sorted_runs]
    global_min = int(np.floor(min(q_mins))) - 1
    global_max = int(np.ceil(max(q_maxs))) + 1

    int_bins  = np.arange(global_min - 0.5, global_max + 1.5, 1)   # periodic
    fine_bins = np.linspace(global_min - 0.5, global_max + 0.5, 40)  # open

    for idx, run_data in enumerate(sorted_runs):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row, col]

        Q_cont = run_data.alpha * run_data.Q_raw

        if boundary == 'periodic':
            Q_rounded = run_data.Q_rescaled
            rounded_label = r'$Q_{\rm rounded}$' if idx == 0 else None
            counts, bin_edges, _ = ax.hist(Q_rounded, bins=int_bins, density=False,
                color='steelblue', edgecolor='steelblue', linewidth=0.5, alpha=0.3,
                label=rounded_label)

            # Fine-binned histogram of unrounded alpha*Q_raw — same style as open
            cont_label = r'$\alpha \hat{Q}$ (unrounded)' if idx == 0 else None
            ax.hist(Q_cont, bins=fine_bins, density=False,
                color='steelblue', edgecolor='steelblue', linewidth=1.2, alpha=0.7,
                label=cont_label)

            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            fit_x_min, fit_x_max = global_min - 1, global_max + 1
            fit_Q_for_mean = Q_rounded
        else:  # open — continuous Q, no rounding
            cont_label = r'$\alpha \hat{Q}$' if idx == 0 else None
            counts, bin_edges, _ = ax.hist(Q_cont, bins=fine_bins, density=False,
                color='steelblue', edgecolor='steelblue', linewidth=1.2, alpha=0.7,
                label=cont_label)

            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            fit_x_min, fit_x_max = global_min - 1, global_max + 1
            fit_Q_for_mean = Q_cont

        try:
            sigma0 = max(np.std(fit_Q_for_mean), 0.5)
            if np.sum(counts > 0) >= 3 and sigma0 > 0.1:
                popt, _ = curve_fit(gaussian, bin_centers, counts,
                    p0=[np.mean(fit_Q_for_mean), sigma0, counts.max()],
                    bounds=([-np.inf, 0.1, 0], [np.inf, np.inf, np.inf]))
                x_fit = np.linspace(fit_x_min, fit_x_max, 200)
                fit_label = 'Gaussian fit' if idx == 0 else None
                ax.plot(x_fit, gaussian(x_fit, *popt), 'r-', linewidth=1.5, label=fit_label)
        except:
            pass

        ax.set_xlim(global_min - 0.5, global_max + 0.5)

        a = lattice_spacing(run_data.beta)
        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel(r'$Q$', fontsize=9)
        ax.set_ylabel(r'$N$', fontsize=9)
        ax.tick_params(labelsize=8)

        # Add panel label (A, B, C...)
        panel_label = chr(ord('A') + idx)
        ax.text(0.02, 0.98, panel_label, transform=ax.transAxes, fontsize=11,
                va='top', ha='left')
    
    for idx in range(n, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row, col].set_visible(False)
    
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.78, 0.98), frameon=False, fontsize=10)
    
    plt.tight_layout(rect=[0, 0, 0.77, 1])
    filepath = os.path.join(output_dir, f"histograms_{gauge_group}_{boundary}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def plot_susceptibility_combined(results_periodic: List[AnalysisResult], 
                                  results_open: List[AnalysisResult],
                                  output_dir: str, gauge_group: str):
    """Plot susceptibility vs lattice spacing."""
    fig, ax = plt.subplots(figsize=(8, 5))
    
    ref_value = CHI_T_FOURTH_ROOT_SU2 if gauge_group == "su2" else CHI_T_FOURTH_ROOT_SU3
    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    
    point_counter = 1
    for results, label, marker, color in [
        (results_periodic, 'Periodic', 'o', 'steelblue'),
        (results_open, 'Open', 's', 'darkorange')
    ]:
        if len(results) == 0:
            continue
        
        betas = np.array([r.beta for r in results])
        chi_fourth_a = np.array([r.chi_t_fourth_root_a for r in results])
        a_values = np.array([lattice_spacing(b) for b in betas])
        
        idx = np.argsort(a_values)[::-1]
        a_values, chi_fourth_a = a_values[idx], chi_fourth_a[idx]
        
        chi_fourth_MeV = chi_fourth_a / a_values
        
        ax.plot(a_values, chi_fourth_MeV, marker, markersize=7, linestyle='none',
                color=color, label=label)
        
        # Add numeric labels to each point
        for i, (x, y) in enumerate(zip(a_values, chi_fourth_MeV)):
            ax.annotate(str(point_counter), (x, y), textcoords='offset points',
                        xytext=(3, 3), fontsize=7, color=color)
            point_counter += 1
    
    ax.axhline(y=ref_value, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    ax.set_xlabel(r'$a$ (fm)')
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
    """Plot tau_int vs lattice spacing."""
    fig, ax = plt.subplots(figsize=(8, 5))
    
    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    
    point_counter = 1
    for results, label, marker, color in [
        (results_periodic, 'Periodic', 'o', 'steelblue'),
        (results_open, 'Open', 's', 'darkorange')
    ]:
        if len(results) == 0:
            continue
        
        betas = np.array([r.beta for r in results])
        tau_int = np.array([r.tau_int for r in results])
        dtau_int = np.array([r.dtau_int for r in results])
        a_values = np.array([lattice_spacing(b) for b in betas])
        
        idx = np.argsort(a_values)[::-1]
        a_values, tau_int, dtau_int = a_values[idx], tau_int[idx], dtau_int[idx]
        
        ax.errorbar(a_values, tau_int, yerr=dtau_int, fmt=marker, markersize=7,
                    linestyle='none', capsize=3, color=color, label=label)
        
        # Add numeric labels to each point
        for i, (x, y) in enumerate(zip(a_values, tau_int)):
            ax.annotate(str(point_counter), (x, y), textcoords='offset points',
                        xytext=(3, 3), fontsize=7, color=color)
            point_counter += 1
    
    ax.set_xlabel(r'$a$ (fm)')
    ax.set_ylabel(r'$\tau_{\rm int}(Q^2)$')
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    filepath = os.path.join(output_dir, f"tau_int_{gauge_group}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def plot_thermalization_comparison(runs: List[RunData], output_dir: str, gauge_group: str):
    """Plot all thermalization curves (plaquette vs MC time) with detected thermalization points."""
    runs_with_plaq = [r for r in runs if len(r.plaq_data) > 0]
    if len(runs_with_plaq) == 0:
        return
    
    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    
    fig, ax = plt.subplots(figsize=(10, 5))
    
    cmap = plt.cm.viridis
    colors = [cmap(i / max(1, len(runs_with_plaq) - 1)) for i in range(len(runs_with_plaq))]
    
    runs_sorted = sorted(runs_with_plaq, key=lambda r: r.beta)
    
    for idx, run in enumerate(runs_sorted):
        color = colors[idx]
        a = lattice_spacing(run.beta)
        label = f'$a = {a:.4f}$ fm ({run.boundary[:4]})'
        
        # Plot plaquette vs MC time
        ax.plot(run.mc_time, run.plaq_data, '-', color=color, linewidth=0.8, alpha=0.8, label=label)
        
        # Mark thermalization point
        if run.start_conf > 0 and run.start_conf < run.mc_time[-1]:
            ax.axvline(x=run.start_conf, color=color, linestyle='--', alpha=0.8, linewidth=1.5)
            ax.scatter([run.start_conf], [run.plaq_data[np.searchsorted(run.mc_time, run.start_conf)]], 
                       color=color, s=50, zorder=5, marker='o', edgecolors='black', linewidth=0.5)
    
    ax.set_xlabel('MC Sweeps')
    ax.set_ylabel(r'$\langle P_{\mu\nu} \rangle$')
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    filepath = os.path.join(output_dir, f"thermalization_comparison_{gauge_group}.png")
    plt.savefig(filepath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")
