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
    CHI_T_FOURTH_ROOT_SU2, CHI_T_FOURTH_ROOT_SU3
)

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'figure.dpi': 150,
})


def plot_Q_vs_mctime_grid(runs: List[RunData], output_dir: str, gauge_group: str, boundary: str):
    """Plot Q vs MC time grid for one boundary type."""
    n = len(runs)
    if n == 0:
        return
    
    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    
    n_cols = min(3, n)
    n_rows = (n + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows))
    if n == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    sorted_runs = sorted(runs, key=lambda x: x.beta)
    
    for idx, run_data in enumerate(sorted_runs):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row, col]
        
        Q = run_data.Q_rescaled
        mc_time = np.arange(len(Q))
        
        ax.scatter(mc_time, Q, s=3, color='steelblue', alpha=0.7)
        
        # Add dashed lines at integer Q values
        Q_min_int, Q_max_int = int(np.floor(Q.min())), int(np.ceil(Q.max()))
        for q_int in range(Q_min_int, Q_max_int + 1):
            ax.axhline(y=q_int, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        
        a = lattice_spacing(run_data.beta)
        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel('MC time', fontsize=9)
        ax.set_ylabel(r'$Q_{\rm re}$', fontsize=9)
        ax.tick_params(labelsize=8)
        
        # Add panel label (A, B, C...)
        panel_label = chr(ord('A') + idx)
        ax.text(0.02, 0.98, panel_label, transform=ax.transAxes, fontsize=11,
                va='top', ha='left')
    
    for idx in range(n, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row, col].set_visible(False)
    
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
    
    for idx, run_data in enumerate(sorted_runs):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row, col]
        
        Q = run_data.Q_rescaled
        Q_min, Q_max = int(Q.min()) - 1, int(Q.max()) + 1
        bins = np.arange(Q_min - 0.5, Q_max + 1.5, 1)
        
        hist_label = r'$Q_{\rm re}$' if idx == 0 else None
        counts, bin_edges, _ = ax.hist(Q, bins=bins, density=False,
            color='steelblue', edgecolor='black', linewidth=0.5, alpha=0.9, rwidth=0.15,
            label=hist_label)
        
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        try:
            sigma0 = max(np.std(Q), 0.5)
            if np.sum(counts > 0) >= 3 and sigma0 > 0.1:
                popt, _ = curve_fit(gaussian, bin_centers, counts,
                    p0=[np.mean(Q), sigma0, counts.max()],
                    bounds=([-np.inf, 0.1, 0], [np.inf, np.inf, np.inf]))
                x_fit = np.linspace(Q_min - 1, Q_max + 1, 100)
                fit_label = 'Gaussian fit' if idx == 0 else None
                ax.plot(x_fit, gaussian(x_fit, *popt), 'r-', linewidth=1.5, label=fit_label)
        except:
            pass
        
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
    
    ax.axhline(y=ref_value, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax.invert_xaxis()
    
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
    
    ax.invert_xaxis()
    
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
