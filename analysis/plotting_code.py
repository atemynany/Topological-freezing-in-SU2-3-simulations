#!/usr/bin/env python3
"""
Plotting functions for SU(2)/SU(3) topological charge analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
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


def _apply_axis_style(ax, x_integer: bool = False, y_integer: bool = False,
                      x_nbins: int = 6, y_nbins: int = 6):
    """Apply a uniform axis style across all plots."""
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.tick_params(labelsize=8)
    ax.minorticks_off()
    if x_integer:
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=x_nbins, integer=True, min_n_ticks=4))
    else:
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=x_nbins, min_n_ticks=4))
    if y_integer:
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=y_nbins, integer=True, min_n_ticks=4))
    else:
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=y_nbins, min_n_ticks=4))


def _legend_outside(ax, handles=None, labels=None, ncol: int = 1, fontsize: int = 8):
    """Place a minimal legend outside the plotting area."""
    if handles is None or labels is None:
        handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None
    return ax.legend(
        handles,
        labels,
        loc='upper left',
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0,
        frameon=False,
        fontsize=fontsize,
        ncol=ncol,
        handlelength=1.4,
        labelspacing=0.3,
    )


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

    # Global y-range across all runs so panels are comparable
    global_min = int(min(r.Q_rescaled.min() for r in sorted_runs))
    global_max = int(max(r.Q_rescaled.max() for r in sorted_runs))
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
        _apply_axis_style(ax, x_integer=True, y_integer=True)
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

    # Global x-range across all runs so panels are comparable
    global_min = int(min(r.Q_rescaled.min() for r in sorted_runs)) - 1
    global_max = int(max(r.Q_rescaled.max() for r in sorted_runs)) + 1
    bins = np.arange(global_min - 0.5, global_max + 1.5, 1)

    for idx, run_data in enumerate(sorted_runs):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row, col]

        Q = run_data.Q_rescaled

        hist_label = 'Data' if idx == 0 else None
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
                x_fit = np.linspace(global_min - 1, global_max + 1, 200)
                fit_label = 'Gaussian' if idx == 0 else None
                ax.plot(x_fit, gaussian(x_fit, *popt), 'r-', linewidth=1.5, label=fit_label)
        except:
            pass

        ax.set_xlim(global_min - 0.5, global_max + 0.5)

        a = lattice_spacing(run_data.beta)
        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel(r'$Q$', fontsize=9)
        ax.set_ylabel(r'$N$', fontsize=9)
        _apply_axis_style(ax, x_integer=True, y_integer=True)

        # Add panel label (A, B, C...)
        panel_label = chr(ord('A') + idx)
        ax.text(0.02, 0.98, panel_label, transform=ax.transAxes, fontsize=11,
                va='top', ha='left')
    
    for idx in range(n, n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        axes[row, col].set_visible(False)
    
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.82, 0.98), frameon=False,
                   fontsize=8, handlelength=1.4, labelspacing=0.3)
    
    plt.tight_layout(rect=[0, 0, 0.80, 1])
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
    ax.axhline(y=220, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    ax.set_xlabel(r'$a$ (fm)')
    ax.set_ylabel(r'$\chi_t^{1/4}$ (MeV)')
    _apply_axis_style(ax)
    _legend_outside(ax, fontsize=8)
    
    plt.tight_layout(rect=[0, 0, 0.84, 1])
    filepath = os.path.join(output_dir, f"susceptibility_{gauge_group}.png")
    plt.savefig(filepath, dpi=200, bbox_inches='tight')
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
    _apply_axis_style(ax)
    _legend_outside(ax, fontsize=8)
    
    plt.tight_layout(rect=[0, 0, 0.84, 1])
    filepath = os.path.join(output_dir, f"tau_int_{gauge_group}.png")
    plt.savefig(filepath, dpi=200, bbox_inches='tight')
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
        # Plot plaquette vs MC time
        ax.plot(run.mc_time, run.plaq_data, '-', color=color, linewidth=0.8, alpha=0.8)
        
        # Mark thermalization point
        if run.start_conf > 0 and run.start_conf < run.mc_time[-1]:
            ax.axvline(x=run.start_conf, color=color, linestyle='--', alpha=0.8, linewidth=1.5)
            ax.scatter([run.start_conf], [run.plaq_data[np.searchsorted(run.mc_time, run.start_conf)]], 
                       color=color, s=50, zorder=5, marker='o', edgecolors='black', linewidth=0.5)
    
    ax.set_xlabel('MC Sweeps')
    ax.set_ylabel(r'$\langle P_{\mu\nu} \rangle$')
    _apply_axis_style(ax, x_integer=True)
    legend_handles = [
        Line2D([0], [0], color='steelblue', lw=1.2, label='Plaquette'),
        Line2D([0], [0], color='black', lw=1.0, ls='--', marker='o', markersize=4,
               label='start conf'),
    ]
    _legend_outside(ax, handles=legend_handles,
                    labels=[h.get_label() for h in legend_handles], fontsize=8)
    
    plt.tight_layout(rect=[0, 0, 0.84, 1])
    filepath = os.path.join(output_dir, f"thermalization_comparison_{gauge_group}.png")
    plt.savefig(filepath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def plot_therm_topcharge_grid(runs: List[RunData], output_dir: str, gauge_group: str, boundary: str):
    """Grid of thermalization Q vs MC sweep, one panel per run. A,B,C labels + a title."""
    therm_fname = "therm_topcharge_su3.dat" if gauge_group == "su3" else "therm_topcharge.dat"
    runs_with_data = [r for r in sorted(runs, key=lambda x: x.beta)
                      if os.path.isfile(os.path.join(r.run_dir, "output", therm_fname))]
    n = len(runs_with_data)
    if n == 0:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    # Global y-range across all runs
    all_q = []
    for run_data in runs_with_data:
        fp = os.path.join(run_data.run_dir, "output", therm_fname)
        d = np.loadtxt(fp, comments='#')
        all_q.append(d[:, 1])
    global_min = int(np.floor(min(q.min() for q in all_q)))
    global_max = int(np.ceil(max(q.max() for q in all_q)))

    fig, axes = plt.subplots(n, 1, figsize=(7, 3 * n))
    axes = np.array(axes).flatten()

    for idx, (run_data, q_raw) in enumerate(zip(runs_with_data, all_q)):
        ax = axes[idx]

        q_data = np.loadtxt(os.path.join(run_data.run_dir, "output", therm_fname), comments='#')
        q_sweeps = q_data[:, 0].astype(int)

        ax.plot(q_sweeps, q_raw, color='steelblue', lw=0.8, alpha=0.25, zorder=2)
        ax.scatter(q_sweeps, q_raw, s=3, color='steelblue', alpha=0.8, zorder=3)

        for q_int in range(global_min, global_max + 1):
            ax.axhline(y=q_int, color='gray', ls='--', lw=0.5, alpha=0.5)

        ax.set_ylim(global_min - 0.5, global_max + 0.5)
        a = lattice_spacing(run_data.beta)
        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel('MC Sweep', fontsize=9)
        ax.set_ylabel(r'$Q$', fontsize=9)
        _apply_axis_style(ax, x_integer=True, y_integer=True)
        ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                fontsize=11, va='top', ha='left')

    plt.tight_layout()
    filepath = os.path.join(output_dir, f"therm_topcharge_{gauge_group}_{boundary}.png")
    plt.savefig(filepath, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(filepath)}")


def plot_therm_topcharge(plaq_file: str, therm_q_file: str, output_dir: str,
                         run_name: str = "", gauge_group: str = ""):
    """Thermalization plot: rescaled Q_re vs MC sweep for a single run.

    Plaquette thermalization is handled by plot_thermalization_comparison.
    Saves one file: therm_topcharge[_<run_name>].png
    """
    os.makedirs(output_dir, exist_ok=True)
    gg_suffix = f"_{gauge_group}" if gauge_group else ""
    suffix = f"{gg_suffix}_{run_name}" if run_name else gg_suffix

    # --- load raw thermalization topcharge ---
    try:
        q_data   = np.loadtxt(therm_q_file, comments='#')
        q_sweeps = q_data[:, 0].astype(int)
        q_raw    = q_data[:, 1]
    except Exception:
        return  # nothing to plot

    if len(q_raw) == 0:
        return

    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(q_sweeps, q_raw, color='steelblue', lw=0.8, alpha=0.25, zorder=2)
    ax.scatter(q_sweeps, q_raw, s=3, color='steelblue', alpha=0.8, zorder=3)
    q_min_int = int(np.floor(q_raw.min()))
    q_max_int = int(np.ceil(q_raw.max()))
    for q_int in range(q_min_int, q_max_int + 1):
        ax.axhline(y=q_int, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.set_xlabel('MC Sweep')
    ax.set_ylabel(r'$Q$ (unsmeared)')
    _apply_axis_style(ax, x_integer=True, y_integer=True)
    plt.tight_layout()
    path = os.path.join(output_dir, f'therm_topcharge{suffix}.png')
    plt.savefig(path, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(path)}")
