#!/usr/bin/env python3
"""
Topological charge density analysis per time slice.

Loads topcharge_timeslice.dat produced by meas_topcharge_su2,
optionally bins over n_bin consecutive time slices to reduce noise
(useful when the spatial volume is small), computes errors via the
Gamma method (pyerrors), and plots q(t) vs physical time.

File format expected (written by meas_topcharge_su2.cc):
  # smear_steps  config_number  t  q_t
  40  100  0  0.00123456
  ...
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pyerrors as pe
from pathlib import Path
from typing import Optional
from autocorrelation import autocorrelation


def _ensure_endpoint_ticks(ax):
    """Ensure the x-axis endpoints are always labeled."""
    import matplotlib.ticker as mticker
    xlim = ax.get_xlim()
    ticks = list(ax.get_xticks())
    for ep in [xlim[0], xlim[1]]:
        rounded = round(ep)
        if not any(abs(t - rounded) < 0.1 for t in ticks):
            ticks.append(rounded)
    ticks = sorted(set(int(round(t)) for t in ticks if xlim[0] - 0.1 <= t <= xlim[1] + 0.1))
    ax.set_xticks(ticks)


def _apply_axis_style(ax, x_integer: bool = False, y_integer: bool = False,
                      x_nbins: int = 6, y_nbins: int = 6):
    """Apply a uniform axis style across plots."""
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


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_timeslice_data(filepath: str, smear: int = None) -> tuple[np.ndarray, np.ndarray]:
    """Load topcharge_timeslice.dat.

    Returns
    -------
    t_values : 1-D int array of time slice indices (unique, sorted)
    q_matrix : 2-D float array, shape (n_configs, n_timeslices)
               q_matrix[i, j] = q(t_j) for config i
    """
    smear_col, conf_col, t_col, q_col = [], [], [], []

    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                smear_col.append(int(parts[0]))
                conf_col.append(int(parts[1]))
                t_col.append(int(parts[2]))
                q_col.append(float(parts[3]))
            except ValueError:
                continue

    smear_col = np.array(smear_col)
    conf_col  = np.array(conf_col)
    t_col     = np.array(t_col)
    q_col     = np.array(q_col)

    # Select smearing level (default: maximum)
    if smear is None:
        smear = int(smear_col.max())
    mask = smear_col == smear
    conf_col, t_col, q_col = conf_col[mask], t_col[mask], q_col[mask]

    configs   = np.unique(conf_col)
    t_values  = np.unique(t_col)
    n_configs = len(configs)
    n_t       = len(t_values)

    q_matrix = np.zeros((n_configs, n_t))
    for i, c in enumerate(configs):
        for j, t in enumerate(t_values):
            sel = (conf_col == c) & (t_col == t)
            if sel.any():
                q_matrix[i, j] = q_col[sel][0]

    return t_values, q_matrix


# ---------------------------------------------------------------------------
# Binning
# ---------------------------------------------------------------------------

def bin_timeslices(t_values: np.ndarray, q_matrix: np.ndarray,
                   n_bin: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Sum q over consecutive windows of n_bin time slices.

    Useful when individual slices are too noisy (small spatial volume).
    Returns bin-centre t values and summed q per bin.

    Parameters
    ----------
    t_values : (n_t,) array of time indices
    q_matrix : (n_configs, n_t) array
    n_bin    : number of consecutive slices to sum (1 = no binning)
    """
    if n_bin <= 1:
        return t_values.astype(float), q_matrix

    n_t       = len(t_values)
    n_bins    = n_t // n_bin
    t_centres = np.array([t_values[i*n_bin : i*n_bin + n_bin].mean()
                          for i in range(n_bins)])
    q_binned  = np.array([q_matrix[:, i*n_bin : i*n_bin + n_bin].sum(axis=1)
                          for i in range(n_bins)]).T   # shape (n_configs, n_bins)
    return t_centres, q_binned


# ---------------------------------------------------------------------------
# Error estimation via Gamma method
# ---------------------------------------------------------------------------

def gamma_method_timeslices(t_centres: np.ndarray, q_binned: np.ndarray,
                             ensemble: str = "ens", S: float = 1.5
                             ) -> tuple[np.ndarray, np.ndarray]:
    """Apply pyerrors Gamma method independently to each time bin.

    Parameters
    ----------
    t_centres : (n_bins,) bin-centre time values
    q_binned  : (n_configs, n_bins) charge per bin per config
    ensemble  : pyerrors ensemble name
    S         : Gamma method window parameter (higher = more conservative)

    Returns
    -------
    means  : (n_bins,) mean q per bin
    errors : (n_bins,) Gamma-method error on the mean
    """
    n_bins = q_binned.shape[1]
    means  = np.zeros(n_bins)
    errors = np.zeros(n_bins)

    for j in range(n_bins):
        obs = pe.Obs([q_binned[:, j]], [ensemble])
        obs.gamma_method(S=S)
        means[j]  = obs.value
        errors[j] = obs.dvalue

    return means, errors


# ---------------------------------------------------------------------------
# Full pipeline
# ---------------------------------------------------------------------------

def analyse_timeslices(filepath: str,
                       lattice_spacing_fm: float,
                       n_bin: int = 1,
                       smear: int = None,
                       ensemble: str = "ens",
                       S: float = 1.5,
                       open_bc: bool = False,
                       n_exclude: int = 2) -> dict:
    """Load, bin, and compute Gamma-method errors for timeslice charge density.

    Parameters
    ----------
    filepath           : path to topcharge_timeslice.dat
    lattice_spacing_fm : a in fm (to convert t_lattice -> t_physical)
    n_bin              : number of adjacent time slices to sum per bin
    smear              : smearing level to select (None = max)
    ensemble           : pyerrors ensemble label
    S                  : Gamma method S parameter
    open_bc            : if True, shade excluded boundary regions in plots
    n_exclude          : number of slices excluded per boundary (for shading only)

    Returns
    -------
    dict with keys: t_phys, means, errors, t_centres_latt
    """
    t_values, q_matrix = load_timeslice_data(filepath, smear=smear)
    t_centres, q_binned = bin_timeslices(t_values, q_matrix, n_bin=n_bin)
    means, errors = gamma_method_timeslices(t_centres, q_binned, ensemble=ensemble, S=S)

    t_phys = t_centres * lattice_spacing_fm  # physical time in fm

    # Per-timeslice susceptibility: chi_tilde(t) = <q(t)^2> / V_tilde
    # V_tilde = L^3 * n_bin (spatial volume times number of binned slices)
    n_bins = q_binned.shape[1]
    chi_tilde = np.zeros(n_bins)
    chi_tilde_err = np.zeros(n_bins)
    for j in range(n_bins):
        q2_obs = pe.Obs([q_binned[:, j]**2], [ensemble])
        q2_obs.gamma_method(S=S)
        chi_tilde[j] = q2_obs.value
        chi_tilde_err[j] = q2_obs.dvalue

    # Per-timeslice integrated autocorrelation time
    tau_int_ts = np.zeros(n_bins)
    dtau_int_ts = np.zeros(n_bins)
    for j in range(n_bins):
        try:
            ac = autocorrelation(q_binned[:, j], name=f"{ensemble}_t{j}", S=S)
            tau_int_ts[j] = ac.tau_int
            dtau_int_ts[j] = ac.dtau_int
        except Exception:
            tau_int_ts[j] = 0.5
            dtau_int_ts[j] = 0.0

    return {
        "t_phys":         t_phys,
        "t_centres_latt": t_centres,
        "means":          means,
        "errors":         errors,
        "chi_tilde":      chi_tilde,
        "chi_tilde_err":  chi_tilde_err,
        "tau_int_ts":     tau_int_ts,
        "dtau_int_ts":    dtau_int_ts,
        "n_configs":      q_matrix.shape[0],
        "n_bin":          n_bin,
        "a_fm":           lattice_spacing_fm,
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_timeslice_density(results: dict,
                           ax: Optional[plt.Axes] = None,
                           label: str = "",
                           open_bc: bool = False,
                           n_exclude: int = 2,
                           color: str = None) -> plt.Axes:
    """Plot q(t) with Gamma-method error bars.

    Parameters
    ----------
    results  : output of analyse_timeslices()
    ax       : matplotlib Axes (created if None)
    label    : legend label
    open_bc  : if True, shade excluded boundary regions
    n_exclude: slices excluded per boundary
    color    : line/marker color
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4))

    t    = results["t_phys"]
    mean = results["means"]
    err  = results["errors"]
    a    = results["a_fm"]

    ax.errorbar(t, mean, yerr=err, fmt='o', capsize=3,
                label=label or f"a={a:.3f} fm, n_bin={results['n_bin']}",
                color=color)
    ax.axhline(0, color='grey', lw=0.8, ls='--')

    if open_bc:
        t_excl = n_exclude * a
        t_max  = t.max()
        ax.axvspan(0, t_excl, alpha=0.10, color='grey', label='boundary excluded')
        ax.axvspan(t_max, t_max + t_excl, alpha=0.10, color='grey')

    ax.set_xlabel("t [fm]")
    ax.set_ylabel(r"$q(t)$")
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0,
              frameon=False, fontsize=8)
    _apply_axis_style(ax)
    return ax


def _smooth(y: np.ndarray, window: int = None) -> np.ndarray:
    """Rolling mean for a guiding line. Window defaults to ~10% of length, min 3."""
    w = max(3, len(y) // 10) if window is None else window
    kernel = np.ones(w) / w
    return np.convolve(y, kernel, mode='same')


def plot_timeslice_mctime(filepath: str, n_bin: int = 1, smear: int = None) -> plt.Figure:
    """Plot q(t) vs MC time for all (binned) time slices on one figure.

    Each binned slice is a separate coloured line. A low-alpha smoothed line
    is drawn as a guide to the eye below the scatter points.

    Parameters
    ----------
    filepath : path to topcharge_timeslice.dat
    n_bin    : number of adjacent time slices to sum (must match timeslice_density)
    smear    : smearing level to select (None = max)
    """
    t_values, q_matrix = load_timeslice_data(filepath, smear=smear)
    t_centres, q_binned = bin_timeslices(t_values, q_matrix, n_bin=n_bin)
    n_configs, n_bins = q_binned.shape
    mc_time = np.arange(n_configs)

    cmap = plt.cm.tab20 if n_bins <= 20 else plt.cm.turbo
    colors = [cmap(i / max(1, n_bins - 1)) for i in range(n_bins)]

    fig, ax = plt.subplots(figsize=(10, 4))

    for j, tc in enumerate(t_centres):
        col = colors[j]
        label = f"$t={int(tc)}$" if n_bin == 1 else f"$t\\in[{t_values[j*n_bin]},{t_values[min(j*n_bin+n_bin-1, len(t_values)-1)]}]$"
        y = q_binned[:, j]
        ax.scatter(mc_time, y, s=4, color=col, alpha=0.8, zorder=3)
        ax.plot(mc_time, y, lw=0.8, color=col, alpha=0.3, zorder=2, label=label)

    ax.axhline(0, color='grey', lw=0.6, ls='--')
    ax.set_xlabel("MC time")
    ax.set_ylabel(r"$q(t)$")

    legend = ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0,
                       fontsize=8, frameon=False,
                       ncol=max(1, n_bins // 8))
    # Make legend lines fully opaque
    for handle in legend.legend_handles:
        handle.set_alpha(1.0)
        handle.set_linewidth(2.0)

    _apply_axis_style(ax, x_integer=True)

    fig.tight_layout()
    fig.subplots_adjust(right=0.82)
    return fig


# ---------------------------------------------------------------------------
# Grid plots (one panel per run, matching histogram / Q_vs_mctime style)
# ---------------------------------------------------------------------------

def plot_timeslice_density_grid(ts_results: list, runs: list, output_dir: str,
                                 gauge_group: str, boundary: str):
    """Grid of timeslice density plots, one panel per run.
    Single-column layout. Shared y-scale (largest measurement). x-axis in lattice units (Xa)."""
    from calculations import lattice_spacing_su2, lattice_spacing_su3
    import os

    pairs = [(res, run) for res, run in zip(ts_results, runs) if res is not None]
    n = len(pairs)
    if n == 0:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    # Global y-range: largest |mean ± err| across all panels
    y_max = max(np.max(np.abs(res["means"]) + res["errors"]) for res, _ in pairs)
    y_lim = (-y_max * 1.15, y_max * 1.15)

    fig, axes = plt.subplots(n, 1, figsize=(7, 3 * n))
    axes = np.array(axes).flatten()

    for idx, (res, run_data) in enumerate(pairs):
        ax = axes[idx]
        open_bc = run_data.boundary == "open"
        a = lattice_spacing(run_data.beta)

        t_latt = res["t_centres_latt"]
        means  = res["means"]
        errors = res["errors"]

        x_plot = t_latt

        h_data, = ax.plot([], [], 'o', color='steelblue', label=r'$q(t)$ mean $\pm\sigma$')
        ax.errorbar(x_plot, means, yerr=errors, fmt='o', capsize=3,
                    color='steelblue', zorder=3)
        ax.axhline(0, color='grey', lw=0.8, ls='--')

        h_excl = None
        if open_bc:
            ax.axvspan(-0.5, x_plot[0], alpha=0.10, color='grey')
            ax.axvspan(x_plot[-1], run_data.T - 0.5, alpha=0.10, color='grey')
            if idx == 0:
                import matplotlib.patches as mpatches
                h_excl = mpatches.Patch(color='grey', alpha=0.3, label='boundary excluded')
        ax.set_xlim(-0.5, run_data.T - 0.5)
        _ensure_endpoint_ticks(ax)

        ax.set_ylim(y_lim)
        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel('$t/a$', fontsize=9)
        ax.set_ylabel(r'$q(t)$', fontsize=9)
        _apply_axis_style(ax, x_integer=True)
        ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                fontsize=11, va='top', ha='left')

        if idx == 0:
            legend_handles = [h_data] + ([h_excl] if h_excl is not None else [])

    # Shared legend outside top-right of first panel
    axes[0].legend(legend_handles, [h.get_label() for h in legend_handles],
                   loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0,
                   fontsize=7, frameon=False, handlelength=1.4, labelspacing=0.3)

    plt.tight_layout()
    out = os.path.join(output_dir, f"timeslice_density_{gauge_group}_{boundary}.png")
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"Saved: {os.path.basename(out)}")


def plot_timeslice_density_comparison(
        ts_results_periodic: list, runs_periodic: list,
        ts_results_open: list, runs_open: list,
        output_dir: str, gauge_group: str):
    """Overlay open (square) and periodic (triangle) boundary density in one panel per beta.

    One panel per unique beta value. Open BC uses squares, periodic uses triangles.
    Points are offset slightly in x so error bars don't overlap.
    Grey excluded regions shown for open BC only. Shared y-scale across panels.
    """
    from calculations import lattice_spacing_su2, lattice_spacing_su3
    import matplotlib.patches as mpatches
    import os

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    # Index by beta for each BC
    periodic_by_beta = {r.beta: (res, r) for res, r in zip(ts_results_periodic, runs_periodic)
                        if res is not None}
    open_by_beta     = {r.beta: (res, r) for res, r in zip(ts_results_open, runs_open)
                        if res is not None}

    all_betas = sorted(set(periodic_by_beta) | set(open_by_beta))
    if not all_betas:
        return

    # Global y-range across all datasets
    all_vals = []
    for d in (periodic_by_beta, open_by_beta):
        for res, _ in d.values():
            all_vals.append(np.max(np.abs(res["means"]) + res["errors"]))
    y_max = max(all_vals)
    y_lim = (-y_max * 1.15, y_max * 1.15)

    n = len(all_betas)
    fig, axes = plt.subplots(n, 1, figsize=(7, 3 * n))
    axes = np.array(axes).flatten()

    x_offset = 0.25  # lattice units offset between BC types

    for idx, beta in enumerate(all_betas):
        ax = axes[idx]
        a = lattice_spacing(beta)
        ax.axhline(0, color='grey', lw=0.8, ls='--')

        legend_handles = []

        _, ref_run = (open_by_beta if beta in open_by_beta else periodic_by_beta)[beta]

        # Open BC — squares, offset right
        if beta in open_by_beta:
            res, _ = open_by_beta[beta]
            t_latt = res["t_centres_latt"]
            x_plot = t_latt + x_offset
            h, = ax.plot([], [], 's', color='steelblue', label='open BC')
            ax.errorbar(x_plot, res["means"], yerr=res["errors"],
                        fmt='s', capsize=3, color='steelblue', zorder=3)
            ax.axvspan(-0.5, t_latt[0], alpha=0.10, color='grey')
            ax.axvspan(t_latt[-1], ref_run.T - 0.5, alpha=0.10, color='grey')
            legend_handles.append(h)

        # Periodic BC — triangles, offset left
        if beta in periodic_by_beta:
            res, _ = periodic_by_beta[beta]
            t_latt = res["t_centres_latt"]
            x_plot = t_latt - x_offset
            h, = ax.plot([], [], '^', color='tomato', label='periodic BC')
            ax.errorbar(x_plot, res["means"], yerr=res["errors"],
                        fmt='^', capsize=3, color='tomato', zorder=3)
            legend_handles.append(h)

        if beta in open_by_beta:
            h_excl = mpatches.Patch(color='grey', alpha=0.3, label='boundary excluded')
            legend_handles.append(h_excl)

        ax.set_xlim(-0.5, ref_run.T - 0.5)
        _ensure_endpoint_ticks(ax)
        ax.set_ylim(y_lim)
        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel('$t/a$', fontsize=9)
        ax.set_ylabel(r'$q(t)$', fontsize=9)
        ax.tick_params(labelsize=8)
        ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                fontsize=11, va='top', ha='left')

        if idx == 0:
            first_handles = legend_handles

    axes[0].legend(first_handles, [h.get_label() for h in first_handles],
                   loc='upper left', bbox_to_anchor=(1.02, 1.0),
                   fontsize=7, frameon=False)

    plt.tight_layout()
    out = os.path.join(output_dir, f"timeslice_density_comparison_{gauge_group}.png")
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.basename(out)}")


def plot_timeslice_mctime_grid(ts_files: list, runs: list, n_bin: int,
                                output_dir: str, gauge_group: str, boundary: str):
    """Grid of q(t) vs MC time, one panel per run. A,B,C labels + a title.
    Each binned time slice is a separate coloured line; one shared legend."""
    from calculations import lattice_spacing_su2, lattice_spacing_su3
    import os

    valid = [(f, r) for f, r in zip(ts_files, runs) if f is not None]
    n = len(valid)
    if n == 0:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    # Pre-compute consistent colormap from first file
    t_values_ref, _ = load_timeslice_data(valid[0][0])
    t_centres_ref, _ = bin_timeslices(t_values_ref, np.zeros((1, len(t_values_ref))), n_bin)
    n_bins = len(t_centres_ref)
    cmap = plt.cm.tab20 if n_bins <= 20 else plt.cm.turbo
    colors = [cmap(i / max(1, n_bins - 1)) for i in range(n_bins)]

    def bin_label(j):
        lo = t_values_ref[j * n_bin]
        hi = t_values_ref[min(j * n_bin + n_bin - 1, len(t_values_ref) - 1)]
        return f"t={int(lo)}" if n_bin == 1 else f"t={int(lo)}-{int(hi)}"

    # Pre-load all data to compute global y-range
    all_q_binned = []
    for ts_file, _ in valid:
        t_v, q_m = load_timeslice_data(ts_file)
        _, q_b = bin_timeslices(t_v, q_m, n_bin)
        all_q_binned.append(q_b)
    global_abs = max(np.abs(q_b).max() for q_b in all_q_binned)
    y_lim = (-global_abs * 1.15, global_abs * 1.15)

    fig, axes = plt.subplots(n, 1, figsize=(7, 3 * n))
    axes = np.array(axes).flatten()

    legend_handles = []
    max_legend_entries = 8
    selected = np.linspace(0, n_bins - 1, min(n_bins, max_legend_entries), dtype=int)
    selected_idx = set(int(x) for x in selected)
    for idx, ((ts_file, run_data), q_binned) in enumerate(zip(valid, all_q_binned)):
        ax = axes[idx]
        mc_time = np.arange(q_binned.shape[0])

        for j in range(n_bins):
            y = q_binned[:, j]
            col_color = colors[j]
            h, = ax.plot(mc_time, y, lw=0.8, color=col_color, alpha=0.3, zorder=2,
                         label=bin_label(j))
            ax.scatter(mc_time, y, s=3, color=col_color, alpha=0.8, zorder=3)
            if idx == 0 and j in selected_idx:
                legend_handles.append(h)

        ax.axhline(0, color='grey', lw=0.6, ls='--')
        ax.set_ylim(y_lim)
        a = lattice_spacing(run_data.beta)
        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel('MC time', fontsize=9)
        ax.set_ylabel(r'$q(t)$', fontsize=9)
        _apply_axis_style(ax, x_integer=True)
        ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                fontsize=11, va='top', ha='left')

    # Shared legend outside the top axes panel, no gap
    leg = axes[0].legend(legend_handles, [h.get_label() for h in legend_handles],
                         loc='upper left', bbox_to_anchor=(1.02, 1.0),
                         borderaxespad=0, fontsize=7, frameon=False,
                         handlelength=1.4, labelspacing=0.3)
    for h in leg.legend_handles:
        h.set_alpha(1.0)
        h.set_linewidth(2.0)

    plt.tight_layout()
    out = os.path.join(output_dir, f"timeslice_mctime_{gauge_group}_{boundary}.png")
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.basename(out)}")


def plot_timeslice_susceptibility_grid(ts_results: list, runs: list, output_dir: str,
                                       gauge_group: str, boundary: str):
    """Grid of per-timeslice susceptibility chi_tilde(t) = <q(t)^2> / V_tilde."""
    from calculations import lattice_spacing_su2, lattice_spacing_su3, HBAR_C
    import os

    pairs = [(res, run) for res, run in zip(ts_results, runs) if res is not None]
    n = len(pairs)
    if n == 0:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    fig, axes = plt.subplots(n, 1, figsize=(7, 3 * n))
    axes = np.array(axes).flatten()

    for idx, (res, run_data) in enumerate(pairs):
        ax = axes[idx]
        a = lattice_spacing(run_data.beta)
        V_tilde = run_data.L ** 3 * res["n_bin"]

        t_latt = res["t_centres_latt"]
        chi_latt = res["chi_tilde"] / V_tilde
        chi_latt_err = res["chi_tilde_err"] / V_tilde

        # Convert to MeV via fourth root: chi^(1/4) [MeV] = (chi_latt / a^4)^(1/4) * hbar_c
        chi_phys = chi_latt / a**4
        chi_phys_err = chi_latt_err / a**4
        chi_fourth = np.sign(chi_phys) * np.abs(chi_phys)**0.25 * HBAR_C
        chi_fourth_err = 0.25 * np.abs(chi_phys)**(-0.75) * chi_phys_err * HBAR_C
        # Avoid NaN for zero values
        chi_fourth_err = np.nan_to_num(chi_fourth_err, nan=0.0)

        x_plot = t_latt

        ref_val = 200.0 if gauge_group == "su2" else 191.0
        ax.axhline(ref_val, color='grey', lw=1.0, ls='--', label=f'ref {ref_val} MeV')

        ax.errorbar(x_plot, chi_fourth, yerr=chi_fourth_err, fmt='o', capsize=3,
                    color='steelblue', zorder=3, label=r'$\tilde{\chi}(t)^{1/4}$')

        open_bc = run_data.boundary == "open"
        if open_bc:
            ax.axvspan(-0.5, x_plot[0], alpha=0.10, color='grey')
            ax.axvspan(x_plot[-1], run_data.T - 0.5, alpha=0.10, color='grey')
        ax.set_xlim(-0.5, run_data.T - 0.5)
        _ensure_endpoint_ticks(ax)

        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel('$t/a$', fontsize=9)
        ax.set_ylabel(r'$\tilde{\chi}(t)^{1/4}$ [MeV]', fontsize=9)
        _apply_axis_style(ax, x_integer=True)
        ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                fontsize=11, va='top', ha='left')
        if idx == 0:
            ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0,
                  fontsize=7, frameon=False, handlelength=1.4, labelspacing=0.3)

    plt.tight_layout()
    out = os.path.join(output_dir, f"timeslice_susceptibility_{gauge_group}_{boundary}.png")
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.basename(out)}")


def plot_timeslice_tauint_grid(ts_results: list, runs: list, output_dir: str,
                                gauge_group: str, boundary: str):
    """Grid of per-timeslice integrated autocorrelation time tau_int(t)."""
    from calculations import lattice_spacing_su2, lattice_spacing_su3
    import os

    pairs = [(res, run) for res, run in zip(ts_results, runs) if res is not None]
    n = len(pairs)
    if n == 0:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    fig, axes = plt.subplots(n, 1, figsize=(7, 3 * n))
    axes = np.array(axes).flatten()

    for idx, (res, run_data) in enumerate(pairs):
        ax = axes[idx]
        a = lattice_spacing(run_data.beta)

        t_latt = res["t_centres_latt"]
        tau = res["tau_int_ts"]
        dtau = res["dtau_int_ts"]

        x_plot = t_latt

        ax.errorbar(x_plot, tau, yerr=dtau, fmt='s', capsize=3,
                    color='tomato', zorder=3, label=r'$\tau_{\mathrm{int}}$')

        open_bc = run_data.boundary == "open"
        if open_bc:
            ax.axvspan(-0.5, x_plot[0], alpha=0.10, color='grey')
            ax.axvspan(x_plot[-1], run_data.T - 0.5, alpha=0.10, color='grey')
        ax.set_xlim(-0.5, run_data.T - 0.5)
        _ensure_endpoint_ticks(ax)

        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel('$t/a$', fontsize=9)
        ax.set_ylabel(r'$\tau_{\mathrm{int}}(t)$', fontsize=9)
        _apply_axis_style(ax, x_integer=True)
        ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                fontsize=11, va='top', ha='left')
        if idx == 0:
            ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0,
                  fontsize=7, frameon=False, handlelength=1.4, labelspacing=0.3)

    plt.tight_layout()
    out = os.path.join(output_dir, f"timeslice_tauint_{gauge_group}_{boundary}.png")
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.basename(out)}")


def plot_timeslice_susceptibility_cropped(ts_results: list, runs: list, output_dir: str,
                                          gauge_group: str, boundary: str):
    """Cropped susceptibility plot: only evaluated timeslices, re-indexed from 0."""
    from calculations import lattice_spacing_su2, lattice_spacing_su3, HBAR_C
    import os

    pairs = [(res, run) for res, run in zip(ts_results, runs) if res is not None]
    n = len(pairs)
    if n == 0:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    fig, axes = plt.subplots(n, 1, figsize=(7, 3 * n))
    axes = np.array(axes).flatten()

    for idx, (res, run_data) in enumerate(pairs):
        ax = axes[idx]
        a = lattice_spacing(run_data.beta)
        V_tilde = run_data.L ** 3 * res["n_bin"]

        chi_latt = res["chi_tilde"] / V_tilde
        chi_latt_err = res["chi_tilde_err"] / V_tilde

        chi_phys = chi_latt / a**4
        chi_phys_err = chi_latt_err / a**4
        chi_fourth = np.sign(chi_phys) * np.abs(chi_phys)**0.25 * HBAR_C
        chi_fourth_err = 0.25 * np.abs(chi_phys)**(-0.75) * chi_phys_err * HBAR_C
        chi_fourth_err = np.nan_to_num(chi_fourth_err, nan=0.0)

        x_plot = np.arange(len(chi_fourth))

        ref_val = 200.0 if gauge_group == "su2" else 191.0
        ax.axhline(ref_val, color='grey', lw=1.0, ls='--', label=f'ref {ref_val} MeV')

        ax.errorbar(x_plot, chi_fourth, yerr=chi_fourth_err, fmt='o', capsize=3,
                    color='steelblue', zorder=3, label=r'$\tilde{\chi}(t)^{1/4}$')

        ax.set_xlim(-0.5, len(chi_fourth) - 0.5)
        _ensure_endpoint_ticks(ax)
        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel('$t/a$ (evaluated only)', fontsize=9)
        ax.set_ylabel(r'$\tilde{\chi}(t)^{1/4}$ [MeV]', fontsize=9)
        _apply_axis_style(ax, x_integer=True)
        ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                fontsize=11, va='top', ha='left')
        if idx == 0:
            ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0,
                  fontsize=7, frameon=False, handlelength=1.4, labelspacing=0.3)

    plt.tight_layout()
    out = os.path.join(output_dir, f"timeslice_susceptibility_cropped_{gauge_group}_{boundary}.png")
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.basename(out)}")


def plot_timeslice_tauint_cropped(ts_results: list, runs: list, output_dir: str,
                                   gauge_group: str, boundary: str):
    """Cropped tauint plot: only evaluated timeslices, re-indexed from 0."""
    from calculations import lattice_spacing_su2, lattice_spacing_su3
    import os

    pairs = [(res, run) for res, run in zip(ts_results, runs) if res is not None]
    n = len(pairs)
    if n == 0:
        return

    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    fig, axes = plt.subplots(n, 1, figsize=(7, 3 * n))
    axes = np.array(axes).flatten()

    for idx, (res, run_data) in enumerate(pairs):
        ax = axes[idx]
        a = lattice_spacing(run_data.beta)

        tau = res["tau_int_ts"]
        dtau = res["dtau_int_ts"]

        x_plot = np.arange(len(tau))

        ax.errorbar(x_plot, tau, yerr=dtau, fmt='s', capsize=3,
                    color='tomato', zorder=3, label=r'$\tau_{\mathrm{int}}$')

        ax.set_xlim(-0.5, len(tau) - 0.5)
        _ensure_endpoint_ticks(ax)
        ax.set_title(f'$a = {a:.4f}$ fm', fontsize=10)
        ax.set_xlabel('$t/a$ (evaluated only)', fontsize=9)
        ax.set_ylabel(r'$\tau_{\mathrm{int}}(t)$', fontsize=9)
        _apply_axis_style(ax, x_integer=True)
        ax.text(0.02, 0.98, chr(ord('A') + idx), transform=ax.transAxes,
                fontsize=11, va='top', ha='left')
        if idx == 0:
            ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), borderaxespad=0,
                  fontsize=7, frameon=False, handlelength=1.4, labelspacing=0.3)

    plt.tight_layout()
    out = os.path.join(output_dir, f"timeslice_tauint_cropped_{gauge_group}_{boundary}.png")
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {os.path.basename(out)}")


# ---------------------------------------------------------------------------
# CLI convenience
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Timeslice topological charge analysis")
    parser.add_argument("filepath", help="Path to topcharge_timeslice.dat")
    parser.add_argument("--a", type=float, required=True, help="Lattice spacing in fm")
    parser.add_argument("--n-bin", type=int, default=2, help="Time slices per bin (default 2)")
    parser.add_argument("--smear", type=int, default=None, help="Smearing steps to select")
    parser.add_argument("--open-bc", action="store_true")
    parser.add_argument("--n-exclude", type=int, default=2)
    parser.add_argument("--S", type=float, default=1.5, help="Gamma method S parameter")
    parser.add_argument("--output", default=None, help="Save figure to this path")
    args = parser.parse_args()

    res = analyse_timeslices(
        args.filepath,
        lattice_spacing_fm=args.a,
        n_bin=args.n_bin,
        smear=args.smear,
        open_bc=args.open_bc,
        n_exclude=args.n_exclude,
        S=args.S,
    )

    fig, ax = plt.subplots(figsize=(8, 4))
    plot_timeslice_density(res, ax=ax, open_bc=args.open_bc, n_exclude=args.n_exclude)
    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=150)
        print(f"Saved to {args.output}")
    else:
        plt.show()
