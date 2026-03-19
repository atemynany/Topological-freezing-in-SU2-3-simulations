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
import pyerrors as pe
from pathlib import Path
from typing import Optional


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
                       n_bin: int = 2,
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

    return {
        "t_phys":         t_phys,
        "t_centres_latt": t_centres,
        "means":          means,
        "errors":         errors,
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

    ax.errorbar(t, mean, yerr=err, fmt='o-', capsize=3,
                label=label or f"a={a:.3f} fm, n_bin={results['n_bin']}",
                color=color)
    ax.axhline(0, color='grey', lw=0.8, ls='--')

    if open_bc:
        t_excl = n_exclude * a
        t_max  = t.max()
        ax.axvspan(0, t_excl, alpha=0.10, color='grey', label='boundary excluded')
        ax.axvspan(t_max, t_max + t_excl, alpha=0.10, color='grey')

    ax.set_xlabel("t [fm]")
    ax.set_ylabel(r"$q(t) = Q(t)/4\pi^2$")
    ax.set_title("Topological charge density per time slice")
    ax.legend()
    return ax


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
