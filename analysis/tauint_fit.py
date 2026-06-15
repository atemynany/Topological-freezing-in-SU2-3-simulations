#!/usr/bin/env python3
"""Power-law fit  tau_int(a) = c * a^(-z)  for hand-picked ensembles.

Workflow:
    1. Run the regular analysis once to get the AnalysisResult lists.
    2. Look at output/figures_analysis/.../tau_int_su2.png and pick the
       point numbers you want included (the plot's annotation labels).
    3. Call fit_selected(results_periodic, results_open, indices=[...]).

Indexing matches plot_tau_int_combined exactly: periodic block first
(sorted by a descending), then open block (sorted by a descending).

The fit minimizes weighted chi^2 with weights 1/sigma_i^2.  Parameter
errors come from a parametric bootstrap (resample tau_i from
N(tau_i, sigma_i), refit n_boot times) so the reported sigma_z reflects
the empirical, possibly non-Gaussian distribution of z_hat.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, List, Tuple

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2_dist


def power_law(a: np.ndarray, c: float, z: float) -> np.ndarray:
    return c * np.power(a, -z)


@dataclass
class FitResult:
    c: float
    z: float
    sigma_c: float
    sigma_z: float
    chi2: float
    dof: int
    chi2_red: float
    p_value: float
    a: np.ndarray
    tau: np.ndarray
    sigma: np.ndarray
    bootstrap_c: np.ndarray
    bootstrap_z: np.ndarray

    def summary(self) -> str:
        lines = [
            f"  N points       : {len(self.a)}",
            f"  c              : {self.c:.4g} +/- {self.sigma_c:.2g}",
            f"  z              : {self.z:.4g} +/- {self.sigma_z:.2g}",
            f"  chi^2 / dof    : {self.chi2:.3f} / {self.dof}  =  {self.chi2_red:.3f}",
            f"  p-value        : {self.p_value:.3g}",
        ]
        return "\n".join(lines)


def _initial_guess(a: np.ndarray, tau: np.ndarray) -> Tuple[float, float]:
    """Log-log linear regression to seed the nonlinear fit."""
    mask = (a > 0) & (tau > 0)
    if mask.sum() < 2:
        return 1.0, 2.0
    slope, intercept = np.polyfit(np.log(a[mask]), np.log(tau[mask]), 1)
    return float(np.exp(intercept)), float(-slope)


def _single_fit(a: np.ndarray, tau: np.ndarray, sigma: np.ndarray
                ) -> Tuple[float, float]:
    p0 = _initial_guess(a, tau)
    popt, _ = curve_fit(
        power_law, a, tau,
        sigma=sigma, absolute_sigma=True,
        p0=p0, maxfev=20000,
    )
    return float(popt[0]), float(popt[1])


def fit_power_law(
    a: Sequence[float],
    tau: Sequence[float],
    sigma: Sequence[float],
    n_boot: int = 4000,
    seed: int = 42,
    sigma_floor_frac: float = 0.05,
) -> FitResult:
    """Weighted chi^2 fit with parametric bootstrap.

    sigma_floor_frac: if a point reports sigma_i = 0 (e.g. frozen-topology
    ensembles where pyerrors returns dtau=0), replace it with
    max(sigma_floor_frac * tau_i, 1e-6) so that point does not get
    infinite weight.  Set to 0 to disable.
    """
    a = np.asarray(a, dtype=float)
    tau = np.asarray(tau, dtype=float)
    sigma = np.asarray(sigma, dtype=float)

    if a.size < 3:
        raise ValueError(f"need at least 3 points to fit (got {a.size})")

    sigma_eff = sigma.copy()
    bad = sigma_eff <= 0
    if bad.any():
        sigma_eff[bad] = np.maximum(sigma_floor_frac * tau[bad], 1e-6)

    c_hat, z_hat = _single_fit(a, tau, sigma_eff)

    resid = (tau - power_law(a, c_hat, z_hat)) / sigma_eff
    chi2_val = float(np.sum(resid ** 2))
    dof = len(a) - 2
    chi2_red = chi2_val / dof if dof > 0 else float("nan")
    p_val = float(1.0 - chi2_dist.cdf(chi2_val, dof)) if dof > 0 else float("nan")

    rng = np.random.default_rng(seed)
    boot_c = np.empty(n_boot)
    boot_z = np.empty(n_boot)
    fail = 0
    for k in range(n_boot):
        tau_b = rng.normal(tau, sigma_eff)
        if np.any(tau_b <= 0):
            tau_b = np.maximum(tau_b, 1e-6)
        try:
            c_b, z_b = _single_fit(a, tau_b, sigma_eff)
            boot_c[k] = c_b
            boot_z[k] = z_b
        except Exception:
            boot_c[k] = np.nan
            boot_z[k] = np.nan
            fail += 1

    if fail:
        print(f"  [bootstrap] {fail}/{n_boot} resamples failed to converge")

    boot_c = boot_c[~np.isnan(boot_c)]
    boot_z = boot_z[~np.isnan(boot_z)]
    sigma_c = float(np.std(boot_c, ddof=1)) if boot_c.size > 1 else float("nan")
    sigma_z = float(np.std(boot_z, ddof=1)) if boot_z.size > 1 else float("nan")

    return FitResult(
        c=c_hat, z=z_hat,
        sigma_c=sigma_c, sigma_z=sigma_z,
        chi2=chi2_val, dof=dof,
        chi2_red=chi2_red, p_value=p_val,
        a=a, tau=tau, sigma=sigma_eff,
        bootstrap_c=boot_c, bootstrap_z=boot_z,
    )


def _ordered_results(results, lattice_spacing):
    """Mirror plot_tau_int_combined: sort by a descending."""
    results = list(results)
    if not results:
        return results, np.array([])
    a_vals = np.array([lattice_spacing(r.beta) for r in results])
    order = np.argsort(a_vals)[::-1]
    return [results[i] for i in order], a_vals[order]


def select_indexed_points(
    results_periodic,
    results_open,
    indices: Sequence[int],
    lattice_spacing,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List]:
    """Pick AnalysisResults by their plot label number (1-indexed).

    Returns (a, tau_int, dtau_int, [AnalysisResult, ...]) for the selection.
    """
    per_sorted, _ = _ordered_results(results_periodic, lattice_spacing)
    open_sorted, _ = _ordered_results(results_open, lattice_spacing)
    flat = per_sorted + open_sorted

    chosen = []
    for idx in indices:
        if idx < 1 or idx > len(flat):
            raise IndexError(f"point index {idx} out of range 1..{len(flat)}")
        chosen.append(flat[idx - 1])

    a = np.array([lattice_spacing(r.beta) for r in chosen])
    tau = np.array([r.tau_int for r in chosen])
    dtau = np.array([r.dtau_int for r in chosen])
    return a, tau, dtau, chosen


def fit_selected(
    results_periodic,
    results_open,
    indices: Sequence[int],
    gauge_group: str = "su2",
    n_boot: int = 4000,
    seed: int = 42,
) -> FitResult:
    """End-to-end: pick points by plot label, return FitResult."""
    from calculations import lattice_spacing_su2, lattice_spacing_su3
    spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3

    a, tau, dtau, chosen = select_indexed_points(
        results_periodic, results_open, indices, spacing,
    )

    print(f"Selected {len(chosen)} points for {gauge_group}:")
    print(f"  {'idx':>4} {'beta':>6} {'a(fm)':>9} {'BC':>9}"
          f" {'tau_int':>10} {'dtau':>10}")
    for i, (j, r) in enumerate(zip(indices, chosen)):
        print(f"  {j:>4} {r.beta:>6.3f} {spacing(r.beta):>9.4f} {r.boundary:>9}"
              f" {r.tau_int:>10.3f} {r.dtau_int:>10.3f}")

    fit = fit_power_law(a, tau, dtau, n_boot=n_boot, seed=seed)
    print("\nFit  tau = c * a^(-z):")
    print(fit.summary())
    return fit


def plot_fit(fit: FitResult, output_path: str, title: str = "") -> None:
    """Quick diagnostic plot: data + fitted curve + bootstrap band."""
    import matplotlib.pyplot as plt

    a_grid = np.geomspace(fit.a.min() * 0.9, fit.a.max() * 1.1, 200)
    tau_grid = power_law(a_grid, fit.c, fit.z)

    boot_curves = np.array([
        power_law(a_grid, c_b, z_b)
        for c_b, z_b in zip(fit.bootstrap_c, fit.bootstrap_z)
    ])
    lo = np.percentile(boot_curves, 16, axis=0)
    hi = np.percentile(boot_curves, 84, axis=0)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.fill_between(a_grid, lo, hi, alpha=0.25, color="tab:blue",
                    label="bootstrap 68%")
    ax.plot(a_grid, tau_grid, "-", color="tab:blue",
            label=fr"$c\,a^{{-z}}$, $z={fit.z:.2f}\pm{fit.sigma_z:.2f}$")
    ax.errorbar(fit.a, fit.tau, yerr=fit.sigma, fmt="o", color="black",
                capsize=3, label="data")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$a$ (fm)")
    ax.set_ylabel(r"$\tau_{\rm int}(Q^2)$")
    if title:
        ax.set_title(title)
    ax.legend(frameon=False)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {output_path}")


def fit_and_plot_boundaries(
    results_periodic,
    results_open,
    output_dir: str,
    gauge_group: str = "su2",
    min_a: float = 0.012,
    exclude_largest_a: int = 0,
    n_boot: int = 4000,
    seed: int = 42,
) -> dict[str, FitResult]:
    """Fit periodic/open tau_int(a) separately and save diagnostic plots.

    min_a excludes the a ~= 0.01 fm points from the scaling fit.
    exclude_largest_a removes the coarsest distinct lattice spacings.
    Returns the fits that had enough points to run.
    """
    import os
    from calculations import lattice_spacing_su2, lattice_spacing_su3

    spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    fit_dir = os.path.join(output_dir, "fits")
    os.makedirs(fit_dir, exist_ok=True)
    suffix = "_fine_a" if exclude_largest_a > 0 else "_no_a001"

    fits = {}
    for name, results in [("periodic", results_periodic), ("open", results_open)]:
        rows = []
        excluded = []
        results_sorted = sorted(results, key=lambda r: spacing(r.beta), reverse=True)
        unique_a_desc = np.array(sorted({spacing(r.beta) for r in results_sorted}, reverse=True))
        excluded_large_a = unique_a_desc[:exclude_largest_a] if exclude_largest_a > 0 else np.array([])
        for result in results_sorted:
            a_value = spacing(result.beta)
            row = (a_value, result.tau_int, result.dtau_int)
            if a_value >= min_a and not np.isin(a_value, excluded_large_a):
                rows.append(row)
            else:
                excluded.append(row)

        if len(rows) < 3:
            print(
                f"  [warn] skipping {name} tau_int fit: "
                f"need at least 3 points after excluding a < {min_a:g} fm"
            )
            continue

        fit = fit_power_law(
            [row[0] for row in rows],
            [row[1] for row in rows],
            [row[2] for row in rows],
            n_boot=n_boot,
            seed=seed,
        )
        fits[name] = fit

        output_path = os.path.join(
            fit_dir, f"tau_int_fit_{name}{suffix}_{gauge_group}.png"
        )
        title_extra = ", exclude coarsest 2 a" if exclude_largest_a == 2 else ""
        plot_fit(
            fit,
            output_path,
            title=fr"{gauge_group.upper()} {name}: exclude $a \simeq 0.01$ fm{title_extra}",
        )
        print(
            f"  {name} fit: tau_int(a) = ({fit.c:.4g} +/- {fit.sigma_c:.2g}) "
            f"* a^(-{fit.z:.4g} +/- {fit.sigma_z:.2g}); "
            f"N={len(rows)}, excluded={len(excluded)}"
        )

    if fits:
        summary_path = os.path.join(fit_dir, f"tau_int_fit_summary{suffix}_{gauge_group}.txt")
        with open(summary_path, "w") as f:
            f.write("# Fit form: tau_int(a) = c * a^(-z)\n")
            f.write(f"# Excluded points with a < {min_a:g} fm\n")
            f.write(f"# Excluded coarsest distinct lattice spacings: {exclude_largest_a}\n")
            f.write("# boundary n_points c sigma_c z sigma_z chi2 dof chi2_red p_value\n")
            for name, fit in fits.items():
                f.write(
                    f"{name} {len(fit.a)} {fit.c:.12g} {fit.sigma_c:.12g} "
                    f"{fit.z:.12g} {fit.sigma_z:.12g} {fit.chi2:.12g} "
                    f"{fit.dof:d} {fit.chi2_red:.12g} {fit.p_value:.12g}\n"
                )
        print(f"Saved: {summary_path}")

    return fits


def plot_tau_int_with_fits(
    results_periodic,
    results_open,
    output_path: str,
    gauge_group: str = "su2",
    min_a: float = 0.012,
    exclude_largest_a: int = 0,
    n_boot: int = 4000,
    seed: int = 42,
) -> dict[str, FitResult]:
    """Save tau_int(a) with periodic/open data and fitted scaling curves."""
    import os
    import matplotlib.pyplot as plt
    from calculations import is_t160_l80, lattice_spacing_su2, lattice_spacing_su3

    spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    fig, ax = plt.subplots(figsize=(8, 5))

    fits = {}
    have_highlight_label = False
    all_a = np.array([
        spacing(r.beta)
        for results in (results_periodic, results_open)
        for r in results
    ])
    fit_styles = {
        "Periodic": {"color": "firebrick", "linestyle": "-", "linewidth": 1.1},
        "Open": {"color": "seagreen", "linestyle": "--", "linewidth": 1.1},
    }
    for results, label, marker, color, line_style in [
        (results_periodic, "Periodic", "o", "steelblue", "-"),
        (results_open, "Open", "s", "darkorange", "--"),
    ]:
        if not results:
            continue

        a_values = np.array([spacing(r.beta) for r in results])
        tau = np.array([r.tau_int for r in results])
        dtau = np.array([r.dtau_int for r in results])
        order = np.argsort(a_values)[::-1]
        a_values, tau, dtau = a_values[order], tau[order], dtau[order]
        results_sorted = [results[i] for i in order]

        highlight = np.array([is_t160_l80(r) for r in results_sorted])
        normal = ~highlight
        if np.any(normal):
            ax.errorbar(a_values[normal], tau[normal], yerr=dtau[normal],
                        fmt=marker, markersize=7, linestyle="none",
                        capsize=3, color=color, label=label)
        if np.any(highlight):
            highlight_label = None if have_highlight_label else "T160_L80"
            have_highlight_label = True
            ax.errorbar(a_values[highlight], tau[highlight], yerr=dtau[highlight],
                        fmt="^", markersize=8, linestyle="none",
                        capsize=3, color="crimson", label=highlight_label,
                        zorder=5)

        fit_mask = a_values >= min_a
        if exclude_largest_a > 0:
            unique_a_desc = np.array(sorted(set(a_values), reverse=True))
            excluded_large_a = unique_a_desc[:exclude_largest_a]
            fit_mask &= ~np.isin(a_values, excluded_large_a)
        excluded_mask = ~fit_mask
        if np.any(excluded_mask):
            normal_excluded = excluded_mask & normal
            highlight_excluded = excluded_mask & highlight
            if np.any(normal_excluded):
                ax.scatter(a_values[normal_excluded], tau[normal_excluded],
                       marker=marker, s=90, facecolors="none", edgecolors=color,
                       linewidths=1.4, alpha=0.8)
            if np.any(highlight_excluded):
                ax.scatter(a_values[highlight_excluded], tau[highlight_excluded],
                           marker="^", s=90, facecolors="none", edgecolors="crimson",
                           linewidths=1.4, alpha=0.9, zorder=5.5)

        if np.sum(fit_mask) >= 3:
            fit = fit_power_law(
                a_values[fit_mask], tau[fit_mask], dtau[fit_mask],
                n_boot=n_boot, seed=seed,
            )
            fits[label.lower()] = fit

            if all_a.size:
                a_grid = np.linspace(all_a.min(), all_a.max(), 400)
            else:
                a_grid = np.linspace(a_values[fit_mask].min(), a_values[fit_mask].max(), 400)
            style = fit_styles[label]
            ax.plot(a_grid, power_law(a_grid, fit.c, fit.z), style["linestyle"],
                    color=style["color"], linewidth=style["linewidth"], alpha=0.95,
                    zorder=2.5,
                    label=fr"{label} fit: $z={fit.z:.2f}\pm{fit.sigma_z:.2f}$")

    ax.axhline(0.5, linestyle="--", color="gray", linewidth=1, alpha=0.7,
               zorder=1, label=r"$\tau_{\rm int}=0.5$")
    ax.set_xlabel(r"$a$ (fm)")
    ax.set_ylabel(r"$\tau_{\rm int}(Q^2)$")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {output_path}")
    return fits


def _collect_results(gauge_group: str, results_dirs):
    """Mirror analysis.main()'s loader so this script runs standalone."""
    from calculations import analyze_run, load_run_data
    from result_paths import find_run_dirs

    run_dirs = find_run_dirs(results_dirs, gauge_group)

    results_periodic, results_open = [], []
    for run_dir in run_dirs:
        run_data = load_run_data(run_dir, gauge_group)
        if run_data is None:
            continue
        if gauge_group == "su2" and run_data.beta > 4:
            continue
        if gauge_group == "su3" and run_data.beta < 4:
            continue
        result = analyze_run(run_data)
        if run_data.boundary == "periodic":
            results_periodic.append(result)
        else:
            results_open.append(result)
    return results_periodic, results_open


if __name__ == "__main__":
    import argparse
    import os
    import sys

    script_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, script_dir)

    parser = argparse.ArgumentParser(
        description="Fit tau_int(a) = c * a^(-z) on hand-picked ensembles.")
    parser.add_argument("--su2", action="store_true")
    parser.add_argument("--su3", action="store_true")
    parser.add_argument(
        "--points", required=True,
        help="comma-separated list of point indices from tau_int_<gauge>.png "
             "(e.g. '7,8,9,10')")
    parser.add_argument("--results-dir", action="append", default=None)
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--plot", default=None,
                        help="optional path to save diagnostic plot")
    args = parser.parse_args()

    if args.su2 == args.su3:
        parser.error("choose exactly one of --su2 / --su3")

    gauge = "su2" if args.su2 else "su3"
    base_dir = os.path.dirname(script_dir)
    from result_paths import resolve_results_dirs
    results_dirs = resolve_results_dirs(base_dir, args.results_dir)

    results_periodic, results_open = _collect_results(gauge, results_dirs)
    indices = [int(x) for x in args.points.split(",") if x.strip()]

    fit = fit_selected(results_periodic, results_open, indices,
                       gauge_group=gauge, n_boot=args.n_boot, seed=args.seed)
    if args.plot:
        plot_fit(fit, args.plot,
                 title=fr"{gauge.upper()}: $\tau_{{\rm int}}(a) = c\,a^{{-z}}$")
