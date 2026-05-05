#!/usr/bin/env python3
"""
Per-smear-level analysis driver.

Reuses the existing analysis/plotting functions, but runs the pipeline
once per APE smearing level present in the data and writes each set of
figures to its own subdirectory:  output/figures_analysis/smear<N>/

Smear levels are auto-detected from topcharge.dat (so this still works
on legacy single-level data).
"""

import os
import sys
import glob
import argparse
from typing import List, Optional

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

import timeslice_analysis as tsa
from calculations import (
    RunData, AnalysisResult,
    load_run_data, analyze_run, load_topcharge,
    find_optimal_alpha, QEstimator,
    lattice_spacing_su2, lattice_spacing_su3,
)
from plotting_code import (
    plot_Q_vs_mctime_grid, plot_histograms_grid,
    plot_susceptibility_combined, plot_tau_int_combined,
    plot_thermalization_comparison,
)
from timeslice_analysis import (
    analyse_timeslices,
    plot_timeslice_density_grid, plot_timeslice_mctime_grid,
    plot_timeslice_susceptibility_grid, plot_timeslice_tauint_grid,
    plot_timeslice_susceptibility_cropped, plot_timeslice_tauint_cropped,
)


def detect_smear_levels(topcharge_files: List[str]) -> List[int]:
    """Sorted unique smear levels found across the given topcharge.dat files."""
    smears = set()
    for fp in topcharge_files:
        with open(fp) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                try:
                    if len(parts) == 3:
                        smears.add(int(parts[1]))      # SU(3): conf smear Q
                    elif len(parts) >= 4:
                        smears.add(int(parts[0]))      # SU(2): smear conf Q plaq
                except ValueError:
                    continue
    return sorted(smears)


def _topcharge_path(run_dir: str, gauge_group: str) -> Optional[str]:
    for cand in (
        os.path.join(run_dir, "output", "topcharge.dat"),
        os.path.join(run_dir, "output", f"topcharge_{gauge_group}.dat"),
    ):
        if os.path.exists(cand):
            return cand
    return None


def rebuild_for_smear(run_data: RunData, smear: int,
                      gauge_group: str) -> Optional[RunData]:
    """Return a copy of run_data with Q_raw/Q_rescaled/alpha at this smear level.

    Mirrors the per-boundary logic in load_run_data so periodic runs get
    integer-rounded Q via QEstimator, while open BC keeps continuous α·Q.
    """
    fp = _topcharge_path(run_data.run_dir, gauge_group)
    if fp is None:
        return None
    Q_raw = load_topcharge(fp, smear=smear)
    if len(Q_raw) == 0:
        return None
    if run_data.boundary == "open":
        alpha = find_optimal_alpha(Q_raw)
        Q_rescaled = alpha * Q_raw
    else:
        est = QEstimator(Q_raw)
        Q_rescaled = est.Q_rescaled
        alpha = est.alpha
    return RunData(
        beta=run_data.beta, seed=run_data.seed, run_dir=run_data.run_dir,
        T=run_data.T, L=run_data.L, boundary=run_data.boundary,
        Q_raw=Q_raw, Q_rescaled=Q_rescaled, alpha=alpha,
        plaq_data=run_data.plaq_data, mc_time=run_data.mc_time,
        start_conf=run_data.start_conf,
    )


def print_table(results: List[AnalysisResult], gauge_group: str, smear: int):
    spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    print(f"\n{'='*95}")
    print(f"  {gauge_group.upper()} Analysis Summary  |  smear = {smear}")
    print(f"{'='*95}")
    print(f"{'β':>6} {'a(fm)':>8} {'BC':>10} {'<Q²>':>8} {'χ_t^1/4':>10} "
          f"{'τ_int':>8} {'N_conf':>8} {'α':>6}")
    print(f"{'':>6} {'':>8} {'':>10} {'':>8} {'(MeV)':>10} {'':>8} {'':>8} {'':>6}")
    print(f"{'-'*95}")
    for r in sorted(results, key=lambda x: (x.boundary, x.beta)):
        a = spacing(r.beta)
        chi_MeV = r.chi_t_fourth_root_a / a
        print(f"{r.beta:>6.2f} {a:>8.4f} {r.boundary:>10} {r.Q2_mean:>8.1f} "
              f"{chi_MeV:>10.1f} {r.tau_int:>8.1f} {r.n_conf:>8} {r.alpha:>6.3f}")
    print(f"{'='*95}\n")


def run_for_smear(runs_master: List[RunData], gauge_group: str,
                  output_root: str, smear: int, bin_ref: int = 4):
    out_dir = os.path.join(output_root, f"smear{smear}")
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n--- smear = {smear}  ->  {out_dir}")

    runs_periodic, runs_open = [], []
    res_periodic, res_open = [], []
    for rd in runs_master:
        rebuilt = rebuild_for_smear(rd, smear, gauge_group)
        if rebuilt is None:
            continue
        result = analyze_run(rebuilt)
        if rebuilt.boundary == "periodic":
            runs_periodic.append(rebuilt)
            res_periodic.append(result)
        else:
            runs_open.append(rebuilt)
            res_open.append(result)

    print(f"  {len(runs_periodic)} periodic, {len(runs_open)} open")
    print_table(res_periodic + res_open, gauge_group, smear)

    if runs_periodic:
        plot_Q_vs_mctime_grid(runs_periodic, out_dir, gauge_group, "periodic")
        plot_histograms_grid(runs_periodic, out_dir, gauge_group, "periodic")
    if runs_open:
        plot_Q_vs_mctime_grid(runs_open, out_dir, gauge_group, "open")
        plot_histograms_grid(runs_open, out_dir, gauge_group, "open")
    plot_susceptibility_combined(res_periodic, res_open, out_dir, gauge_group)
    plot_tau_int_combined(res_periodic, res_open, out_dir, gauge_group)

    all_runs = runs_periodic + runs_open
    if all_runs:
        plot_thermalization_comparison(all_runs, out_dir, gauge_group)

    # Constant-physical-subvolume binning. Two scales:
    #   analysis n_bin = round(a_ref / a)               → density / susc / tauint plots
    #   mctime  n_bin = round(bin_ref * a_ref / a)      → mctime line plot only
    # Reference is the smallest-β run (coarsest a).
    spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    all_runs_for_ref = runs_periodic + runs_open
    ref_beta = min(rd.beta for rd in all_runs_for_ref) if all_runs_for_ref else None
    a_ref = spacing(ref_beta) if ref_beta is not None else None
    if a_ref is not None:
        print(f"  subvolume reference: smallest β = {ref_beta}, a_ref = {a_ref:.4f} fm  "
              f"(mctime bin_ref = {bin_ref} → ≈ {bin_ref * a_ref:.4f} fm)")

    def _n_bin_for(beta: float) -> int:
        if a_ref is None:
            return 1
        return max(1, int(round(a_ref / spacing(beta))))

    def _n_bin_mctime_for(beta: float) -> int:
        if a_ref is None:
            return 1
        return max(1, int(round(bin_ref * a_ref / spacing(beta))))

    # plot_timeslice_mctime_grid loads timeslice data internally without a
    # smear arg, so pin the loader to this smear level for the duration.
    orig_loader = tsa.load_timeslice_data
    tsa.load_timeslice_data = lambda fp, smear=None, _orig=orig_loader, _s=smear: _orig(fp, smear=_s)
    try:
        for runs_grp, boundary in [(runs_periodic, "periodic"), (runs_open, "open")]:
            ts_files, ts_results, ts_runs, ts_nbins, ts_nbins_mctime = [], [], [], [], []
            for rd in sorted(runs_grp, key=lambda x: x.beta):
                ts_file = os.path.join(rd.run_dir, "output", "topcharge_timeslice.dat")
                if not os.path.isfile(ts_file):
                    continue
                a_fm = spacing(rd.beta)
                n_bin_run = _n_bin_for(rd.beta)
                n_bin_mctime = _n_bin_mctime_for(rd.beta)
                try:
                    res = analyse_timeslices(
                        ts_file,
                        lattice_spacing_fm=a_fm,
                        n_bin=n_bin_run,
                        smear=smear,
                        open_bc=(rd.boundary == "open"),
                        ensemble=f"b{rd.beta}_{rd.boundary}",
                    )
                except Exception as e:
                    print(f"  [warn] timeslice analysis failed for {rd.run_dir}: {e}")
                    continue
                ts_files.append(ts_file)
                ts_results.append(res)
                ts_runs.append(rd)
                ts_nbins.append(n_bin_run)
                ts_nbins_mctime.append(n_bin_mctime)

            if ts_runs:
                print(f"  [{boundary}] analysis n_bin: " +
                      ", ".join(f"β={r.beta}->{nb}" for r, nb in zip(ts_runs, ts_nbins)))
                print(f"  [{boundary}] mctime   n_bin: " +
                      ", ".join(f"β={r.beta}->{nb}" for r, nb in zip(ts_runs, ts_nbins_mctime)))
                plot_timeslice_density_grid(ts_results, ts_runs, out_dir, gauge_group, boundary)
                plot_timeslice_mctime_grid(ts_files, ts_runs, ts_nbins_mctime,
                                           out_dir, gauge_group, boundary)
                plot_timeslice_susceptibility_grid(ts_results, ts_runs, out_dir, gauge_group, boundary)
                plot_timeslice_tauint_grid(ts_results, ts_runs, out_dir, gauge_group, boundary)
                plot_timeslice_susceptibility_cropped(ts_results, ts_runs, out_dir, gauge_group, boundary)
                plot_timeslice_tauint_cropped(ts_results, ts_runs, out_dir, gauge_group, boundary)
    finally:
        tsa.load_timeslice_data = orig_loader


def main():
    parser = argparse.ArgumentParser(
        description="Run analysis once per APE smearing level present in the data.")
    parser.add_argument("--su2", action="store_true")
    parser.add_argument("--su3", action="store_true")
    parser.add_argument("--results-dir", type=str, default=None)
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Root dir; per-smear subdirs are created inside it.")
    parser.add_argument("--smear", type=int, action="append",
                        help="Restrict to a specific smear level (repeatable). "
                             "Default: every level detected in the data.")
    parser.add_argument("--bin-ref", type=int, default=4,
                        help="Bin width (in raw timeslices) for the smallest-β run. "
                             "Other runs are binned to match the same physical thickness. "
                             "Default: 4.")
    args = parser.parse_args()

    gauge_group = "su3" if args.su3 else "su2"
    base_dir = os.path.dirname(script_dir)
    results_dir = args.results_dir or os.path.join(base_dir, "data", "results")
    output_root = args.output_dir or os.path.join(base_dir, "output", "figures_analysis")
    os.makedirs(output_root, exist_ok=True)

    print(f"\n{gauge_group.upper()} per-smear analysis")
    print(f"Results: {results_dir}")
    print(f"Output root: {output_root}")

    run_dirs = sorted(glob.glob(
        os.path.join(results_dir, f"T*_L*_b*_*_seed*_{gauge_group}")))

    runs_master: List[RunData] = []
    topcharge_files: List[str] = []
    for run_dir in run_dirs:
        rd = load_run_data(run_dir, gauge_group)
        if rd is None:
            continue
        if gauge_group == "su2" and rd.beta > 4:
            continue
        if gauge_group == "su3" and rd.beta < 4:
            continue
        runs_master.append(rd)
        fp = _topcharge_path(run_dir, gauge_group)
        if fp:
            topcharge_files.append(fp)

    if not runs_master:
        print("No runs found.")
        return

    detected = detect_smear_levels(topcharge_files)
    smear_levels = args.smear or detected
    print(f"Detected smear levels: {detected}")
    print(f"Processing: {smear_levels}\n")

    for s in smear_levels:
        run_for_smear(runs_master, gauge_group, output_root, s, bin_ref=args.bin_ref)

    print(f"\nDone. Per-smear plots under: {output_root}/smear<level>/")


if __name__ == "__main__":
    main()
