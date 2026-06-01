#!/usr/bin/env python3
"""
Main analysis script for SU(2)/SU(3) with periodic and open boundaries.
Plots susceptibility and tau_int with both boundary types in same plot.
"""

import os
import sys
import argparse
from typing import List

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from calculations import (
    RunData, AnalysisResult,
    load_run_data, analyze_run,
    lattice_spacing_su2, lattice_spacing_su3,
    topcharge_timeslice_file_path,
)
from plotting_code import (
    plot_Q_vs_mctime_grid, plot_histograms_grid,
    plot_susceptibility_combined, plot_tau_int_combined,
    plot_thermalization_comparison
)
from timeslice_analysis import (
    analyse_timeslices,
    plot_timeslice_density_grid, plot_timeslice_mctime_grid,
    plot_timeslice_edge_bulk_mctime_grid,
    plot_timeslice_susceptibility_grid, plot_timeslice_tauint_grid,
    plot_timeslice_susceptibility_cropped, plot_timeslice_tauint_cropped,
)
from action_density_analysis import analyze_action_density
from result_paths import find_run_dirs, resolve_results_dirs


def print_table(results: List[AnalysisResult], gauge_group: str):
    """Print summary table."""
    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    
    print(f"\n{'='*95}")
    print(f"  {gauge_group.upper()} Analysis Summary")
    print(f"{'='*95}")
    print(f"{'β':>6} {'a(fm)':>8} {'BC':>10} {'<Q²>':>8} {'χ_t^1/4':>10} {'τ_int':>8} {'N_conf':>8} {'α':>6}")
    print(f"{'':>6} {'':>8} {'':>10} {'':>8} {'(MeV)':>10} {'':>8} {'':>8} {'':>6}")
    print(f"{'-'*95}")
    
    for r in sorted(results, key=lambda x: (x.boundary, x.beta)):
        a = lattice_spacing(r.beta)
        chi_MeV = r.chi_t_fourth_root_a / a
        print(f"{r.beta:>6.2f} {a:>8.4f} {r.boundary:>10} {r.Q2_mean:>8.1f} {chi_MeV:>10.1f} "
              f"{r.tau_int:>8.1f} {r.n_conf:>8} {r.alpha:>6.3f}")
    print(f"{'='*95}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--su2", action="store_true")
    parser.add_argument("--su3", action="store_true")
    parser.add_argument("--results-dir", type=str, action="append", default=None)
    parser.add_argument("--output-dir", type=str, default=None)
    args = parser.parse_args()
    
    gauge_group = "su3" if args.su3 else "su2"
    
    base_dir = os.path.dirname(script_dir)
    results_dirs = resolve_results_dirs(base_dir, args.results_dir)
    output_dir = args.output_dir or os.path.join(base_dir, "output", "figures_analysis")
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n{gauge_group.upper()} Analysis")
    print("Results:")
    for results_dir in results_dirs:
        print(f"  {results_dir}")
    print(f"Output: {output_dir}\n")
    
    run_dirs = find_run_dirs(results_dirs, gauge_group)
    
    runs_periodic, runs_open = [], []
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
            runs_periodic.append(run_data)
            results_periodic.append(result)
        else:
            runs_open.append(run_data)
            results_open.append(result)
    
    print(f"Found: {len(runs_periodic)} periodic, {len(runs_open)} open\n")
    
    print_table(results_periodic + results_open, gauge_group)
    
    if runs_periodic:
        plot_Q_vs_mctime_grid(runs_periodic, output_dir, gauge_group, "periodic")
        plot_histograms_grid(runs_periodic, output_dir, gauge_group, "periodic")
    if runs_open:
        plot_Q_vs_mctime_grid(runs_open, output_dir, gauge_group, "open")
        plot_histograms_grid(runs_open, output_dir, gauge_group, "open")
    
    plot_susceptibility_combined(results_periodic, results_open, output_dir, gauge_group)
    plot_tau_int_combined(results_periodic, results_open, output_dir, gauge_group)
    
    # Thermalization comparison plot (plaquette, all runs)
    all_runs = runs_periodic + runs_open
    if all_runs:
        plot_thermalization_comparison(all_runs, output_dir, gauge_group)
        analyze_action_density(all_runs, output_dir, gauge_group)

    # Timeslice grid plots — collect results per boundary group
    lattice_spacing = lattice_spacing_su2 if gauge_group == "su2" else lattice_spacing_su3
    n_bin = 2  # group 2 consecutive slices per evaluation point (subvolume = 2*L^3)
    for runs_group, boundary in [(runs_periodic, "periodic"), (runs_open, "open")]:
        ts_files, ts_results, ts_runs = [], [], []
        for run_data in sorted(runs_group, key=lambda x: x.beta):
            ts_file = topcharge_timeslice_file_path(run_data.run_dir, gauge_group)
            if ts_file is None:
                continue
            a_fm = lattice_spacing(run_data.beta)
            open_bc = run_data.boundary == "open"
            try:
                res = analyse_timeslices(
                    ts_file,
                    lattice_spacing_fm=a_fm,
                    n_bin=n_bin,
                    open_bc=open_bc,
                    alpha=run_data.alpha,
                    ensemble=f"b{run_data.beta}_{run_data.boundary}",
                )
            except Exception as e:
                print(f"  [warn] timeslice analysis failed for {run_data.run_dir}: {e}")
                continue
            ts_files.append(ts_file)
            ts_results.append(res)
            ts_runs.append(run_data)

        if ts_runs:
            plot_timeslice_density_grid(ts_results, ts_runs, output_dir, gauge_group, boundary)
            plot_timeslice_mctime_grid(ts_files, ts_runs, n_bin, output_dir, gauge_group, boundary)
            plot_timeslice_edge_bulk_mctime_grid(ts_files, ts_runs, output_dir, gauge_group, boundary)
            plot_timeslice_susceptibility_grid(ts_results, ts_runs, output_dir, gauge_group, boundary)
            plot_timeslice_tauint_grid(ts_results, ts_runs, output_dir, gauge_group, boundary)
            plot_timeslice_susceptibility_cropped(ts_results, ts_runs, output_dir, gauge_group, boundary)
            plot_timeslice_tauint_cropped(ts_results, ts_runs, output_dir, gauge_group, boundary)

    print(f"\nDone! Plots in: {output_dir}")


if __name__ == "__main__":
    main()
