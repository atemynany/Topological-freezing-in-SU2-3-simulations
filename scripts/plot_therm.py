#!/usr/bin/env python3
"""
Plot thermalization for a single run directory.
Produces two separate figures (plaquette and Q_re vs MC sweep).

Usage:
    python3 scripts/plot_therm.py <run_dir> [--su3] [--output-dir <dir>]

Reads:
    <run_dir>/output/plaquette.dat           (SU2)
    <run_dir>/output/therm_topcharge.dat     (SU2)
  or
    <run_dir>/output/plaquette_su3.dat       (SU3)
    <run_dir>/output/therm_topcharge_su3.dat (SU3)

Saves (default: output/figures_analysis/):
    therm_plaquette_<run_name>.png
    therm_topcharge_<run_name>.png
"""

import sys
import os
import argparse

script_dir   = os.path.dirname(os.path.abspath(__file__))
project_dir  = os.path.dirname(script_dir)
analysis_dir = os.path.join(project_dir, "analysis")
sys.path.insert(0, analysis_dir)

from plotting_code import plot_therm_topcharge


def main():
    parser = argparse.ArgumentParser(description="Plot thermalization for a run directory.")
    parser.add_argument("run_dir", help="Path to run directory")
    parser.add_argument("--su3", action="store_true", help="SU(3) run (default: SU(2))")
    parser.add_argument("--output-dir", default=None,
                        help="Where to save plots (default: output/figures_analysis/)")
    args = parser.parse_args()

    run_dir      = args.run_dir
    run_name     = os.path.basename(os.path.normpath(run_dir))
    output_subdir = os.path.join(run_dir, "output")

    output_dir = args.output_dir or os.path.join(project_dir, "output", "figures_analysis")

    if args.su3:
        plaq_file    = os.path.join(output_subdir, "plaquette_su3.dat")
        therm_q_file = os.path.join(output_subdir, "therm_topcharge_su3.dat")
    else:
        plaq_file    = os.path.join(output_subdir, "plaquette.dat")
        therm_q_file = os.path.join(output_subdir, "therm_topcharge.dat")

    if not os.path.isfile(plaq_file) and not os.path.isfile(therm_q_file):
        print(f"Error: neither plaquette nor therm_topcharge file found in {output_subdir}")
        sys.exit(1)

    gauge_group = "su3" if args.su3 else "su2"
    plot_therm_topcharge(plaq_file, therm_q_file, output_dir, run_name=run_name, gauge_group=gauge_group)


if __name__ == "__main__":
    main()
