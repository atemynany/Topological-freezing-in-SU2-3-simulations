#!/usr/bin/env python3
"""Plot thermalization progress from running simulations."""

import argparse
import glob
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_dat(path):
    """Load a .dat file, skipping comment lines."""
    try:
        data = np.loadtxt(path, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except Exception:
        return None


def plot_run(run_dir, output_dir):
    """Plot plaquette and therm_topcharge for a single run."""
    plaq_path = os.path.join(run_dir, 'plaquette.dat')
    qtherm_path = os.path.join(run_dir, 'therm_topcharge.dat')

    plaq = load_dat(plaq_path)
    qtherm = load_dat(qtherm_path)

    if plaq is None or len(plaq) < 2:
        return

    run_name = os.path.basename(run_dir)
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    axes[0].plot(plaq[:, 0], plaq[:, 1], 'b-', linewidth=0.5)
    axes[0].set_ylabel('Plaquette')
    axes[0].set_title(run_name)

    if qtherm is not None and len(qtherm) > 1:
        axes[1].plot(qtherm[:, 0], qtherm[:, 1], 'r-', linewidth=0.5)
    axes[1].set_ylabel('Q (unsmeared)')
    axes[1].set_xlabel('Sweep')

    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f'{run_name}.png')
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f'Saved: {out_path}')


def main():
    parser = argparse.ArgumentParser(description='Monitor running simulations')
    parser.add_argument('--run-dir', type=str, help='Specific run directory')
    parser.add_argument('--beta', type=str, help='Filter by beta value')
    parser.add_argument('--data-dir', type=str, default='data/results',
                        help='Base data directory')
    parser.add_argument('--output-dir', type=str, default='output/monitoring',
                        help='Output directory for plots')
    args = parser.parse_args()

    if args.run_dir:
        run_dirs = [args.run_dir]
    else:
        pattern = os.path.join(args.data_dir, 'T*_L*_b*_*_seed*')
        run_dirs = sorted(glob.glob(pattern))

    if args.beta:
        run_dirs = [d for d in run_dirs if f'_b{args.beta}_' in d]

    if not run_dirs:
        print('No runs found.')
        return

    for run_dir in run_dirs:
        if os.path.exists(os.path.join(run_dir, 'plaquette.dat')):
            plot_run(run_dir, args.output_dir)


if __name__ == '__main__':
    main()
