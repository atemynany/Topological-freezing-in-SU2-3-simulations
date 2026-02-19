#!/usr/bin/env python3
"""
Plot thermalization for SU(3) lattice gauge simulation.
Reads plaquette_su3.dat and generates thermalization plot.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

# Thermalization detection parameters
SMOOTHING_WINDOW = 10
DERIVATIVE_THRESHOLD = 1e-5

def find_thermalization_point(mc_time, plaq_data, window=SMOOTHING_WINDOW, threshold=DERIVATIVE_THRESHOLD):
    """Find thermalization point using numerical derivative."""
    kernel = np.ones(window) / window
    plaq_smoothed = np.convolve(plaq_data, kernel, mode='valid')
    mc_smoothed = mc_time[window//2 : window//2 + len(plaq_smoothed)]
    
    dt = np.diff(mc_smoothed)
    dP = np.diff(plaq_smoothed)
    derivative = dP / dt
    
    deriv_smoothed = np.convolve(derivative, kernel, mode='valid') if len(derivative) > window else derivative
    
    zero_crossings = np.where(np.abs(deriv_smoothed) < threshold)[0]
    
    if len(zero_crossings) < 2:
        median_deriv = np.median(np.abs(deriv_smoothed))
        zero_crossings = np.where(np.abs(deriv_smoothed) < median_deriv)[0]
    
    if len(zero_crossings) >= 2:
        therm_idx_in_smoothed = zero_crossings[1]
        therm_idx = therm_idx_in_smoothed + window
    else:
        therm_idx = len(plaq_data) // 5
    
    therm_idx = min(therm_idx, len(plaq_data) - 1)
    return therm_idx, derivative, deriv_smoothed


def main():
    parser = argparse.ArgumentParser(description='Plot SU(3) thermalization')
    parser.add_argument('plaq_file', nargs='?', default='output/plaquette_su3.dat',
                        help='Plaquette data file (default: output/plaquette_su3.dat)')
    parser.add_argument('-o', '--output', default='output/figures',
                        help='Output directory for figures')
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    # Load plaquette data
    print(f"Loading {args.plaq_file}...")
    mc_time = []
    plaq_data = []
    with open(args.plaq_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mc_time.append(int(parts[0]))
                    plaq_data.append(float(parts[1]))
                except ValueError:
                    continue
    
    plaq_data = np.array(plaq_data)
    mc_time = np.array(mc_time)
    
    print(f"Loaded {len(plaq_data)} measurements")
    print(f"Plaquette range: {plaq_data.min():.6f} to {plaq_data.max():.6f}")
    
    # Find thermalization
    therm_idx, derivative, deriv_smoothed = find_thermalization_point(mc_time, plaq_data)
    
    equilibrium_mean = np.mean(plaq_data[therm_idx:])
    equilibrium_std = np.std(plaq_data[therm_idx:])
    
    print(f"Thermalization at sweep {mc_time[therm_idx]} (index {therm_idx})")
    print(f"Equilibrium: {equilibrium_mean:.6f} ± {equilibrium_std:.6f}")
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), height_ratios=[2, 1], sharex=True)
    
    # Plaquette
    ax1.plot(mc_time, plaq_data, 'b-', linewidth=0.8, label=r'$\langle P \rangle$')
    ax1.axhline(y=equilibrium_mean, color='r', linestyle='--', linewidth=1.5,
               label=f'Equilibrium: {equilibrium_mean:.4f} ± {equilibrium_std:.4f}')
    ax1.fill_between(mc_time, equilibrium_mean - equilibrium_std, equilibrium_mean + equilibrium_std,
                    color='red', alpha=0.2)
    ax1.axvline(x=mc_time[therm_idx], color='g', linestyle='-', linewidth=2, alpha=0.8,
               label=f'Thermalization: {mc_time[therm_idx]} sweeps')
    ax1.set_ylabel(r'$\langle P \rangle$')
    ax1.set_title('SU(3) Thermalization')
    ax1.legend(loc='lower right')
    ax1.grid(True, alpha=0.3)
    
    # Derivative
    deriv_mc_time = mc_time[SMOOTHING_WINDOW:SMOOTHING_WINDOW + len(deriv_smoothed)]
    ax2.plot(deriv_mc_time, deriv_smoothed, 'b-', linewidth=0.8)
    ax2.axhline(y=0, color='k', linestyle='-', linewidth=1)
    ax2.axhline(y=DERIVATIVE_THRESHOLD, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax2.axhline(y=-DERIVATIVE_THRESHOLD, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax2.axvline(x=mc_time[therm_idx], color='g', linestyle='-', linewidth=2, alpha=0.8)
    ax2.set_xlabel('MC Sweeps')
    ax2.set_ylabel(r'$d\langle P \rangle / dt$')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    outfile = os.path.join(args.output, 'thermalization_su3.png')
    plt.savefig(outfile, dpi=200)
    print(f"Saved: {outfile}")
    plt.close()


if __name__ == '__main__':
    main()
