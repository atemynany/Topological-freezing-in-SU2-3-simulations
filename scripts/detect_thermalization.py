#!/usr/bin/env python3
"""
Detect thermalization point from plaquette data and output conservatively rounded start_conf.

Usage:
    python3 detect_thermalization.py <plaquette_file> <save_interval>

Output:
    Prints only the rounded start_conf value (for shell script capture)

The start_conf is rounded UP to the next save_interval boundary to ensure
no correlated (pre-thermalization) data is used.

Example:
    - Thermalization detected at sweep 47, save_interval=10 → outputs 50
    - Thermalization detected at sweep 50, save_interval=10 → outputs 50
    - Thermalization detected at sweep 51, save_interval=10 → outputs 60
"""

import numpy as np
import sys
import os

# ==============================================================================
# Configuration
# ==============================================================================
SMOOTHING_WINDOW = 10
DERIVATIVE_THRESHOLD = 1e-5

# ==============================================================================
# Thermalization Detection (reusing existing logic)
# ==============================================================================

def find_thermalization_point(mc_time, plaq_data, window=SMOOTHING_WINDOW, threshold=DERIVATIVE_THRESHOLD):
    """Find thermalization point using numerical derivative."""
    
    # Handle edge case: not enough data
    if len(plaq_data) < 2 * window:
        return len(plaq_data) // 5, np.array([]), np.array([])
    
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


def round_up_to_interval(value, interval):
    """Round up to the next interval boundary.
    
    Examples with interval=10:
        47 → 50
        50 → 50
        51 → 60
    """
    return ((value + interval - 1) // interval) * interval


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 detect_thermalization.py <plaquette_file> <save_interval>", file=sys.stderr)
        sys.exit(1)
    
    plaq_file = sys.argv[1]
    save_interval = int(sys.argv[2])
    
    if not os.path.exists(plaq_file):
        print(f"Error: Plaquette file not found: {plaq_file}", file=sys.stderr)
        sys.exit(1)
    
    # Load plaquette data
    mc_time = []
    plaq_data = []
    
    with open(plaq_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    sweep = int(parts[0])
                    val = float(parts[1])
                    mc_time.append(sweep)
                    plaq_data.append(val)
                except ValueError:
                    continue
    
    if len(plaq_data) == 0:
        print("Error: No valid plaquette data found", file=sys.stderr)
        sys.exit(1)
    
    mc_time = np.array(mc_time)
    plaq_data = np.array(plaq_data)
    
    # Find thermalization point
    therm_idx, _, _ = find_thermalization_point(mc_time, plaq_data)
    therm_sweep = mc_time[therm_idx]
    
    # Round up conservatively to next save_interval boundary
    start_conf = round_up_to_interval(int(therm_sweep), save_interval)
    
    # Ensure start_conf is at least one save_interval
    start_conf = max(start_conf, save_interval)
    
    # Print to stderr for debugging
    print(f"Thermalization detected at sweep {therm_sweep}, using start_conf={start_conf}", file=sys.stderr)
    
    # Print only the value to stdout (for shell capture)
    print(start_conf)


if __name__ == "__main__":
    main()
