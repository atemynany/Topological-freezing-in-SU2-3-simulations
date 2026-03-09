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
RELATIVE_THRESHOLD = 0.05  # Thermalized when |deriv| < 5% of initial |deriv|

# ==============================================================================
# Thermalization Detection (reusing existing logic)
# ==============================================================================

def find_thermalization_point(mc_time, plaq_data, window=SMOOTHING_WINDOW, rel_threshold=RELATIVE_THRESHOLD):
    """Find thermalization point using numerical derivative.
    
    Thermalization is detected when the smoothed derivative drops below
    rel_threshold * max(|derivative|), indicating equilibrium reached.
    """
    
    # Handle edge case: not enough data
    if len(plaq_data) < 2 * window:
        return len(plaq_data) // 5, np.array([]), np.array([])
    
    kernel = np.ones(window) / window
    plaq_smoothed = np.convolve(plaq_data, kernel, mode='valid')
    
    dP = np.diff(plaq_smoothed)
    derivative = dP  # dt=1 for consecutive sweeps
    
    deriv_smoothed = np.convolve(np.abs(derivative), kernel, mode='valid') if len(derivative) > window else np.abs(derivative)
    
    # Use relative threshold: thermalized when derivative < rel_threshold * max
    max_deriv = np.max(deriv_smoothed[:min(20, len(deriv_smoothed))])  # Max from early sweeps
    threshold = rel_threshold * max_deriv
    
    # Find first point where derivative stays below threshold
    below_threshold = deriv_smoothed < threshold
    
    # Look for sustained period below threshold (at least 5 consecutive points)
    therm_idx_smoothed = None
    consecutive = 0
    for i, below in enumerate(below_threshold):
        if below:
            consecutive += 1
            if consecutive >= 5 and therm_idx_smoothed is None:
                therm_idx_smoothed = i - 4  # Start of the sustained period
                break
        else:
            consecutive = 0
    
    if therm_idx_smoothed is None:
        # Fallback: use first crossing
        crossings = np.where(below_threshold)[0]
        therm_idx_smoothed = crossings[0] if len(crossings) > 0 else len(plaq_data) // 5
    
    # Account for offset from two convolutions
    therm_idx = therm_idx_smoothed + window
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
