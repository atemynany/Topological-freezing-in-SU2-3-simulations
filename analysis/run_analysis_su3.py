#!/usr/bin/env python3
"""
Plot thermalization and topological charge analysis for SU(3) lattice gauge simulation.

Equivalent to run_analysis.py but for SU(3) data formats:
- plaquette_su3.dat: sweep plaquette
- topcharge_su3.dat: conf smear_step Q
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from Q_estimator import QEstimator
from topological_susceptibility import topological_susceptibility

# ==============================================================================
# Configuration Parameters
# ==============================================================================

# Paths - use script directory to find project root
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(script_dir)  # Parent of analysis/ is project root
output_dir = os.path.join(base_dir, "output")
plaq_file = os.path.join(output_dir, "plaquette_su3.dat")
topcharge_file = os.path.join(output_dir, "topcharge_su3.dat")

# Allow command line override
if len(sys.argv) >= 2:
    plaq_file = sys.argv[1]
if len(sys.argv) >= 3:
    topcharge_file = sys.argv[2]

# Histogram binning
HISTOGRAM_BINS = 35

# Thermalization detection parameters
SMOOTHING_WINDOW = 10
DERIVATIVE_THRESHOLD = 1e-5

# ==============================================================================

def find_thermalization_point(mc_time, plaq_data, window=SMOOTHING_WINDOW, threshold=DERIVATIVE_THRESHOLD):
    """
    Find thermalization point by computing numerical derivative and finding
    where it becomes approximately zero (equilibrium reached).
    """
    if len(plaq_data) < 2 * window:
        return len(plaq_data) // 5, np.array([]), np.array([])
    
    # Moving average smoothing
    kernel = np.ones(window) / window
    plaq_smoothed = np.convolve(plaq_data, kernel, mode='valid')
    mc_smoothed = mc_time[window//2 : window//2 + len(plaq_smoothed)]
    
    # Compute numerical derivative: dP/dt
    dt = np.diff(mc_smoothed)
    dP = np.diff(plaq_smoothed)
    derivative = dP / dt
    
    # Smooth the derivative as well
    deriv_smoothed = np.convolve(derivative, kernel, mode='valid') if len(derivative) > window else derivative
    
    # Find points where derivative is approximately zero
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

# ==============================================================================

# Create figures directory for SU(3)
fig_dir = os.path.join(output_dir, "figures_su3")
os.makedirs(fig_dir, exist_ok=True)

print("="*60)
print("SU(3) Lattice Gauge Theory Analysis")
print("="*60)

# ============================================================================
# 1. Plot Thermalization (Plaquette vs MC time)
# ============================================================================
print("\n[1/5] Loading plaquette data...")

mc_time = []
plaq_data = []

if os.path.exists(plaq_file):
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

    plaq_data = np.array(plaq_data)
    mc_time = np.array(mc_time)

    print(f"      Loaded {len(plaq_data)} plaquette measurements")
    print(f"      Mean plaquette: {np.mean(plaq_data):.6f}")
    print(f"      Std plaquette: {np.std(plaq_data):.6f}")

    # Find thermalization point
    therm_idx, derivative, deriv_smoothed = find_thermalization_point(mc_time, plaq_data)

    # Calculate equilibrium statistics after thermalization
    equilibrium_mean = np.mean(plaq_data[therm_idx:])
    equilibrium_std = np.std(plaq_data[therm_idx:])

    print(f"      Thermalization at sweep {mc_time[therm_idx]} (index {therm_idx})")
    print(f"      Equilibrium value: {equilibrium_mean:.6f} ± {equilibrium_std:.6f}")

    # Plot thermalization with derivative subplot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), height_ratios=[2, 1], sharex=True)

    ax1.plot(mc_time, plaq_data, 'b-', linewidth=0.8, label=r'$\langle P_{\mu\nu} \rangle$')
    ax1.axhline(y=equilibrium_mean, color='r', linestyle='--', linewidth=1.5, 
               label=f'Equilibrium: {equilibrium_mean:.4f} ± {equilibrium_std:.4f}')
    ax1.fill_between(mc_time, equilibrium_mean - equilibrium_std, equilibrium_mean + equilibrium_std,
                    color='red', alpha=0.2)
    ax1.axvline(x=mc_time[therm_idx], color='g', linestyle='-', linewidth=2, alpha=0.8,
               label=f'Thermalization: {mc_time[therm_idx]} sweeps')
    ax1.set_ylabel(r'$\langle P_{\mu\nu} \rangle$')
    ax1.set_title('SU(3) Thermalization')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)

    if len(deriv_smoothed) > 0:
        deriv_mc_time = mc_time[SMOOTHING_WINDOW:SMOOTHING_WINDOW + len(deriv_smoothed)]
        ax2.plot(deriv_mc_time, deriv_smoothed, 'b-', linewidth=0.8, label=r'$d\langle P \rangle / dt$ (smoothed)')
        ax2.axhline(y=0, color='k', linestyle='-', linewidth=1)
        ax2.axhline(y=DERIVATIVE_THRESHOLD, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        ax2.axhline(y=-DERIVATIVE_THRESHOLD, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        ax2.axvline(x=mc_time[therm_idx], color='g', linestyle='-', linewidth=2, alpha=0.8)
    ax2.set_xlabel('Monte Carlo Sweeps')
    ax2.set_ylabel(r'$d\langle P \rangle / dt$')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right')   

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'thermalization.png'), dpi=200)
    print(f"      Saved: {fig_dir}/thermalization.png")
    plt.close()
else:
    print(f"      Warning: {plaq_file} not found, skipping thermalization plot")

# ============================================================================
# 2. Load and process topological charge data
# ============================================================================
print("\n[2/5] Loading topological charge data...")

config_nums = []
smear_steps = []
Q_values = []

if not os.path.exists(topcharge_file):
    print(f"      Error: {topcharge_file} not found")
    sys.exit(1)

with open(topcharge_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split()
        # SU(3) format: conf smear_step Q
        if len(parts) >= 3:
            conf = int(parts[0])
            smear = int(parts[1])
            Q = float(parts[2])
            config_nums.append(conf)
            smear_steps.append(smear)
            Q_values.append(Q)

config_nums = np.array(config_nums)
smear_steps = np.array(smear_steps)
Q_values = np.array(Q_values)

print(f"      Loaded {len(Q_values)} measurements")

unique_smear = np.unique(smear_steps)
unique_configs = np.unique(config_nums)

print(f"      Smearing steps: {unique_smear}")
print(f"      Configurations: {unique_configs[0]} to {unique_configs[-1]} ({len(unique_configs)} configs)")

# ============================================================================
# 3. Plot Q vs Monte Carlo Time (at fixed smearing)
# ============================================================================
print("\n[3/5] Plotting Q vs MC time...")

# Use maximum smearing available
optimal_smear = unique_smear[-1]
mask = smear_steps == optimal_smear
configs = config_nums[mask]
Q = Q_values[mask]

# Sort by config number
sort_idx = np.argsort(configs)
configs = configs[sort_idx]
Q = Q[sort_idx]

# Calculate statistics
n = len(Q)
Q_mean = np.mean(Q)
Q_std = np.std(Q, ddof=1)

# Jackknife error on the mean
theta_jack = np.array([np.mean(np.delete(Q, i)) for i in range(n)])
Q_jack_error = np.sqrt((n - 1) / n * np.sum((theta_jack - np.mean(theta_jack))**2))

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(configs, Q, marker='o', linestyle='', markersize=4, color='steelblue',
        label=rf'$Q_{{top}}$ (n={n} configs)')

# Integer Q reference lines
for q_int in [-3, -2, -1, 0, 1, 2, 3]:
    ax.axhline(q_int, color='gray', ls='--', lw=1, alpha=0.5)

ax.set_xlabel('Monte Carlo Time (Sweeps)')
ax.set_ylabel(r'$Q_{top}$')
ax.set_title(f'SU(3) Topological Charge (smearing={optimal_smear})')
ax.legend(loc='best')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_vs_MCtime.png'), dpi=200)
print(f"      Saved: {fig_dir}/Q_vs_MCtime.png")
plt.close()

# ============================================================================
# 4. Q Estimator: Rescaled charge with optimal rescaling
# ============================================================================
print("\n[4/5] Computing Q estimator with optimal rescaling...")

# Use the Q estimator module
est = QEstimator(Q)
alpha_opt = est.alpha
Q_scaled = est.Q_scaled
Q_rescaled = est.Q_rescaled

print(f"      Optimal α = {alpha_opt:.6f}")
print(f"      RMS deviation from integer: {est.result.rms_deviation:.6f}")

# Plot Q_scaled (continuous, no rounding) vs MC time
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(configs, Q_scaled, marker='o', linestyle='', markersize=5, color='steelblue',
        label=rf'$\alpha \hat{{Q}}$, $\alpha={alpha_opt:.4f}$')

for q_int_val in [-3, -2, -1, 0, 1, 2, 3]:
    ax.axhline(q_int_val, color='gray', ls='--', lw=1, alpha=0.5)

ax.set_xlabel('Monte Carlo Time (Sweeps)')
ax.set_ylabel(r'$\alpha \cdot Q$')
ax.set_title('SU(3) Scaled Topological Charge')
ax.legend(loc='best')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_scaled_vs_MCtime.png'), dpi=200)
print(f"      Saved: {fig_dir}/Q_scaled_vs_MCtime.png")
plt.close()

# Plot Q_rescaled (rounded to integers) vs MC time
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(configs, Q_rescaled, marker='o', linestyle='', markersize=5, color='steelblue',
        label=rf'$Q_{{rescaled}} = \mathrm{{round}}(\alpha \hat{{Q}})$, $\alpha={alpha_opt:.4f}$')

for q_int_val in [-3, -2, -1, 0, 1, 2, 3]:
    ax.axhline(q_int_val, color='gray', ls='--', lw=1, alpha=0.5)

ax.set_xlabel('Monte Carlo Time (Sweeps)')
ax.set_ylabel(r'$Q_{rescaled}$')
ax.set_title('SU(3) Rescaled Topological Charge')
ax.legend(loc='best')
ax.grid(True, alpha=0.3)
ax.set_yticks(range(-4, 5))
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_rescaled_vs_MCtime.png'), dpi=200)
print(f"      Saved: {fig_dir}/Q_rescaled_vs_MCtime.png")
plt.close()

# ============================================================================
# 5. Histograms of Q values
# ============================================================================
print("\n[5/5] Plotting Q histograms...")

# Raw Q histogram
fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(Q, bins=HISTOGRAM_BINS, alpha=0.7, color='steelblue', edgecolor='black')
ax.axvline(x=0, color='k', linestyle='-', alpha=0.5)
for q_int in [-2, -1, 1, 2]:
    ax.axvline(x=q_int, color='gray', ls='--', lw=1, alpha=0.5)
ax.set_xlabel(r'$Q_{top}$')
ax.set_ylabel(r'$N_{conf}$')
ax.set_title('SU(3) Topological Charge Distribution')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_histogram.png'), dpi=200)
print(f"      Saved: {fig_dir}/Q_histogram.png")
plt.close()

# Scaled Q histogram
fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(Q_scaled, bins=HISTOGRAM_BINS, alpha=0.7, color='steelblue', edgecolor='black')
ax.axvline(x=0, color='k', linestyle='-', alpha=0.5)
for q_int in [-2, -1, 1, 2]:
    ax.axvline(x=q_int, color='gray', ls='--', lw=1, alpha=0.5)
ax.set_xlabel(r'$\alpha \cdot Q$')
ax.set_ylabel(r'$N_{conf}$')
ax.set_title('SU(3) Scaled Topological Charge Distribution')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_scaled_histogram.png'), dpi=200)
print(f"      Saved: {fig_dir}/Q_scaled_histogram.png")
plt.close()

# Rescaled Q histogram
fig, ax = plt.subplots(figsize=(7, 5))
unique_Qi, counts_Qi = np.unique(Q_rescaled, return_counts=True)
ax.hist(Q_rescaled, bins=HISTOGRAM_BINS, alpha=0.7, color='steelblue', edgecolor='black')
ax.axvline(x=0, color='k', linestyle='-', alpha=0.5)
for q_int in [-2, -1, 1, 2]:
    ax.axvline(x=q_int, color='gray', ls='--', lw=1, alpha=0.5)
ax.set_xlabel(r'$Q_{res}$')
ax.set_ylabel(r'$N_{conf}$')
ax.set_title('SU(3) Rescaled Topological Charge Distribution')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_rescaled_histogram.png'), dpi=200)
print(f"      Saved: {fig_dir}/Q_rescaled_histogram.png")
plt.close()

# ============================================================================
# 6. Summary Statistics
# ============================================================================
print("\n" + "="*60)
print("Summary Statistics")
print("="*60)

for smear in unique_smear:
    mask_s = smear_steps == smear
    Q_s = Q_values[mask_s]
    print(f"Smearing = {smear:2d}:  <Q> = {np.mean(Q_s):8.4f} ± {np.std(Q_s):6.4f}  "
          f"(min={np.min(Q_s):7.4f}, max={np.max(Q_s):7.4f})")

print("-"*60)
print("Q Estimator (rescaled charge with optimal α):")
print(f"  α_optimal = {alpha_opt:.6f}")
print(f"  <Q_rescaled>  = {np.mean(Q_rescaled):.4f}")
print(f"  <Q_rescaled²> = {np.mean(Q_rescaled**2):.4f}")
print(f"  Distribution: ", end="")
for qi, c in zip(unique_Qi, counts_Qi):
    print(f"Q={int(qi):+d}:{c} ", end="")
print()

# Topological susceptibility (requires lattice spacing)
# TODO: Set your lattice parameters here
T_lattice = 8   # temporal extent (adjust for your simulation!)
L_lattice = 8   # spatial extent  (adjust for your simulation!)
a_fm = 0.1      # lattice spacing in fm (adjust for your β!)

Q2_mean = np.mean(Q_rescaled**2)
chi = topological_susceptibility(Q2_mean, T_lattice, L_lattice, a_fm)

print("-"*60)
print("Topological Susceptibility:")
print(f"  Lattice: {T_lattice} x {L_lattice}^3, a = {a_fm} fm")
print(f"  <Q²> = {Q2_mean:.4f}")
print(f"  χ_t (lattice) = {chi['chi_t_lattice']:.6e}")
print(f"  χ_t^(1/4) = {chi['chi_t_fourth_root_MeV']:.1f} MeV")
print(f"  Reference SU(3): ~175 MeV (Sommer r0 scale)")

print("="*60)
print(f"\nFigures saved to: {fig_dir}")
print("Done!")
