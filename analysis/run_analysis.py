#!/usr/bin/env python3
"""
Plot thermalization and topological charge analysis for SU(2) lattice gauge simulation.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

from Q_estimator import QEstimator
from topological_susceptibility import topological_susceptibility

# ==============================================================================
# Configuration Parameters
# ==============================================================================

# Paths
base_dir = "/home/alex/Desktop/workspace/master_thesis/lattice_qcd_topolgical_charge/Topological_charge_my_implementation"
output_dir = os.path.join(base_dir, "output")
plaq_file = os.path.join(output_dir, "plaquette.dat")
topcharge_file = os.path.join(output_dir, "topcharge.dat")

# Histogram binning
HISTOGRAM_BINS = 35  # Number of bins for Q histogram (can also use bin edges array)

# Thermalization detection parameters
SMOOTHING_WINDOW = 10  # Window size for smoothing the derivative
DERIVATIVE_THRESHOLD = 1e-5  # Threshold for "zero" derivative

# ==============================================================================

def find_thermalization_point(mc_time, plaq_data, window=SMOOTHING_WINDOW, threshold=DERIVATIVE_THRESHOLD):
    """
    Find thermalization point by computing numerical derivative and finding
    where it becomes approximately zero (equilibrium reached).
    
    Returns the second point where |derivative| < threshold (to skip initial noise).
    """
    
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
        # Fall back: find where derivative magnitude is below median
        median_deriv = np.median(np.abs(deriv_smoothed))
        zero_crossings = np.where(np.abs(deriv_smoothed) < median_deriv)[0]
    
    if len(zero_crossings) >= 2:
        # Take the second point where derivative is ~0 (skip first to avoid initial fluctuations)
        therm_idx_in_smoothed = zero_crossings[1]
        # Map back to original array index
        therm_idx = therm_idx_in_smoothed + window
    else:
        # Fallback to 20% if no clear thermalization
        therm_idx = len(plaq_data) // 5
    
    # Ensure we don't exceed array bounds
    therm_idx = min(therm_idx, len(plaq_data) - 1)
    
    return therm_idx, derivative, deriv_smoothed

# ==============================================================================

# Create figures directory
fig_dir = os.path.join(output_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

# ============================================================================
# 1. Plot Thermalization (Plaquette vs MC time)
# ============================================================================
print("Loading plaquette data...")

mc_time = []
plaq_data = []
with open(plaq_file, 'r') as f:
    for line in f:
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

print(f"Loaded {len(plaq_data)} plaquette measurements")
print(f"Mean plaquette: {np.mean(plaq_data):.6f}")
print(f"Std plaquette: {np.std(plaq_data):.6f}")

# Find thermalization point using numerical derivative
therm_result = find_thermalization_point(mc_time, plaq_data)
therm_idx = therm_result[0]
derivative = therm_result[1]
deriv_smoothed = therm_result[2]

# Calculate equilibrium statistics after thermalization
equilibrium_mean = np.mean(plaq_data[therm_idx:])
equilibrium_std = np.std(plaq_data[therm_idx:])

print(f"Thermalization point detected at sweep {mc_time[therm_idx]} (index {therm_idx})")
print(f"Equilibrium value: {equilibrium_mean:.6f} ± {equilibrium_std:.6f}")

# Plot thermalization with derivative subplot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), height_ratios=[2, 1], sharex=True)

# Top panel: Plaquette vs MC time
ax1.plot(mc_time, plaq_data, 'b-', linewidth=0.8, label=r'$\langle P_{\mu\nu} \rangle$')
ax1.axhline(y=equilibrium_mean, color='r', linestyle='--', linewidth=1.5, 
           label=f'Equilibrium: {equilibrium_mean:.4f} ± {equilibrium_std:.4f}')
ax1.fill_between(mc_time, equilibrium_mean - equilibrium_std, equilibrium_mean + equilibrium_std,
                color='red', alpha=0.2)
ax1.axvline(x=mc_time[therm_idx], color='g', linestyle='-', linewidth=2, alpha=0.8,
           label=f'Thermalization: {mc_time[therm_idx]} sweeps')
ax1.set_ylabel(r'$\langle P_{\mu\nu} \rangle$')
ax1.legend(loc='upper right')
ax1.grid(True, alpha=0.3)

# Bottom panel: Derivative of plaquette
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
print(f"Saved thermalization plot to {fig_dir}/thermalization.png")
plt.close()

# ============================================================================
# 2. Load and process topological charge data
# ============================================================================
print("\nLoading topological charge data...")

smear_steps = []
config_nums = []
Q_values = []
plaq_values = []

with open(topcharge_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 4:
            smear = int(parts[0])
            conf = int(parts[1])
            Q = float(parts[2])
            plaq = float(parts[3])
            smear_steps.append(smear)
            config_nums.append(conf)
            Q_values.append(Q)
            plaq_values.append(plaq)

smear_steps = np.array(smear_steps)
config_nums = np.array(config_nums)
Q_values = np.array(Q_values)
plaq_values = np.array(plaq_values)

print(f"Loaded {len(Q_values)} measurements")

unique_smear = np.unique(smear_steps)
unique_configs = np.unique(config_nums)

print(f"Smearing steps: {unique_smear}")
print(f"Configurations: {unique_configs[0]} to {unique_configs[-1]} ({len(unique_configs)} configs)")

# ============================================================================
# 3. Plot Q vs Monte Carlo Time (at fixed smearing)
# ============================================================================
print("\nPlotting Q vs MC time...")

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

# Show mean with jackknife error as horizontal band
# ax.axhline(Q_mean, color='r', linestyle='-', linewidth=1.5)
# ax.fill_between(configs, Q_mean - Q_jack_error, Q_mean + Q_jack_error, 
#                color='red', alpha=0.3, label=rf'$\langle Q \rangle = {Q_mean:.2f} \pm {Q_jack_error:.2f}$ (JK)')

# Integer Q reference lines
for q_int in [-3, -2, -1, 0, 1, 2, 3]:
    ax.axhline(q_int, color='gray', ls='--', lw=1, alpha=0.5)

ax.set_xlabel('Monte Carlo Time (Sweeps)')
ax.set_ylabel(r'$Q_{top}$')
ax.legend(loc='best')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_vs_MCtime.png'), dpi=200)
print(f"Saved Q vs MC time plot to {fig_dir}/Q_vs_MCtime.png")
plt.close()

# ============================================================================
# 4. Q Estimator: Rescaled charge with optimal rescaling
# ============================================================================
print("\nComputing Q estimator with optimal rescaling...")

# Use the Q estimator module
est = QEstimator(Q)
alpha_opt = est.alpha
Q_scaled = est.Q_scaled
Q_rescaled = est.Q_rescaled

print(f"Optimal α = {alpha_opt:.6f}")
print(f"RMS deviation from integer: {est.result.rms_deviation:.6f}")

# Plot Q_scaled (continuous, no rounding) vs MC time
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(configs, Q_scaled, marker='o', linestyle='', markersize=5, color='steelblue',
        label=rf'$\alpha \hat{{Q}}$, $\alpha={alpha_opt:.4f}$')

# Integer Q reference lines
for q_int_val in [-3, -2, -1, 0, 1, 2, 3]:
    ax.axhline(q_int_val, color='gray', ls='--', lw=1, alpha=0.5)

ax.set_xlabel('Monte Carlo Time (Sweeps)')
ax.set_ylabel(r'$\alpha \cdot Q$')
ax.legend(loc='best')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_scaled_vs_MCtime.png'), dpi=200)
print(f"Saved Q scaled (continuous) plot to {fig_dir}/Q_scaled_vs_MCtime.png")
plt.close()

# Plot Q_rescaled (rounded to integers) vs MC time
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(configs, Q_rescaled, marker='o', linestyle='', markersize=5, color='steelblue',
        label=rf'$Q_{{rescaled}} = \mathrm{{round}}(\alpha \hat{{Q}})$, $\alpha={alpha_opt:.4f}$')

# Integer Q reference lines
for q_int_val in [-3, -2, -1, 0, 1, 2, 3]:
    ax.axhline(q_int_val, color='gray', ls='--', lw=1, alpha=0.5)

ax.set_xlabel('Monte Carlo Time (Sweeps)')
ax.set_ylabel(r'$Q_{rescaled}$')
ax.legend(loc='best')
ax.grid(True, alpha=0.3)
ax.set_yticks(range(-4, 5))
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_rescaled_vs_MCtime.png'), dpi=200)
print(f"Saved Q rescaled plot to {fig_dir}/Q_rescaled_vs_MCtime.png")
plt.close()

# ============================================================================
# 5. Histogram of Q values at optimal smearing
# ============================================================================
print("\nPlotting Q histogram...")

fig, ax = plt.subplots(figsize=(7, 5))

Q_optimal = Q_values[smear_steps == optimal_smear]

ax.hist(Q_optimal, bins=HISTOGRAM_BINS, alpha=0.7, 
        color='steelblue', edgecolor='black')
ax.axvline(x=0, color='k', linestyle='-', alpha=0.5)

# Mark integer values
for q_int in [-2, -1, 1, 2]:
    ax.axvline(x=q_int, color='gray', ls='--', lw=1, alpha=0.5)

ax.set_xlabel(r'$Q_{top}$')
ax.set_ylabel(r'$N_{conf}$')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_histogram.png'), dpi=200)
print(f"Saved Q histogram to {fig_dir}/Q_histogram.png")
plt.close()

# Histogram of scaled Q (continuous, no rounding)
fig, ax = plt.subplots(figsize=(7, 5))

ax.hist(Q_scaled, bins=HISTOGRAM_BINS, alpha=0.7, 
        color='steelblue', edgecolor='black')
ax.axvline(x=0, color='k', linestyle='-', alpha=0.5)

for q_int in [-2, -1, 1, 2]:
    ax.axvline(x=q_int, color='gray', ls='--', lw=1, alpha=0.5)

ax.set_xlabel(r'$\alpha \cdot Q$')
ax.set_ylabel(r'$N_{conf}$')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_scaled_histogram.png'), dpi=200)
print(f"Saved Q scaled histogram to {fig_dir}/Q_scaled_histogram.png")
plt.close()

# Histogram of rescaled Q estimator (rounded to integers)
fig, ax = plt.subplots(figsize=(7, 5))
unique_Qi, counts_Qi = np.unique(Q_rescaled, return_counts=True)

ax.hist(Q_rescaled, bins=HISTOGRAM_BINS, alpha=0.7, color='steelblue', edgecolor='black')
ax.axvline(x=0, color='k', linestyle='-', alpha=0.5)

for q_int in [-2, -1, 1, 2]:
    ax.axvline(x=q_int, color='gray', ls='--', lw=1, alpha=0.5)

ax.set_xlabel(r'$Q_{res}$')
ax.set_ylabel(r'$N_{conf}$')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'Q_rescaled_histogram.png'), dpi=200)
print(f"Saved Q rescaled histogram to {fig_dir}/Q_rescaled_histogram.png")
plt.close()

# ============================================================================
# 6. Summary Statistics
# ============================================================================
print("\n" + "="*50)
print("Summary Statistics")
print("="*50)

for smear in unique_smear:
    mask = smear_steps == smear
    Q = Q_values[mask]
    print(f"Smearing = {smear:2d}:  <Q> = {np.mean(Q):8.4f} ± {np.std(Q):6.4f}  "
          f"(min={np.min(Q):7.4f}, max={np.max(Q):7.4f})")

print("-"*50)
print("Q Estimator (rescaled charge with optimal α):")
print(f"  α_optimal = {alpha_opt:.6f}")
print(f"  <Q_rescaled>  = {np.mean(Q_rescaled):.4f}")
print(f"  <Q_rescaled²> = {np.mean(Q_rescaled**2):.4f}")
print(f"  Distribution: ", end="")
for qi, c in zip(unique_Qi, counts_Qi):
    print(f"Q={qi:+d}:{c} ", end="")
print()

# Topological susceptibility (requires lattice spacing)
# TODO: Set your lattice spacing here (in fm)
T_lattice = 16  # temporal extent
L_lattice = 16  # spatial extent  
a_fm = 0.1      # lattice spacing in fm (adjust for your simulation!)

Q2_mean = np.mean(Q_rescaled**2)
chi = topological_susceptibility(Q2_mean, T_lattice, L_lattice, a_fm)

print("-"*50)
print("Topological Susceptibility:")
print(f"  Lattice: {T_lattice} x {L_lattice}^3, a = {a_fm} fm")
print(f"  <Q²> = {Q2_mean:.4f}")
print(f"  χ_t (lattice) = {chi['chi_t_lattice']:.6e}")
print(f"  χ_t^(1/4) = {chi['chi_t_fourth_root_MeV']:.1f} MeV")
print(f"  Reference SU(2): 200(15) MeV")

print("="*50)
print(f"\nFigures saved to: {fig_dir}")
print("Done!")
