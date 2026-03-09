# SU(2)/SU(3) Topological Charge Freezing Study

Monte Carlo simulation of topological charge in SU(2) and SU(3) lattice gauge theory to study topological freezing at fine lattice spacings.

## Build

```bash
./scripts/build.sh release
```

## Main Workflow: Beta Scan

Run simulations across multiple β values and analyze:

```bash
# 1. Edit beta values to scan
nano input/beta_scan_su2.txt

# 2. Run scan (generates configs + measures Q)
./scripts/run_beta_scan.sh --su2

# 3. Analyze results
conda activate master_thesis
python analysis/analysis.py --su2
```

Results saved to `data/results/` and plots to `output/figures_analysis/`.

## Key Executables

| Binary | Purpose |
|--------|---------|
| `mc_heatbath` | SU(2) heatbath config generation |
| `mc_heatbath_su3` | SU(3) heatbath config generation |
| `meas_topcharge` | SU(2) topological charge measurement |
| `meas_topcharge_su3` | SU(3) topological charge measurement |

## Input Files

- `input/base_params_su2.txt` - Simulation parameters (T, L, sweeps, smearing)
- `input/beta_scan_su2.txt` - List of β values to scan

## Output

- `data/results/T{T}_L{L}_b{beta}_seed{seed}/` - Per-run results
  - `configs/` - Gauge configurations
  - `output/topcharge.dat` - Q measurements
  - `run_info.txt` - Run metadata

## Analysis

Key outputs from `analysis.py`:
- **χ_t^(1/4) vs β** - Susceptibility (should match ~200 MeV for SU(2))
- **τ_int(Q²) vs β** - Autocorrelation time (increases with β → freezing)
- **Q histograms** - Distribution narrowing indicates freezing

## Requirements

- CMake ≥ 3.16, C++17 compiler
- Python: numpy, matplotlib, pyerrors (conda env: `master_thesis`)
