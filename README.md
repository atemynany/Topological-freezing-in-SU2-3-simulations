# SU(2)/SU(3) Topological Charge Freezing Study

Monte Carlo study of topological charge freezing in SU(2) and SU(3) pure-gauge lattice field theory. Generates gauge configurations via heatbath, measures the topological charge with APE smearing, and extracts χ_t and τ_int across multiple β values and boundary conditions.

Works on **Linux** and **macOS** (including Apple Silicon).

## Requirements

- CMake ≥ 3.16
- C++17 compiler (GCC, Clang, or Apple Clang)
- Python 3: numpy, matplotlib, scipy, pyerrors
- Rust + maturin *(optional — fast analysis extension)*

### macOS extras

```bash
xcode-select --install   # clang + make
brew install cmake
brew install libomp      # optional OpenMP
```

## Build

```bash
./scripts/build.sh release
```

## Main Workflow

```bash
# run full beta scan (generates configs, measures Q, runs analysis)
./run_simulation.sh --su2
./run_simulation.sh --su3
./run_simulation.sh --su2 --open    # open boundary conditions
./run_simulation.sh --dry-run       # preview without running

# or step by step
./scripts/run_beta_scan.sh --su2
python3 analysis/analysis.py --su2
```

Results go to `data/results/` and plots to `output/figures_analysis/`.

## Key Executables

| Binary | Purpose |
|--------|---------|
| `mc_heatbath` | SU(2) heatbath + thermalization monitoring |
| `mc_heatbath_su3` | SU(3) heatbath + thermalization monitoring |
| `meas_topcharge` | SU(2) topological charge with APE smearing |
| `meas_topcharge_su3` | SU(3) topological charge with APE smearing |
| `compute_instanton_Q` | Instanton topological charge test |

## Input Files

- `input/base_params_su2.txt` — SU(2) simulation parameters
- `input/base_params_su3.txt` — SU(3) simulation parameters
- `input/beta_scan_su2.txt` — β values to scan
- `input/beta_scan_su3.txt` — SU(3) β values

## Data Pipeline

```
mc_heatbath[_su3] -i input.txt
   ├── every sweep:  raw Q (unsmeared) → output/therm_topcharge[_su3].dat
   ├── every sweep:  plaquette         → output/plaquette[_su3].dat
   └── every N:      save config       → output/configs[_su3]/conf*.

meas_topcharge[_su3] -i input.txt
   └── load config → APE smear N steps → Q → output/topcharge[_su3].dat

analysis/analysis.py --su2/--su3
   └── load topcharge files → alpha rescaling → χ_t, τ_int, plots
```

## Thermalization Monitoring

The heatbath records the **unsmeared** topological charge at every MC sweep. This is a qualitative indicator — without smearing the raw Q values are not integer-quantized (UV noise dominates for SU(3)), but the time series reveals whether the system is thermalizing. The actual physics (histograms, χ_t, τ_int) is computed from the **smeared** Q values in the production measurement step.

Plots are saved to `output/figures_analysis/`:
- `therm_topcharge_su2_<run>.png` / `therm_topcharge_su3_<run>.png` — raw Q vs MC sweep
- `thermalization_comparison_su2/3.png` — plaquette vs MC sweep (all runs overlaid)

Use these to visually confirm thermalization before relying on automated detection.

## Analysis Outputs

All plots go to `output/figures_analysis/`:

| Plot | Description |
|------|-------------|
| `Q_vs_mctime_su2/3_*.png` | Q vs MC time — all β share the same y-axis scale |
| `histograms_su2/3_*.png` | Q distributions — all β share the same x-axis scale |
| `susceptibility_su2/3.png` | χ_t^(1/4) vs lattice spacing |
| `tau_int_su2/3.png` | τ_int(Q²) vs lattice spacing — shows freezing |
| `thermalization_comparison_*.png` | Plaquette thermalization curves |

## Run Directory Structure

```
data/results/T{T}_L{L}_b{beta}_{boundary}_seed{seed}/
├── input.txt
├── output/
│   ├── plaquette[_su3].dat          # plaquette vs sweep
│   ├── therm_topcharge[_su3].dat    # raw Q vs sweep (thermalization)
│   ├── topcharge[_su3].dat          # smeared Q (production)
│   └── configs[_su3]/               # saved gauge configurations
└── run_info.txt
```

## Rust Analysis Extension (Optional)

`rust_analysis/` accelerates analysis hot-paths with rayon parallelism. Falls back to pure Python automatically if not installed.

| Function | Description |
|---|---|
| `find_optimal_alpha(q)` | Parallel golden-section search for α |
| `jackknife_error(data)` | O(N) delete-one jackknife |
| `bootstrap_alpha(q, n_boot)` | Parallel bootstrap → `(α, σ_α)` |
| `bootstrap_Q2(q, n_boot)` | Parallel bootstrap → `(⟨Q²⟩, σ)` |

```bash
cd rust_analysis && maturin develop --release
```

## Project Structure

```
├── src/              # C++ simulation sources
├── include/          # Headers (topcharge, heatbath, smearing, progressbar)
├── _Utility/         # Core SU(2) library (fields, RNG, smearing, I/O)
├── analysis/         # Python analysis (calculations, plotting, autocorrelation)
├── rust_analysis/    # Optional Rust extension
├── scripts/          # Workflow shell scripts
├── input/            # Parameter files
├── tests/            # Catch2 unit tests
└── run_simulation.sh # Main entry point
```
