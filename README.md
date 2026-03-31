# SU(2)/SU(3) Topological Charge Freezing Study

Monte Carlo simulation of topological charge in SU(2) and SU(3) lattice gauge theory to study topological freezing at fine lattice spacings.

Works on **Linux** and **macOS** (including Apple Silicon).

## Requirements

- CMake ≥ 3.16
- C++17 compiler (GCC, Clang, or Apple Clang)
- Python 3: numpy, matplotlib, pyerrors

### macOS extras

```bash
# install Xcode command line tools (gives you clang + make)
xcode-select --install

# install cmake (via Homebrew)
brew install cmake

# optional: OpenMP support for parallel loops
brew install libomp
```

On macOS without `libomp`, everything still compiles and runs — just without OpenMP parallelization. The SU(3) heatbath uses checkerboard decomposition for thread-safe parallel sweeps; set `OMP_NUM_THREADS` to control parallelism.

## Build

```bash
./scripts/build.sh release
```

Or manually:

```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --parallel 4
```

CMake auto-detects macOS and sets up the right flags. If you installed `libomp` via Homebrew, OpenMP is picked up automatically.

## Run Tests

```bash
./scripts/run_tests.sh
```

## Main Workflow: Beta Scan

Run simulations across multiple β values and analyze:

```bash
# 1. Edit beta values to scan
nano input/beta_scan_su2.txt

# 2. Run scan (generates configs + measures Q)
./scripts/run_beta_scan.sh --su2

# 3. Analyze results
python3 analysis/analysis.py --su2
```

Or use the all-in-one entry point:

```bash
./run_simulation.sh --su2          # SU(2) with default betas
./run_simulation.sh --su3          # SU(3) with default betas
./run_simulation.sh --su2 --open   # SU(2) with open boundaries
./run_simulation.sh --dry-run      # preview without running
```

Results saved to `data/results/` and plots to `output/figures_analysis/`.

## Key Executables

| Binary | Purpose |
|--------|---------|
| `mc_heatbath` | SU(2) heatbath config generation (+ smeared therm Q) |
| `mc_heatbath_su3` | SU(3) heatbath config generation (OpenMP checkerboard, smeared therm Q) |
| `meas_topcharge` | SU(2) topological charge measurement |
| `meas_topcharge_su3` | SU(3) topological charge measurement |
| `compute_instanton_Q` | Instanton topological charge computation |

## Input Files

- `input/base_params_su2.txt` — Simulation parameters (T, L, sweeps, smearing)
- `input/base_params_su3.txt` — SU(3) simulation parameters
- `input/beta_scan_su2.txt` — List of β values to scan
- `input/beta_scan_su3.txt` — SU(3) β values

## Output

- `data/results/T{T}_L{L}_b{beta}_{boundary}_seed{seed}/` — Per-run results
  - `configs/` — Gauge configurations
  - `output/topcharge.dat` — Q measurements (APE-smeared)
  - `output/plaquette[_su3].dat` — Plaquette history (every sweep)
  - `output/therm_topcharge.dat` — Q per sweep for thermalization monitoring
  - `run_info.txt` — Run metadata

## Analysis

Key outputs from `analysis.py`:
- **χ_t^(1/4) vs β** — Susceptibility (should match ~200 MeV for SU(2))
- **τ_int(Q²) vs β** — Autocorrelation time (increases with β → freezing)
- **Q histograms** — Distribution narrowing indicates freezing

## Project Structure

```
├── src/                  # C/C++ simulation sources
├── include/              # Header files
├── _Utility/             # Core library (fields, RNG, smearing)
├── tests/                # Catch2 unit tests
├── scripts/              # Build, run, and analysis scripts
├── analysis/             # Python analysis and plotting
├── input/                # Parameter and beta scan files
└── run_simulation.sh     # Main entry point
```
