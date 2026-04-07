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

## Main Workflow

The simulation runs in two phases: heatbath (config generation) then topological charge measurement.

```bash
./run_heatbath_scan.sh --su2                # Phase 1: generate configs
./run_heatbath_scan.sh --su3                # SU(3) configs
./run_heatbath_scan.sh --su2 --dry-run      # preview without running
./run_heatbath_scan.sh --su2 --jobs 4       # parallel (4 jobs)

./run_topcharge_scan.sh --su2               # Phase 2: measure Q
./run_topcharge_scan.sh --su3 --jobs 4      # parallel
```

Setup files in `input/setup_*_su2.txt` / `input/setup_*_su3.txt` define lattice parameters and beta values. Each file contains one lattice configuration; multiple `beta` lines produce multiple runs. `--su2`/`--su3` discovers all matching setup files automatically.

Results saved to `data/results/` and plots to `output/figures_analysis/`.

## Cluster (Fuchs HPC)

```bash
# On the cluster login node — build once before submitting
rm -rf build && ./scripts/build.sh release

# Submit jobs
sbatch cluster/fuchs_heatbath.sh     # heatbath scan
sbatch cluster/fuchs_topcharge.sh    # topcharge measurement

# Monitor
tail -f logs/scan_<jobid>.out
```

Sync results to local (run from local machine):

```bash
./cluster/sync_results.sh
```

## Key Executables

| Binary | Purpose |
|--------|---------|
| `mc_heatbath` | SU(2) heatbath config generation (+ smeared therm Q) |
| `mc_heatbath_su3` | SU(3) heatbath config generation (OpenMP checkerboard, smeared therm Q) |
| `meas_topcharge` | SU(2) topological charge measurement |
| `meas_topcharge_su3` | SU(3) topological charge measurement |
| `compute_instanton_Q` | Instanton topological charge computation |

## Input Files

- `input/setup_*_su2.txt` — SU(2) setup files (lattice params + beta values)
- `input/setup_*_su3.txt` — SU(3) setup files
- `input/topcharge_params_su2.txt` — Topcharge measurement params (smearing, exclusion)
- `input/topcharge_params_su3.txt` — SU(3) topcharge params

Example setup file:
```
T                   8
L                   8
start_type          cold
boundary            periodic
num_sweeps          1000
save_interval       10

beta 2.3
beta 2.5
beta 2.7
```

## Output

- `data/results/T{T}_L{L}_b{beta}_{boundary}_seed{seed}_{su2|su3}/` — Per-run results
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
├── input/                # Setup files and topcharge params
├── run_heatbath_scan.sh  # Phase 1: config generation
└── run_topcharge_scan.sh # Phase 2: topcharge measurement
```
