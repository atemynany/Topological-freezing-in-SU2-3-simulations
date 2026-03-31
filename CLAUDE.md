# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
./scripts/build.sh release          # Release build (O3, fast-math)
./scripts/build.sh debug            # Debug build
```

Requires CMake ≥ 3.16 and a C++17 compiler. On macOS, OpenMP via `brew install libomp` is optional — builds fine without it.

Manual build: `mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && cmake --build . --parallel`

Executables go to `build/bin/`, library to `build/lib/`.

## Tests

```bash
./scripts/run_tests.sh              # All tests (Catch2 via CTest)
./scripts/run_tests.sh test_topcharge  # Single test suite by name (-R pattern)
```

Auto-builds if `build/` is missing. Five test suites: `test_heatbath`, `test_linear_algebra`, `test_topcharge`, `test_smearing`, `test_plaquette`.

## Simulation Workflow

```bash
./run_simulation.sh --su2              # SU(2) periodic BC
./run_simulation.sh --su3              # SU(3) periodic BC
./run_simulation.sh --su2 --open       # Open boundary conditions
./run_simulation.sh --dry-run          # Preview without running
```

This delegates to `scripts/run_beta_scan.sh`, which for each beta: runs heatbath → detects thermalization → measures topological charge with APE smearing.

Parameters: `input/base_params_su2.txt` / `input/base_params_su3.txt`. Beta values: `input/beta_scan_su2.txt` / `input/beta_scan_su3.txt`. Copy from `input/example_*.txt` to create these.

## Analysis

```bash
python3 analysis/analysis.py --su2     # or --su3
```

Requires Python 3 with `numpy`, `matplotlib`, `pyerrors` (Conda env `master_thesis`). Output goes to `output/figures_analysis/`.

## Architecture

### Gauge group duality

SU(2) and SU(3) are implemented as parallel but separate codepaths. Each has its own:
- Heatbath MC: `src/MC_heatbath.cc` (Kennedy-Pendleton) / `src/MC_heatbath_su3.cc` (Cabibbo-Marinari + overrelaxation + OpenMP checkerboard)
- Topcharge measurement: `src/meas_topcharge_su2.cc` / `src/meas_topcharge_su3.cc`
- Linear algebra: `_Utility/include/linear_algebra.hh` (SU(2), 8 doubles/link) / `include/su3_linear_algebra.hh` (SU(3), 18 doubles/link)
- Topcharge computation: `include/topcharge_su2.hh` / `include/topcharge_su3.hh`

### Core library (`_Utility/`)

Shared infrastructure for both gauge groups:
- `fields.hh/cc` — Gauge field allocation (32-byte aligned via `posix_memalign`)
- `geometry.hh` — Lattice site/link indexing; `ggi(site, mu)` for link offsets
- `ranlux.hh/cc`, `ranlxd.c`, `ranlxs.c` — RANLUX RNG
- `smearing_techniques.hh/cc` — APE smearing (+ HYP, FAT variants in `smearing_techniques_all.cc`)
- `io.hh/cc` — Binary config I/O and parameter parsing

### Global state pattern

Both SU(2) and SU(3) simulations use `extern` globals for lattice dimensions (`T`, `L`), the gauge field pointer, boundary condition flags, and precomputed neighbor tables. These are set once at startup from input files.

### OpenMP parallelization

The SU(3) heatbath (`MC_heatbath_su3.cc`) uses checkerboard (even/odd) site decomposition for thread-safe OpenMP parallel sweeps. Both heatbath and overrelaxation loops are parallelized. Per-thread RNG seeding via `rlxd_init(2, seed + tid * 997)`. The SU(2) heatbath does not yet use checkerboard parallelization.

### Thermalization topcharge

Both heatbath programs output `therm_topcharge.dat` (Q measured every sweep). SU(2) measures unsmeared Q; SU(3) applies APE smearing (`smear_steps`, `smear_alpha` params) before measuring Q for better signal.

### Data layout

- SU(2): link at `(site, mu)` → 8 contiguous doubles at offset `(4*site + mu) * 8`
- SU(3): link at `(site, mu)` → 18 contiguous doubles at offset `(4*site + mu) * 18`
- All linear algebra and topcharge functions are inline for performance

### Boundary conditions

Open BC (Lüscher-Schaefer): temporal plaquettes at `t=0` and `t=T-1` weighted by 0.5 in staple sum. Topcharge measurement excludes configurable boundary slices.

### Data flow

`mc_heatbath[_su3]` → `configs/conf_*.bin` + `plaquette.dat` + `therm_topcharge.dat` → `detect_thermalization.py` → `meas_topcharge[_su3]` → `topcharge[_su3].dat` → `analysis.py`

Results per run: `data/results/T{T}_L{L}_b{beta}_{boundary}_seed{seed}/`

### Analysis pipeline (`analysis/`)

- `analysis.py` — Entry point, scans `data/results/` for runs
- `calculations.py` — Core: `RunData`/`AnalysisResult` containers, `find_optimal_alpha()` for Q rescaling, `analyze_run()` for χ_t and τ_int
- `plotting_code.py` — Publication figures (Q vs MC time, histograms, susceptibility, τ_int)
- `timeslice_analysis.py` — Spatial topcharge density analysis
- `autocorrelation.py` — Wraps `pyerrors` Gamma method for integrated autocorrelation time

Key analysis formula: `Q_rescaled = round(α * Q_raw)` where α ≈ 0.9–1.1 found via golden-section search, then `χ_t = ⟨Q²⟩/V`.
