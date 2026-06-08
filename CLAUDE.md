# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
./scripts/build.sh release          # Release build (O3, fast-math)
./scripts/build.sh release --march=x86-64  # Optional architecture override
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

Auto-builds if `build/` is missing. Seven test suites: `test_heatbath`, `test_linear_algebra`, `test_topcharge`, `test_topcharge_su3`, `test_smearing`, `test_plaquette`, `test_gauge_invariance`.

Parity check for the fused SU(2) pipeline (not part of CTest):

```bash
./scripts/test_parity_heatbath_topcharge.sh   # legacy vs fused, OMP=1 forced
```

Forces `OMP_NUM_THREADS=1` because multi-thread heatbath is nondeterministic (checkerboard scheduling affects RNG draw order). Expected: plaquette bit-identical, max |ΔQ|=0, max |Δq(t)|=0.

## Simulation Workflow

Two-phase workflow: heatbath (config generation) then topological charge measurement.

```bash
./run_heatbath_scan.sh --su2               # Phase 1: generate configs
./run_heatbath_scan.sh --su2 --jobs 4      # parallel (4 beta values at once)
./run_heatbath_scan.sh --su2 --dry-run     # preview without running

./run_topcharge_scan.sh --su2              # Phase 2: measure Q
./run_topcharge_scan.sh --su2 --jobs 4     # parallel
```

Input files live in three subdirs of `input/`:

- `input/heatbath_input/` — active two-phase setup files must match `setup_*_su2.txt` or `setup_*_su3.txt`; `_finished.txt` files are archived/parked and are not discovered by the scan glob. Topcharge params live here as `topcharge_params_su{2,3}.txt`, with `_finished.txt` fallback.
- `input/heatbath_topcharge_input/` — active fused SU(2) setup files must match `setup_*_su2.txt`; `_finished.txt` files are parked
- `input/tests_input/` — smoke-test setups for the fused binary

The topcharge scan auto-detects periodic vs open BC from the run directory name,
but preserves an explicit `exclude_boundary_slices` in `input.txt`. A periodic
run with `exclude_boundary_slices > 0` is a temporal-subvolume measurement and
must be analyzed with continuous αQ, not integer-rounded Q.

### Fused SU(2) pipeline (in-memory, no configs on disk)

`heatbath_topcharge_su2` runs heatbath and measures Q on an APE-smeared scratch copy at every `save_interval`. It writes `plaquette.dat`, `action_density.dat`, `topcharge.dat`, and `topcharge_timeslice.dat` — the gauge field never hits disk. Merges `setup_*_su2.txt` with `topcharge_params_su2.txt` automatically.

```bash
./run_heatbath_topcharge_scan.sh                 # all setup_*_su2.txt
./run_heatbath_topcharge_scan.sh --jobs 4
./run_heatbath_topcharge_scan.sh --exclude "T65\|T81"
./run_heatbath_topcharge_scan.sh --dry-run
```

Cluster: `cluster/fuchs_heatbath_topcharge.sh`, `cluster/hhlr_heatbath_topcharge.sh`. Both run single multi-threaded jobs (not parallel single-core) because smearing + topcharge are OMP-parallel.

Source: `src/heatbath_topcharge_su2.cc`. SU(3) fused variant not implemented.

## Cluster

Cluster: `fuchs.hhlr-gu.de`, user `barros`, project dir `/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations/`

Fuchs wrappers: `cluster/fuchs_heatbath.sh`, `cluster/fuchs_topcharge.sh`,
`cluster/fuchs_heatbath_topcharge.sh`, `cluster/fuchs_continue.sh`,
`cluster/fuchs_action_density_T160_L80.sh`, and thread-test helpers.

HHLR-Goethe wrappers: `cluster/hhlr_heatbath.sh`, `cluster/hhlr_topcharge.sh`,
`cluster/hhlr_heatbath_topcharge.sh`, `cluster/hhlr_continue.sh`, plus
specialized loop/qtarget scripts. HHLR jobs run binaries from `/home` and write
large outputs to `/work` through `HEATBATH_RESULTS_DIR`.

```bash
sbatch cluster/fuchs_heatbath.sh
sbatch cluster/fuchs_topcharge.sh
sbatch cluster/hhlr_heatbath.sh
sbatch cluster/hhlr_topcharge.sh
tail -1 data/results/T*_su2/output/plaquette.dat  # Check sweep progress
squeue -u barros                                   # Job status
```

**Before first submit:** build on the cluster login node:
```bash
rm -rf build && ./scripts/build.sh release
```

**Sync results to local** (run from local machine, skips binary configs):
```bash
./cluster/sync_results.sh
```

## Analysis

```bash
python3 analysis/analysis.py --su2 --results-dir data/results
python3 analysis/analysis.py --su3 --results-dir data/results
python3 analysis/analysis_per_smear.py --su2 --results-dir data/results
```

Requires Python 3 with `numpy`, `matplotlib`, `pyerrors` (Conda env `master_thesis`). Without `--results-dir`, analysis defaults to synced cluster roots `data/results_home` and `data/results_work` and writes separate output folders below `output/figures_analysis/`. Output root can be changed with `--output-dir`.

## Architecture

### Gauge group duality

SU(2) and SU(3) are implemented as parallel but separate codepaths. Each has its own:
- Heatbath MC: `src/MC_heatbath.cc` (Kennedy-Pendleton + OpenMP checkerboard) / `src/MC_heatbath_su3.cc` (Cabibbo-Marinari + overrelaxation + OpenMP checkerboard)
- Topcharge measurement: `src/meas_topcharge_su2.cc` / `src/meas_topcharge_su3.cc`
- Linear algebra: `_Utility/include/linear_algebra.hh` (SU(2), 8 doubles/link) / `include/su3_linear_algebra.hh` (SU(3), 18 doubles/link)
- Topcharge computation: `include/topcharge_su2.hh` / `include/topcharge_su3.hh`
- Optional SU(2) action-density recomputation: `src/meas_action_density_su2.cc`

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

Both SU(2) and SU(3) heatbath use checkerboard (even/odd) site decomposition for thread-safe OpenMP parallel sweeps. Per-thread RNG seeding via `rlxd_init(2, seed + tid * 997)`. The SU(3) heatbath also parallelizes overrelaxation sweeps. For small lattices (below ~32^4), `OMP_NUM_THREADS=1` is faster due to thread synchronization overhead.

### Data layout

- SU(2): link at `(site, mu)` → 8 contiguous doubles at offset `(4*site + mu) * 8`
- SU(3): link at `(site, mu)` → 18 contiguous doubles at offset `(4*site + mu) * 18`
- All linear algebra and topcharge functions are inline for performance

### Boundary conditions

Open BC (Lüscher-Schaefer): temporal plaquettes at `t=0` and `t=T-1` weighted by 0.5 in staple sum. Topcharge measurement excludes configurable boundary slices.

### Data flow

Two-phase: `mc_heatbath[_su3]` → `configs/conf.*` + `plaquette[_su3].dat` + `action_density[_su3].dat` → `detect_thermalization.py` → `meas_topcharge[_su3]` → `topcharge[_su3].dat` + `topcharge_timeslice.dat` → `analysis.py`

Fused SU(2): `heatbath_topcharge_su2` → `plaquette.dat` + `action_density.dat` + `topcharge.dat` + `topcharge_timeslice.dat` → `analysis.py` (no `configs/`)

Results per run: `data/results/T{T}_L{L}_b{beta}_{boundary}_seed{seed}/`

### Analysis pipeline (`analysis/`)

- `analysis.py` — Entry point; scans explicit `--results-dir` roots or synced defaults
- `result_paths.py` — Resolves explicit results roots or synced defaults (`data/results_home`, `data/results_work`) and supports qtarget aggregate layouts
- `calculations.py` — Core: `RunData`/`AnalysisResult` containers, `find_optimal_alpha()` for Q rescaling, `analyze_run()` for χ_t and τ_int
- `plotting_code.py` — Publication figures (Q vs MC time, histograms, susceptibility, τ_int)
- `timeslice_analysis.py` — Temporal topcharge density, susceptibility-contribution, and tau_int(t) analysis
- `action_density_analysis.py` — Action-density plots and summaries
- `analysis_per_smear.py` — Full plotting pass per detected APE smearing level
- `autocorrelation.py` — Wraps `pyerrors` Gamma method for integrated autocorrelation time

Key analysis formula: full-lattice periodic runs use `Q_rescaled = round(α * Q_raw)` where α is found by golden-section search, while open-boundary or temporal-subvolume runs use continuous `α * Q_raw`; then `χ_t = ⟨Q²⟩/V`.
