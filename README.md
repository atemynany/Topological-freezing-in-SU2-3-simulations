# SU(2)/SU(3) Topological Charge Freezing Study

Monte Carlo simulation of topological charge in SU(2) and SU(3) lattice gauge theory to study topological freezing at fine lattice spacings.

Works on **Linux** and **macOS** (including Apple Silicon).

## Requirements

- CMake ≥ 3.16
- C++17 compiler (GCC, Clang, or Apple Clang)
- Python 3: numpy, scipy, matplotlib, pyerrors
- Optional: Rust toolchain + `maturin` (only if you want to use the `rust_analysis/` accelerator)

### macOS extras

```bash
# install Xcode command line tools (gives you clang + make)
xcode-select --install

# install cmake (via Homebrew)
brew install cmake

# optional: OpenMP support for parallel loops
brew install libomp
```

On macOS without `libomp`, everything still compiles and runs — just without OpenMP parallelization. Both SU(2) and SU(3) heatbath use checkerboard (even/odd) decomposition for thread-safe OpenMP parallel sweeps; set `OMP_NUM_THREADS` to control parallelism. For small lattices (below ~32^4), single-threaded is faster due to thread overhead.

## Build

```bash
./scripts/build.sh release
./scripts/build.sh release --march=x86-64   # optional: override -march on non-Apple platforms
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

Seven test suites: `test_heatbath`, `test_linear_algebra`, `test_topcharge`, `test_topcharge_su3`, `test_smearing`, `test_plaquette`, `test_gauge_invariance`.

Parity check for the fused SU(2) pipeline (outside CTest):

```bash
./scripts/test_parity_heatbath_topcharge.sh      # legacy vs fused, OMP=1
```

Quick smoke-test input files for the fused binary:

```bash
./build/bin/heatbath_topcharge_su2 -i input/tests_input/setup_test_16x16_periodic_su2.txt
./build/bin/heatbath_topcharge_su2 -i input/tests_input/setup_test_16x16_open_su2.txt
```

500-sweep runs, β=2.5, 16⁴. Outputs land in `data/results/test_16x16_{periodic,open}_b2.5/`.

## Main Workflow

### Two-phase pipeline (default, SU(2) and SU(3))

Heatbath (config generation) then topological charge measurement.

```bash
./run_heatbath_scan.sh --su2                # Phase 1: generate configs
./run_heatbath_scan.sh --su3                # SU(3) configs
./run_heatbath_scan.sh --su2 --dry-run      # preview without running
./run_heatbath_scan.sh --su2 --jobs 4       # parallel (4 jobs)

./run_topcharge_scan.sh --su2               # Phase 2: measure Q
./run_topcharge_scan.sh --su3 --jobs 4      # parallel
./run_topcharge_scan.sh --su2 --exclude "T65\|T81"  # skip specific ensembles
```

`run_heatbath_scan.sh` detects thermalization from the plaquette history and writes
the resulting `start_conf` to both `input.txt` and `run_info.txt`. The topcharge
scan preserves that per-run value and only falls back to
`input/heatbath_input/topcharge_params_su*.txt` for legacy runs without
`start_conf`. Set `HEATBATH_RESULTS_DIR=/path/to/results` to redirect both scan
phases to a non-default results directory.

### Fused pipeline (SU(2), in-memory, no gauge configs written)

For SU(2) scans where you don't need persisted gauge configs, use the fused binary `heatbath_topcharge_su2`. It runs heatbath and measures Q on an APE-smeared scratch copy at every `save_interval`, writing `plaquette.dat`, `action_density.dat`, `topcharge.dat`, and `topcharge_timeslice.dat`.

```bash
./run_heatbath_topcharge_scan.sh            # discovers input/heatbath_topcharge_input/setup_*_su2.txt
./run_heatbath_topcharge_scan.sh --jobs 4
./run_heatbath_topcharge_scan.sh --dry-run
./run_heatbath_topcharge_scan.sh --exclude "T65\|T81"
```

Parity with the two-phase pipeline is verified at single-thread determinism (`scripts/test_parity_heatbath_topcharge.sh`): plaquette bit-identical, max |ΔQ| = 0. Multi-thread heatbath trajectories are nondeterministic by design.

Setup files live in three subdirs of `input/`:

- `input/heatbath_input/` — active setups for the two-phase pipeline (`setup_*_su2.txt` or `setup_*_su3.txt`) plus `topcharge_params_su{2,3}.txt` (with `_finished.txt` fallback for parameter files); setup files ending in `_finished.txt` are archived and are not picked up by the scan globs
- `input/heatbath_topcharge_input/` — active setups for the fused SU(2) pipeline (`setup_*_su2.txt`); `_finished.txt` files are parked
- `input/tests_input/` — small smoke-test setups (e.g. `setup_test_16x16_periodic_su2.txt`)

Each setup file contains one lattice configuration; multiple `beta` lines produce multiple runs. `--su2`/`--su3` discovers all matching setup files automatically in the relevant subdir.

Results saved to `data/results/` and plots to `output/figures_analysis/`.

## Cluster (Fuchs / HHLR-Goethe)

Two cluster targets are supported, each with its own `sbatch` wrapper. Build once on the login node before submitting:

```bash
rm -rf build && ./scripts/build.sh release
```

### Fuchs HPC

```bash
sbatch cluster/fuchs_heatbath.sh             # two-phase: heatbath scan
sbatch cluster/fuchs_topcharge.sh            # two-phase: topcharge measurement
sbatch cluster/fuchs_heatbath_topcharge.sh   # fused SU(2) (no disk configs)
sbatch cluster/fuchs_heatbath_thread_test.sh # OpenMP thread-scaling test
sbatch cluster/fuchs_continue.sh             # resume specific runs from existing input.txt
```

### HHLR-Goethe

The HHLR scripts run from `/home` (binary) but write outputs to `/work` via `HEATBATH_RESULTS_DIR`:

```bash
sbatch cluster/hhlr_heatbath.sh             # two-phase: heatbath scan (40 cores, 11 days)
sbatch cluster/hhlr_topcharge.sh            # two-phase: topcharge measurement
sbatch cluster/hhlr_heatbath_topcharge.sh   # fused SU(2)
sbatch cluster/hhlr_continue.sh             # resume specific runs (edit the input.txt list inside)
```

The `*_continue.sh` scripts re-launch `mc_heatbath` against existing `data/results/.../input.txt` files — useful when a job hits the wall-time and you want to extend a run without re-thermalising. Edit the script to point at the runs you want to resume.

### Monitoring

```bash
tail -1 data/results/T*_su2/output/plaquette.dat   # check sweep number
squeue -u barros                                   # job status
sacct                                              # job history
```

### Syncing results to local

```bash
./cluster/sync_results.sh
```

## Key Executables

| Binary | Purpose |
|--------|---------|
| `mc_heatbath` | SU(2) heatbath config generation (OpenMP checkerboard) |
| `mc_heatbath_su3` | SU(3) heatbath config generation (OpenMP checkerboard + overrelaxation) |
| `meas_topcharge` | SU(2) topological charge measurement |
| `meas_topcharge_su3` | SU(3) topological charge measurement |
| `meas_action_density_su2` | Recompute SU(2) Wilson action density from saved configs, optionally on a temporal bulk window |
| `heatbath_topcharge_su2` | Fused SU(2) heatbath + in-memory topcharge (no gauge configs written) |
| `compute_instanton_Q` | Instanton topological charge computation |

## Input Files

- `input/heatbath_input/setup_*_su2.txt` / `setup_*_su3.txt` — active setups for the two-phase pipeline (lattice params + beta values)
- `input/heatbath_input/topcharge_params_su2.txt` / `_su3.txt` — topcharge measurement params (smearing, boundary exclusion), with `_finished.txt` fallback; also read by the fused SU(2) pipeline
- `input/heatbath_topcharge_input/setup_*_su2.txt` — active setups for the fused SU(2) pipeline
- `input/tests_input/setup_test_*_su2.txt` — smoke-test setups for the fused binary

The topcharge scan script auto-detects periodic vs open boundary conditions from
the run directory name (`_periodic_` or `_open_`). It preserves an explicit
`exclude_boundary_slices` already present in a run input; otherwise open BC runs
use `exclude_boundary_slices_open` from the param file and periodic runs use 0.
A periodic run with `exclude_boundary_slices > 0` is a temporal-subvolume
measurement, so Q is not integer-quantized.

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

`data/results/T{T}_L{L}_b{beta}_{boundary}_seed{seed}_{su2|su3}/`:
- `output/plaquette[_su3].dat` — Plaquette and Wilson action density history (every sweep)
- `output/action_density[_su3].dat` — Wilson action density history only
- `output/topcharge.dat` — Q measurements (APE-smeared)
- `output/topcharge_su3.dat` — SU(3) Q measurements (two-phase SU(3) runs)
- `output/topcharge_timeslice.dat` — Per-timeslice q(t)
- `configs/` — Gauge configurations (two-phase scan only; absent in fused SU(2) runs)
- `input.txt`, `run_info.txt`, `*.log` — Input, metadata, stdout

For the heatbath plaquette history, the columns are:

```text
sweep  average_plaquette  action_density
```

The action density is computed from the normalized average plaquette using the
4D Wilson gauge-action relation:

```text
s = 6 * beta * (1 - <P>)
```

The shared C++ helper is `avg_action_density(beta, avg_plaquette)` in
`include/action_density.hh`.

## Analysis

Run the standard analysis:

```bash
python3 analysis/analysis.py --su2 --results-dir data/results
python3 analysis/analysis.py --su3 --results-dir data/results
```

Without `--results-dir`, the analysis looks for synced cluster roots
`data/results_home` and `data/results_work` and analyzes them separately into
`output/figures_analysis/results_home/` and
`output/figures_analysis/results_work/`. Pass `--results-dir` repeatedly to
intentionally combine several roots, and use `--output-dir` to redirect plots
and summaries.

Or generate one figure set per detected APE smearing level:

```bash
python3 analysis/analysis_per_smear.py --su2 --results-dir data/results
python3 analysis/analysis_per_smear.py --su3 --results-dir data/results --smear 20
```

Key outputs:
- **χ_t^(1/4) vs β** — Susceptibility (should match ~200 MeV for SU(2))
- **τ_int(Q²) vs β** — Autocorrelation time (increases with β → freezing)
- **Q histograms** — Distribution narrowing indicates freezing
- **q(t) timeslice density** — Temporal topcharge density, susceptibility contribution, and τ_int(t) grids via `timeslice_analysis.py`
- **Action density** — `action_density_{su2,su3}.png` plus `action_density_summary_{su2,su3}.txt` when action-density files are present

`analysis_per_smear.py` auto-detects every APE smearing level present in `topcharge.dat` or `topcharge_su3.dat` and writes each set of figures into `smear<N>/` below each dataset output folder, for example `output/figures_analysis/results_home/smear15/`. Useful when comparing how the susceptibility plateau and τ_int depend on smearing depth. Its timeslice binning scales with the coarsest lattice spacing so subvolumes remain approximately fixed in physical units.

Standalone helpers (do not need the full pipeline):
- `analysis/plot_topcharge.py` — Q vs smearing steps / MC time from a raw `topcharge.dat`
- `analysis/analyze_topcharge.py` — jackknife mean / error / bias on Q
- `analysis/generate_project_overview.py` — creates `output/project_overview.html` from synced results, figures, and run metadata
- `scale_set.ipynb` — SU(2) and SU(3) lattice-spacing scale-setting formulas (arXiv:1811.02800 and hep-lat/0108008)

### Optional: Rust analysis extension

`rust_analysis/` is a PyO3 module (`qcd_analysis`) that provides parallelised, hot-path replacements for the slow Python routines: `find_optimal_alpha`, `jackknife_error`, `bootstrap_alpha`, `bootstrap_q2`. Build with maturin:

```bash
cd rust_analysis
maturin develop --release         # installs qcd_analysis into the active Python env
```

The Python analysis falls back to its pure-Python implementations if the module is not installed.

### Supplementary docs

- `context.md` — physics goal, observables, and pipeline summary
- `additional_infos.md` — bug-scan notes and design decisions worth knowing before editing the analysis
- `CLAUDE.md` — concise operational notes for coding agents
- `_Utility/README.md` — utility-library contents and how it is built through the root CMake project
- `su3_static_potential_check/README.md` — minimal cl2qcd open-boundary patch notes

The analysis uses the Gamma method via `pyerrors` for autocorrelation-aware
errors. Only full-lattice periodic Q is integer-rounded to `Q_re`; open-boundary
or temporal-subvolume measurements use continuous `alpha * Q` because those
charges are not globally integer-quantized.

## Git and SSH Remote

The repository is configured for SSH pushes to GitHub:

```bash
git remote -v
# origin  git@github.com:atemynany/Topological-freezing-in-SU2-3-simulations.git
```

If VS Code prints a line like `refs/heads/master refs/remotes/master`, that is
only a ref query. The tracked upstream should be `master -> origin/master`; check
or repair it with:

```bash
git branch -vv
git branch --set-upstream-to=origin/master master
git remote set-url origin git@github.com:atemynany/Topological-freezing-in-SU2-3-simulations.git
ssh -T git@github.com
```

SSH keys themselves live outside the repository in `~/.ssh/` and must be added to
the GitHub account. This repo should not store private keys.

## Project Structure

```
├── src/                  # C/C++ simulation sources
├── include/              # Header files
├── _Utility/             # Core library (fields, RNG, smearing)
├── tests/                # Catch2 unit tests
├── scripts/              # Build, run, and analysis scripts
├── analysis/             # Python analysis and plotting
├── rust_analysis/        # Optional Rust+PyO3 accelerator (qcd_analysis)
├── cluster/              # SLURM scripts for Fuchs and HHLR-Goethe
├── input/                # Setup files and topcharge params
├── scale_set.ipynb                # Lattice-spacing scale-setting notebook
├── context.md                     # Physics + pipeline summary
├── additional_infos.md            # Bug-scan notes and design decisions
├── run_heatbath_scan.sh           # Phase 1: config generation
├── run_topcharge_scan.sh          # Phase 2: topcharge measurement
└── run_heatbath_topcharge_scan.sh # SU(2) fused pipeline (no disk configs)
```
