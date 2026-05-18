# SU(2)/SU(3) Topological Charge Freezing Study

Monte Carlo simulation of topological charge in SU(2) and SU(3) lattice gauge theory to study topological freezing at fine lattice spacings.

Works on **Linux** and **macOS** (including Apple Silicon).

## Requirements

- CMake ≥ 3.16
- C++17 compiler (GCC, Clang, or Apple Clang)
- Python 3: numpy, scipy, matplotlib, pyerrors

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

For SU(2) scans where you don't need persisted gauge configs, use the fused binary `heatbath_topcharge_su2`. It runs heatbath and measures Q on an APE-smeared scratch copy at every `save_interval`, writing only `plaquette.dat`, `topcharge.dat`, and `topcharge_timeslice.dat`.

```bash
./run_heatbath_topcharge_scan.sh            # discovers input/heatbath_topcharge_input/setup_*_su2.txt
./run_heatbath_topcharge_scan.sh --jobs 4
./run_heatbath_topcharge_scan.sh --dry-run
./run_heatbath_topcharge_scan.sh --exclude "T65\|T81"
```

Parity with the two-phase pipeline is verified at single-thread determinism (`scripts/test_parity_heatbath_topcharge.sh`): plaquette bit-identical, max |ΔQ| = 0. Multi-thread heatbath trajectories are nondeterministic by design.

Setup files live in three subdirs of `input/`:

- `input/heatbath_input/` — production setups for the two-phase pipeline (`setup_*_su{2,3}_finished.txt`) plus `topcharge_params_su{2,3}.txt`
- `input/heatbath_topcharge_input/` — setups for the fused SU(2) pipeline (drop `setup_*_su2.txt` files here)
- `input/tests_input/` — small smoke-test setups (e.g. `setup_test_16x16_periodic_su2.txt`)

Each setup file contains one lattice configuration; multiple `beta` lines produce multiple runs. `--su2`/`--su3` discovers all matching setup files automatically in the relevant subdir.

Results saved to `data/results/` and plots to `output/figures_analysis/`.

## Cluster (Fuchs HPC)

```bash
# On the cluster login node — build once before submitting
rm -rf build && ./scripts/build.sh release

# Submit jobs
sbatch cluster/fuchs_heatbath.sh            # two-phase: heatbath scan
sbatch cluster/fuchs_topcharge.sh           # two-phase: topcharge measurement
sbatch cluster/fuchs_heatbath_topcharge.sh  # fused SU(2) (no disk configs)

# Monitor progress
tail -1 data/results/T*_su2/output/plaquette.dat   # check sweep number
squeue -u barros                                     # job status
sacct                                                # job history
```

Sync results to local (run from local machine):

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
| `heatbath_topcharge_su2` | Fused SU(2) heatbath + in-memory topcharge (no gauge configs written) |
| `compute_instanton_Q` | Instanton topological charge computation |

## Input Files

- `input/heatbath_input/setup_*_su2.txt` / `setup_*_su3.txt` — setups for the two-phase pipeline (lattice params + beta values)
- `input/heatbath_input/topcharge_params_su2.txt` / `_su3.txt` — topcharge measurement params (smearing, boundary exclusion); also read by the fused SU(2) pipeline
- `input/heatbath_topcharge_input/setup_*_su2.txt` — setups for the fused SU(2) pipeline
- `input/tests_input/setup_test_*_su2.txt` — smoke-test setups for the fused binary

The topcharge scan script auto-detects periodic vs open boundary conditions from the run directory name (`_periodic_` or `_open_`). For open BC runs, `exclude_boundary_slices_open` from the param file is applied; periodic runs use 0.

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
- `output/topcharge.dat` — Q measurements (APE-smeared)
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
python3 analysis/analysis.py --su2
python3 analysis/analysis.py --su3
```

Or generate one figure set per detected APE smearing level:

```bash
python3 analysis/analysis_per_smear.py --su2
python3 analysis/analysis_per_smear.py --su3 --smear 20
```

Key outputs:
- **χ_t^(1/4) vs β** — Susceptibility (should match ~200 MeV for SU(2))
- **τ_int(Q²) vs β** — Autocorrelation time (increases with β → freezing)
- **Q histograms** — Distribution narrowing indicates freezing

The analysis uses the Gamma method via `pyerrors` for autocorrelation-aware
errors. Periodic-boundary runs use integer-rounded `Q_re`; open-boundary runs use
the continuous `alpha * Q` because global Q is not quantized with open temporal
boundaries.

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
├── input/                # Setup files and topcharge params
├── run_heatbath_scan.sh           # Phase 1: config generation
├── run_topcharge_scan.sh          # Phase 2: topcharge measurement
└── run_heatbath_topcharge_scan.sh # SU(2) fused pipeline (no disk configs)
```
