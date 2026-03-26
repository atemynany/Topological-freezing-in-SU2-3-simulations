# SU(2) Parallelization Design

## Goal

Speed up SU(2) heatbath configuration generation for large lattices (up to 90^4) by:
1. Running independent (beta, seed) combinations simultaneously on a cluster
2. Parallelizing each individual heatbath sweep across CPU cores via checkerboard decomposition

SU(2) only. SU(3) is out of scope.

## Deliverables

### 1. SLURM Job Array Script (`scripts/slurm_beta_scan.sh`)

Submits the full beta scan as a SLURM job array. One array index per (beta, seed) pair.

**Behavior per job:**
- Read beta and seed from a generated parameter list (indexed by `$SLURM_ARRAY_TASK_ID`)
- Run `mc_heatbath` with input parameters
- Run `detect_thermalization.py` on the plaquette output
- Run `meas_topcharge` on post-thermalization configs
- Each step checks the exit code of the previous step. On failure: log the error to `job_errors.log` in the run directory and exit that job. Other array jobs are unaffected.

**Submission:**
```bash
sbatch --array=0-N scripts/slurm_beta_scan.sh
```

where N = (number of betas × number of seeds) - 1.

**Configuration:** The script reads from `input/base_params_su2.txt` and `input/beta_scan_su2.txt` plus a seeds list. SLURM resource requests (nodes, cores, walltime) are set at the top of the script and should be adjusted per cluster.

### 2. Checkerboard OpenMP in `MC_heatbath.cc`

**Algorithm:**

Every site gets a parity: `parity = (t + x + y + z) % 2`. Even sites have only odd neighbors and vice versa on a 4D hypercubic lattice. This allows updating all sites of one parity in parallel without data races.

**Changes to `src/MC_heatbath.cc`:**

1. At initialization, precompute two arrays:
   - `even_sites[]` — indices where `(t+x+y+z) % 2 == 0`
   - `odd_sites[]` — indices where `(t+x+y+z) % 2 == 1`
   Each has `volume/2` entries.

2. Replace the single sweep loop:
   ```cpp
   // BEFORE (sequential)
   for (int site = 0; site < volume; site++) {
       for (int mu = 0; mu < 4; mu++) {
           // compute staple, heatbath update
       }
   }
   ```
   with two half-sweeps:
   ```cpp
   // AFTER (parallel, checkerboard)
   // even half-sweep
   #pragma omp parallel for schedule(static)
   for (int i = 0; i < n_even; i++) {
       int site = even_sites[i];
       for (int mu = 0; mu < 4; mu++) {
           // compute staple, heatbath update (uses per-thread RNG)
       }
   }
   // odd half-sweep
   #pragma omp parallel for schedule(static)
   for (int i = 0; i < n_odd; i++) {
       int site = odd_sites[i];
       for (int mu = 0; mu < 4; mu++) {
           // compute staple, heatbath update (uses per-thread RNG)
       }
   }
   ```

3. **Per-thread RNG:** The current RANLUX implementation (`_Utility/src/ranlxd.c`) stores its state in file-scope `static` variables — a single global state. Calling `ranlxd()`/`DRand()` from multiple OpenMP threads is a data race.

   **Fix:** Add `_Thread_local` (C11) to the static state variables in `ranlxd.c`. Each thread then gets its own independent RANLUX stream. At the start of each `#pragma omp parallel` region, each thread calls `rlxd_init(2, seed + omp_get_thread_num())` to seed its stream. This is minimally invasive — same RNG algorithm, just thread-local state.

   Files affected: `_Utility/src/ranlxd.c` (add `_Thread_local` to statics), `_Utility/src/ranlxs.c` (same treatment for consistency).

4. **Open boundary conditions:** The checkerboard split is BC-agnostic — even/odd parity doesn't depend on boundary type. The existing boundary weight factors (half-weight staples at `t=0` and `t=T-1`) are preserved since they live inside the per-site staple computation, not in the loop structure.

5. **Implicit barrier:** The implicit OpenMP barrier at the end of each `#pragma omp parallel for` ensures all even-site updates complete before odd-site updates begin.

6. **Thread count:** The SLURM script sets `OMP_NUM_THREADS` explicitly to match the requested cores per node.

7. **Comments:** Minimal `//` comments on changed lines only. No bloat.

**Physics correctness:** The checkerboard sweep satisfies detailed balance identically to the sequential sweep — same heatbath accept/reject, same Boltzmann weight. Only the traversal order changes. This is a standard technique in lattice QCD (used by openQCD, MILC, etc.).

### 3. Monitoring Script (`scripts/monitor_runs.py`)

Reads live output files from running jobs and plots thermalization progress.

**Usage:**
```bash
python3 scripts/monitor_runs.py                    # all active runs
python3 scripts/monitor_runs.py --beta 2.5         # specific beta
python3 scripts/monitor_runs.py --run-dir data/results/T16_L16_b2.500_periodic_seed42/
```

**Behavior:**
- Scans `data/results/` for directories with `plaquette.dat`
- For each run, reads `plaquette.dat` and `therm_topcharge.dat` (which grow as the heatbath runs)
- Produces a 2-panel plot per run: plaquette vs sweep (top), unsmeared Q vs sweep (bottom)
- Saves to `output/monitoring/` as PNG files
- Works over SSH without X forwarding (saves files, no interactive display required)

**Dependencies:** numpy, matplotlib (already required).

### 4. Correctness Tests

#### a. New Catch2 test: `tests/test_checkerboard.cc`

- Verify even/odd site classification: every even site's neighbors are odd and vice versa
- Verify that after a checkerboard sweep, all links remain valid SU(2) matrices (determinant = 1, unitarity)
- Verify link update count: every link updated exactly once per sweep

#### b. Validation comparison

A short script or test that:
- Runs 500 sweeps on 8^4 at a fixed beta with the sequential code (before checkerboard changes)
- Runs 500 sweeps on 8^4 at the same beta with checkerboard code
- Compares: average plaquette should agree within statistical uncertainty
- Note: exact sweep-by-sweep match is not expected because thread scheduling changes RNG consumption order. Statistical agreement is the correct criterion.

## Files Modified

- `src/MC_heatbath.cc` — checkerboard sweep, per-thread RNG
- `_Utility/src/ranlxd.c` — `_Thread_local` on static state variables
- `_Utility/src/ranlxs.c` — same `_Thread_local` treatment
- `CMakeLists.txt` — add `test_checkerboard` target

## Files Created

- `scripts/slurm_beta_scan.sh` — SLURM array job script
- `scripts/monitor_runs.py` — live monitoring/plotting
- `tests/test_checkerboard.cc` — checkerboard correctness tests

## Out of Scope

- SU(3) (will be done separately later)
- MPI / multi-node parallelism within a single run
- GPU acceleration
- Changes to the analysis pipeline
- Changes to topcharge measurement code
