# Additional Notes: Bug Findings & Implementation Details

---

## 1. Bug Scan Results

### Fixed

**APE smearing index guard used `index_mt_3` as a boolean** (`_Utility/src/smearing_techniques.cc`)

Three spatial-link smearing branches checked `index_mt_3` instead of
`index_mt_3 >= 0`. In periodic runs that could skip a valid link when the
computed offset was exactly zero; in open-boundary runs it could treat `-1` as
true. Fixed all three checks. After rebuilding, the full CTest suite passes.

**Detected thermalization cut was not preserved for topcharge scans** (`run_heatbath_scan.sh`, `run_topcharge_scan.sh`)

`run_heatbath_scan.sh` detected `START_CONF` from plaquette thermalization and
wrote it to `run_info.txt`, but the generated `input.txt` did not contain
`start_conf`. Then `run_topcharge_scan.sh` re-patched `input.txt` using the
shared topcharge parameter file, so the per-run thermalization cut could be lost.
Fixed by writing `start_conf` into the heatbath-generated `input.txt` and making
the topcharge scan prefer that per-run value, falling back to the shared params
file only for legacy runs.

**Parallel scan failures are reported without aborting completed work** (`run_heatbath_scan.sh`, `run_topcharge_scan.sh`, `run_heatbath_topcharge_scan.sh`)

Parallel mode counts failed child jobs and prints a warning, but deliberately
does not exit with status 1 at the end. This keeps long production scans from
wasting resources when most jobs completed successfully; inspect the per-run log
files for the failed subset and rerun only those ensembles.

**Build metadata still described an SU(2)-only project** (`CMakeLists.txt`, `scripts/build.sh`)

The build already produced SU(3) executables, but the project name, install
targets, and build-script executable list were stale. Updated them to include
SU(2) and SU(3) binaries consistently.

**chi_t used `Q_raw²` instead of `Q_rescaled²`** (`analysis/calculations.py`)

`Q_raw` is the smeared but not-yet-rescaled float value (e.g. −2.71). `Q_rescaled = round(α · Q_raw)` is the integer-rounded charge (e.g. −2). The susceptibility, tau_int, and histograms must all use the same Q. Since `Q_raw ≈ Q_rescaled / α`, using `Q_raw²` inflated chi_t by `1/α²` (≈ +24% for α = 0.898) and chi_t^(1/4) by ≈ +6%. Fixed to use `Q_rescaled²` everywhere in `analyze_run()`.

---

### Not bugs (investigated and ruled out)

**Topcharge normalization `1/(4π²)`** — correct for SU(N) in the fundamental representation with `Tr(T^a T^b) = δ^{ab}/2`. The local density formula sums over the three independent dual pairs (01,23), (02,13), (03,12), each with a factor of 2 to account for the (μν) and (νμ) contributions. The overall sign is a convention choice (the code defines the local density as `−ε Tr(FF)`).

**`chi_t_fourth_root_a` dimensional chain** — the division by `a` happens in `plotting_code.py` (`chi_fourth_MeV = chi_fourth_a / a_values`). The stored intermediate `chi_t_fourth_root_a = (chi_t_lattice)^(1/4) · ℏc` has units MeV·fm; dividing by `a` [fm] gives MeV. Correct.

**`start_conf` not filtering in analysis** — the measurement binaries process only the configuration range written into each run's `input.txt`. Analysis-side `start_conf` is metadata only and does not need to filter the already-produced `topcharge[_su3].dat` files.

**SU(3) link projection (`su3_proj`)** — the Gram-Schmidt approach (normalise row 0, orthogonalise and normalise row 1, construct row 2 as the conjugate cross-product) is the standard correct method for projecting onto SU(3). Verified mathematically sound.

**Progress bar** — `conf_count` increments to `total_confs` on the last iteration, giving progress = 1.0 (100%). Correct.

**Alpha rescaling range [0.8, 1.2]** — all current ensembles have α within this window (observed range: 0.80–1.13). For significantly coarser or finer lattices in the future, this window should be widened. No action needed now.

---

### Design decisions (not bugs, but worth knowing)

**Jackknife ignores autocorrelation** — `jackknife_error()` is documented to ignore autocorrelations. It is only used as a rough secondary estimate alongside the Gamma method. The Gamma method (`autocorrelation.py` via `pyerrors`) is the primary error estimator and correctly accounts for autocorrelations.

**Gamma method `S = 1.5`** — the windowing parameter S = 1.5 is the `pyerrors` default and appropriate for most ensembles. It can be adjusted if the automatic window selection behaves poorly for very long autocorrelations.

**Thermalization detection rounds up** — `detect_thermalization.py` rounds the detected thermalization point up to the nearest `save_interval`. This conservatively discards a few sweeps of already-thermalised data but guarantees no contamination.

---

## 2. Open Boundary Conditions

### What type of BC?

The implementation is the **Lüscher-Schaefer open boundary condition** (Lüscher & Schaefer, JHEP 2011, arXiv:1105.4749). It is **neither Dirichlet nor von Neumann** in the classical PDE sense.

- **Dirichlet** would fix link values at the boundary to a constant (e.g. U = 1). The boundary links here are fully free to update.
- **von Neumann** sets the normal derivative to zero (no flux). The Lüscher-Schaefer BC is effectively the opposite — it removes the action barrier to topology flowing through the boundary.

### How it is implemented

Temporal plaquettes at `t = 0` and `t = T−1` enter the action with weight **1/2** instead of 1. In the heatbath this is done by multiplying each staple at the boundary by 0.5:

```cpp
// MC_heatbath_su3.cc  (identical logic in MC_heatbath.cc)
const bool at_boundary = open_bc && (it == 0 || it == T_size - 1);
...
if (at_boundary) su3_ti_eq_re(T2, 0.5);   // half-weight boundary staple
```

The factor of 1/2 is physically correct: a bulk plaquette is shared between two adjacent time slices and contributes to both their local actions; a boundary plaquette has only one neighbour, so its natural weight is half.

### Why it cures topological freezing

With **periodic BC**, the topological charge Q is strictly an integer and the sectors Q = 0, ±1, ±2, … are separated by large action barriers that diverge as a → 0. Local updates cannot tunnel between sectors — τ_int diverges.

With **open BC**, the boundary removes the topological sector structure. Topological charge can flow continuously through the temporal boundaries, so Q is no longer quantized globally and the barriers disappear. τ_int remains finite even at fine lattice spacings.

### In the topological charge measurement

The boundary time slices are excluded from the topological charge sum via the `exclude_boundary_slices` parameter. The gauge field near the boundary is affected by the open BC and does not represent bulk physics, so those slices are discarded to avoid contaminating the measurement.

```
# base_params_su3.txt
exclude_boundary_slices 0    # set > 0 for open BC runs (typically 1–2)
```

---

## 4. Scale Setting Utilities (`scale_set.ipynb`)

Three utility functions added to the notebook for interactive use:

### `beta_from_a_su2(a_fm)`
Analytically inverts `lattice_spacing_su2` by solving the quadratic in `delta_beta = beta - 2.6`:

```
ln(t0/a²) = 1.285 + 6.409·x - 0.7411·x²   →   quadratic formula
```

Picks the root closest to the physical regime (beta ~ 2–3).

### `beta_from_a_su3(a_fm)`
Inverts `lattice_spacing_su3` via `numpy.roots` on the cubic in `delta_beta = beta - 6.0`:

```
ln(a/r0) = -1.6804 - 1.7331·x + 0.7849·x² - 0.4428·x³
```

Picks the real root closest to x = 0 (physical regime beta ~ 5.5–6.5).

### `adjust_lattice(T_ref, L_ref, a_ref_fm, a_new_fm, open_bc, n_exclude)`
Given a reference volume at `a_ref`, computes integer `(T_new, L_new)` at `a_new` that preserve physical volume:

- **Spatial**: `L_new = round(L_ref * a_ref / a_new)` (no boundary effect)
- **Periodic temporal**: `T_new = round(T_ref * a_ref / a_new)`
- **Open BC temporal**: physical extent is `T_eff = (T - 2*n_exclude) * a`, so `T_new = round(T_eff_ref * a_ref/a_new) + 2*n_exclude`. `n_exclude` must match `exclude_boundary_slices` in the input file.

---

## 5. Timeslice Topological Charge Density

### C++ output (`meas_topcharge_su2.cc`)

After computing the total Q for each config, `meas_topcharge` now also writes `topcharge_timeslice.dat` in the same output directory:

```
# smear_steps  config_number  t  q_t
```

`q_t` is the charge density summed over all spatial sites at timeslice `t`, normalized by `4π²`:

```
q(t) = (1/4π²) · Σ_{x,y,z} q_local(t,x,y,z)
```

For open BC, only timeslices `t ∈ [exclude_boundary_slices, T - exclude_boundary_slices)` are written. **No new input parameters required** — the existing `boundary` and `exclude_boundary_slices` fields control this automatically.

### Python analysis (`analysis/timeslice_analysis.py`)

Full pipeline: load → bin → Gamma method errors → plot.

**Key design choice — temporal binning**: because the spatial volume is small, individual timeslice measurements are very noisy. The `n_bin` parameter sums over `n_bin` consecutive slices per bin. Typical values: `n_bin=2` or `3`. Match `n_exclude` in the call to `analyse_timeslices` with `exclude_boundary_slices` in the input file.

**Error estimation**: `pyerrors` Gamma method applied independently to each time bin across MC configurations. Properly accounts for MC autocorrelation (large `tau_int` near freezing). Same `S=1.5` default as the rest of the analysis.

**CLI usage**:
```bash
python analysis/timeslice_analysis.py \
    data/results/T16_L4_b2.50_open_seed12345/output/topcharge_timeslice.dat \
    --a 0.093 --n-bin 2 --open-bc --n-exclude 2 --output output/figures_analysis/timeslice.pdf
```

**Programmatic usage**:
```python
from analysis.timeslice_analysis import analyse_timeslices, plot_timeslice_density

res = analyse_timeslices("topcharge_timeslice.dat", lattice_spacing_fm=0.093,
                         n_bin=2, open_bc=True, n_exclude=2)
plot_timeslice_density(res, open_bc=True)
```
