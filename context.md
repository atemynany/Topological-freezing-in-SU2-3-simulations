# Project Context: Topological Charge Freezing in SU(2)/SU(3) Lattice Gauge Theory

This document summarises the theoretical background, implementation, and analysis pipeline for the master thesis project studying topological freezing in pure-gauge SU(2) and SU(3) lattice field theory.

---

## 1. Physics Goal

The project measures the **topological susceptibility** χ_t and the **integrated autocorrelation time** τ_int(Q²) across a range of lattice spacings for both SU(2) and SU(3) pure-gauge theories. The central question is how topological freezing (the inability of the Monte Carlo to tunnel between topological sectors) sets in as the continuum limit is approached, and whether open boundary conditions mitigate it.

---

## 2. Topological Charge

The topological charge Q is an integer-valued observable that characterises the global winding number of a gauge configuration:

$$Q = \frac{1}{32\pi^2} \int d^4x \, \varepsilon_{\mu\nu\rho\sigma} \, \mathrm{Tr}(F_{\mu\nu} F_{\rho\sigma})$$

On the lattice, this is approximated using the **clover (symmetric) discretisation** of F_μν:

$$F_{\mu\nu}^{\text{clover}} = \frac{1}{4}(P_{++} + P_{-+} + P_{--} + P_{+-})$$

where each P is a 1×1 plaquette oriented differently around the site. The field strength is then extracted as the anti-Hermitian part: `(C - C†)/2`. The local charge density sums over dual pairs of field strength components, and the total Q is normalised by 1/(4π²).

Without smearing, the clover Q is dominated by UV lattice noise and not close to an integer. **APE smearing** is applied to the gauge field before measuring Q in production.

---

## 3. APE Smearing

APE smearing iteratively replaces each link with a weighted average of itself and its staples, then projects back to SU(N):

$$U_\mu'(x) = \mathcal{P}_{SU(N)}\!\left[(1-\alpha)\,U_\mu(x) + \frac{\alpha}{6}\sum_{\nu\neq\mu} V_{\mu\nu}(x)\right]$$

where V_μν is the sum of the two staples in the ν direction. The SU(3) projection uses iterative polar decomposition (via Cayley-Hamilton). Typically 20 steps at α = 0.3 (SU3) or α = 0.45 (SU2) are used. Smearing removes UV noise, making Q converge to near-integer values.

**Important**: smearing is applied to the gauge field. You cannot smear Q values directly. The pipeline is always: gauge field → smear → compute Q.

---

## 4. Topological Freezing and Critical Slowing Down

As β → ∞ (continuum limit, a → 0), the action barriers between topological sectors diverge. Local update algorithms (heatbath, overrelaxation) cannot tunnel between sectors, so Q becomes frozen. The autocorrelation time scales as:

$$\tau_{\text{int}}(Q) \propto a^{-z}, \quad z \approx 5 \text{ (SU(3) pure gauge)}$$

**Open boundary conditions (OBC)** in the temporal direction allow topological charge to flow through the boundaries, breaking the topological sector structure and avoiding the divergence of τ_int. The tradeoff is reduced translational invariance and boundary effects that must be excluded from measurements.

---

## 5. Topological Susceptibility

$$\chi_t = \frac{\langle Q^2 \rangle}{V}, \quad \chi_t^{1/4} \approx 200 \text{ MeV (SU(2))}, \quad \approx 180 \text{ MeV (SU(3))}$$

where V = T·L³ is the lattice volume. In code, χ_t is computed from the rescaled integer charges Q_re (see §6).

---

## 6. Q Rescaling (Alpha Estimator)

Raw smeared Q values are close to but not exactly integers due to residual discretisation errors. A rescaling factor α is found that minimises the deviation from the nearest integer:

$$\alpha^* = \arg\min_\alpha \langle (\alpha Q - \text{round}(\alpha Q))^2 \rangle$$

This is solved via **golden-section search** on [0.8, 1.2]. The rescaled integer charge is then:

$$Q_{\text{re}} = \text{round}(\alpha^* Q)$$

α ≈ 1 for well-smeared configurations. This is only meaningful when Q is close to integers (i.e., after sufficient smearing). Applying it to unsmeared thermalization Q is not physically meaningful.

The Rust extension (`rust_analysis/`) provides a parallel implementation of this search via rayon.

---

## 7. Scale Setting

### SU(3) — Necco-Sommer formula (Sommer scale r₀ = 0.5 fm)

$$\ln\!\left(\frac{a}{r_0}\right) = -1.6804 - 1.7331(\beta-6) + 0.7849(\beta-6)^2 - 0.4428(\beta-6)^3$$

### SU(2) — Gradient flow scale t₀

$$\ln\!\left(\frac{t_0}{a^2}\right) = 1.285 + 6.409(\beta - 2.600) - 0.7411(\beta - 2.600)^2, \quad t_0^{1/2} \approx 0.1065 \text{ fm}$$

Both are implemented in `analysis/calculations.py` as `lattice_spacing_su2/su3(beta)`.

---

## 8. Autocorrelation Analysis

τ_int is estimated using the **Gamma method** (via the `pyerrors` library), which provides an unbiased estimate of the integrated autocorrelation time including the tail contribution. The input observable is Q²_re (the squared rescaled topological charge). τ_int > 0.5 indicates autocorrelations; large τ_int signals freezing.

A Rust-accelerated jackknife error estimate (O(N) via the delete-one identity) is also available in `rust_analysis/`.

---

## 9. Monte Carlo Algorithm

### SU(2) Heatbath
Implemented in `src/MC_heatbath.cc`. Uses the **Kennedy-Pendleton** algorithm to exactly sample from the heatbath distribution for each SU(2) link. Optionally uses open boundary conditions (temporal links at t=0 and t=T−1 get half-weight staples).

### SU(3) Heatbath
Implemented in `src/MC_heatbath_su3.cc`. Uses the **Cabibbo-Marinari** decomposition: the SU(3) matrix is updated by cycling through three embedded SU(2) subgroups, each updated with the SU(2) heatbath. Followed by overrelaxation steps (reflect link through the staple sum, stays exactly on SU(3) manifold).

**OpenMP parallelization**: The SU(3) heatbath uses checkerboard (even/odd) site decomposition for thread-safe parallel sweeps. Sites are classified by parity `(t+x+y+z) % 2` at startup. Each half-sweep updates only even or odd sites, so no two threads touch neighboring links simultaneously. Both the heatbath and overrelaxation loops are parallelized this way. Per-thread RNG seeding uses `rlxd_init(2, seed + tid * 997)`.

Both programs:
- Write the plaquette and Wilson action density at every sweep to `plaquette[_su3].dat`
- Save gauge configurations every `save_interval` sweeps
- Validate all input parameters at startup

---

## 10. Production Topological Charge Measurement

`src/meas_topcharge_su3.cc` and `src/meas_topcharge_su2.cc`:
- Load saved configurations from disk
- Apply N APE smearing steps to a temporary copy of the gauge field
- Measure Q on the smeared copy
- Write `conf smear_step Q` to `topcharge_su3.dat` (SU3) or `smear_step conf Q plaquette` to `topcharge.dat` (SU2)
- Display a progress bar during the run

The smearing does **not** modify the stored configurations — it operates on an in-memory copy.

---

## 11. Thermalization

Thermalization is detected automatically by `scripts/detect_thermalization.py` using a sliding-window test on the plaquette. The detected `start_conf` is written to `run_info.txt` and used as the starting configuration for the topological charge measurement.

Thermalization plot per gauge group (saved to `output/figures_analysis/`):
- **`thermalization_comparison_su2/3.png`**: plaquette vs MC sweep for all β values overlaid, with detected thermalization points marked.

---

## 12. Analysis Pipeline

`analysis/analysis.py --su2/--su3`:
1. Scans `data/results/T*_L*_b*_*_seed*/` for completed runs
2. Loads `topcharge[_su3].dat` — auto-detects SU2 (4-col) or SU3 (3-col) format
3. Applies alpha rescaling to get Q_re
4. Computes χ_t = ⟨Q²_re⟩/V and τ_int via Gamma method
5. Generates plots — histogram and Q-vs-mctime grids use a **shared axis range** across all β panels so narrowing of distributions with increasing β is visually apparent

`analysis/calculations.py` — data loading, QEstimator, jackknife, autocorrelation
`analysis/plotting_code.py` — all matplotlib figures
`analysis/autocorrelation.py` — Gamma method wrapper

---

## 13. File Formats

| File | Columns | Notes |
|------|---------|-------|
| `topcharge.dat` (SU2) | `smear_step conf Q plaquette` | 4 columns |
| `topcharge_su3.dat` | `conf smear_step Q` | 3 columns |
| `topcharge_timeslice.dat` | `smear_step conf t q_t` | per-timeslice topological charge density |
| `plaquette[_su3].dat` | `sweep plaquette action_density` | every heatbath sweep; `s = 6 beta (1 - <P>)` |

---

## 14. Rust Extension (`rust_analysis/`)

Optional Python extension module built with maturin/PyO3. Provides parallel versions of analysis hot-paths using rayon. Automatically used by `calculations.py` when installed.

Functions: `find_optimal_alpha`, `jackknife_error`, `bootstrap_alpha`, `bootstrap_Q2`

Build: `cd rust_analysis && maturin develop --release`

---

## 15. Key References

1. **Schaefer, Sommer, Virotta (2011)** — Critical slowing down and error analysis. Nucl. Phys. B 845. [arXiv:1009.5228]
2. **Lüscher, Schaefer (2011)** — Lattice QCD without topology barriers. JHEP 07. [arXiv:1105.4749]
3. **Necco, Sommer (2001)** — SU(3) scale setting. Nucl. Phys. B 622. [arXiv:hep-lat/0108008]
4. **Kennedy, Pendleton (1985)** — Improved heat bath method. Phys. Lett. B 156.
5. **Cabibbo, Marinari (1982)** — SU(3) heatbath via SU(2) subgroups. Phys. Lett. B 119.
