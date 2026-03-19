#!/usr/bin/env python3
"""
Data loading and analysis calculations for SU(2)/SU(3) topological charge.

Includes:
- Q Estimator: Rescaled topological charge (Q_rescaled = round(α * Q̂))
- Q Statistics: Statistical analysis with gamma method
- Topological Susceptibility: χ_t = ⟨Q²⟩ / V
- Data loading and run analysis
"""

import numpy as np
import pyerrors as pe
import os
import re
from dataclasses import dataclass
from typing import Optional, Tuple

from analysis.autocorrelation import autocorrelation

# ---------------------------------------------------------------------------
# Optional Rust extension (built with maturin in rust_analysis/).
# Falls back to pure Python if not installed.
# To build: cd rust_analysis && maturin develop --release
# ---------------------------------------------------------------------------
try:
    import qcd_analysis as _rust
    _HAS_RUST = True
except ImportError:
    _rust = None
    _HAS_RUST = False

# Physical constants
HBAR_C = 197.3  # MeV * fm

# Reference values for topological susceptibility (fourth root, in MeV)
CHI_T_FOURTH_ROOT_SU2 = 200.0
CHI_T_FOURTH_ROOT_SU3 = 180.0


# =============================================================================
# Q Estimator
# =============================================================================

def find_optimal_alpha(Q_raw: np.ndarray, alpha_min=0.8, alpha_max=1.2, tol=1e-6) -> float:
    """Find α that minimizes ⟨(αQ̂ - round(αQ̂))²⟩.

    Uses the Rust extension (parallel rayon) when available, pure Python otherwise.
    """
    if _HAS_RUST:
        return _rust.find_optimal_alpha(Q_raw.astype(np.float64), alpha_min, alpha_max, tol)

    # Pure-Python fallback
    def objective(alpha):
        scaled = alpha * Q_raw
        return np.mean((scaled - np.round(scaled))**2)

    phi = (np.sqrt(5) - 1) / 2
    a, b = alpha_min, alpha_max
    c, d = b - phi*(b-a), a + phi*(b-a)

    while abs(b - a) > tol:
        if objective(c) < objective(d):
            b, d = d, c
            c = b - phi*(b-a)
        else:
            a, c = c, d
            d = a + phi*(b-a)

    return (a + b) / 2


def bootstrap_alpha(Q_raw: np.ndarray, n_boot: int = 500,
                    alpha_min=0.8, alpha_max=1.2, tol=1e-6) -> Tuple[float, float]:
    """Bootstrap standard error on the optimal α.

    Returns (alpha_central, alpha_std).
    Requires the Rust extension; raises RuntimeError if not installed.
    """
    if not _HAS_RUST:
        raise RuntimeError(
            "bootstrap_alpha requires the Rust extension. "
            "Build it with: cd rust_analysis && maturin develop --release"
        )
    return _rust.bootstrap_alpha(Q_raw.astype(np.float64), n_boot, alpha_min, alpha_max, tol)


def bootstrap_Q2(Q_raw: np.ndarray, n_boot: int = 500,
                 alpha_min=0.8, alpha_max=1.2, tol=1e-6) -> Tuple[float, float]:
    """Bootstrap standard error on ⟨Q_re²⟩.

    Returns (Q2_central, Q2_std).
    Requires the Rust extension; raises RuntimeError if not installed.
    """
    if not _HAS_RUST:
        raise RuntimeError(
            "bootstrap_Q2 requires the Rust extension. "
            "Build it with: cd rust_analysis && maturin develop --release"
        )
    return _rust.bootstrap_Q2(Q_raw.astype(np.float64), n_boot, alpha_min, alpha_max, tol)


@dataclass
class QEstimatorResult:
    alpha: float
    Q_scaled: np.ndarray
    Q_rescaled: np.ndarray
    mean_deviation_squared: float
    
    @property
    def rms_deviation(self) -> float:
        return np.sqrt(self.mean_deviation_squared)
    
    @property
    def Q2_mean(self) -> float:
        return np.mean(self.Q_rescaled ** 2)


class QEstimator:
    """Rescaled topological charge estimator."""
    
    def __init__(self, Q_raw: np.ndarray):
        alpha = find_optimal_alpha(Q_raw)
        Q_scaled = alpha * Q_raw
        Q_rescaled = np.round(Q_scaled).astype(int)
        
        self._result = QEstimatorResult(
            alpha=alpha,
            Q_scaled=Q_scaled,
            Q_rescaled=Q_rescaled,
            mean_deviation_squared=np.mean((Q_scaled - Q_rescaled)**2)
        )
    
    @property
    def alpha(self) -> float:
        return self._result.alpha
    
    @property
    def Q_scaled(self) -> np.ndarray:
        return self._result.Q_scaled
    
    @property
    def Q_rescaled(self) -> np.ndarray:
        return self._result.Q_rescaled
    
    @property
    def result(self) -> QEstimatorResult:
        return self._result


# =============================================================================
# Q Statistics
# =============================================================================

@dataclass
class QStats:
    Q_mean: float
    Q_error: float
    Q_error_jack: float
    Q2_mean: float
    Q2_error: float
    Q2_error_jack: float
    tau_int: float
    n_eff: float
    alpha: float


def jackknife_error(data: np.ndarray) -> float:
    """Compute jackknife error estimate (ignores autocorrelation).

    O(N) implementation via the delete-one mean identity.
    Uses Rust when available (parallel); otherwise pure Python O(N).
    """
    if _HAS_RUST:
        return _rust.jackknife_error(data.astype(np.float64))

    # Pure-Python O(N) fallback (avoids the original O(N²) loop)
    n = len(data)
    if n < 2:
        return 0.0
    total = data.sum()
    jack = (total - data) / (n - 1)
    return float(np.sqrt((n - 1) / n * np.sum((jack - jack.mean()) ** 2)))


def compute_Q_statistics(Q_raw: np.ndarray, ensemble: str = "ens", S: float = 1.5) -> QStats:
    """Compute Q statistics using gamma method for error estimation."""
    est = QEstimator(Q_raw)
    Q = est.Q_rescaled.astype(float)
    
    q_obs = pe.Obs([Q], [ensemble])
    q_obs.gamma_method(S=S)
    
    q2_obs = pe.Obs([Q**2], [ensemble])
    q2_obs.gamma_method(S=S)
    
    tau = q_obs.e_tauint[ensemble]
    
    Q_jack_err = jackknife_error(Q)
    Q2_jack_err = jackknife_error(Q**2)
    
    return QStats(
        Q_mean=q_obs.value,
        Q_error=q_obs.dvalue,
        Q_error_jack=Q_jack_err,
        Q2_mean=q2_obs.value,
        Q2_error=q2_obs.dvalue,
        Q2_error_jack=Q2_jack_err,
        tau_int=tau,
        n_eff=len(Q) / (2 * tau) if tau > 0 else len(Q),
        alpha=est.alpha
    )


def load_topcharge(filepath: str, smear: int = None) -> np.ndarray:
    """Load topological charge from file.

    Handles two column formats automatically:
      SU2 (4 cols): smear_step  conf  Q  plaquette
      SU3 (3 cols): conf  smear_step  Q
    """
    configs, Q_vals, smears = [], [], []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split()
            try:
                if len(parts) == 3:
                    # SU3 format: conf smear_step Q
                    c, s, q = int(parts[0]), int(parts[1]), float(parts[2])
                elif len(parts) >= 4:
                    # SU2 format: smear_step conf Q [plaquette ...]
                    s, c, q = int(parts[0]), int(parts[1]), float(parts[2])
                else:
                    continue
                if smear is None or s == smear:
                    smears.append(s)
                    configs.append(c)
                    Q_vals.append(q)
            except: pass
    if smear is None and smears:
        max_smear = max(smears)
        configs = [c for s, c in zip(smears, configs) if s == max_smear]
        Q_vals = [q for s, q in zip(smears, Q_vals) if s == max_smear]
    idx = np.argsort(configs)
    return np.array(Q_vals)[idx]


# =============================================================================
# Topological Susceptibility
# =============================================================================

def topological_susceptibility(Q2_mean: float, T: int, L: int, lattice_spacing: float) -> dict:
    """Compute topological susceptibility χ_t = ⟨Q²⟩ / V."""
    volume = T * L ** 3
    chi_t_lattice = Q2_mean / volume
    chi_t_fm4 = chi_t_lattice / (lattice_spacing ** 4)
    chi_t_fourth_root_MeV = (chi_t_fm4 ** 0.25) * HBAR_C
    
    return {
        'chi_t_lattice': chi_t_lattice,
        'chi_t_fm4': chi_t_fm4,
        'chi_t_fourth_root_MeV': chi_t_fourth_root_MeV,
        'reference_SU2_MeV': 200.0
    }


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class RunData:
    beta: float
    seed: int
    run_dir: str
    T: int
    L: int
    boundary: str
    Q_raw: np.ndarray
    Q_rescaled: np.ndarray
    alpha: float
    plaq_data: np.ndarray
    mc_time: np.ndarray
    start_conf: int


@dataclass
class AnalysisResult:
    beta: float
    boundary: str
    Q_mean: float
    Q2_mean: float
    chi_t_lattice: float
    chi_t_fourth_root_a: float
    tau_int: float
    dtau_int: float
    alpha: float
    n_conf: int


def gaussian(x, mu, sigma, A):
    """Gaussian function for histogram fitting."""
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def parse_input_file(run_dir: str) -> dict:
    """Parse input.txt file from a run directory."""
    info = {}
    input_file = os.path.join(run_dir, "input.txt")
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    info[parts[0]] = parts[1]
    return info


def load_run_data(run_dir: str, gauge_group: str = "su2", min_sweeps: int = 1000) -> Optional[RunData]:
    """Load all data from a run directory."""
    info = parse_input_file(run_dir)
    
    num_sweeps = int(info.get('num_sweeps', 0))
    if num_sweeps < min_sweeps:
        return None
    
    # Parse boundary from directory name (e.g., T16_L16_b2.50_open_seed12345)
    dirname = os.path.basename(run_dir)
    if '_open_' in dirname or dirname.endswith('_open'):
        boundary = 'open'
    elif '_periodic_' in dirname or dirname.endswith('_periodic'):
        boundary = 'periodic'
    else:
        boundary = info.get('boundary', 'periodic')
    
    candidates = [
        os.path.join(run_dir, "output", "topcharge.dat"),
        os.path.join(run_dir, "output", f"topcharge_{gauge_group}.dat"),
    ]
    topcharge_file = next((f for f in candidates if os.path.exists(f)), None)
    
    if not topcharge_file:
        return None
    
    Q_raw = load_topcharge(topcharge_file)
    if len(Q_raw) < 50:
        return None
    
    estimator = QEstimator(Q_raw)
    
    match = re.search(r'_b([\d.]+)_', os.path.basename(run_dir))
    beta_from_dir = float(match.group(1)) if match else 0
    
    beta = float(info.get('beta', beta_from_dir))
    seed = int(info.get('seed', 0))
    T = int(info.get('T', 16))
    L = int(info.get('L', 16))
    start_conf = int(info.get('start_conf', 50))
    
    # Load plaquette data
    plaq_candidates = [
        os.path.join(run_dir, "output", "plaquette.dat"),
        os.path.join(run_dir, "output", f"plaquette_{gauge_group}.dat"),
    ]
    plaq_file = next((f for f in plaq_candidates if os.path.exists(f)), None)
    
    mc_time, plaq_data = np.array([]), np.array([])
    if plaq_file:
        mc_list, plaq_list = [], []
        with open(plaq_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mc_list.append(int(parts[0]))
                        plaq_list.append(float(parts[1]))
                    except ValueError:
                        continue
        mc_time = np.array(mc_list)
        plaq_data = np.array(plaq_list)
    
    return RunData(
        beta=beta, seed=seed, run_dir=run_dir, T=T, L=L,
        boundary=boundary, Q_raw=Q_raw, Q_rescaled=estimator.Q_rescaled, alpha=estimator.alpha,
        plaq_data=plaq_data, mc_time=mc_time, start_conf=start_conf
    )


def analyze_run(run_data: RunData) -> AnalysisResult:
    """Perform full analysis on a single run."""
    volume = run_data.T * run_data.L ** 3
    Q = run_data.Q_rescaled.astype(float)   # integer-rounded charge, consistent everywhere
    Q2 = Q ** 2
    Q2_mean = np.mean(Q2)
    chi_t_lattice = Q2_mean / volume

    chi_t_fourth_root_a = (chi_t_lattice ** 0.25) * HBAR_C if chi_t_lattice > 0 else 0.0
    
    try:
        autocorr = autocorrelation(Q2, name=f"Q2_b{run_data.beta}")
        tau_int = autocorr.tau_int
        dtau_int = autocorr.dtau_int
    except:
        tau_int, dtau_int = 0.5, 0.0
    
    return AnalysisResult(
        beta=run_data.beta, boundary=run_data.boundary,
        Q_mean=np.mean(Q), Q2_mean=Q2_mean,
        chi_t_lattice=chi_t_lattice, chi_t_fourth_root_a=chi_t_fourth_root_a,
        tau_int=tau_int, dtau_int=dtau_int,
        alpha=run_data.alpha, n_conf=len(Q_raw)
    )


def lattice_spacing_su2(beta: float) -> float:
    """SU(2) lattice spacing from arXiv:1811.02800."""
    t0_phys = 0.01133  # fm^2
    delta_beta = beta - 2.600
    ln_t0_over_a2 = 1.285 + 6.409 * delta_beta - 0.7411 * delta_beta**2
    t0_over_a2 = np.exp(ln_t0_over_a2)
    return np.sqrt(t0_phys / t0_over_a2)


def lattice_spacing_su3(beta: float) -> float:
    """SU(3) lattice spacing from arXiv:hep-lat/0108008."""
    r0 = 0.5  # fm
    delta_beta = beta - 6.0
    ln_a_over_r0 = -1.6804 - 1.7331 * delta_beta + 0.7849 * delta_beta**2 - 0.4428 * delta_beta**3
    return np.exp(ln_a_over_r0) * r0
