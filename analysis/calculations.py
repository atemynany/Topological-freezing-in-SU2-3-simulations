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
import glob
import hashlib
from dataclasses import dataclass
from typing import Optional, Tuple

from autocorrelation import autocorrelation

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
CHI_T_FOURTH_ROOT_SU3 = 191.0


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


def load_topcharge_from_timeslices(filepath: str, T: int, exclude_boundary_slices: int,
                                   smear: int = None) -> np.ndarray:
    """Load Q by summing q(t) from topcharge_timeslice.dat.

    This is the authoritative path for temporal-subvolume measurements. It also
    repairs legacy periodic-subvolume runs where topcharge.dat contains the
    full-lattice Q but topcharge_timeslice.dat still has enough information to
    reconstruct the intended cropped Q.
    """
    rows = []
    smears = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                s = int(parts[0])
                c = int(parts[1])
                t = int(parts[2])
                q_t = float(parts[3])
            except ValueError:
                continue
            rows.append((s, c, t, q_t))
            smears.append(s)

    if not rows:
        return np.array([])
    if smear is None:
        smear = max(smears)

    t_min = max(0, exclude_boundary_slices)
    t_max = T - max(0, exclude_boundary_slices)
    by_conf = {}
    for s, c, t, q_t in rows:
        if s != smear:
            continue
        if t_min <= t < t_max:
            by_conf[c] = by_conf.get(c, 0.0) + q_t

    return np.array([by_conf[c] for c in sorted(by_conf)])


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
    start_type: str
    exclude_boundary_slices: int


@dataclass
class AnalysisResult:
    beta: float
    boundary: str
    run_dir: str
    T: int
    L: int
    Q_mean: float
    Q2_mean: float
    chi_t_lattice: float
    chi_t_lattice_err: float
    chi_t_fourth_root_a: float
    chi_t_fourth_root_a_err: float
    tau_int: float
    dtau_int: float
    alpha: float
    n_conf: int


def is_t160_l80(obj) -> bool:
    """Return True for the T160_L80 ensemble metadata used in highlight plots."""
    run_dir = getattr(obj, "run_dir", "")
    if "results_work" in os.path.normpath(run_dir).split(os.sep):
        return False
    return int(getattr(obj, "T", -1)) == 160 and int(getattr(obj, "L", -1)) == 80


def gaussian(x, mu, sigma, A):
    """Gaussian function for histogram fitting."""
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def parse_input_file(run_dir: str) -> dict:
    """Parse run metadata from local and CL2QCD run files."""
    info = {}
    metadata_files = [
        "input.txt",
        "run_info.txt",
        "input_su3_topcharge.txt",
        "input_su3heatbath.txt",
    ]

    for filename in metadata_files:
        input_file = os.path.join(run_dir, filename)
        if not os.path.exists(input_file):
            continue
        with open(input_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                if '=' in line:
                    key, value = line.split('=', 1)
                    value = value.strip().split()[0] if value.strip() else ""
                    info[key.strip()] = value
                else:
                    parts = line.split()
                    if len(parts) >= 2:
                        info[parts[0]] = parts[1]
    return info


def normalize_boundary(boundary: str) -> str:
    """Normalize boundary labels used by the local code and CL2QCD scripts."""
    return "open" if boundary in {"open", "obc"} else boundary


def parse_run_dir_metadata(run_dir: str) -> dict:
    """Parse run metadata encoded in directory names."""
    dirname = os.path.basename(run_dir)
    match = re.search(
        r'T(?P<T>\d+)_L(?P<L>\d+)_b(?P<beta>\d+(?:\.\d+)?)_(?P<boundary>open|obc|periodic)_seed(?P<seed>\d+)',
        dirname,
    )
    if match:
        info = match.groupdict()
        info["boundary"] = normalize_boundary(info["boundary"])
        return info

    parts = run_dir.split(os.sep)
    qtarget_dir = next((part for part in reversed(parts) if part.startswith("qtarget_")), None)
    if not qtarget_dir:
        return {}

    match = re.search(
        r'qtarget_T(?P<T>\d+)_L(?P<L>\d+)_b(?P<beta>\d+(?:\.\d+)?)_(?P<boundary>open|obc|periodic)_',
        qtarget_dir,
    )
    if not match:
        return {}

    info = match.groupdict()
    info["boundary"] = normalize_boundary(info["boundary"])
    seed_match = re.search(r'_seed(?P<seed>\d+)$', dirname)
    if seed_match:
        info["seed"] = seed_match.group("seed")
    return info


def topcharge_file_path(run_dir: str, gauge_group: str = "su2") -> Optional[str]:
    """Return the topcharge file for standard and qtarget run layouts."""
    candidates = [
        os.path.join(run_dir, "output", "topcharge.dat"),
        os.path.join(run_dir, "output", f"topcharge_{gauge_group}.dat"),
        os.path.join(run_dir, "topcharge", "topcharge.dat"),
        os.path.join(run_dir, "topcharge", f"topcharge_{gauge_group}.dat"),
    ]
    found = next((f for f in candidates if os.path.exists(f)), None)
    if found:
        return found
    found = qtarget_hot_list_file_path(run_dir, "topcharge.dat")
    if found:
        return found
    return qtarget_aggregate_file_path(run_dir, gauge_group, "topcharge.dat")


def topcharge_timeslice_file_path(run_dir: str, gauge_group: str = "su2") -> Optional[str]:
    """Return the topcharge-timeslice file for standard and qtarget layouts."""
    candidates = [
        os.path.join(run_dir, "output", "topcharge_timeslice.dat"),
        os.path.join(run_dir, "output", f"topcharge_timeslice_{gauge_group}.dat"),
        os.path.join(run_dir, "topcharge", "topcharge_timeslice.dat"),
        os.path.join(run_dir, "topcharge", f"topcharge_timeslice_{gauge_group}.dat"),
    ]
    found = next((f for f in candidates if os.path.exists(f)), None)
    if found:
        return found
    found = qtarget_hot_list_file_path(run_dir, "topcharge_timeslice.dat")
    if found:
        return found
    return qtarget_aggregate_file_path(run_dir, gauge_group, "topcharge_timeslice.dat")


def qtarget_hot_list_file_path(run_dir: str, filename: str) -> Optional[str]:
    """Return a qtarget hot-list output file, preferring it over hot candidates."""
    if not os.path.basename(run_dir).startswith("qtarget_"):
        return None
    pattern = os.path.join(run_dir, "hot_list_*", "output", filename)
    matches = sorted(glob.glob(pattern))
    return matches[-1] if matches else None


def qtarget_candidate_dirs(run_dir: str) -> list[str]:
    """Return qtarget hot-candidate directories, if this is a qtarget run."""
    if not os.path.basename(run_dir).startswith("qtarget_"):
        return []
    pattern = os.path.join(run_dir, "hot_candidates", "try_*_seed*")
    return sorted(path for path in glob.glob(pattern) if os.path.isdir(path))


def qtarget_aggregate_file_path(
    run_dir: str,
    gauge_group: str = "su2",
    filename: str = "topcharge.dat",
) -> Optional[str]:
    """Build a temporary aggregate .dat file for a qtarget run."""
    candidate_dirs = qtarget_candidate_dirs(run_dir)
    if not candidate_dirs:
        return None

    digest = hashlib.sha1(os.path.abspath(run_dir).encode("utf-8")).hexdigest()[:12]
    aggregate_dir = os.path.join("/tmp", "su23_qtarget_aggregates", digest)
    os.makedirs(aggregate_dir, exist_ok=True)
    aggregate_file = os.path.join(aggregate_dir, filename)

    rows_written = 0
    with open(aggregate_file, "w") as out:
        if filename == "topcharge_timeslice.dat":
            out.write("# smear_steps  config_number  t  q_t\n")
        else:
            out.write("# smear_steps  config_number  Q  plaquette\n")

        for conf_idx, candidate_dir in enumerate(candidate_dirs):
            if filename == "topcharge_timeslice.dat":
                source = topcharge_timeslice_file_path(candidate_dir, gauge_group)
            else:
                source = topcharge_file_path(candidate_dir, gauge_group)
            if source is None:
                continue

            with open(source) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.split()
                    if filename == "topcharge_timeslice.dat" and len(parts) >= 4:
                        out.write(f"{parts[0]} {conf_idx} {parts[2]} {parts[3]}\n")
                        rows_written += 1
                    elif filename == "topcharge.dat" and len(parts) >= 4:
                        out.write(f"{parts[0]} {conf_idx} {parts[2]} {parts[3]}\n")
                        rows_written += 1

    return aggregate_file if rows_written else None


def load_qtarget_topcharge(run_dir: str, gauge_group: str = "su2") -> np.ndarray:
    """Load one topcharge value from each qtarget hot candidate."""
    values = []
    for candidate_dir in qtarget_candidate_dirs(run_dir):
        topcharge_file = topcharge_file_path(candidate_dir, gauge_group)
        if topcharge_file is None:
            continue
        q_values = load_topcharge(topcharge_file)
        if len(q_values) > 0:
            values.append(q_values[-1])
    return np.array(values)


def load_run_data(run_dir: str, gauge_group: str = "su2") -> Optional[RunData]:
    """Load all data from a run directory."""
    info = parse_input_file(run_dir)
    dir_info = parse_run_dir_metadata(run_dir)

    if 'nTime' in info and 'T' not in info:
        info['T'] = info['nTime']
    if 'nSpace' in info and 'L' not in info:
        info['L'] = info['nSpace']
    if 'hostSeed' in info and 'seed' not in info:
        info['seed'] = info['hostSeed']
    if info.get('useOpenTemporalGaugeBoundary', '').lower() in {'true', '1'}:
        info['boundary'] = 'open'

    dirname = os.path.basename(run_dir)
    if '_open_' in dirname or '_obc_' in dirname:
        boundary = 'open'
    elif '_periodic_' in dirname:
        boundary = 'periodic'
    else:
        boundary = normalize_boundary(info.get('boundary', dir_info.get('boundary', 'periodic')))

    beta = float(info.get('beta', dir_info.get('beta', 0)))
    seed = int(info.get('seed', dir_info.get('seed', 0)))
    T = int(info.get('T', dir_info.get('T', 16)))
    L = int(info.get('L', dir_info.get('L', 16)))
    start_conf = int(info.get('start_conf', 50))
    start_type = info.get('start_type', dir_info.get('start_type', 'unknown'))
    exclude_boundary_slices = int(info.get('exclude_boundary_slices', 0))

    timeslice_file = topcharge_timeslice_file_path(run_dir, gauge_group)
    topcharge_file = topcharge_file_path(run_dir, gauge_group)
    if exclude_boundary_slices > 0 and timeslice_file:
        Q_raw = load_topcharge_from_timeslices(timeslice_file, T, exclude_boundary_slices)
    elif topcharge_file:
        Q_raw = load_topcharge(topcharge_file)
    else:
        Q_raw = load_qtarget_topcharge(run_dir, gauge_group)

    if len(Q_raw) == 0:
        return None

    continuous_q = (boundary == 'open') or (exclude_boundary_slices > 0)

    if continuous_q:
        # Open-boundary or temporal-subvolume charges are not globally
        # integer-quantized, so keep the alpha-rescaled continuous value.
        alpha = find_optimal_alpha(Q_raw)
        q_rescaled = alpha * Q_raw
    else:
        estimator = QEstimator(Q_raw)
        q_rescaled = estimator.Q_rescaled
        alpha = estimator.alpha
    
    # Load plaquette data
    plaq_candidates = [
        os.path.join(run_dir, "output", "plaquette.dat"),
        os.path.join(run_dir, "output", f"plaquette_{gauge_group}.dat"),
    ]
    plaq_file = next((f for f in plaq_candidates if os.path.exists(f)), None)
    if plaq_file is None:
        plaq_file = qtarget_hot_list_file_path(run_dir, "plaquette.dat")
    
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
    else:
        mc_list, plaq_list = [], []
        for idx, candidate_dir in enumerate(qtarget_candidate_dirs(run_dir)):
            candidate_plaq = os.path.join(candidate_dir, "output", "plaquette.dat")
            if not os.path.isfile(candidate_plaq):
                continue
            with open(candidate_plaq, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mc_list.append(len(mc_list))
                            plaq_list.append(float(parts[1]))
                        except ValueError:
                            continue
        mc_time = np.array(mc_list)
        plaq_data = np.array(plaq_list)
    
    return RunData(
        beta=beta, seed=seed, run_dir=run_dir, T=T, L=L,
        boundary=boundary, Q_raw=Q_raw, Q_rescaled=q_rescaled, alpha=alpha,
        plaq_data=plaq_data, mc_time=mc_time, start_conf=start_conf,
        start_type=start_type,
        exclude_boundary_slices=exclude_boundary_slices
    )


def analyze_run(run_data: RunData) -> AnalysisResult:
    """Perform full analysis on a single run."""
    volume = run_data.T * run_data.L ** 3
    Q = run_data.Q_rescaled.astype(float)
    Q2 = Q ** 2

    # ⟨Q²⟩ and its error via the Gamma method (pyerrors), so χ_t inherits a
    # proper autocorrelation-aware error bar.
    autocorr = autocorrelation(Q2, name=f"Q2_b{run_data.beta}")
    Q2_mean = autocorr.value          # == np.mean(Q2)
    Q2_err  = autocorr.error          # Gamma-method error on ⟨Q²⟩
    tau_int = autocorr.tau_int
    dtau_int = autocorr.dtau_int

    chi_t_lattice = Q2_mean / volume
    chi_t_lattice_err = Q2_err / volume

    if chi_t_lattice > 0:
        chi_t_fourth_root_a = (chi_t_lattice ** 0.25) * HBAR_C
        chi_t_fourth_root_a_err = 0.25 * chi_t_lattice ** (-0.75) * chi_t_lattice_err * HBAR_C
    else:
        chi_t_fourth_root_a = 0.0
        chi_t_fourth_root_a_err = 0.0

    return AnalysisResult(
        beta=run_data.beta, boundary=run_data.boundary, run_dir=run_data.run_dir,
        T=run_data.T, L=run_data.L,
        Q_mean=np.mean(Q), Q2_mean=Q2_mean,
        chi_t_lattice=chi_t_lattice, chi_t_lattice_err=chi_t_lattice_err,
        chi_t_fourth_root_a=chi_t_fourth_root_a,
        chi_t_fourth_root_a_err=chi_t_fourth_root_a_err,
        tau_int=tau_int, dtau_int=dtau_int,
        alpha=run_data.alpha, n_conf=len(Q)
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
