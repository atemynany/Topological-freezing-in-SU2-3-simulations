#!/usr/bin/env python3
import numpy as np
import pyerrors as pe
from dataclasses import dataclass
from Q_estimator import QEstimator

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
    """Compute jackknife error estimate (ignores autocorrelation)."""
    n = len(data)
    theta_jack = np.array([np.mean(np.delete(data, i)) for i in range(n)])
    return np.sqrt((n - 1) / n * np.sum((theta_jack - theta_jack.mean())**2))

def compute_Q_statistics(Q_raw: np.ndarray, ensemble: str = "ens", S: float = 1.5) -> QStats:
    est = QEstimator(Q_raw)
    Q = est.Q_rescaled.astype(float)
    
    # Gamma method (handles autocorrelation)
    q_obs = pe.Obs([Q], [ensemble])
    q_obs.gamma_method(S=S)
    
    q2_obs = pe.Obs([Q**2], [ensemble])
    q2_obs.gamma_method(S=S)
    
    tau = q_obs.e_tauint[ensemble]
    
    # Jackknife errors (for comparison)
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
    configs, Q_vals, smears = [], [], []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    s, c, q = int(parts[0]), int(parts[1]), float(parts[2])
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

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Topcharge file")
    parser.add_argument("--smear", type=int, help="Smearing step (default: max)")
    args = parser.parse_args()
    
    Q_raw = load_topcharge(args.file, args.smear)
    stats = compute_Q_statistics(Q_raw)
    
    print(f"α = {stats.alpha:.4f}")
    print(f"<Q>  = {stats.Q_mean:.4f} ± {stats.Q_error:.4f} (Gamma) ± {stats.Q_error_jack:.4f} (JK)")
    print(f"<Q²> = {stats.Q2_mean:.4f} ± {stats.Q2_error:.4f} (Gamma) ± {stats.Q2_error_jack:.4f} (JK)")
    print(f"τ_int = {stats.tau_int:.2f}, N_eff = {stats.n_eff:.1f}")
