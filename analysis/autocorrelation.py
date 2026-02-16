#!/usr/bin/env python3
import numpy as np
import pyerrors as pe
from dataclasses import dataclass

@dataclass
class AutocorrResult:
    tau_int: float
    dtau_int: float
    thinning: int
    n_eff: float
    value: float
    error: float

def autocorrelation(data: np.ndarray, name: str = "obs", S: float = 1.5) -> AutocorrResult:
    obs = pe.Obs([np.asarray(data).flatten()], [name])
    obs.gamma_method(S=S)
    tau = obs.e_tauint[name]
    return AutocorrResult(
        tau_int=tau,
        dtau_int=obs.e_dtauint[name],
        thinning=max(1, int(np.ceil(2 * tau))),
        n_eff=len(data) / (2 * tau) if tau > 0 else len(data),
        value=obs.value,
        error=obs.dvalue
    )

def load_plaquette(filepath: str) -> np.ndarray:
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split()
            if len(parts) >= 2:
                try: data.append(float(parts[1]))
                except: pass
    return np.array(data)

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
        return np.array([q for s, q in zip(smears, Q_vals) if s == max_smear])
    idx = np.argsort(configs)
    return np.array(Q_vals)[idx]

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--plaq", help="Plaquette file")
    parser.add_argument("--topcharge", help="Topcharge file")
    parser.add_argument("--smear", type=int, help="Smearing step")
    parser.add_argument("--therm", type=int, default=0, help="Thermalization cut index")
    args = parser.parse_args()
    
    if args.plaq:
        data = load_plaquette(args.plaq)[args.therm:]
        r = autocorrelation(data, "plaq")
        print(f"Plaquette: tau_int={r.tau_int:.2f}±{r.dtau_int:.2f}, thin={r.thinning}, N_eff={r.n_eff:.1f}")
    
    if args.topcharge:
        data = load_topcharge(args.topcharge, args.smear)
        r = autocorrelation(data, "Q")
        print(f"Topcharge: tau_int={r.tau_int:.2f}±{r.dtau_int:.2f}, thin={r.thinning}, N_eff={r.n_eff:.1f}")
