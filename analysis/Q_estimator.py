"""
Q Estimator: Rescaled topological charge.

    Q_rescaled = round(α * Q̂)

where α minimizes ⟨(α * Q̂ - Q_rescaled)²⟩
"""

import numpy as np
from dataclasses import dataclass


def find_optimal_alpha(Q_raw: np.ndarray, alpha_min=0.8, alpha_max=1.2, tol=1e-6) -> float:
    """Find α that minimizes ⟨(αQ̂ - round(αQ̂))²⟩."""
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
