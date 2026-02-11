"""
Topological Susceptibility for SU(2) Lattice Gauge Theory

    χ_t = ⟨Q²⟩ / V

Reference for SU(2): χ_t^(1/4) = 200(15) MeV
Conversion: ℏc ≈ 197.3 MeV·fm
"""

import numpy as np

HBAR_C = 197.3  # MeV * fm


def topological_susceptibility(Q2_mean: float, T: int, L: int, lattice_spacing: float) -> dict:

    volume = T * L ** 3
    
    # χ_t in lattice units
    chi_t_lattice = Q2_mean / volume
    
    # χ_t in fm⁻⁴
    chi_t_fm4 = chi_t_lattice / (lattice_spacing ** 4)
    
    # χ_t^(1/4) in MeV
    chi_t_fourth_root_MeV = (chi_t_fm4 ** 0.25) * HBAR_C
    
    return {
        'chi_t_lattice': chi_t_lattice,
        'chi_t_fm4': chi_t_fm4,
        'chi_t_fourth_root_MeV': chi_t_fourth_root_MeV,
        'reference_SU2_MeV': 200.0
    }
