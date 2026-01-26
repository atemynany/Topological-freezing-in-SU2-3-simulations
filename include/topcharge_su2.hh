// ==============================================================================
// topcharge_su2.hh
// ==============================================================================
// Topological charge measurement for SU(2) lattice gauge theory using the
// clover definition of the field strength tensor.
//
// The topological charge is computed using:
//   Q = (1 / 32 * pi^2) * sum_x epsilon_{mu nu rho sigma} Tr(F_{mu nu} F_{rho sigma})
//
// where F_{mu nu} is approximated by the clover (4-leaf) plaquette average.
//
// Author: Alexander de Barros Noll
// Original implementation: Carolin Riehl (2020)
// Date: January 2026
// ==============================================================================

#ifndef TOPCHARGE_SU2_HH
#define TOPCHARGE_SU2_HH

#include <cmath>
#include <cstdlib>

#include "linear_algebra.hh"
#include "geometry.hh"
#include "fields.hh"

// ==============================================================================
// Clover (Field Strength) Computation
// ==============================================================================

/**
 * @brief Computes the plaquette U_mu(x) U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x)
 *        in the positive-positive direction.
 * 
 * @param result Output SU(2) matrix
 * @param gauge_field Gauge field array
 * @param idx Lattice site index
 * @param mu First direction
 * @param nu Second direction
 * @param T Temporal extent
 * @param L Spatial extent
 */
inline void compute_plaq_pp(double *result, double *gauge_field,
                            int it, int ix, int iy, int iz,
                            int mu, int nu, int T, int L) {
    double M1[8], M2[8];
    int idx[4] = {it, ix, iy, iz};
    
    // U_mu(x)
    int i1 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    // U_nu(x + mu)
    idx[mu] += 1;
    int i2 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    // U_mu^dag(x + nu)
    idx[mu] -= 1;
    idx[nu] += 1;
    int i3 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    // U_nu^dag(x)
    idx[nu] -= 1;
    int i4 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0) {
        cm_eq_zero(result);
        return;
    }
    
    cm_eq_cm_ti_cm(M1, gauge_field + i1, gauge_field + i2);
    cm_eq_cm_ti_cm_dag(M2, M1, gauge_field + i3);
    cm_eq_cm_ti_cm_dag(result, M2, gauge_field + i4);
}

/**
 * @brief Computes plaquette in negative mu, positive nu direction.
 */
inline void compute_plaq_mp(double *result, double *gauge_field,
                            int it, int ix, int iy, int iz,
                            int mu, int nu, int T, int L) {
    double M1[8], M2[8];
    int idx[4] = {it, ix, iy, iz};
    
    // U_mu^dag(x - mu)
    idx[mu] -= 1;
    int i1 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    // U_nu(x - mu)
    int i2 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    // U_mu(x - mu + nu)
    idx[nu] += 1;
    int i3 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    // U_nu^dag(x)
    idx[mu] += 1;
    idx[nu] -= 1;
    int i4 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0) {
        cm_eq_zero(result);
        return;
    }
    
    cm_eq_cm_dag_ti_cm(M1, gauge_field + i1, gauge_field + i2);
    cm_eq_cm_ti_cm(M2, M1, gauge_field + i3);
    cm_eq_cm_ti_cm_dag(result, M2, gauge_field + i4);
}

/**
 * @brief Computes plaquette in negative mu, negative nu direction.
 */
inline void compute_plaq_mm(double *result, double *gauge_field,
                            int it, int ix, int iy, int iz,
                            int mu, int nu, int T, int L) {
    double M1[8], M2[8];
    int idx[4] = {it, ix, iy, iz};
    
    // U_mu^dag(x - mu)
    idx[mu] -= 1;
    int i1 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    // U_nu^dag(x - mu - nu)
    idx[nu] -= 1;
    int i2 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    // U_mu(x - mu - nu)
    int i3 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    // U_nu(x - nu)
    idx[mu] += 1;
    int i4 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0) {
        cm_eq_zero(result);
        return;
    }
    
    cm_eq_cm_dag_ti_cm_dag(M1, gauge_field + i1, gauge_field + i2);
    cm_eq_cm_ti_cm(M2, M1, gauge_field + i3);
    cm_eq_cm_ti_cm(result, M2, gauge_field + i4);
}

/**
 * @brief Computes plaquette in positive mu, negative nu direction.
 */
inline void compute_plaq_pm(double *result, double *gauge_field,
                            int it, int ix, int iy, int iz,
                            int mu, int nu, int T, int L) {
    double M1[8], M2[8];
    int idx[4] = {it, ix, iy, iz};
    
    // U_mu(x)
    int i1 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    // U_nu^dag(x + mu - nu)
    idx[mu] += 1;
    idx[nu] -= 1;
    int i2 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    // U_mu^dag(x - nu)
    idx[mu] -= 1;
    int i3 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    // U_nu(x - nu)
    int i4 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0) {
        cm_eq_zero(result);
        return;
    }
    
    cm_eq_cm_ti_cm_dag(M1, gauge_field + i1, gauge_field + i2);
    cm_eq_cm_ti_cm_dag(M2, M1, gauge_field + i3);
    cm_eq_cm_ti_cm(result, M2, gauge_field + i4);
}

/**
 * @brief Computes the clover (4-leaf) sum for the field strength tensor.
 * 
 * C_{mu nu}(x) = P_{mu,nu}(x) + P_{-mu,nu}(x) + P_{-mu,-nu}(x) + P_{mu,-nu}(x)
 * 
 * The imaginary (anti-Hermitian) part gives the lattice field strength tensor:
 *   F_{mu nu} = (C - C^dag) / (8i)
 * 
 * @param clover Output: clover sum matrix
 * @param gauge_field Gauge field array
 * @param it, ix, iy, iz Lattice site coordinates
 * @param mu First direction
 * @param nu Second direction
 * @param T Temporal extent
 * @param L Spatial extent
 */
inline void compute_clover(double *clover, double *gauge_field,
                           int it, int ix, int iy, int iz,
                           int mu, int nu, int T, int L) {
    double P[8], sum[8];
    cm_eq_zero(sum);
    
    // Sum of four plaquettes around the point
    compute_plaq_pp(P, gauge_field, it, ix, iy, iz, mu, nu, T, L);
    cm_pl_eq_cm(sum, P);
    
    compute_plaq_mp(P, gauge_field, it, ix, iy, iz, nu, mu, T, L);
    cm_pl_eq_cm(sum, P);
    
    compute_plaq_mm(P, gauge_field, it, ix, iy, iz, mu, nu, T, L);
    cm_pl_eq_cm(sum, P);
    
    compute_plaq_pm(P, gauge_field, it, ix, iy, iz, nu, mu, T, L);
    cm_pl_eq_cm(sum, P);
    
    // Normalize
    cm_eq_cm_ti_re(clover, sum, 0.25);
}

/**
 * @brief Computes the imaginary part of the clover (field strength approximation).
 * 
 * @param F_mu_nu Output: imaginary part of clover = (C - C^dag) / 2
 * @param gauge_field Gauge field array
 * @param it, ix, iy, iz Lattice site coordinates
 * @param mu First direction
 * @param nu Second direction
 * @param T Temporal extent
 * @param L Spatial extent
 */
inline void compute_field_strength(double *F_mu_nu, double *gauge_field,
                                    int it, int ix, int iy, int iz,
                                    int mu, int nu, int T, int L) {
    double C[8], C_dag[8];
    
    compute_clover(C, gauge_field, it, ix, iy, iz, mu, nu, T, L);
    cm_eq_cm_dag(C_dag, C);
    
    // F = (C - C^dag) / 2 (gives the anti-Hermitian part)
    cm_mi_eq_cm(C, C_dag);
    cm_eq_cm_ti_re(F_mu_nu, C, 0.5);
}

// ==============================================================================
// Topological Charge Density and Total Charge
// ==============================================================================

/**
 * @brief Computes the product Tr(F_{d1,d2} * F_{d3,d4}) at a given site.
 * 
 * This is the core computation for the topological charge density.
 * 
 * @param gauge_field Gauge field array
 * @param it, ix, iy, iz Lattice site coordinates
 * @param d1, d2 Directions for first field strength
 * @param d3, d4 Directions for second field strength
 * @param T Temporal extent
 * @param L Spatial extent
 * @return Tr(Im(C_{d1,d2}) * Im(C_{d3,d4}))
 */
inline double compute_clover_product(double *gauge_field,
                                     int it, int ix, int iy, int iz,
                                     int d1, int d2, int d3, int d4,
                                     int T, int L) {
    double F1[8], F2[8], prod[8];
    
    compute_field_strength(F1, gauge_field, it, ix, iy, iz, d1, d2, T, L);
    compute_field_strength(F2, gauge_field, it, ix, iy, iz, d3, d4, T, L);
    
    cm_eq_cm_ti_cm(prod, F1, F2);
    
    complex tr;
    co_eq_tr_cm(&tr, prod);
    
    return tr.re;
}

/**
 * @brief Computes the total topological charge Q of the gauge configuration.
 * 
 * Q = (1 / 4*pi^2) * sum_x epsilon_{mu nu rho sigma} Tr(F_{mu nu} F_{rho sigma})
 * 
 * Using the reduced formula (faster):
 *   Q = -(1 / 4*pi^2) * sum_x [ Tr(F_01 * F_32) + Tr(F_12 * F_30) + Tr(F_20 * F_31) ]
 * 
 * For smooth configurations, Q should be close to an integer value.
 * 
 * @param gauge_field Gauge field array
 * @param T Temporal extent
 * @param L Spatial extent
 * @return Topological charge Q
 */
inline double compute_topological_charge(double *gauge_field, int T, int L) {
    double Q = 0.0;
    
    // Using the reduced formula (equivalent to full epsilon contraction)
    // Q = -sum_x [ Tr(F_01 * F_32) + Tr(F_12 * F_30) + Tr(F_20 * F_31) ]
    
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    Q -= compute_clover_product(gauge_field, it, ix, iy, iz, 0, 1, 3, 2, T, L);
                    Q -= compute_clover_product(gauge_field, it, ix, iy, iz, 1, 2, 3, 0, T, L);
                    Q -= compute_clover_product(gauge_field, it, ix, iy, iz, 2, 0, 3, 1, T, L);
                }
            }
        }
    }
    
    // Normalize by 4*pi^2
    return Q / (4.0 * M_PI * M_PI);
}

/**
 * @brief Computes the topological charge density at each lattice site.
 * 
 * Useful for visualizing instanton positions and studying local fluctuations.
 * 
 * @param q_density Output array of size T*L*L*L containing charge density at each site
 * @param gauge_field Gauge field array
 * @param T Temporal extent
 * @param L Spatial extent
 */
inline void compute_topological_charge_density(double *q_density, double *gauge_field, int T, int L) {
    const double norm = 1.0 / (4.0 * M_PI * M_PI);
    
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    int idx = ((it * L + ix) * L + iy) * L + iz;
                    
                    double q_local = 0.0;
                    q_local -= compute_clover_product(gauge_field, it, ix, iy, iz, 0, 1, 3, 2, T, L);
                    q_local -= compute_clover_product(gauge_field, it, ix, iy, iz, 1, 2, 3, 0, T, L);
                    q_local -= compute_clover_product(gauge_field, it, ix, iy, iz, 2, 0, 3, 1, T, L);
                    
                    q_density[idx] = q_local * norm;
                }
            }
        }
    }
}

#endif // TOPCHARGE_SU2_HH
