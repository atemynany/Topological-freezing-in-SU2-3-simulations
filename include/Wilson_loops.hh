// ==============================================================================
// Wilson_loops.hh
// ==============================================================================
// Wilson loop and plaquette calculations for SU(2) lattice gauge theory.
// Provides average plaquette measurement for thermalization monitoring.
//
// Author: Alexander de Barros Noll
// Date: January 2026
// ==============================================================================

#ifndef WILSON_LOOPS_HH
#define WILSON_LOOPS_HH

#include "linear_algebra.hh"
#include "geometry.hh"

// ==============================================================================
// External declarations for boundary conditions
// ==============================================================================

extern bool open_boundary_conditions;

// ==============================================================================
// Plaquette Calculation
// ==============================================================================

/**
 * @brief Computes a single plaquette at position (it, ix, iy, iz) in the mu-nu plane.
 * 
 * The plaquette is defined as:
 *   P_{mu,nu}(x) = U_mu(x) * U_nu(x+mu) * U_mu^dag(x+nu) * U_nu^dag(x)
 * 
 * @param plaq Output: 8-component SU(2) matrix representing the plaquette
 * @param gauge_field Pointer to the gauge field array
 * @param it Temporal coordinate
 * @param ix Spatial x-coordinate
 * @param iy Spatial y-coordinate
 * @param iz Spatial z-coordinate
 * @param mu First direction index (0=t, 1=x, 2=y, 3=z)
 * @param nu Second direction index
 * @param T Temporal lattice extent
 * @param L Spatial lattice extent
 */
inline void compute_plaquette(double *plaq, double *gauge_field,
                               int it, int ix, int iy, int iz,
                               int mu, int nu, int T, int L) {
    double M1[8], M2[8];
    
    int idx[4] = {it, ix, iy, iz};
    
    // Get link indices
    int idx_U_mu = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    idx[mu] += 1;
    int idx_U_nu_shifted = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    idx[mu] -= 1;
    idx[nu] += 1;
    int idx_U_mu_shifted = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    idx[nu] -= 1;
    int idx_U_nu = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    // Handle boundary conditions
    if (idx_U_mu < 0 || idx_U_nu_shifted < 0 || idx_U_mu_shifted < 0 || idx_U_nu < 0) {
        cm_eq_id(plaq);
        return;
    }
    
    // P = U_mu(x) * U_nu(x+mu) * U_mu^dag(x+nu) * U_nu^dag(x)
    cm_eq_cm_ti_cm(M1, gauge_field + idx_U_mu, gauge_field + idx_U_nu_shifted);
    cm_eq_cm_ti_cm_dag(M2, M1, gauge_field + idx_U_mu_shifted);
    cm_eq_cm_ti_cm_dag(plaq, M2, gauge_field + idx_U_nu);
}

/**
 * @brief Computes the real part of the trace of a plaquette (normalized).
 * 
 * For SU(2): returns (1/2) * Re(Tr(P))
 * This value is 1 for the identity (cold start) and ~0 for random (hot start).
 * 
 * @param gauge_field Pointer to the gauge field array
 * @param it Temporal coordinate
 * @param ix Spatial x-coordinate
 * @param iy Spatial y-coordinate
 * @param iz Spatial z-coordinate
 * @param mu First direction
 * @param nu Second direction
 * @param T Temporal lattice extent
 * @param L Spatial lattice extent
 * @return Normalized plaquette value in [0, 1]
 */
inline double plaquette_trace(double *gauge_field,
                               int it, int ix, int iy, int iz,
                               int mu, int nu, int T, int L) {
    double plaq[8];
    compute_plaquette(plaq, gauge_field, it, ix, iy, iz, mu, nu, T, L);
    
    complex tr;
    co_eq_tr_cm(&tr, plaq);
    
    // Normalize by N=2 for SU(2)
    return 0.5 * tr.re;
}

// ==============================================================================
// Average Plaquette
// ==============================================================================

/**
 * @brief Computes the spatially and temporally averaged plaquette.
 * 
 * The average plaquette is an order parameter for the confinement-deconfinement
 * transition and serves as a measure of thermalization during Monte Carlo updates.
 * 
 * For SU(2): <P> = 1 for cold (ordered) configurations
 *            <P> ~ 0 for hot (random) configurations
 *            <P> approaches beta-dependent value at equilibrium
 * 
 * @param gauge_field Pointer to the gauge field array
 * @param T Temporal lattice extent
 * @param L Spatial lattice extent
 * @return Average plaquette value
 */
inline double Average_Plaquette(double *gauge_field, int T, int L) {
    double sum = 0.0;
    int count = 0;
    
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    // Sum over all 6 plaquette orientations
                    for (int mu = 0; mu < 4; mu++) {
                        for (int nu = mu + 1; nu < 4; nu++) {
                            // Check boundary conditions for open BC
                            if (open_boundary_conditions) {
                                // Skip temporal plaquettes at boundaries
                                if (mu == 0 && (it == T - 1)) continue;
                                if (nu == 0 && (it == T - 1)) continue;
                            }
                            
                            sum += plaquette_trace(gauge_field, it, ix, iy, iz, mu, nu, T, L);
                            count++;
                        }
                    }
                }
            }
        }
    }
    
    return (count > 0) ? sum / count : 0.0;
}

/**
 * @brief Computes only the spatial average plaquette (excludes temporal plaquettes).
 * 
 * Useful for analyzing spatial correlations and comparing with temporal behavior.
 * 
 * @param gauge_field Pointer to the gauge field array
 * @param T Temporal lattice extent
 * @param L Spatial lattice extent
 * @return Spatial average plaquette value
 */
inline double Average_Spatial_Plaquette(double *gauge_field, int T, int L) {
    double sum = 0.0;
    int count = 0;
    
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    // Sum over spatial plaquettes only (mu, nu > 0)
                    for (int mu = 1; mu < 4; mu++) {
                        for (int nu = mu + 1; nu < 4; nu++) {
                            sum += plaquette_trace(gauge_field, it, ix, iy, iz, mu, nu, T, L);
                            count++;
                        }
                    }
                }
            }
        }
    }
    
    return (count > 0) ? sum / count : 0.0;
}

/**
 * @brief Computes the temporal average plaquette (only time-space plaquettes).
 * 
 * @param gauge_field Pointer to the gauge field array
 * @param T Temporal lattice extent
 * @param L Spatial lattice extent
 * @return Temporal average plaquette value
 */
inline double Average_Temporal_Plaquette(double *gauge_field, int T, int L) {
    double sum = 0.0;
    int count = 0;
    
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    // Sum over temporal-spatial plaquettes (mu=0, nu=1,2,3)
                    for (int nu = 1; nu < 4; nu++) {
                        if (open_boundary_conditions && it == T - 1) continue;
                        
                        sum += plaquette_trace(gauge_field, it, ix, iy, iz, 0, nu, T, L);
                        count++;
                    }
                }
            }
        }
    }
    
    return (count > 0) ? sum / count : 0.0;
}

#endif // WILSON_LOOPS_HH
