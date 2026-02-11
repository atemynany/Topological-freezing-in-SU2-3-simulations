#ifndef PLAQUETTE_HH
#define PLAQUETTE_HH

#include "linear_algebra.hh"
#include "geometry.hh"

extern bool open_boundary_conditions;

inline void compute_plaquette(double *plaq, double *gauge_field,
                               int it, int ix, int iy, int iz,
                               int mu, int nu, int T, int L) {
    double M1[8], M2[8];
    
    int idx[4] = {it, ix, iy, iz};
    
    int idx_U_mu = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    idx[mu] += 1;
    int idx_U_nu_shifted = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    idx[mu] -= 1;
    idx[nu] += 1;
    int idx_U_mu_shifted = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    idx[nu] -= 1;
    int idx_U_nu = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    if (idx_U_mu < 0 || idx_U_nu_shifted < 0 || idx_U_mu_shifted < 0 || idx_U_nu < 0) {
        cm_eq_id(plaq);
        return;
    }
    
    cm_eq_cm_ti_cm(M1, gauge_field + idx_U_mu, gauge_field + idx_U_nu_shifted);
    cm_eq_cm_ti_cm_dag(M2, M1, gauge_field + idx_U_mu_shifted);
    cm_eq_cm_ti_cm_dag(plaq, M2, gauge_field + idx_U_nu);
}

inline double plaquette_trace(double *gauge_field,
                               int it, int ix, int iy, int iz,
                               int mu, int nu, int T, int L) {
    double plaq[8];
    compute_plaquette(plaq, gauge_field, it, ix, iy, iz, mu, nu, T, L);
    
    complex tr;
    co_eq_tr_cm(&tr, plaq);
    
    return 0.5 * tr.re;
}

inline double Average_Plaquette(double *gauge_field, int T, int L) {
    double sum = 0.0;
    int count = 0;
    
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    for (int mu = 0; mu < 4; mu++) {
                        for (int nu = mu + 1; nu < 4; nu++) {
                            if (open_boundary_conditions) {
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

inline double Average_Spatial_Plaquette(double *gauge_field, int T, int L) {
    double sum = 0.0;
    int count = 0;
    
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
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

inline double Average_Temporal_Plaquette(double *gauge_field, int T, int L) {
    double sum = 0.0;
    int count = 0;
    
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
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

#endif
