#ifndef TOPCHARGE_SU2_HH
#define TOPCHARGE_SU2_HH

#include <cmath>
#include <cstdlib>

#include "linear_algebra.hh"
#include "geometry.hh"
#include "fields.hh"

inline void compute_plaq_pp(double *result, double *gauge_field,
                            int it, int ix, int iy, int iz,
                            int mu, int nu, int T, int L) {
    double M1[8], M2[8];
    int idx[4] = {it, ix, iy, iz};
    
    int i1 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    idx[mu] += 1;
    int i2 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    idx[mu] -= 1;
    idx[nu] += 1;
    int i3 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
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

inline void compute_plaq_mp(double *result, double *gauge_field,
                            int it, int ix, int iy, int iz,
                            int mu, int nu, int T, int L) {
    double M1[8], M2[8];
    int idx[4] = {it, ix, iy, iz};
    
    idx[mu] -= 1;
    int i1 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    int i2 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    idx[nu] += 1;
    int i3 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
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

inline void compute_plaq_mm(double *result, double *gauge_field,
                            int it, int ix, int iy, int iz,
                            int mu, int nu, int T, int L) {
    double M1[8], M2[8];
    int idx[4] = {it, ix, iy, iz};
    
    idx[mu] -= 1;
    int i1 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    idx[nu] -= 1;
    int i2 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    int i3 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
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

inline void compute_plaq_pm(double *result, double *gauge_field,
                            int it, int ix, int iy, int iz,
                            int mu, int nu, int T, int L) {
    double M1[8], M2[8];
    int idx[4] = {it, ix, iy, iz};
    
    int i1 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    idx[mu] += 1;
    idx[nu] -= 1;
    int i2 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    idx[mu] -= 1;
    int i3 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), mu);
    
    int i4 = ggi(get_index(idx[0], idx[1], idx[2], idx[3], T, L), nu);
    
    if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0) {
        cm_eq_zero(result);
        return;
    }
    
    cm_eq_cm_ti_cm_dag(M1, gauge_field + i1, gauge_field + i2);
    cm_eq_cm_ti_cm_dag(M2, M1, gauge_field + i3);
    cm_eq_cm_ti_cm(result, M2, gauge_field + i4);
}

inline void compute_clover(double *clover, double *gauge_field,
                           int it, int ix, int iy, int iz,
                           int mu, int nu, int T, int L) {
    double P[8], sum[8];
    cm_eq_zero(sum);
    
    compute_plaq_pp(P, gauge_field, it, ix, iy, iz, mu, nu, T, L);
    cm_pl_eq_cm(sum, P);
    
    compute_plaq_mp(P, gauge_field, it, ix, iy, iz, nu, mu, T, L);
    cm_pl_eq_cm(sum, P);
    
    compute_plaq_mm(P, gauge_field, it, ix, iy, iz, mu, nu, T, L);
    cm_pl_eq_cm(sum, P);
    
    compute_plaq_pm(P, gauge_field, it, ix, iy, iz, nu, mu, T, L);
    cm_pl_eq_cm(sum, P);
    
    cm_eq_cm_ti_re(clover, sum, 0.25);
}

inline void compute_field_strength(double *F_mu_nu, double *gauge_field,
                                    int it, int ix, int iy, int iz,
                                    int mu, int nu, int T, int L) {
    double C[8], C_dag[8];
    
    compute_clover(C, gauge_field, it, ix, iy, iz, mu, nu, T, L);
    cm_eq_cm_dag(C_dag, C);
    
    cm_mi_eq_cm(C, C_dag);
    cm_eq_cm_ti_re(F_mu_nu, C, 0.5);
}

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

inline double compute_topological_charge(double *gauge_field, int T, int L) {
    double Q = 0.0;
    
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
    
    return Q / (4.0 * M_PI * M_PI);
}

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

#endif
