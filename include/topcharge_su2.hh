#ifndef TOPCHARGE_SU2_HH
#define TOPCHARGE_SU2_HH

#include <cmath>
#include <cstdlib>

#include "linear_algebra.hh"
#include "geometry.hh"
#include "fields.hh"

// ---------------------------------------------------------------------------
// Plaquette staples — four orientations of the 1x1 plaquette in the (mu,nu)
// plane, used to build the clover-leaf field-strength tensor.
//
// Naming convention: sign of the mu/nu step taken before the first link.
//   pp = (+mu, +nu)   mp = (-mu, +nu)   mm = (-mu, -nu)   pm = (+mu, -nu)
//
// Each function returns the ordered product of four SU(2) links forming
// one corner of the clover leaf at site (it,ix,iy,iz).
// ---------------------------------------------------------------------------

// P_{++}: U_mu(x) U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x)
inline void compute_plaq_pp(double *result, const double *gauge_field, int it, int ix, int iy, int iz, int mu, int nu, int T, int L) {
    alignas(32) double M1[8], M2[8];
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

// P_{-+}: U_mu^dag(x-mu) U_nu(x-mu) U_mu(x-mu+nu) U_nu^dag(x)
inline void compute_plaq_mp(double *result, const double *gauge_field, int it, int ix, int iy, int iz, int mu, int nu, int T, int L) {
    alignas(32) double M1[8], M2[8];
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

// P_{--}: U_mu^dag(x-mu) U_nu^dag(x-mu-nu) U_mu(x-mu-nu) U_nu(x-nu)
inline void compute_plaq_mm(double *result, const double *gauge_field, int it, int ix, int iy, int iz, int mu, int nu, int T, int L) {
    alignas(32) double M1[8], M2[8];
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

// P_{+-}: U_mu(x) U_nu^dag(x+mu-nu) U_mu^dag(x-nu) U_nu(x-nu)
inline void compute_plaq_pm(double *result, const double *gauge_field, int it, int ix, int iy, int iz, int mu, int nu, int T, int L) {
    alignas(32) double M1[8], M2[8];
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

// ---------------------------------------------------------------------------
// Clover leaf: C_μν = (1/4)(P_{++} + P_{-+} + P_{--} + P_{+-})
// Symmetrised sum of all four corner plaquettes in the (mu,nu) plane.
// ---------------------------------------------------------------------------
inline void compute_clover(double *clover, const double *gauge_field, int it, int ix, int iy, int iz, int mu, int nu, int T, int L) {
    alignas(32) double P[8], sum[8];
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

// ---------------------------------------------------------------------------
// Field-strength tensor: F_μν = (C_μν - C_μν^†) / 2  (anti-Hermitian)
// Takes the anti-Hermitian part of the clover to extract the lattice
// approximation to the continuum F_μν.
// ---------------------------------------------------------------------------
inline void compute_field_strength(double *F_mu_nu, const double *gauge_field, int it, int ix, int iy, int iz, int mu, int nu, int T, int L) {
    alignas(32) double C[8], C_dag[8];

    compute_clover(C, gauge_field, it, ix, iy, iz, mu, nu, T, L);
    cm_eq_cm_dag(C_dag, C);

    cm_mi_eq_cm(C, C_dag);
    cm_eq_cm_ti_re(F_mu_nu, C, 0.5);
}

// Re Tr(F_{d1 d2} F_{d3 d4}) at a single site — helper used in the
// local density below.
inline double compute_clover_product(const double *gauge_field, int it, int ix, int iy, int iz, int d1, int d2, int d3, int d4, int T, int L) {
    alignas(32) double F1[8], F2[8], prod[8];

    compute_field_strength(F1, gauge_field, it, ix, iy, iz, d1, d2, T, L);
    compute_field_strength(F2, gauge_field, it, ix, iy, iz, d3, d4, T, L);

    cm_eq_cm_ti_cm(prod, F1, F2);

    complex tr;
    co_eq_tr_cm(&tr, prod);

    return tr.re;
}

// ---------------------------------------------------------------------------
// Local topological charge density at site (it,ix,iy,iz):
//   q(x) = -(1/8π²) ε_μνρσ Tr(F_μν F_ρσ)
//         = -Tr(F_01 F_32) - Tr(F_12 F_30) - Tr(F_20 F_31)
// The 1/(4π²) normalisation is applied in the calling functions.
// ---------------------------------------------------------------------------
inline double compute_local_topcharge_density(const double *gauge_field, int it, int ix, int iy, int iz, int T, int L) {
    alignas(32) double F01[8], F32[8];
    alignas(32) double F12[8], F30[8];
    alignas(32) double F20[8], F31[8];
    alignas(32) double prod[8];

    compute_field_strength(F01, gauge_field, it, ix, iy, iz, 0, 1, T, L);
    compute_field_strength(F32, gauge_field, it, ix, iy, iz, 3, 2, T, L);
    compute_field_strength(F12, gauge_field, it, ix, iy, iz, 1, 2, T, L);
    compute_field_strength(F30, gauge_field, it, ix, iy, iz, 3, 0, T, L);
    compute_field_strength(F20, gauge_field, it, ix, iy, iz, 2, 0, T, L);
    compute_field_strength(F31, gauge_field, it, ix, iy, iz, 3, 1, T, L);

    double q_local = 0.0;
    complex tr;

    cm_eq_cm_ti_cm(prod, F01, F32);
    co_eq_tr_cm(&tr, prod);
    q_local -= tr.re;

    cm_eq_cm_ti_cm(prod, F12, F30);
    co_eq_tr_cm(&tr, prod);
    q_local -= tr.re;

    cm_eq_cm_ti_cm(prod, F20, F31);
    co_eq_tr_cm(&tr, prod);
    q_local -= tr.re;

    return q_local;
}

// ---------------------------------------------------------------------------
// Total topological charge: Q = (1/4π²) Σ_x q(x), parallelised over sites.
// ---------------------------------------------------------------------------
inline double compute_topological_charge(double * __restrict__ gauge_field, int T, int L) {
    double Q = 0.0;
    const int total_sites = T * L * L * L;

    #pragma omp parallel for reduction(+:Q) schedule(static)
    for (int site = 0; site < total_sites; site++) {
        const int iz = site % L;
        const int iy = (site / L) % L;
        const int ix = (site / (L * L)) % L;
        const int it = site / (L * L * L);

        Q += compute_local_topcharge_density(gauge_field, it, ix, iy, iz, T, L);
    }

    return Q / (4.0 * M_PI * M_PI);
}

// Topological charge with boundary exclusion for open BC:
// sums only over t ∈ [exclude_boundary, T - exclude_boundary).
inline double compute_topological_charge_open(double * __restrict__ gauge_field, int T, int L, int exclude_boundary) {
    double Q = 0.0;

    const int t_min = exclude_boundary;
    const int t_max = T - exclude_boundary;

    if (t_min >= t_max) {
        return compute_topological_charge(gauge_field, T, L);
    }

    #pragma omp parallel for reduction(+:Q) schedule(static)
    for (int it = t_min; it < t_max; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    Q += compute_local_topcharge_density(gauge_field, it, ix, iy, iz, T, L);
                }
            }
        }
    }

    return Q / (4.0 * M_PI * M_PI);
}

#endif // TOPCHARGE_SU2_HH
