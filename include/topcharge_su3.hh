// topcharge_su3.hh - SU(3) topological charge using clover discretization
//
// Q = (1/32π²) Σ_x ε_μνρσ Tr(F_μν(x) F_ρσ(x))
// F_μν from clover (average of 4 plaquettes)
//
#ifndef TOPCHARGE_SU3_HH
#define TOPCHARGE_SU3_HH

#include <cmath>
#include <cstdlib>
#include "su3_linear_algebra.hh"

// Geometry: site index and link index for SU(3)
// These must be set by the calling code
extern int T_size, L_size;
extern std::vector<int> neighbor_plus[4];
extern std::vector<int> neighbor_minus[4];
extern std::vector<int> link_index_su3;

inline int su3_get_site(int t, int x, int y, int z) {
    int tt = (t + T_size) % T_size;
    int xx = (x + L_size) % L_size;
    int yy = (y + L_size) % L_size;
    int zz = (z + L_size) % L_size;
    return ((tt * L_size + xx) * L_size + yy) * L_size + zz;
}

inline int su3_link_idx(int site, int mu) {
    return (4 * site + mu) * 18;
}

// P_{++}: U_mu(x) U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x)
inline void su3_plaq_pp(double * __restrict__ result,
                        const double * __restrict__ gf,
                        int it, int ix, int iy, int iz,
                        int mu, int nu) {
    alignas(32) double M1[18], M2[18];
    int idx[4] = {it, ix, iy, iz};
    int site = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    
    int i1 = su3_link_idx(site, mu);
    
    idx[mu] += 1;
    int site_mu = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    int i2 = su3_link_idx(site_mu, nu);
    
    idx[mu] -= 1; idx[nu] += 1;
    int site_nu = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    int i3 = su3_link_idx(site_nu, mu);
    
    idx[nu] -= 1;
    int i4 = su3_link_idx(site, nu);
    
    su3_eq_su3_ti_su3(M1, gf + i1, gf + i2);
    su3_eq_su3_ti_su3_dag(M2, M1, gf + i3);
    su3_eq_su3_ti_su3_dag(result, M2, gf + i4);
}

// P_{-+}: U_mu^dag(x-mu) U_nu(x-mu) U_mu(x-mu+nu) U_nu^dag(x)
inline void su3_plaq_mp(double * __restrict__ result,
                        const double * __restrict__ gf,
                        int it, int ix, int iy, int iz,
                        int mu, int nu) {
    alignas(32) double M1[18], M2[18];
    int idx[4] = {it, ix, iy, iz};
    
    idx[mu] -= 1;
    int site_mmu = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    int i1 = su3_link_idx(site_mmu, mu);
    int i2 = su3_link_idx(site_mmu, nu);
    
    idx[nu] += 1;
    int site_mmu_pnu = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    int i3 = su3_link_idx(site_mmu_pnu, mu);
    
    idx[mu] += 1; idx[nu] -= 1;
    int site0 = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    int i4 = su3_link_idx(site0, nu);
    
    su3_eq_su3_dag_ti_su3(M1, gf + i1, gf + i2);
    su3_eq_su3_ti_su3(M2, M1, gf + i3);
    su3_eq_su3_ti_su3_dag(result, M2, gf + i4);
}

// P_{--}: U_mu^dag(x-mu) U_nu^dag(x-mu-nu) U_mu(x-mu-nu) U_nu(x-nu)
inline void su3_plaq_mm(double * __restrict__ result,
                        const double * __restrict__ gf,
                        int it, int ix, int iy, int iz,
                        int mu, int nu) {
    alignas(32) double M1[18], M2[18];
    int idx[4] = {it, ix, iy, iz};
    
    idx[mu] -= 1;
    int site_mmu = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    int i1 = su3_link_idx(site_mmu, mu);
    
    idx[nu] -= 1;
    int site_mmu_mnu = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    int i2 = su3_link_idx(site_mmu_mnu, nu);
    int i3 = su3_link_idx(site_mmu_mnu, mu);
    
    idx[mu] += 1;
    int site_mnu = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    int i4 = su3_link_idx(site_mnu, nu);
    
    su3_eq_su3_dag_ti_su3_dag(M1, gf + i1, gf + i2);
    su3_eq_su3_ti_su3(M2, M1, gf + i3);
    su3_eq_su3_ti_su3(result, M2, gf + i4);
}

// P_{+-}: U_mu(x) U_nu^dag(x+mu-nu) U_mu^dag(x-nu) U_nu(x-nu)
inline void su3_plaq_pm(double * __restrict__ result,
                        const double * __restrict__ gf,
                        int it, int ix, int iy, int iz,
                        int mu, int nu) {
    alignas(32) double M1[18], M2[18];
    int idx[4] = {it, ix, iy, iz};
    int site0 = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    int i1 = su3_link_idx(site0, mu);
    
    idx[mu] += 1; idx[nu] -= 1;
    int site_pmu_mnu = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    int i2 = su3_link_idx(site_pmu_mnu, nu);
    
    idx[mu] -= 1;
    int site_mnu = su3_get_site(idx[0], idx[1], idx[2], idx[3]);
    int i3 = su3_link_idx(site_mnu, mu);
    int i4 = su3_link_idx(site_mnu, nu);
    
    su3_eq_su3_ti_su3_dag(M1, gf + i1, gf + i2);
    su3_eq_su3_ti_su3_dag(M2, M1, gf + i3);
    su3_eq_su3_ti_su3(result, M2, gf + i4);
}

// Clover: C_μν = (1/4)(P_++ + P_-+ + P_-- + P_+-)
inline void su3_clover(double * __restrict__ clover,
                       const double * __restrict__ gf,
                       int it, int ix, int iy, int iz,
                       int mu, int nu) {
    alignas(32) double P[18], sum[18];
    su3_eq_zero(sum);
    
    su3_plaq_pp(P, gf, it, ix, iy, iz, mu, nu);
    su3_pl_eq_su3(sum, P);
    
    su3_plaq_mp(P, gf, it, ix, iy, iz, nu, mu);
    su3_pl_eq_su3(sum, P);
    
    su3_plaq_mm(P, gf, it, ix, iy, iz, mu, nu);
    su3_pl_eq_su3(sum, P);
    
    su3_plaq_pm(P, gf, it, ix, iy, iz, nu, mu);
    su3_pl_eq_su3(sum, P);
    
    su3_ti_eq_re(sum, 0.25);
    su3_eq_su3(clover, sum);
}

// Field strength: F_μν = (C - C†) / 2  (anti-hermitian)
inline void su3_field_strength(double * __restrict__ F,
                               const double * __restrict__ gf,
                               int it, int ix, int iy, int iz,
                               int mu, int nu) {
    alignas(32) double C[18], C_dag[18];
    
    su3_clover(C, gf, it, ix, iy, iz, mu, nu);
    su3_eq_su3_dag(C_dag, C);
    
    // F = (C - C†) / 2
    for (int i = 0; i < 18; i++) {
        F[i] = 0.5 * (C[i] - C_dag[i]);
    }
}

// Local topological charge density (un-normalized)
// q = -ε_μνρσ Tr(F_μν F_ρσ) = -2[Tr(F_01 F_23) - Tr(F_02 F_13) + Tr(F_03 F_12)]
inline double su3_local_topcharge_density(const double * __restrict__ gf,
                                          int it, int ix, int iy, int iz) {
    alignas(32) double F01[18], F23[18];
    alignas(32) double F02[18], F13[18];
    alignas(32) double F03[18], F12[18];
    alignas(32) double prod[18];
    
    su3_field_strength(F01, gf, it, ix, iy, iz, 0, 1);
    su3_field_strength(F23, gf, it, ix, iy, iz, 2, 3);
    su3_field_strength(F02, gf, it, ix, iy, iz, 0, 2);
    su3_field_strength(F13, gf, it, ix, iy, iz, 1, 3);
    su3_field_strength(F03, gf, it, ix, iy, iz, 0, 3);
    su3_field_strength(F12, gf, it, ix, iy, iz, 1, 2);
    
    double q = 0.0;
    
    // Tr(F_01 F_23)
    su3_eq_su3_ti_su3(prod, F01, F23);
    q -= su3_re_tr(prod);
    
    // -Tr(F_02 F_13)
    su3_eq_su3_ti_su3(prod, F02, F13);
    q += su3_re_tr(prod);
    
    // Tr(F_03 F_12)
    su3_eq_su3_ti_su3(prod, F03, F12);
    q -= su3_re_tr(prod);
    
    return 2.0 * q;  // factor of 2 from ε-tensor
}

// Total topological charge: Q = (1/4π²) Σ_x q(x)
inline double su3_topological_charge(const double * __restrict__ gf, int T, int L) {
    double Q = 0.0;
    const int volume = T * L * L * L;
    
    #pragma omp parallel for reduction(+:Q) schedule(static)
    for (int site = 0; site < volume; site++) {
        const int iz = site % L;
        const int iy = (site / L) % L;
        const int ix = (site / (L * L)) % L;
        const int it = site / (L * L * L);
        
        Q += su3_local_topcharge_density(gf, it, ix, iy, iz);
    }
    
    // Normalization: 1/(4π²)
    return Q / (4.0 * M_PI * M_PI);
}

#endif // TOPCHARGE_SU3_HH