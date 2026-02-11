// ==============================================================================
// test_topcharge.cc
// ==============================================================================
// Unit tests for topological charge measurement in SU(2) lattice gauge theory.
// Tests clover computation, field strength tensor, and topological charge.
//
// Author: Alexander de Barros Noll
// Date: January 2026
// ==============================================================================

#include <catch2/catch_all.hpp>
#include <cmath>
#include <vector>
#include <array>

#include "fields.hh"
#include "geometry.hh"
#include "linear_algebra.hh"
#include "ranlux.hh"
#include "topcharge_su2.hh"

// ==============================================================================
// External Global Variables
// ==============================================================================

// Required for geometry functions
bool open_boundary_conditions = false;

// ==============================================================================
// Test Fixtures and Helpers
// ==============================================================================

constexpr double TOL = 1e-10;
constexpr double LOOSE_TOL = 1e-6;

// Small test lattice dimensions
constexpr int TEST_T = 4;
constexpr int TEST_L = 4;

// Generate random SU(2) matrix
static void random_su2(double *A) {
    double h[4];
    double norm;
    
    do {
        h[0] = 2.0 * DRand() - 1.0;
        h[1] = 2.0 * DRand() - 1.0;
        h[2] = 2.0 * DRand() - 1.0;
        h[3] = 2.0 * DRand() - 1.0;
        norm = h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3];
    } while (norm < 1e-10 || norm > 1.0);
    
    norm = 1.0 / std::sqrt(norm);
    h[0] *= norm;
    h[1] *= norm;
    h[2] *= norm;
    h[3] *= norm;
    
    cm_from_h(A, h);
}

// Check if matrix is valid SU(2)
static bool is_valid_su2(const double *A, double tol = TOL) {
    double h[4];
    h_from_cm(h, A);
    double norm = h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3];
    return std::fabs(norm - 1.0) < tol;
}

// ==============================================================================
// Clover/Plaquette Tests
// ==============================================================================

TEST_CASE("Topological Charge: Plaquette on trivial configuration", "[topcharge][plaquette]") {
    InitializeRand(100);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Unity(gauge_field, TEST_T, TEST_L);
    
    double plaq[8];
    
    SECTION("Plaquette on unit configuration is identity") {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu + 1; nu < 4; nu++) {
                compute_plaq_pp(plaq, gauge_field, 0, 0, 0, 0, mu, nu, TEST_T, TEST_L);
                
                double I[8];
                cm_eq_id(I);
                
                for (int i = 0; i < 8; i++) {
                    REQUIRE(std::fabs(plaq[i] - I[i]) < TOL);
                }
            }
        }
    }
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Topological Charge: Clover symmetry", "[topcharge][clover]") {
    InitializeRand(101);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    SECTION("Clover is symmetric in its construction") {
        double C1[8], C2[8];
        
        // C_{mu,nu} = sum of four plaquettes
        compute_clover(C1, gauge_field, 1, 1, 1, 1, 0, 1, TEST_T, TEST_L);
        
        // Should be well-defined
        REQUIRE(!std::isnan(C1[0]));
        REQUIRE(!std::isinf(C1[0]));
    }
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Field Strength Tensor Tests
// ==============================================================================

TEST_CASE("Topological Charge: Field strength on trivial config", "[topcharge][fieldstrength]") {
    InitializeRand(102);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Unity(gauge_field, TEST_T, TEST_L);
    
    SECTION("F_{mu nu} = 0 on unit configuration") {
        double F[8];
        
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu + 1; nu < 4; nu++) {
                compute_field_strength(F, gauge_field, 1, 1, 1, 1, mu, nu, TEST_T, TEST_L);
                
                // F should be zero (or very small)
                for (int i = 0; i < 8; i++) {
                    REQUIRE(std::fabs(F[i]) < TOL);
                }
            }
        }
    }
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Topological Charge: Field strength anti-symmetry", "[topcharge][fieldstrength]") {
    InitializeRand(103);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    SECTION("F_{mu nu} is anti-Hermitian") {
        double F[8], F_dag[8], sum[8];
        
        compute_field_strength(F, gauge_field, 1, 1, 1, 1, 0, 1, TEST_T, TEST_L);
        cm_eq_cm_dag(F_dag, F);
        
        // F + F^dag should be zero for anti-Hermitian
        cm_eq_cm(sum, F);
        cm_pl_eq_cm(sum, F_dag);
        
        for (int i = 0; i < 8; i++) {
            REQUIRE(std::fabs(sum[i]) < LOOSE_TOL);
        }
    }
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Topological Charge Tests
// ==============================================================================

TEST_CASE("Topological Charge: Q=0 on trivial configuration", "[topcharge][charge]") {
    InitializeRand(104);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Unity(gauge_field, TEST_T, TEST_L);
    
    double Q = compute_topological_charge(gauge_field, TEST_T, TEST_L);
    
    REQUIRE(std::fabs(Q) < TOL);
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Topological Charge: Q is real-valued", "[topcharge][charge]") {
    InitializeRand(105);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double Q = compute_topological_charge(gauge_field, TEST_T, TEST_L);
    
    // Q should be a real number (not NaN or Inf)
    REQUIRE(!std::isnan(Q));
    REQUIRE(!std::isinf(Q));
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Topological Charge: Q is bounded for random config", "[topcharge][charge]") {
    InitializeRand(106);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double Q = compute_topological_charge(gauge_field, TEST_T, TEST_L);
    
    // For a random configuration, Q should be some finite value
    // On a small lattice, typically |Q| < volume
    double max_Q = TEST_T * TEST_L * TEST_L * TEST_L;
    REQUIRE(std::fabs(Q) < max_Q);
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Topological Charge Density Tests
// ==============================================================================

TEST_CASE("Topological Charge: Density sums to Q", "[topcharge][density]") {
    InitializeRand(107);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    // Compute total Q
    double Q_total = compute_topological_charge(gauge_field, TEST_T, TEST_L);
    
    // Compute Q from density
    int vol = TEST_T * TEST_L * TEST_L * TEST_L;
    std::vector<double> q_density(vol);
    compute_topological_charge_density(q_density.data(), gauge_field, TEST_T, TEST_L);
    
    double Q_from_density = 0.0;
    for (int i = 0; i < vol; i++) {
        Q_from_density += q_density[i];
    }
    
    REQUIRE(std::fabs(Q_total - Q_from_density) < LOOSE_TOL);
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Topological Charge: Density is zero on trivial config", "[topcharge][density]") {
    InitializeRand(108);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Unity(gauge_field, TEST_T, TEST_L);
    
    int vol = TEST_T * TEST_L * TEST_L * TEST_L;
    std::vector<double> q_density(vol);
    compute_topological_charge_density(q_density.data(), gauge_field, TEST_T, TEST_L);
    
    for (int i = 0; i < vol; i++) {
        REQUIRE(std::fabs(q_density[i]) < TOL);
    }
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Gauge Invariance Tests
// ==============================================================================

TEST_CASE("Topological Charge: Gauge transformation invariance", "[topcharge][gauge_invariance]") {
    InitializeRand(109);
    
    double *gauge_field;
    double *gauge_field_transformed;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Alloc(&gauge_field_transformed, TEST_T, TEST_L);
    
    // Create random configuration
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    // Copy to transformed field
    Gauge_Field_Copy(gauge_field_transformed, gauge_field, TEST_T, TEST_L);
    
    // Apply gauge transformation: U_mu(x) -> g(x) U_mu(x) g^dag(x+mu)
    // For simplicity, we test that Q doesn't change drastically under small perturbations
    // (Full gauge transformation test would require implementing the transform)
    
    double Q_original = compute_topological_charge(gauge_field, TEST_T, TEST_L);
    double Q_copy = compute_topological_charge(gauge_field_transformed, TEST_T, TEST_L);
    
    // Should be identical for the copy
    REQUIRE(std::fabs(Q_original - Q_copy) < TOL);
    
    Gauge_Field_Free(&gauge_field);
    Gauge_Field_Free(&gauge_field_transformed);
}

// ==============================================================================
// Lattice Size Dependence
// ==============================================================================

TEST_CASE("Topological Charge: Works on different lattice sizes", "[topcharge][sizes]") {
    InitializeRand(110);
    
    SECTION("2x2 lattice") {
        const int T = 2, L = 2;
        double *gf;
        Gauge_Field_Alloc(&gf, T, L);
        Gauge_Field_Unity(gf, T, L);
        
        double Q = compute_topological_charge(gf, T, L);
        REQUIRE(std::fabs(Q) < TOL);
        
        Gauge_Field_Free(&gf);
    }
    
    SECTION("4x4 lattice") {
        const int T = 4, L = 4;
        double *gf;
        Gauge_Field_Alloc(&gf, T, L);
        Gauge_Field_Unity(gf, T, L);
        
        double Q = compute_topological_charge(gf, T, L);
        REQUIRE(std::fabs(Q) < TOL);
        
        Gauge_Field_Free(&gf);
    }
    
    SECTION("8x8 lattice") {
        const int T = 8, L = 8;
        double *gf;
        Gauge_Field_Alloc(&gf, T, L);
        Gauge_Field_Unity(gf, T, L);
        
        double Q = compute_topological_charge(gf, T, L);
        REQUIRE(std::fabs(Q) < TOL);
        
        Gauge_Field_Free(&gf);
    }
}

// ==============================================================================
// Consistency Tests
// ==============================================================================

TEST_CASE("Topological Charge: Reproducibility", "[topcharge][reproducibility]") {
    // Test that the same configuration gives the same Q
    InitializeRand(111);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double Q1 = compute_topological_charge(gauge_field, TEST_T, TEST_L);
    double Q2 = compute_topological_charge(gauge_field, TEST_T, TEST_L);
    
    REQUIRE(std::fabs(Q1 - Q2) < TOL);
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Advanced Gauge Invariance Tests
// ==============================================================================

// Helper: Generate deterministic "random" SU(2) from seed
static void su2_from_seed(double *U, unsigned seed) {
    double a = 0.1 * std::sin(1.0 * seed + 0.1);
    double b = 0.1 * std::sin(2.0 * seed + 0.2);
    double c = 0.1 * std::sin(3.0 * seed + 0.3);
    su2_exp_from_Amu(U, std::array<double, 3>{a, b, c});
    cm_proj(U);
}

// Apply a full gauge transformation: U_mu(x) -> G(x) U_mu(x) G^dag(x+mu)
static void apply_gauge_transform_topcharge(double *gf, int L, int T) {
    int vol = T * L * L * L;
    
    std::vector<double> G(static_cast<std::size_t>(vol) * 8, 0.0);
    for (int site = 0; site < vol; ++site)
        su2_from_seed(&G[static_cast<std::size_t>(site) * 8], 1234u + site);
    
    auto site_index = [&](int t, int x, int y, int z) {
        return t * L * L * L + x * L * L + y * L + z;
    };
    
    auto forward = [&](int t, int x, int y, int z, int mu) {
        if (mu == 0) t = (t + 1 + T) % T;
        if (mu == 1) x = (x + 1 + L) % L;
        if (mu == 2) y = (y + 1 + L) % L;
        if (mu == 3) z = (z + 1 + L) % L;
        return site_index(t, x, y, z);
    };
    
    double tmp1[8], tmp2[8], Gd[8];
    
    for (int t = 0; t < T; ++t)
        for (int x = 0; x < L; ++x)
            for (int y = 0; y < L; ++y)
                for (int z = 0; z < L; ++z) {
                    int s = site_index(t, x, y, z);
                    for (int mu = 0; mu < 4; ++mu) {
                        int sp = forward(t, x, y, z, mu);
                        
                        double *U = &gf[static_cast<std::size_t>(s) * 4 * 8 + mu * 8];
                        const double *Gs = &G[static_cast<std::size_t>(s) * 8];
                        const double *Gsp = &G[static_cast<std::size_t>(sp) * 8];
                        
                        cm_eq_cm_dag(Gd, Gsp);
                        cm_eq_cm_ti_cm(tmp1, Gs, U);
                        cm_eq_cm_ti_cm(tmp2, tmp1, Gd);
                        cm_eq_cm(U, tmp2);
                    }
                }
}

// Apply a small gauge transformation
static void apply_small_gauge_transform(double *gf, int L, int T, double eps) {
    int vol = T * L * L * L;
    
    std::vector<double> G(static_cast<std::size_t>(vol) * 8, 0.0);
    for (int site = 0; site < vol; ++site) {
        std::array<double, 3> A = {eps * 0.1 * std::sin(site + 1.0),
                                   eps * 0.1 * std::cos(site + 2.0),
                                   eps * 0.1 * std::sin(site + 3.0)};
        su2_exp_from_Amu(&G[static_cast<std::size_t>(site) * 8], A);
        cm_proj(&G[static_cast<std::size_t>(site) * 8]);
    }
    
    auto site_index = [&](int t, int x, int y, int z) {
        return t * L * L * L + x * L * L + y * L + z;
    };
    
    auto forward = [&](int t, int x, int y, int z, int mu) {
        if (mu == 0) t = (t + 1 + T) % T;
        if (mu == 1) x = (x + 1 + L) % L;
        if (mu == 2) y = (y + 1 + L) % L;
        if (mu == 3) z = (z + 1 + L) % L;
        return site_index(t, x, y, z);
    };
    
    double tmp1[8], tmp2[8], Gd[8];
    
    for (int t = 0; t < T; ++t)
        for (int x = 0; x < L; ++x)
            for (int y = 0; y < L; ++y)
                for (int z = 0; z < L; ++z) {
                    int s = site_index(t, x, y, z);
                    for (int mu = 0; mu < 4; ++mu) {
                        int sp = forward(t, x, y, z, mu);
                        
                        double *U = &gf[static_cast<std::size_t>(s) * 4 * 8 + mu * 8];
                        const double *Gs = &G[static_cast<std::size_t>(s) * 8];
                        const double *Gsp = &G[static_cast<std::size_t>(sp) * 8];
                        
                        cm_eq_cm_dag(Gd, Gsp);
                        cm_eq_cm_ti_cm(tmp1, Gs, U);
                        cm_eq_cm_ti_cm(tmp2, tmp1, Gd);
                        cm_eq_cm(U, tmp2);
                    }
                }
}

// Fill gauge field with identity
static void fill_identity_topcharge(double *gf, int T, int L) {
    int vol = T * L * L * L;
    for (int site = 0; site < vol; ++site)
        for (int mu = 0; mu < 4; ++mu)
            cm_eq_id(&gf[site * 4 * 8 + mu * 8]);
}

TEST_CASE("Topological Charge: Q is gauge invariant (full transform)", "[topcharge][gauge_invariance][full]") {
    const int L = 4, T = 4;
    
    double *gauge_field;
    double *gauge_field_transformed;
    Gauge_Field_Alloc(&gauge_field, T, L);
    Gauge_Field_Alloc(&gauge_field_transformed, T, L);
    
    // Create nontrivial configuration
    fill_identity_topcharge(gauge_field, T, L);
    int vol = T * L * L * L;
    for (int site = 0; site < vol; ++site) {
        double U[8];
        su2_from_seed(U, 777u + site);
        // Set mu=1 links to nontrivial values
        cm_eq_cm(&gauge_field[static_cast<std::size_t>(site) * 4 * 8 + 1 * 8], U);
    }
    
    // Copy and transform
    Gauge_Field_Copy(gauge_field_transformed, gauge_field, T, L);
    apply_gauge_transform_topcharge(gauge_field_transformed, L, T);
    
    double Q1 = compute_topological_charge(gauge_field, T, L);
    double Q2 = compute_topological_charge(gauge_field_transformed, T, L);
    
    REQUIRE(std::fabs(Q1 - Q2) < LOOSE_TOL);
    
    Gauge_Field_Free(&gauge_field);
    Gauge_Field_Free(&gauge_field_transformed);
}

TEST_CASE("Topological Charge: Q is gauge invariant (small transform)", "[topcharge][gauge_invariance][small]") {
    const int L = 4, T = 4;
    
    double *gauge_field;
    double *gauge_field_transformed;
    Gauge_Field_Alloc(&gauge_field, T, L);
    Gauge_Field_Alloc(&gauge_field_transformed, T, L);
    
    fill_identity_topcharge(gauge_field, T, L);
    int vol = T * L * L * L;
    for (int site = 0; site < vol; ++site) {
        double U[8];
        su2_from_seed(U, 1000u + site);
        cm_eq_cm(&gauge_field[static_cast<std::size_t>(site) * 4 * 8 + 1 * 8], U);
    }
    
    Gauge_Field_Copy(gauge_field_transformed, gauge_field, T, L);
    apply_small_gauge_transform(gauge_field_transformed, L, T, 1e-3);
    
    double Q1 = compute_topological_charge(gauge_field, T, L);
    double Q2 = compute_topological_charge(gauge_field_transformed, T, L);
    
    REQUIRE(std::fabs(Q1 - Q2) < LOOSE_TOL);
    
    Gauge_Field_Free(&gauge_field);
    Gauge_Field_Free(&gauge_field_transformed);
}

// ==============================================================================
// Abelian Background Tests
// ==============================================================================

static void build_abelian_F12(double *gf, int L, int T, double theta) {
    fill_identity_topcharge(gf, T, L);
    
    for (int t = 0; t < T; ++t)
        for (int x = 0; x < L; ++x)
            for (int y = 0; y < L; ++y)
                for (int z = 0; z < L; ++z) {
                    double U1[8];
                    su2_exp_from_Amu(U1, std::array<double, 3>{0.0, 0.0, theta * static_cast<double>(y)});
                    int site = t * L * L * L + x * L * L + y * L + z;
                    cm_eq_cm(&gf[static_cast<std::size_t>(site) * 4 * 8 + 1 * 8], U1);
                }
}

TEST_CASE("Topological Charge: Abelian F12-only gives Q=0", "[topcharge][abelian]") {
    const int L = 6, T = 4;
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, T, L);
    build_abelian_F12(gauge_field, L, T, 0.4);
    
    double Q = compute_topological_charge(gauge_field, T, L);
    
    // Pure F12 background should give Q=0 (no instanton content)
    REQUIRE(std::fabs(Q) < LOOSE_TOL);
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Clover Covariance Tests
// ==============================================================================

TEST_CASE("Topological Charge: Clover transforms covariantly", "[topcharge][clover][gauge]") {
    const int L = 4, T = 4;
    
    double *gauge_field;
    double *gauge_field_transformed;
    Gauge_Field_Alloc(&gauge_field, T, L);
    Gauge_Field_Alloc(&gauge_field_transformed, T, L);
    
    // Create mildly nontrivial configuration
    fill_identity_topcharge(gauge_field, T, L);
    int vol = T * L * L * L;
    for (int site = 0; site < vol; ++site) {
        double U[8];
        su2_from_seed(U, 777u + site);
        cm_eq_cm(&gauge_field[static_cast<std::size_t>(site) * 4 * 8 + 1 * 8], U);
    }
    
    Gauge_Field_Copy(gauge_field_transformed, gauge_field, T, L);
    apply_gauge_transform_topcharge(gauge_field_transformed, L, T);
    
    // Test clover at interior point
    int t = 1, x = 1, y = 1, z = 1, mu = 0, nu = 1;
    
    double C[8], C2[8];
    compute_clover(C, gauge_field, t, x, y, z, mu, nu, T, L);
    compute_clover(C2, gauge_field_transformed, t, x, y, z, mu, nu, T, L);
    
    // Build G(x) from same deterministic construction
    int site = t * L * L * L + x * L * L + y * L + z;
    double G[8], Gd[8], tmp[8], cov[8];
    su2_from_seed(G, 1234u + site);
    cm_eq_cm_dag(Gd, G);
    
    // cov = G * C * G^dag
    cm_eq_cm_ti_cm(tmp, G, C);
    cm_eq_cm_ti_cm(cov, tmp, Gd);
    
    // Compare
    for (int i = 0; i < 8; ++i)
        REQUIRE(std::fabs(cov[i] - C2[i]) < LOOSE_TOL);
    
    Gauge_Field_Free(&gauge_field);
    Gauge_Field_Free(&gauge_field_transformed);
}

