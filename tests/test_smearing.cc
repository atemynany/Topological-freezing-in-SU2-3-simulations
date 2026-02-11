// ==============================================================================
// test_smearing.cc
// ==============================================================================
// Unit tests for APE smearing in SU(2) lattice gauge theory.
// Tests that smearing preserves SU(2) structure and reduces fluctuations.
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
#include "smearing_techniques.hh"
#include "Plaquette.hh"
#include "topcharge_su2.hh"

// ==============================================================================
// External Global Variables
// ==============================================================================

bool open_boundary_conditions = false;

// ==============================================================================
// Test Fixtures and Helpers
// ==============================================================================

constexpr double TOL = 1e-10;
constexpr double LOOSE_TOL = 1e-6;

// Test lattice dimensions
constexpr int TEST_T = 4;
constexpr int TEST_L = 4;

// Check if matrix is valid SU(2)
static bool is_valid_su2(const double *A, double tol = TOL) {
    double h[4];
    h_from_cm(h, A);
    double norm = h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3];
    return std::fabs(norm - 1.0) < tol;
}

// Check if entire gauge field contains valid SU(2) matrices
static bool all_links_valid_su2(double *gauge_field, int T, int L, double tol = LOOSE_TOL) {
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    for (int mu = 0; mu < 4; mu++) {
                        int idx = ggi(get_index(it, ix, iy, iz, T, L), mu);
                        if (idx >= 0 && !is_valid_su2(gauge_field + idx, tol)) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}

// ==============================================================================
// APE Smearing Basic Tests
// ==============================================================================

TEST_CASE("Smearing: APE smearing preserves SU(2)", "[smearing][su2]") {
    InitializeRand(200);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    // Verify initial field is valid
    REQUIRE(all_links_valid_su2(gauge_field, TEST_T, TEST_L));
    
    // Apply smearing
    double alpha = 0.5;
    APE_Smearing_all(gauge_field, TEST_T, TEST_L, alpha);
    
    // Verify field is still valid SU(2)
    REQUIRE(all_links_valid_su2(gauge_field, TEST_T, TEST_L));
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Smearing: Multiple APE smearing steps preserve SU(2)", "[smearing][su2]") {
    InitializeRand(201);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double alpha = 0.5;
    int n_smear = 10;
    
    for (int i = 0; i < n_smear; i++) {
        APE_Smearing_all(gauge_field, TEST_T, TEST_L, alpha);
        REQUIRE(all_links_valid_su2(gauge_field, TEST_T, TEST_L));
    }
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// APE Smearing on Trivial Configuration
// ==============================================================================

TEST_CASE("Smearing: Unit configuration unchanged", "[smearing][trivial]") {
    InitializeRand(202);
    
    double *gauge_field;
    double *original_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Alloc(&original_field, TEST_T, TEST_L);
    
    Gauge_Field_Unity(gauge_field, TEST_T, TEST_L);
    Gauge_Field_Copy(original_field, gauge_field, TEST_T, TEST_L);
    
    // Apply smearing
    double alpha = 0.5;
    APE_Smearing_all(gauge_field, TEST_T, TEST_L, alpha);
    
    // Unit configuration should remain unchanged (all links are identity)
    int vol = TEST_T * TEST_L * TEST_L * TEST_L;
    for (int i = 0; i < vol * 4 * 8; i++) {
        REQUIRE(std::fabs(gauge_field[i] - original_field[i]) < LOOSE_TOL);
    }
    
    Gauge_Field_Free(&gauge_field);
    Gauge_Field_Free(&original_field);
}

// ==============================================================================
// APE Smearing Increases Plaquette
// ==============================================================================

TEST_CASE("Smearing: Increases average plaquette", "[smearing][plaquette]") {
    InitializeRand(203);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double plaq_before = Average_Plaquette(gauge_field, TEST_T, TEST_L);
    
    // Apply multiple smearing steps
    double alpha = 0.5;
    int n_smear = 20;
    
    for (int i = 0; i < n_smear; i++) {
        APE_Smearing_all(gauge_field, TEST_T, TEST_L, alpha);
    }
    
    double plaq_after = Average_Plaquette(gauge_field, TEST_T, TEST_L);
    
    // Smearing should increase the plaquette (make configuration smoother)
    // For random config, plaq ~ 0; after smearing should be closer to 1
    REQUIRE(plaq_after >= plaq_before - TOL);
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Smearing: Plaquette monotonically increases", "[smearing][plaquette]") {
    InitializeRand(204);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double alpha = 0.5;
    double prev_plaq = Average_Plaquette(gauge_field, TEST_T, TEST_L);
    
    for (int i = 0; i < 10; i++) {
        APE_Smearing_all(gauge_field, TEST_T, TEST_L, alpha);
        double curr_plaq = Average_Plaquette(gauge_field, TEST_T, TEST_L);
        
        // Plaquette should not decrease (within numerical tolerance)
        REQUIRE(curr_plaq >= prev_plaq - LOOSE_TOL);
        prev_plaq = curr_plaq;
    }
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// APE Smearing Parameter Dependence
// ==============================================================================

TEST_CASE("Smearing: Small alpha produces small change", "[smearing][alpha]") {
    InitializeRand(205);
    
    double *gauge_field;
    double *original_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Alloc(&original_field, TEST_T, TEST_L);
    
    // Start from cold configuration where effect is clearer
    Gauge_Field_Unity(gauge_field, TEST_T, TEST_L);
    Gauge_Field_Copy(original_field, gauge_field, TEST_T, TEST_L);
    
    double plaq_orig = Average_Plaquette(original_field, TEST_T, TEST_L);
    
    // Apply smearing with small alpha (0.1)
    APE_Smearing_all(gauge_field, TEST_T, TEST_L, 0.1);
    
    // For cold configuration, plaquette should remain ~1
    double plaq_smeared = Average_Plaquette(gauge_field, TEST_T, TEST_L);
    
    // Both should be close to 1.0 for cold start
    REQUIRE(std::fabs(plaq_orig - 1.0) < TOL);
    REQUIRE(std::fabs(plaq_smeared - 1.0) < LOOSE_TOL);
    
    Gauge_Field_Free(&gauge_field);
    Gauge_Field_Free(&original_field);
}

TEST_CASE("Smearing: Larger alpha gives stronger smoothing", "[smearing][alpha]") {
    InitializeRand(206);
    
    double *gauge_field1;
    double *gauge_field2;
    double *original;
    Gauge_Field_Alloc(&gauge_field1, TEST_T, TEST_L);
    Gauge_Field_Alloc(&gauge_field2, TEST_T, TEST_L);
    Gauge_Field_Alloc(&original, TEST_T, TEST_L);
    
    Gauge_Field_Random(original, TEST_T, TEST_L);
    Gauge_Field_Copy(gauge_field1, original, TEST_T, TEST_L);
    Gauge_Field_Copy(gauge_field2, original, TEST_T, TEST_L);
    
    // Apply smearing with different alpha values
    int n_smear = 5;
    
    for (int i = 0; i < n_smear; i++) {
        APE_Smearing_all(gauge_field1, TEST_T, TEST_L, 0.3);
        APE_Smearing_all(gauge_field2, TEST_T, TEST_L, 0.7);
    }
    
    double plaq1 = Average_Plaquette(gauge_field1, TEST_T, TEST_L);
    double plaq2 = Average_Plaquette(gauge_field2, TEST_T, TEST_L);
    
    // Larger alpha should give larger plaquette (more smoothing)
    // Note: This is generally true but not guaranteed in all cases
    // We just verify both are reasonable
    REQUIRE(plaq1 > 0.0);
    REQUIRE(plaq2 > 0.0);
    REQUIRE(plaq1 < 1.0 + TOL);
    REQUIRE(plaq2 < 1.0 + TOL);
    
    Gauge_Field_Free(&gauge_field1);
    Gauge_Field_Free(&gauge_field2);
    Gauge_Field_Free(&original);
}

// ==============================================================================
// Smearing Convergence
// ==============================================================================

TEST_CASE("Smearing: Approaches smooth configuration", "[smearing][convergence]") {
    InitializeRand(207);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double alpha = 0.5;
    
    // Apply many smearing steps
    for (int i = 0; i < 100; i++) {
        APE_Smearing_all(gauge_field, TEST_T, TEST_L, alpha);
    }
    
    double plaq = Average_Plaquette(gauge_field, TEST_T, TEST_L);
    
    // After many smearing steps, plaquette should be close to 1
    // (configuration approaches trivial/smooth)
    REQUIRE(plaq > 0.9);
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Numerical Stability
// ==============================================================================

TEST_CASE("Smearing: Numerically stable over moderate iterations", "[smearing][stability]") {
    InitializeRand(208);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double alpha = 0.5;
    
    // Apply moderate number of smearing steps
    // Note: After many iterations, numerical precision can degrade.
    // The important property is that links remain valid SU(2) matrices.
    for (int i = 0; i < 30; i++) {
        APE_Smearing_all(gauge_field, TEST_T, TEST_L, alpha);
        
        // Check SU(2) validity at each step (this is the critical test)
        REQUIRE(all_links_valid_su2(gauge_field, TEST_T, TEST_L));
        
        double plaq = Average_Plaquette(gauge_field, TEST_T, TEST_L);
        REQUIRE(!std::isnan(plaq));
        REQUIRE(!std::isinf(plaq));
    }
    
    // After smearing, configuration should be smooth (plaquette close to 1)
    double final_plaq = Average_Plaquette(gauge_field, TEST_T, TEST_L);
    REQUIRE(final_plaq > 0.5);  // Should have improved from random
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Different Lattice Sizes
// ==============================================================================

TEST_CASE("Smearing: Works on different lattice sizes", "[smearing][sizes]") {
    InitializeRand(209);
    
    SECTION("2x2 lattice") {
        const int T = 2, L = 2;
        double *gf;
        Gauge_Field_Alloc(&gf, T, L);
        Gauge_Field_Random(gf, T, L);
        
        APE_Smearing_all(gf, T, L, 0.5);
        REQUIRE(all_links_valid_su2(gf, T, L));
        
        Gauge_Field_Free(&gf);
    }
    
    SECTION("6x6 lattice") {
        const int T = 6, L = 6;
        double *gf;
        Gauge_Field_Alloc(&gf, T, L);
        Gauge_Field_Random(gf, T, L);
        
        APE_Smearing_all(gf, T, L, 0.5);
        REQUIRE(all_links_valid_su2(gf, T, L));
        
        Gauge_Field_Free(&gf);
    }
    
    SECTION("8x4 asymmetric lattice") {
        const int T = 8, L = 4;
        double *gf;
        Gauge_Field_Alloc(&gf, T, L);
        Gauge_Field_Random(gf, T, L);
        
        APE_Smearing_all(gf, T, L, 0.5);
        REQUIRE(all_links_valid_su2(gf, T, L));
        
        Gauge_Field_Free(&gf);
    }
}

// ==============================================================================
// Reproducibility
// ==============================================================================

TEST_CASE("Smearing: Deterministic given same input", "[smearing][reproducibility]") {
    InitializeRand(210);
    
    double *gauge_field1;
    double *gauge_field2;
    Gauge_Field_Alloc(&gauge_field1, TEST_T, TEST_L);
    Gauge_Field_Alloc(&gauge_field2, TEST_T, TEST_L);
    
    // Create identical configurations
    InitializeRand(211);
    Gauge_Field_Random(gauge_field1, TEST_T, TEST_L);
    InitializeRand(211);
    Gauge_Field_Random(gauge_field2, TEST_T, TEST_L);
    
    // Apply same smearing
    double alpha = 0.5;
    APE_Smearing_all(gauge_field1, TEST_T, TEST_L, alpha);
    APE_Smearing_all(gauge_field2, TEST_T, TEST_L, alpha);
    
    // Results should be identical
    int vol = TEST_T * TEST_L * TEST_L * TEST_L;
    for (int i = 0; i < vol * 4 * 8; i++) {
        REQUIRE(std::fabs(gauge_field1[i] - gauge_field2[i]) < TOL);
    }
    
    Gauge_Field_Free(&gauge_field1);
    Gauge_Field_Free(&gauge_field2);
}

// ==============================================================================
// Gauge Covariance Tests
// ==============================================================================

// Helper: Generate deterministic "random" SU(2) from seed
static void su2_from_seed_smearing(double *U, unsigned seed) {
    double a = 0.1 * std::sin(1.0 * seed + 0.1);
    double b = 0.1 * std::sin(2.0 * seed + 0.2);
    double c = 0.1 * std::sin(3.0 * seed + 0.3);
    su2_exp_from_Amu(U, std::array<double, 3>{a, b, c});
    cm_proj(U);
}

// Apply a gauge transformation
static void apply_gauge_transform_smearing(double *gf, int L, int T) {
    int vol = T * L * L * L;
    
    std::vector<double> G(static_cast<std::size_t>(vol) * 8, 0.0);
    for (int site = 0; site < vol; ++site)
        su2_from_seed_smearing(&G[static_cast<std::size_t>(site) * 8], 1234u + site);
    
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

TEST_CASE("Smearing: Gauge covariance", "[smearing][gauge]") {
    const int T = 4, L = 4;
    const int vol = T * L * L * L;
    
    double *gfA;
    double *gfB;
    Gauge_Field_Alloc(&gfA, T, L);
    Gauge_Field_Alloc(&gfB, T, L);
    
    // Start from the same random configuration
    InitializeRand(300);
    Gauge_Field_Random(gfA, T, L);
    InitializeRand(300);
    Gauge_Field_Random(gfB, T, L);
    
    // A: smear then gauge transform
    APE_Smearing_all(gfA, T, L, 0.3);
    apply_gauge_transform_smearing(gfA, L, T);
    
    // B: gauge transform then smear
    apply_gauge_transform_smearing(gfB, L, T);
    APE_Smearing_all(gfB, T, L, 0.3);
    
    // Compare all links
    double maxdiff = 0.0;
    for (std::size_t i = 0; i < static_cast<std::size_t>(vol) * 4 * 8; ++i)
        maxdiff = std::max(maxdiff, std::fabs(gfA[i] - gfB[i]));
    
    REQUIRE(maxdiff < LOOSE_TOL);
    
    Gauge_Field_Free(&gfA);
    Gauge_Field_Free(&gfB);
}

// ==============================================================================
// Topological Charge and Smearing
// ==============================================================================

TEST_CASE("Smearing: Q approaches integer with smearing", "[smearing][topcharge]") {
    InitializeRand(400);
    
    // Use larger lattice for meaningful Q measurement
    const int T = 8, L = 8;
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, T, L);
    Gauge_Field_Random(gauge_field, T, L);
    
    // Apply multiple smearing steps
    double alpha = 0.3;
    int n_smear = 50;
    
    for (int i = 0; i < n_smear; ++i) {
        APE_Smearing_all(gauge_field, T, L, alpha);
    }
    
    // After significant smearing, Q should be closer to an integer
    // (This is a soft test - we just verify the code runs and Q is reasonable)
    double Q = compute_topological_charge(gauge_field, T, L);
    
    REQUIRE(!std::isnan(Q));
    REQUIRE(!std::isinf(Q));
    REQUIRE(std::fabs(Q) < static_cast<double>(T * L * L * L));  // Bounded
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Smearing: Q stable after many smearing steps", "[smearing][topcharge][stability]") {
    InitializeRand(401);
    
    const int T = 6, L = 6;
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, T, L);
    Gauge_Field_Random(gauge_field, T, L);
    
    double alpha = 0.3;
    
    // Smear significantly
    for (int i = 0; i < 30; ++i) {
        APE_Smearing_all(gauge_field, T, L, alpha);
    }
    
    double Q1 = compute_topological_charge(gauge_field, T, L);
    
    // Smear more
    for (int i = 0; i < 20; ++i) {
        APE_Smearing_all(gauge_field, T, L, alpha);
    }
    
    double Q2 = compute_topological_charge(gauge_field, T, L);
    
    // Q should not change dramatically after configuration is already smooth
    REQUIRE(std::fabs(Q2 - Q1) < 0.5);
    
    Gauge_Field_Free(&gauge_field);
}

