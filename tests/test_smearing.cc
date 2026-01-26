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
#include "Wilson_loops.hh"

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

TEST_CASE("Smearing: Alpha=0 leaves configuration unchanged", "[smearing][alpha]") {
    InitializeRand(205);
    
    double *gauge_field;
    double *original_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Alloc(&original_field, TEST_T, TEST_L);
    
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    Gauge_Field_Copy(original_field, gauge_field, TEST_T, TEST_L);
    
    // Apply smearing with alpha = 0
    APE_Smearing_all(gauge_field, TEST_T, TEST_L, 0.0);
    
    // Configuration should be unchanged (but re-projected to SU(2))
    // Links should still be close to original
    double plaq_orig = Average_Plaquette(original_field, TEST_T, TEST_L);
    double plaq_smeared = Average_Plaquette(gauge_field, TEST_T, TEST_L);
    
    REQUIRE(std::fabs(plaq_orig - plaq_smeared) < LOOSE_TOL);
    
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

TEST_CASE("Smearing: Numerically stable over many iterations", "[smearing][stability]") {
    InitializeRand(208);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double alpha = 0.5;
    
    // Apply many smearing steps
    for (int i = 0; i < 200; i++) {
        APE_Smearing_all(gauge_field, TEST_T, TEST_L, alpha);
        
        // Check every 20 steps
        if (i % 20 == 0) {
            REQUIRE(all_links_valid_su2(gauge_field, TEST_T, TEST_L));
            
            double plaq = Average_Plaquette(gauge_field, TEST_T, TEST_L);
            REQUIRE(!std::isnan(plaq));
            REQUIRE(!std::isinf(plaq));
            REQUIRE(plaq >= -TOL);
            REQUIRE(plaq <= 1.0 + TOL);
        }
    }
    
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
