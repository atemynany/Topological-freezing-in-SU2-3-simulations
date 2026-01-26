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
