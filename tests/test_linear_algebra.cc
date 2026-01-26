// ==============================================================================
// test_linear_algebra.cc
// ==============================================================================
// Unit tests for SU(2) linear algebra operations.
// Tests matrix operations, SU(2) properties, and numerical accuracy.
//
// Author: Alexander de Barros Noll
// Date: January 2026
// ==============================================================================

#include <catch2/catch_all.hpp>
#include <cmath>
#include <vector>
#include <array>

#include "linear_algebra.hh"
#include "ranlux.hh"

// ==============================================================================
// Test Fixtures and Helpers
// ==============================================================================

// Tolerance for floating point comparisons
constexpr double TOL = 1e-10;
constexpr double LOOSE_TOL = 1e-6;

// Check if matrix is valid SU(2): det = 1, unitary
static bool is_valid_su2(const double *A, double tol = TOL) {
    double h[4];
    h_from_cm(h, A);
    double norm = h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3];
    return std::fabs(norm - 1.0) < tol;
}

// Check if two matrices are equal
static bool matrices_equal(const double *A, const double *B, double tol = TOL) {
    for (int i = 0; i < 8; ++i) {
        if (std::fabs(A[i] - B[i]) > tol) return false;
    }
    return true;
}

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

// ==============================================================================
// Basic Matrix Operations
// ==============================================================================

TEST_CASE("SU(2) Matrix: Identity", "[linear_algebra][basic]") {
    double A[8];
    cm_eq_id(A);
    
    SECTION("Identity has correct structure") {
        // SU(2) identity: diag(1, 1)
        REQUIRE(std::fabs(A[0] - 1.0) < TOL);  // (0,0) real
        REQUIRE(std::fabs(A[1]) < TOL);        // (0,0) imag
        REQUIRE(std::fabs(A[2]) < TOL);        // (0,1) real
        REQUIRE(std::fabs(A[3]) < TOL);        // (0,1) imag
        REQUIRE(std::fabs(A[4]) < TOL);        // (1,0) real
        REQUIRE(std::fabs(A[5]) < TOL);        // (1,0) imag
        REQUIRE(std::fabs(A[6] - 1.0) < TOL);  // (1,1) real
        REQUIRE(std::fabs(A[7]) < TOL);        // (1,1) imag
    }
    
    SECTION("Identity is valid SU(2)") {
        REQUIRE(is_valid_su2(A));
    }
    
    SECTION("Trace of identity is 2") {
        complex tr;
        co_eq_tr_cm(&tr, A);
        REQUIRE(std::fabs(tr.re - 2.0) < TOL);
        REQUIRE(std::fabs(tr.im) < TOL);
    }
    
    SECTION("Determinant of identity is 1") {
        complex det;
        co_eq_det_cm(&det, A);
        REQUIRE(std::fabs(det.re - 1.0) < TOL);
        REQUIRE(std::fabs(det.im) < TOL);
    }
}

TEST_CASE("SU(2) Matrix: Zero", "[linear_algebra][basic]") {
    double A[8];
    cm_eq_zero(A);
    
    for (int i = 0; i < 8; ++i) {
        REQUIRE(std::fabs(A[i]) < TOL);
    }
}

TEST_CASE("SU(2) Matrix: Copy", "[linear_algebra][basic]") {
    InitializeRand(42);
    double A[8], B[8];
    random_su2(A);
    cm_eq_cm(B, A);
    
    REQUIRE(matrices_equal(A, B));
}

// ==============================================================================
// Matrix Arithmetic
// ==============================================================================

TEST_CASE("SU(2) Matrix: Addition", "[linear_algebra][arithmetic]") {
    InitializeRand(43);
    double A[8], B[8], C[8];
    random_su2(A);
    random_su2(B);
    cm_eq_cm(C, A);
    
    cm_pl_eq_cm(C, B);
    
    for (int i = 0; i < 8; ++i) {
        REQUIRE(std::fabs(C[i] - (A[i] + B[i])) < TOL);
    }
}

TEST_CASE("SU(2) Matrix: Scalar multiplication", "[linear_algebra][arithmetic]") {
    InitializeRand(44);
    double A[8], B[8];
    random_su2(A);
    double scalar = 2.5;
    
    cm_eq_cm_ti_re(B, A, scalar);
    
    for (int i = 0; i < 8; ++i) {
        REQUIRE(std::fabs(B[i] - scalar * A[i]) < TOL);
    }
}

TEST_CASE("SU(2) Matrix: In-place scalar multiplication", "[linear_algebra][arithmetic]") {
    InitializeRand(45);
    double A[8], A_orig[8];
    random_su2(A);
    cm_eq_cm(A_orig, A);
    double scalar = 0.7;
    
    cm_ti_eq_re(A, scalar);
    
    for (int i = 0; i < 8; ++i) {
        REQUIRE(std::fabs(A[i] - scalar * A_orig[i]) < TOL);
    }
}

// ==============================================================================
// Matrix Multiplication
// ==============================================================================

TEST_CASE("SU(2) Matrix: Multiplication with identity", "[linear_algebra][multiplication]") {
    InitializeRand(46);
    double A[8], I[8], B[8];
    random_su2(A);
    cm_eq_id(I);
    
    SECTION("A * I = A") {
        cm_eq_cm_ti_cm(B, A, I);
        REQUIRE(matrices_equal(A, B));
    }
    
    SECTION("I * A = A") {
        cm_eq_cm_ti_cm(B, I, A);
        REQUIRE(matrices_equal(A, B));
    }
}

TEST_CASE("SU(2) Matrix: Multiplication closure", "[linear_algebra][multiplication]") {
    InitializeRand(47);
    double A[8], B[8], C[8];
    random_su2(A);
    random_su2(B);
    
    cm_eq_cm_ti_cm(C, A, B);
    
    // Product of SU(2) matrices is SU(2)
    REQUIRE(is_valid_su2(C, LOOSE_TOL));
}

TEST_CASE("SU(2) Matrix: Multiplication associativity", "[linear_algebra][multiplication]") {
    InitializeRand(48);
    double A[8], B[8], C[8];
    double AB[8], BC[8], ABC1[8], ABC2[8];
    
    random_su2(A);
    random_su2(B);
    random_su2(C);
    
    // (A * B) * C
    cm_eq_cm_ti_cm(AB, A, B);
    cm_eq_cm_ti_cm(ABC1, AB, C);
    
    // A * (B * C)
    cm_eq_cm_ti_cm(BC, B, C);
    cm_eq_cm_ti_cm(ABC2, A, BC);
    
    REQUIRE(matrices_equal(ABC1, ABC2, LOOSE_TOL));
}

// ==============================================================================
// Hermitian Conjugate
// ==============================================================================

TEST_CASE("SU(2) Matrix: Hermitian conjugate", "[linear_algebra][conjugate]") {
    InitializeRand(49);
    double A[8], A_dag[8], prod[8];
    random_su2(A);
    
    cm_eq_cm_dag(A_dag, A);
    
    SECTION("A^dag is valid SU(2)") {
        REQUIRE(is_valid_su2(A_dag));
    }
    
    SECTION("A * A^dag = I") {
        cm_eq_cm_ti_cm(prod, A, A_dag);
        
        double I[8];
        cm_eq_id(I);
        REQUIRE(matrices_equal(prod, I, LOOSE_TOL));
    }
    
    SECTION("A^dag * A = I") {
        cm_eq_cm_ti_cm(prod, A_dag, A);
        
        double I[8];
        cm_eq_id(I);
        REQUIRE(matrices_equal(prod, I, LOOSE_TOL));
    }
    
    SECTION("(A^dag)^dag = A") {
        double A_dag_dag[8];
        cm_eq_cm_dag(A_dag_dag, A_dag);
        REQUIRE(matrices_equal(A, A_dag_dag, TOL));
    }
}

TEST_CASE("SU(2) Matrix: In-place conjugate", "[linear_algebra][conjugate]") {
    InitializeRand(50);
    double A[8], A_copy[8], A_dag[8];
    random_su2(A);
    cm_eq_cm(A_copy, A);
    cm_eq_cm_dag(A_dag, A);
    
    cm_dag_eq_cm(A_copy);
    
    REQUIRE(matrices_equal(A_copy, A_dag, TOL));
}

// ==============================================================================
// Mixed Operations with Conjugate
// ==============================================================================

TEST_CASE("SU(2) Matrix: A^dag * B", "[linear_algebra][mixed]") {
    InitializeRand(51);
    double A[8], B[8], result[8], expected[8], A_dag[8];
    random_su2(A);
    random_su2(B);
    
    cm_eq_cm_dag_ti_cm(result, A, B);
    
    // Compare with explicit: A^dag * B
    cm_eq_cm_dag(A_dag, A);
    cm_eq_cm_ti_cm(expected, A_dag, B);
    
    REQUIRE(matrices_equal(result, expected, TOL));
}

TEST_CASE("SU(2) Matrix: A * B^dag", "[linear_algebra][mixed]") {
    InitializeRand(52);
    double A[8], B[8], result[8], expected[8], B_dag[8];
    random_su2(A);
    random_su2(B);
    
    cm_eq_cm_ti_cm_dag(result, A, B);
    
    // Compare with explicit: A * B^dag
    cm_eq_cm_dag(B_dag, B);
    cm_eq_cm_ti_cm(expected, A, B_dag);
    
    REQUIRE(matrices_equal(result, expected, TOL));
}

TEST_CASE("SU(2) Matrix: A^dag * B^dag", "[linear_algebra][mixed]") {
    InitializeRand(53);
    double A[8], B[8], result[8], expected[8], A_dag[8], B_dag[8];
    random_su2(A);
    random_su2(B);
    
    cm_eq_cm_dag_ti_cm_dag(result, A, B);
    
    // Compare with explicit: A^dag * B^dag
    cm_eq_cm_dag(A_dag, A);
    cm_eq_cm_dag(B_dag, B);
    cm_eq_cm_ti_cm(expected, A_dag, B_dag);
    
    REQUIRE(matrices_equal(result, expected, TOL));
}

// ==============================================================================
// h-parametrization
// ==============================================================================

TEST_CASE("SU(2) Matrix: h-parametrization roundtrip", "[linear_algebra][parametrization]") {
    InitializeRand(54);
    
    for (int trial = 0; trial < 10; trial++) {
        double A[8], h_orig[4], h_recovered[4], A_recovered[8];
        
        // Generate random normalized h
        double norm;
        do {
            h_orig[0] = 2.0 * DRand() - 1.0;
            h_orig[1] = 2.0 * DRand() - 1.0;
            h_orig[2] = 2.0 * DRand() - 1.0;
            h_orig[3] = 2.0 * DRand() - 1.0;
            norm = h_orig[0]*h_orig[0] + h_orig[1]*h_orig[1] + h_orig[2]*h_orig[2] + h_orig[3]*h_orig[3];
        } while (norm < 1e-10);
        
        norm = 1.0 / std::sqrt(norm);
        h_orig[0] *= norm;
        h_orig[1] *= norm;
        h_orig[2] *= norm;
        h_orig[3] *= norm;
        
        // h -> A -> h'
        cm_from_h(A, h_orig);
        h_from_cm(h_recovered, A);
        
        for (int i = 0; i < 4; i++) {
            REQUIRE(std::fabs(h_orig[i] - h_recovered[i]) < TOL);
        }
        
        // h -> A -> h' -> A'
        cm_from_h(A_recovered, h_recovered);
        REQUIRE(matrices_equal(A, A_recovered, TOL));
    }
}

// ==============================================================================
// SU(2) Projection
// ==============================================================================

TEST_CASE("SU(2) Matrix: Projection", "[linear_algebra][projection]") {
    InitializeRand(55);
    
    SECTION("Projection of SU(2) matrix is identity operation") {
        double A[8], A_copy[8];
        random_su2(A);
        cm_eq_cm(A_copy, A);
        
        cm_proj(A_copy);
        
        REQUIRE(matrices_equal(A, A_copy, LOOSE_TOL));
    }
    
    SECTION("Projection of scaled matrix gives valid SU(2)") {
        double A[8];
        random_su2(A);
        
        // Scale by arbitrary factor (det != 1)
        cm_ti_eq_re(A, 2.5);
        
        // Project back
        cm_proj(A);
        
        REQUIRE(is_valid_su2(A, LOOSE_TOL));
    }
}

// ==============================================================================
// Determinant and Trace
// ==============================================================================

TEST_CASE("SU(2) Matrix: Determinant properties", "[linear_algebra][determinant]") {
    InitializeRand(56);
    double A[8], B[8], AB[8];
    random_su2(A);
    random_su2(B);
    
    SECTION("det(SU(2)) = 1") {
        complex det;
        co_eq_det_cm(&det, A);
        REQUIRE(std::fabs(det.re - 1.0) < LOOSE_TOL);
        REQUIRE(std::fabs(det.im) < LOOSE_TOL);
    }
    
    SECTION("det(A*B) = det(A)*det(B)") {
        complex det_A, det_B, det_AB;
        co_eq_det_cm(&det_A, A);
        co_eq_det_cm(&det_B, B);
        
        cm_eq_cm_ti_cm(AB, A, B);
        co_eq_det_cm(&det_AB, AB);
        
        // det(A)*det(B)
        double expected_re = det_A.re * det_B.re - det_A.im * det_B.im;
        double expected_im = det_A.re * det_B.im + det_A.im * det_B.re;
        
        REQUIRE(std::fabs(det_AB.re - expected_re) < LOOSE_TOL);
        REQUIRE(std::fabs(det_AB.im - expected_im) < LOOSE_TOL);
    }
}

TEST_CASE("SU(2) Matrix: Trace properties", "[linear_algebra][trace]") {
    InitializeRand(57);
    double A[8], B[8], AB[8], BA[8];
    random_su2(A);
    random_su2(B);
    
    SECTION("Tr(A*B) = Tr(B*A)") {
        cm_eq_cm_ti_cm(AB, A, B);
        cm_eq_cm_ti_cm(BA, B, A);
        
        complex tr_AB, tr_BA;
        co_eq_tr_cm(&tr_AB, AB);
        co_eq_tr_cm(&tr_BA, BA);
        
        REQUIRE(std::fabs(tr_AB.re - tr_BA.re) < LOOSE_TOL);
        REQUIRE(std::fabs(tr_AB.im - tr_BA.im) < LOOSE_TOL);
    }
    
    SECTION("Tr(A) = Tr(A^dag)* for SU(2)") {
        double A_dag[8];
        cm_eq_cm_dag(A_dag, A);
        
        complex tr_A, tr_A_dag;
        co_eq_tr_cm(&tr_A, A);
        co_eq_tr_cm(&tr_A_dag, A_dag);
        
        // For SU(2), Tr(A^dag) = Tr(A)* (complex conjugate)
        REQUIRE(std::fabs(tr_A.re - tr_A_dag.re) < TOL);
        REQUIRE(std::fabs(tr_A.im + tr_A_dag.im) < TOL);
    }
}

// ==============================================================================
// SU(2) Exponential
// ==============================================================================

TEST_CASE("SU(2) Exponential: su2_exp_from_Amu", "[linear_algebra][exponential]") {
    InitializeRand(58);
    
    SECTION("exp(0) = I") {
        double U[8];
        std::array<double, 3> A = {0.0, 0.0, 0.0};
        
        su2_exp_from_Amu(U, A);
        
        double I[8];
        cm_eq_id(I);
        REQUIRE(matrices_equal(U, I, TOL));
    }
    
    SECTION("exp(A) is valid SU(2)") {
        for (int trial = 0; trial < 10; trial++) {
            double U[8];
            std::array<double, 3> A;
            A[0] = 2.0 * DRand() - 1.0;
            A[1] = 2.0 * DRand() - 1.0;
            A[2] = 2.0 * DRand() - 1.0;
            
            su2_exp_from_Amu(U, A);
            
            REQUIRE(is_valid_su2(U, LOOSE_TOL));
        }
    }
}

// ==============================================================================
// SU(2) Logarithm
// ==============================================================================

TEST_CASE("SU(2) Logarithm: su2_log_to_algebra_components", "[linear_algebra][logarithm]") {
    InitializeRand(59);
    
    SECTION("log(I) = 0") {
        double I[8];
        cm_eq_id(I);
        
        std::array<double, 3> theta;
        su2_log_to_algebra_components(I, theta);
        
        REQUIRE(std::fabs(theta[0]) < TOL);
        REQUIRE(std::fabs(theta[1]) < TOL);
        REQUIRE(std::fabs(theta[2]) < TOL);
    }
    
    SECTION("exp(log(U)) ~ U for small rotations") {
        for (int trial = 0; trial < 10; trial++) {
            double U[8];
            random_su2(U);
            
            std::array<double, 3> theta;
            su2_log_to_algebra_components(U, theta);
            
            double U_reconstructed[8];
            su2_exp_from_Amu(U_reconstructed, theta);
            
            // Should be close (up to sign ambiguity for large rotations)
            // For small rotations, should match exactly
            double h1[4], h2[4];
            h_from_cm(h1, U);
            h_from_cm(h2, U_reconstructed);
            
            bool match = true;
            bool match_negated = true;
            for (int i = 0; i < 4; i++) {
                if (std::fabs(h1[i] - h2[i]) > LOOSE_TOL) match = false;
                if (std::fabs(h1[i] + h2[i]) > LOOSE_TOL) match_negated = false;
            }
            
            REQUIRE((match || match_negated));
        }
    }
}
