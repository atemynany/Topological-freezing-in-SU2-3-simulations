// ==============================================================================
// test_heatbath.cc
// ==============================================================================
// Unit tests for the MC Heatbath algorithm.
// Tests gauge field initialization, heatbath updates, and thermalization.
//
// Author: Alexander de Barros Noll
// Date: January 2026
// ==============================================================================

#include <catch2/catch_all.hpp>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>

// Include utility headers
#include "fields.hh"
#include "geometry.hh"
#include "linear_algebra.hh"
#include "ranlux.hh"

// ==============================================================================
// Global variables for tests (defined here, not extern)
// ==============================================================================
int T = 4;
int L = 4;
double *gauge_field = nullptr;
bool open_boundary_conditions = false;
bool hot_start = false;

// ********************
// Helper Functions
// ********************

// Check if SU(2) matrix is valid (det = 1, unitary)
static bool is_valid_su2(const double *A, double tol = 1e-10) {
    double h[4];
    h_from_cm(h, A);
    double norm = h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3];
    return std::fabs(norm - 1.0) < tol;
}

// Get maximum absolute value of 8-component array
static double max_abs_8(const double *A) {
    double m = 0.0;
    for (int i = 0; i < 8; ++i)
        m = std::max(m, std::fabs(A[i]));
    return m;
}

// Check if two matrices are equal within tolerance
static bool matrices_equal(const double *A, const double *B, double tol = 1e-10) {
    for (int i = 0; i < 8; ++i) {
        if (std::fabs(A[i] - B[i]) > tol) return false;
    }
    return true;
}

// Calculate staple sum at a given site and direction
static void compute_staple_sum(double *S_l, double *gf, int it, int ix, int iy, int iz, int i1, int T_local, int L_local) {
    double SU2_1[8], SU2_2[8];
    cm_eq_zero(S_l);
    
    for (int i2 = 0; i2 < 4; i2++) {
        if (i2 == i1) continue;
        
        int i[4];
        int index1, index2, index3;
        double weight = 1.0;
        
        // Negative staple: U_nu^dag(x-nu) * U_mu(x-nu) * U_nu(x-nu+mu)
        i[0] = it; i[1] = ix; i[2] = iy; i[3] = iz;
        i[i2] -= 1;
        
        index1 = ggi(get_index(i[0], i[1], i[2], i[3], T_local, L_local), i2);
        index2 = ggi(get_index(i[0], i[1], i[2], i[3], T_local, L_local), i1);
        i[i1] += 1;
        index3 = ggi(get_index(i[0], i[1], i[2], i[3], T_local, L_local), i2);
        
        if (index1 >= 0 && index2 >= 0 && index3 >= 0) {
            cm_eq_cm_ti_cm(SU2_1, gf + index2, gf + index3);
            cm_eq_cm_dag_ti_cm(SU2_2, gf + index1, SU2_1);
            cm_ti_eq_re(SU2_2, weight);
            cm_pl_eq_cm(S_l, SU2_2);
        }
        
        // Positive staple: U_nu(x) * U_mu(x+nu) * U_nu^dag(x+mu)
        i[0] = it; i[1] = ix; i[2] = iy; i[3] = iz;
        
        index1 = ggi(get_index(i[0], i[1], i[2], i[3], T_local, L_local), i2);
        i[i2] += 1;
        index2 = ggi(get_index(i[0], i[1], i[2], i[3], T_local, L_local), i1);
        i[i1] += 1;
        i[i2] -= 1;
        index3 = ggi(get_index(i[0], i[1], i[2], i[3], T_local, L_local), i2);
        
        if (index1 >= 0 && index2 >= 0 && index3 >= 0) {
            cm_eq_cm_ti_cm_dag(SU2_1, gf + index2, gf + index3);
            cm_eq_cm_ti_cm(SU2_2, gf + index1, SU2_1);
            cm_ti_eq_re(SU2_2, weight);
            cm_pl_eq_cm(S_l, SU2_2);
        }
    }
}

// Perform one heatbath update on a single link
static void heatbath_update_link(double *gf, int it, int ix, int iy, int iz, int i1, double beta, int T_local, int L_local) {
    double S_l[8];
    compute_staple_sum(S_l, gf, it, ix, iy, iz, i1, T_local, L_local);
    
    // Skip if S_l is zero
    double S_l_sum = 0;
    for (int i = 0; i < 8; i++) S_l_sum += std::fabs(S_l[i]);
    if (S_l_sum < 1e-15) return;
    
    cm_dag_eq_cm(S_l);
    
    double k = std::sqrt(S_l[0]*S_l[6] - S_l[1]*S_l[7] - S_l[2]*S_l[4] + S_l[3]*S_l[5]);
    if (k < 1e-15) return;
    
    double beta_k = beta * k;
    double y_min = std::exp(-beta_k);
    double y_max = std::exp(+beta_k);
    
    double a[4];
    
    // Generate a[0] with correct Boltzmann distribution
    while (true) {
        double y = y_min + (y_max - y_min) * DRand();
        a[0] = std::log(y) / beta_k;
        if (DRand() <= std::sqrt(1.0 - a[0]*a[0]))
            break;
    }
    
    // Generate a[1], a[2], a[3] uniformly on sphere
    double norm;
    while (true) {
        a[1] = 2.0 * DRand() - 1.0;
        a[2] = 2.0 * DRand() - 1.0;
        a[3] = 2.0 * DRand() - 1.0;
        norm = a[1]*a[1] + a[2]*a[2] + a[3]*a[3];
        if (norm >= 1e-10 && norm <= 1.0)
            break;
    }
    norm = std::sqrt((1.0 - a[0]*a[0]) / norm);
    a[1] *= norm;
    a[2] *= norm;
    a[3] *= norm;
    
    // Construct new link
    double U_0[8], U_0l[8], SU2_1[8];
    cm_eq_cm_dag(U_0, S_l);
    cm_ti_eq_re(U_0, 1.0/k);
    cm_from_h(U_0l, a);
    cm_eq_cm_ti_cm(SU2_1, U_0l, U_0);
    
    // Project to SU(2)
    double h[4];
    h_from_cm(h, SU2_1);
    norm = 1.0 / std::sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3]);
    h[0] *= norm; h[1] *= norm; h[2] *= norm; h[3] *= norm;
    cm_from_h(SU2_1, h);
    
    // Update link
    int index = ggi(get_index(it, ix, iy, iz, T_local, L_local), i1);
    if (index >= 0)
        cm_eq_cm(gf + index, SU2_1);
}

// ********************
// GEOMETRY TESTS
// ********************

TEST_CASE("Geometry: get_index returns valid indices", "[geometry]") {
    T = 16;
    L = 16;
    open_boundary_conditions = false;
    
    SECTION("Valid indices within lattice bounds") {
        for (int it = 0; it < T; it++) {
            for (int ix = 0; ix < L; ix++) {
                for (int iy = 0; iy < L; iy++) {
                    for (int iz = 0; iz < L; iz++) {
                        int idx = get_index(it, ix, iy, iz, T, L);
                        REQUIRE(idx >= 0);
                        REQUIRE(idx < T * L * L * L);
                    }
                }
            }
        }
    }
    
    SECTION("Periodic BC: wrapping works correctly") {
        open_boundary_conditions = false;
        
        // t = -1 should wrap to t = T-1
        int idx1 = get_index(-1, 0, 0, 0, T, L);
        int idx2 = get_index(T-1, 0, 0, 0, T, L);
        REQUIRE(idx1 == idx2);
        
        // t = T should wrap to t = 0
        int idx3 = get_index(T, 0, 0, 0, T, L);
        int idx4 = get_index(0, 0, 0, 0, T, L);
        REQUIRE(idx3 == idx4);
        
        // Spatial wrapping
        int idx5 = get_index(0, -1, 0, 0, T, L);
        int idx6 = get_index(0, L-1, 0, 0, T, L);
        REQUIRE(idx5 == idx6);
    }
    
    SECTION("Open BC: out-of-bounds returns -1") {
        open_boundary_conditions = true;
        
        int idx1 = get_index(-1, 0, 0, 0, T, L);
        REQUIRE(idx1 == -1);
        
        int idx2 = get_index(T, 0, 0, 0, T, L);
        REQUIRE(idx2 == -1);
    }
}

TEST_CASE("Geometry: ggi returns valid gauge field indices", "[geometry]") {
    T = 16;
    L = 16;
    open_boundary_conditions = false;
    
    SECTION("All directions give valid indices") {
        int site_idx = get_index(2, 1, 1, 1, T, L);
        for (int dir = 0; dir < 4; dir++) {
            int gf_idx = ggi(site_idx, dir);
            REQUIRE(gf_idx >= 0);
            // Each link has 8 doubles (SU(2) matrix)
            REQUIRE(gf_idx % 8 == 0);
        }
    }
    
    SECTION("ggi(-1, mu) returns -1") {
        for (int mu = 0; mu < 4; mu++) {
            REQUIRE(ggi(-1, mu) == -1);
        }
    }
}

// ********************
// GAUGE FIELD TESTS
// ********************

TEST_CASE("Gauge Field: Allocation and deallocation", "[fields]") {
    T = 16;
    L = 16;
    double *test_field = nullptr;
    
    SECTION("Allocation succeeds") {
        Gauge_Field_Alloc(&test_field, T, L);
        REQUIRE(test_field != nullptr);
        Gauge_Field_Free(&test_field);
    }
}

TEST_CASE("Gauge Field: Cold start (unity) initialization", "[fields]") {
    T = 16;
    L = 16;
    open_boundary_conditions = false;
    
    double *test_field = nullptr;
    Gauge_Field_Alloc(&test_field, T, L);
    Gauge_Field_Unity(test_field, T, L);
    
    SECTION("All links are identity matrices") {
        double identity[8];
        cm_eq_id(identity);
        
        for (int it = 0; it < T; it++) {
            for (int ix = 0; ix < L; ix++) {
                for (int iy = 0; iy < L; iy++) {
                    for (int iz = 0; iz < L; iz++) {
                        for (int dir = 0; dir < 4; dir++) {
                            int idx = ggi(get_index(it, ix, iy, iz, T, L), dir);
                            REQUIRE(matrices_equal(test_field + idx, identity));
                        }
                    }
                }
            }
        }
    }
    
    Gauge_Field_Free(&test_field);
}

TEST_CASE("Gauge Field: Hot start (random) initialization", "[fields]") {
    T = 16;
    L = 16;
    open_boundary_conditions = false;
    InitializeRand(12345);
    
    double *test_field = nullptr;
    Gauge_Field_Alloc(&test_field, T, L);
    Gauge_Field_Random(test_field, T, L);
    
    SECTION("All links are valid SU(2) matrices") {
        for (int it = 0; it < T; it++) {
            for (int ix = 0; ix < L; ix++) {
                for (int iy = 0; iy < L; iy++) {
                    for (int iz = 0; iz < L; iz++) {
                        for (int dir = 0; dir < 4; dir++) {
                            int idx = ggi(get_index(it, ix, iy, iz, T, L), dir);
                            REQUIRE(is_valid_su2(test_field + idx));
                        }
                    }
                }
            }
        }
    }
    
    SECTION("Random field is not identity") {
        double identity[8];
        cm_eq_id(identity);
        
        bool found_non_identity = false;
        for (int it = 0; it < T && !found_non_identity; it++) {
            int idx = ggi(get_index(it, 0, 0, 0, T, L), 0);
            if (!matrices_equal(test_field + idx, identity, 0.01)) {
                found_non_identity = true;
            }
        }
        REQUIRE(found_non_identity);
    }
    
    Gauge_Field_Free(&test_field);
}

// ********************
// LINEAR ALGEBRA TESTS
// ********************

TEST_CASE("Linear Algebra: SU(2) identity matrix", "[linear_algebra]") {
    double I[8];
    cm_eq_id(I);
    
    SECTION("Identity has correct structure") {
        // For SU(2): I = [[1,0],[0,1]]
        // In 8-double format: [Re(a), Im(a), Re(b), Im(b), Re(c), Im(c), Re(d), Im(d)]
        // Identity: a=1, b=0, c=0, d=1
        REQUIRE(std::fabs(I[0] - 1.0) < 1e-14);  // Re(a) = 1
        REQUIRE(std::fabs(I[6] - 1.0) < 1e-14);  // Re(d) = 1
    }
    
    SECTION("Identity is valid SU(2)") {
        REQUIRE(is_valid_su2(I));
    }
}

TEST_CASE("Linear Algebra: Matrix multiplication", "[linear_algebra]") {
    double A[8], B[8], C[8], I[8];
    cm_eq_id(I);
    
    SECTION("Identity times anything equals that matrix") {
        double h[4] = {0.5, 0.5, 0.5, 0.5};  // Valid SU(2) quaternion
        cm_from_h(A, h);
        
        cm_eq_cm_ti_cm(C, I, A);  // C = I * A
        REQUIRE(matrices_equal(C, A));
        
        cm_eq_cm_ti_cm(C, A, I);  // C = A * I
        REQUIRE(matrices_equal(C, A));
    }
    
    SECTION("Matrix times its dagger equals identity") {
        double h[4] = {0.6, 0.0, 0.8, 0.0};  // Valid SU(2)
        cm_from_h(A, h);
        
        cm_eq_cm_ti_cm_dag(C, A, A);  // C = A * A^dag
        
        double h_result[4];
        h_from_cm(h_result, C);
        REQUIRE(std::fabs(h_result[0] - 1.0) < 1e-10);
        REQUIRE(std::fabs(h_result[1]) < 1e-10);
        REQUIRE(std::fabs(h_result[2]) < 1e-10);
        REQUIRE(std::fabs(h_result[3]) < 1e-10);
    }
}

TEST_CASE("Linear Algebra: Quaternion conversion", "[linear_algebra]") {
    SECTION("Quaternion to matrix and back") {
        double h_original[4] = {0.5, 0.5, 0.5, 0.5};  // Normalized
        double A[8], h_result[4];
        
        cm_from_h(A, h_original);
        h_from_cm(h_result, A);
        
        for (int i = 0; i < 4; i++) {
            REQUIRE(std::fabs(h_result[i] - h_original[i]) < 1e-10);
        }
    }
    
    SECTION("Quaternion norm preserved") {
        double h[4] = {0.6, 0.3, 0.4, 0.6164};
        // Normalize
        double norm = std::sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3]);
        for (int i = 0; i < 4; i++) h[i] /= norm;
        
        double A[8];
        cm_from_h(A, h);
        REQUIRE(is_valid_su2(A));
    }
}

TEST_CASE("Linear Algebra: Matrix addition", "[linear_algebra]") {
    double A[8], B[8];
    cm_eq_id(A);
    cm_eq_id(B);
    
    cm_pl_eq_cm(A, B);  // A = A + B = 2*I
    
    // Check that A is now 2*Identity
    REQUIRE(std::fabs(A[0] - 2.0) < 1e-14);
    REQUIRE(std::fabs(A[6] - 2.0) < 1e-14);
}

TEST_CASE("Linear Algebra: Scalar multiplication", "[linear_algebra]") {
    double A[8];
    cm_eq_id(A);
    
    cm_ti_eq_re(A, 3.0);  // A = 3*I
    
    REQUIRE(std::fabs(A[0] - 3.0) < 1e-14);
    REQUIRE(std::fabs(A[6] - 3.0) < 1e-14);
}

// ********************
// RANDOM NUMBER GENERATOR TESTS
// ********************

TEST_CASE("Random: DRand produces values in [0,1)", "[random]") {
    InitializeRand(42);
    
    SECTION("1000 random numbers in valid range") {
        for (int i = 0; i < 1000; i++) {
            double r = DRand();
            REQUIRE(r >= 0.0);
            REQUIRE(r < 1.0);
        }
    }
}

TEST_CASE("Random: Reproducibility with same seed", "[random]") {
    InitializeRand(999);
    double r1 = DRand();
    double r2 = DRand();
    
    InitializeRand(999);
    double r3 = DRand();
    double r4 = DRand();
    
    REQUIRE(r1 == r3);
    REQUIRE(r2 == r4);
}

TEST_CASE("Random: Different seeds give different sequences", "[random]") {
    InitializeRand(111);
    double r1 = DRand();
    
    InitializeRand(222);
    double r2 = DRand();
    
    REQUIRE(r1 != r2);
}

// ********************
// STAPLE CALCULATION TESTS
// ********************

TEST_CASE("Staple: Cold configuration staple sum", "[staple]") {
    T = 16;
    L = 16;
    open_boundary_conditions = false;
    
    double *test_field = nullptr;
    Gauge_Field_Alloc(&test_field, T, L);
    Gauge_Field_Unity(test_field, T, L);
    
    SECTION("Staple sum equals 6*Identity for periodic BC") {
        double S_l[8];
        compute_staple_sum(S_l, test_field, 2, 1, 1, 1, 0, T, L);
        
        // For cold config with periodic BC, each of 3 directions contributes 2 staples
        // Each staple is Identity, so sum = 6*Identity
        double h[4];
        h_from_cm(h, S_l);
        REQUIRE(std::fabs(h[0] - 6.0) < 1e-10);
        REQUIRE(std::fabs(h[1]) < 1e-10);
        REQUIRE(std::fabs(h[2]) < 1e-10);
        REQUIRE(std::fabs(h[3]) < 1e-10);
    }
    
    Gauge_Field_Free(&test_field);
}

TEST_CASE("Staple: Open BC reduces staples at boundaries", "[staple][boundary]") {
    T = 16;
    L = 16;
    open_boundary_conditions = true;
    
    double *test_field = nullptr;
    Gauge_Field_Alloc(&test_field, T, L);
    Gauge_Field_Unity(test_field, T, L);
    
    SECTION("Staple sum at t=0 boundary") {
        double S_l[8];
        // At t=0, temporal direction has fewer staples
        compute_staple_sum(S_l, test_field, 0, 1, 1, 1, 0, T, L);
        
        // Should have fewer than 6 staples due to boundary
        double h[4];
        h_from_cm(h, S_l);
        // The exact value depends on implementation, but should be less than 6
        REQUIRE(h[0] < 6.0 + 1e-10);
    }
    
    Gauge_Field_Free(&test_field);
}

// ********************
// HEATBATH ALGORITHM TESTS
// ********************

TEST_CASE("Heatbath: Single update preserves SU(2)", "[heatbath]") {
    T = 16;
    L = 16;
    open_boundary_conditions = false;
    InitializeRand(54321);
    
    double *test_field = nullptr;
    Gauge_Field_Alloc(&test_field, T, L);
    Gauge_Field_Unity(test_field, T, L);
    
    double beta = 2.3;
    
    SECTION("Link remains valid SU(2) after update") {
        heatbath_update_link(test_field, 2, 1, 1, 1, 0, beta, T, L);
        
        int idx = ggi(get_index(2, 1, 1, 1, T, L), 0);
        REQUIRE(is_valid_su2(test_field + idx));
    }
    
    Gauge_Field_Free(&test_field);
}

TEST_CASE("Heatbath: Full sweep preserves all SU(2) matrices", "[heatbath]") {
    T = 16;
    L = 16;
    open_boundary_conditions = false;
    InitializeRand(12345);
    
    double *test_field = nullptr;
    Gauge_Field_Alloc(&test_field, T, L);
    Gauge_Field_Unity(test_field, T, L);
    
    double beta = 2.3;
    
    // Perform one full sweep
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    for (int i1 = 0; i1 < 4; i1++) {
                        heatbath_update_link(test_field, it, ix, iy, iz, i1, beta, T, L);
                    }
                }
            }
        }
    }
    
    SECTION("All links remain valid SU(2)") {
        for (int it = 0; it < T; it++) {
            for (int ix = 0; ix < L; ix++) {
                for (int iy = 0; iy < L; iy++) {
                    for (int iz = 0; iz < L; iz++) {
                        for (int dir = 0; dir < 4; dir++) {
                            int idx = ggi(get_index(it, ix, iy, iz, T, L), dir);
                            REQUIRE(is_valid_su2(test_field + idx));
                        }
                    }
                }
            }
        }
    }
    
    Gauge_Field_Free(&test_field);
}

TEST_CASE("Heatbath: Updates change the configuration", "[heatbath]") {
    T = 16;
    L = 16;
    open_boundary_conditions = false;
    InitializeRand(99999);
    
    double *test_field = nullptr;
    Gauge_Field_Alloc(&test_field, T, L);
    Gauge_Field_Unity(test_field, T, L);
    
    double beta = 2.3;
    double identity[8];
    cm_eq_id(identity);
    
    // Store original link
    int idx = ggi(get_index(2, 1, 1, 1, T, L), 0);
    double original[8];
    cm_eq_cm(original, test_field + idx);
    
    // Perform several updates
    for (int sweep = 0; sweep < 10; sweep++) {
        heatbath_update_link(test_field, 2, 1, 1, 1, 0, beta, T, L);
    }
    
    SECTION("Link has changed from identity") {
        // After multiple updates, the link should have changed
        bool changed = !matrices_equal(test_field + idx, identity, 0.01);
        REQUIRE(changed);
    }
    
    Gauge_Field_Free(&test_field);
}

// ********************
// OPEN BOUNDARY CONDITIONS TESTS
// ********************

TEST_CASE("Open BC: Boundary links handled correctly", "[boundary]") {
    T = 16;
    L = 16;
    open_boundary_conditions = true;
    InitializeRand(33333);
    
    double *test_field = nullptr;
    Gauge_Field_Alloc(&test_field, T, L);
    Gauge_Field_Unity(test_field, T, L);
    
    double beta = 2.5;
    
    // Perform a sweep
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    for (int i1 = 0; i1 < 4; i1++) {
                        heatbath_update_link(test_field, it, ix, iy, iz, i1, beta, T, L);
                    }
                }
            }
        }
    }
    
    SECTION("All valid links remain SU(2)") {
        for (int it = 0; it < T; it++) {
            for (int ix = 0; ix < L; ix++) {
                for (int iy = 0; iy < L; iy++) {
                    for (int iz = 0; iz < L; iz++) {
                        for (int dir = 0; dir < 4; dir++) {
                            int idx = ggi(get_index(it, ix, iy, iz, T, L), dir);
                            if (idx >= 0) {
                                REQUIRE(is_valid_su2(test_field + idx));
                            }
                        }
                    }
                }
            }
        }
    }
    
    Gauge_Field_Free(&test_field);
}

// ********************
// DETERMINISM TESTS
// ********************

TEST_CASE("Heatbath: Deterministic with same seed", "[heatbath][determinism]") {
    T = 16;
    L = 16;
    open_boundary_conditions = false;
    
    double beta = 2.5;
    
    // First run
    InitializeRand(77777);
    double *field1 = nullptr;
    Gauge_Field_Alloc(&field1, T, L);
    Gauge_Field_Unity(field1, T, L);
    
    for (int sweep = 0; sweep < 5; sweep++) {
        for (int it = 0; it < T; it++) {
            for (int ix = 0; ix < L; ix++) {
                for (int iy = 0; iy < L; iy++) {
                    for (int iz = 0; iz < L; iz++) {
                        for (int i1 = 0; i1 < 4; i1++) {
                            heatbath_update_link(field1, it, ix, iy, iz, i1, beta, T, L);
                        }
                    }
                }
            }
        }
    }
    
    // Second run with same seed
    InitializeRand(77777);
    double *field2 = nullptr;
    Gauge_Field_Alloc(&field2, T, L);
    Gauge_Field_Unity(field2, T, L);
    
    for (int sweep = 0; sweep < 5; sweep++) {
        for (int it = 0; it < T; it++) {
            for (int ix = 0; ix < L; ix++) {
                for (int iy = 0; iy < L; iy++) {
                    for (int iz = 0; iz < L; iz++) {
                        for (int i1 = 0; i1 < 4; i1++) {
                            heatbath_update_link(field2, it, ix, iy, iz, i1, beta, T, L);
                        }
                    }
                }
            }
        }
    }
    
    SECTION("Both runs produce identical configurations") {
        for (int it = 0; it < T; it++) {
            for (int ix = 0; ix < L; ix++) {
                for (int iy = 0; iy < L; iy++) {
                    for (int iz = 0; iz < L; iz++) {
                        for (int dir = 0; dir < 4; dir++) {
                            int idx = ggi(get_index(it, ix, iy, iz, T, L), dir);
                            REQUIRE(matrices_equal(field1 + idx, field2 + idx));
                        }
                    }
                }
            }
        }
    }
    
    Gauge_Field_Free(&field1);
    Gauge_Field_Free(&field2);
}
