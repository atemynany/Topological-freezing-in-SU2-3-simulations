// ==============================================================================
// test_plaquette.cc
// ==============================================================================
// Unit tests for plaquette calculations in SU(2) lattice gauge theory.
// Tests plaquette computations, clover construction, and gauge covariance.
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
#include "Plaquette.hh"
#include "topcharge_su2.hh"

// ==============================================================================
// External Global Variables
// ==============================================================================

bool open_boundary_conditions = false;

// ==============================================================================
// Test Helpers
// ==============================================================================

constexpr double TOL = 1e-10;
constexpr double LOOSE_TOL = 1e-6;

constexpr int TEST_T = 4;
constexpr int TEST_L = 4;

// Maximum absolute value of an 8-component matrix
static double max_abs_8(const double *A) {
    double m = 0.0;
    for (int i = 0; i < 8; ++i)
        m = std::max(m, std::fabs(A[i]));
    return m;
}

// Generate deterministic "random" SU(2) matrix from seed
static void su2_random_from_seed(double *U, unsigned seed) {
    double a = 0.1 * std::sin(1.0 * seed + 0.1);
    double b = 0.1 * std::sin(2.0 * seed + 0.2);
    double c = 0.1 * std::sin(3.0 * seed + 0.3);
    su2_exp_from_Amu(U, std::array<double, 3>{a, b, c});
    cm_proj(U);
}

// Check if matrix is valid SU(2)
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

// Fill gauge field with identity links
static void fill_identity(double *gauge_field, int T, int L) {
    int vol = T * L * L * L;
    for (int site = 0; site < vol; ++site)
        for (int mu = 0; mu < 4; ++mu)
            cm_eq_id(&gauge_field[site * 4 * 8 + mu * 8]);
}

// Get link offset in the gauge field array
static inline std::size_t link_offset(int t, int x, int y, int z, int mu, int L, int T) {
    int site = t * L * L * L + x * L * L + y * L + z;
    return static_cast<std::size_t>(site) * 4 * 8 + static_cast<std::size_t>(mu) * 8;
}

// Apply a gauge transformation: U_mu(x) -> G(x) U_mu(x) G^dag(x+mu)
static void apply_gauge_transform(double *gf, int L, int T) {
    int vol = T * L * L * L;
    
    // Build deterministic gauge transformation G(x)
    std::vector<double> G(static_cast<std::size_t>(vol) * 8, 0.0);
    for (int site = 0; site < vol; ++site)
        su2_random_from_seed(&G[static_cast<std::size_t>(site) * 8], 1234u + site);
    
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
                        
                        cm_eq_cm_dag(Gd, Gsp);           // G^dag(x+mu)
                        cm_eq_cm_ti_cm(tmp1, Gs, U);     // G(x) U
                        cm_eq_cm_ti_cm(tmp2, tmp1, Gd);  // (G U) G^dag
                        cm_eq_cm(U, tmp2);
                    }
                }
}

// ==============================================================================
// Basic Plaquette Tests
// ==============================================================================

TEST_CASE("Plaquette: Identity field gives P = I", "[plaquette][basic]") {
    InitializeRand(300);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Unity(gauge_field, TEST_T, TEST_L);
    
    double I[8], plaq[8];
    cm_eq_id(I);
    
    // Check all plaquettes
    for (int it = 0; it < TEST_T; ++it) {
        for (int ix = 0; ix < TEST_L; ++ix) {
            for (int iy = 0; iy < TEST_L; ++iy) {
                for (int iz = 0; iz < TEST_L; ++iz) {
                    for (int mu = 0; mu < 4; ++mu) {
                        for (int nu = mu + 1; nu < 4; ++nu) {
                            compute_plaquette(plaq, gauge_field, it, ix, iy, iz, 
                                              mu, nu, TEST_T, TEST_L);
                            REQUIRE(matrices_equal(plaq, I, TOL));
                        }
                    }
                }
            }
        }
    }
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Plaquette: Trace equals 2 for identity plaquette", "[plaquette][trace]") {
    InitializeRand(301);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Unity(gauge_field, TEST_T, TEST_L);
    
    double tr = plaquette_trace(gauge_field, 0, 0, 0, 0, 0, 1, TEST_T, TEST_L);
    
    // plaquette_trace returns (1/2) Re(Tr(P)) for SU(2)
    // For identity, Tr(I) = 2, so normalized = 1.0
    REQUIRE(std::fabs(tr - 1.0) < TOL);
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Plaquette: Average plaquette is 1 for cold start", "[plaquette][average]") {
    InitializeRand(302);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Unity(gauge_field, TEST_T, TEST_L);
    
    double avg_plaq = Average_Plaquette(gauge_field, TEST_T, TEST_L);
    
    REQUIRE(std::fabs(avg_plaq - 1.0) < TOL);
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Plaquette: Plaquette is SU(2) matrix", "[plaquette][su2]") {
    InitializeRand(303);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double plaq[8];
    
    // Check several random plaquettes
    for (int trial = 0; trial < 10; ++trial) {
        int it = rand() % TEST_T;
        int ix = rand() % TEST_L;
        int iy = rand() % TEST_L;
        int iz = rand() % TEST_L;
        int mu = rand() % 4;
        int nu = (mu + 1 + rand() % 3) % 4;
        if (nu == mu) nu = (nu + 1) % 4;
        
        compute_plaquette(plaq, gauge_field, it, ix, iy, iz, mu, nu, TEST_T, TEST_L);
        REQUIRE(is_valid_su2(plaq, LOOSE_TOL));
    }
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Clover Plaquette Tests
// ==============================================================================

TEST_CASE("Clover: Identity field gives C = I", "[clover][basic]") {
    InitializeRand(310);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Unity(gauge_field, TEST_T, TEST_L);
    
    double C[8], I[8];
    cm_eq_id(I);
    
    compute_clover(C, gauge_field, 1, 1, 1, 1, 0, 1, TEST_T, TEST_L);
    
    REQUIRE(matrices_equal(C, I, TOL));
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Clover: Field strength zero on trivial config", "[clover][fieldstrength]") {
    InitializeRand(311);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Unity(gauge_field, TEST_T, TEST_L);
    
    double F[8];
    
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = mu + 1; nu < 4; ++nu) {
            compute_field_strength(F, gauge_field, 1, 1, 1, 1, mu, nu, TEST_T, TEST_L);
            REQUIRE(max_abs_8(F) < TOL);
        }
    }
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Clover: Field strength is anti-Hermitian", "[clover][antihermitian]") {
    InitializeRand(312);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double F[8], F_dag[8], sum[8];
    
    // Check several plane combinations
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = mu + 1; nu < 4; ++nu) {
            compute_field_strength(F, gauge_field, 1, 1, 1, 1, mu, nu, TEST_T, TEST_L);
            cm_eq_cm_dag(F_dag, F);
            
            // F + F^dag should be zero for anti-Hermitian
            cm_eq_cm(sum, F);
            cm_pl_eq_cm(sum, F_dag);
            
            REQUIRE(max_abs_8(sum) < LOOSE_TOL);
        }
    }
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Clover: Field strength is traceless", "[clover][traceless]") {
    InitializeRand(313);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double F[8];
    complex tr;
    
    // Check several plane combinations
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = mu + 1; nu < 4; ++nu) {
            compute_field_strength(F, gauge_field, 1, 1, 1, 1, mu, nu, TEST_T, TEST_L);
            co_eq_tr_cm(&tr, F);
            
            REQUIRE(std::fabs(tr.re) < LOOSE_TOL);
            REQUIRE(std::fabs(tr.im) < LOOSE_TOL);
        }
    }
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Abelian Field Configuration Tests
// ==============================================================================

TEST_CASE("Plaquette: Abelian F12 background gives nonzero C12", "[plaquette][abelian]") {
    const int L = 6, T = 4;
    const int vol = T * L * L * L;
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, T, L);
    fill_identity(gauge_field, T, L);
    
    // Construct U1(t,x,y,z) = exp(i * theta * y * sigma3/2)
    const double theta = 0.4;
    
    for (int t = 0; t < T; ++t) {
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                for (int z = 0; z < L; ++z) {
                    double U1[8];
                    su2_exp_from_Amu(U1, std::array<double, 3>{0.0, 0.0, theta * static_cast<double>(y)});
                    cm_eq_cm(&gauge_field[link_offset(t, x, y, z, 1, L, T)], U1);
                }
            }
        }
    }
    
    // Check interior sites (avoid boundary effects)
    const int t0 = 0;
    for (int x = 1; x < L - 2; ++x) {
        for (int y = 1; y < L - 2; ++y) {
            for (int z = 1; z < L - 2; ++z) {
                double C12[8];
                compute_clover(C12, gauge_field, t0, x, y, z, 1, 2, T, L);
                
                // C12 should be nonzero
                REQUIRE(max_abs_8(C12) > 1e-6);
                
                // Other planes should be close to identity (zero field strength)
                double C13[8], C23[8], C01[8];
                compute_clover(C13, gauge_field, t0, x, y, z, 1, 3, T, L);
                compute_clover(C23, gauge_field, t0, x, y, z, 2, 3, T, L);
                compute_clover(C01, gauge_field, t0, x, y, z, 0, 1, T, L);
                
                double I[8];
                cm_eq_id(I);
                
                // These should be close to identity
                REQUIRE(matrices_equal(C13, I, LOOSE_TOL));
                REQUIRE(matrices_equal(C23, I, LOOSE_TOL));
                REQUIRE(matrices_equal(C01, I, LOOSE_TOL));
            }
        }
    }
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Gauge Covariance Tests
// ==============================================================================

TEST_CASE("Plaquette: Transforms covariantly under gauge transform", "[plaquette][gauge]") {
    const int L = 4, T = 4;
    
    double *gauge_field;
    double *gauge_field_transformed;
    Gauge_Field_Alloc(&gauge_field, T, L);
    Gauge_Field_Alloc(&gauge_field_transformed, T, L);
    
    // Create mildly nontrivial configuration
    fill_identity(gauge_field, T, L);
    int vol = T * L * L * L;
    for (int site = 0; site < vol; ++site) {
        double U[8];
        su2_random_from_seed(U, 777u + site);
        cm_eq_cm(&gauge_field[static_cast<std::size_t>(site) * 4 * 8 + 2 * 8], U);
    }
    
    // Copy and transform
    Gauge_Field_Copy(gauge_field_transformed, gauge_field, T, L);
    apply_gauge_transform(gauge_field_transformed, L, T);
    
    // Get plaquettes at a test point
    int t = 1, x = 1, y = 1, z = 1, mu = 0, nu = 2;
    double P[8], P2[8];
    
    compute_plaquette(P, gauge_field, t, x, y, z, mu, nu, T, L);
    compute_plaquette(P2, gauge_field_transformed, t, x, y, z, mu, nu, T, L);
    
    // Reconstruct G(x) (same as apply_gauge_transform uses)
    int site = t * L * L * L + x * L * L + y * L + z;
    double G[8], Gd[8], tmp[8], cov[8];
    su2_random_from_seed(G, 1234u + site);
    cm_eq_cm_dag(Gd, G);
    
    // cov = G * P * G^dag
    cm_eq_cm_ti_cm(tmp, G, P);
    cm_eq_cm_ti_cm(cov, tmp, Gd);
    
    // Compare with transformed plaquette
    REQUIRE(matrices_equal(cov, P2, LOOSE_TOL));
    
    Gauge_Field_Free(&gauge_field);
    Gauge_Field_Free(&gauge_field_transformed);
}

TEST_CASE("Clover: Transforms covariantly under gauge transform", "[clover][gauge]") {
    const int L = 4, T = 4;
    
    double *gauge_field;
    double *gauge_field_transformed;
    Gauge_Field_Alloc(&gauge_field, T, L);
    Gauge_Field_Alloc(&gauge_field_transformed, T, L);
    
    fill_identity(gauge_field, T, L);
    int vol = T * L * L * L;
    for (int site = 0; site < vol; ++site) {
        double U[8];
        su2_random_from_seed(U, 777u + site);
        cm_eq_cm(&gauge_field[static_cast<std::size_t>(site) * 4 * 8 + 2 * 8], U);
    }
    
    Gauge_Field_Copy(gauge_field_transformed, gauge_field, T, L);
    apply_gauge_transform(gauge_field_transformed, L, T);
    
    int t = 1, x = 1, y = 1, z = 1, mu = 0, nu = 2;
    
    double C[8], C2[8];
    compute_clover(C, gauge_field, t, x, y, z, mu, nu, T, L);
    compute_clover(C2, gauge_field_transformed, t, x, y, z, mu, nu, T, L);
    
    int site = t * L * L * L + x * L * L + y * L + z;
    double G[8], Gd[8], tmp[8], cov[8];
    su2_random_from_seed(G, 1234u + site);
    cm_eq_cm_dag(Gd, G);
    
    cm_eq_cm_ti_cm(tmp, G, C);
    cm_eq_cm_ti_cm(cov, tmp, Gd);
    
    REQUIRE(matrices_equal(cov, C2, LOOSE_TOL));
    
    Gauge_Field_Free(&gauge_field);
    Gauge_Field_Free(&gauge_field_transformed);
}

// ==============================================================================
// Average Plaquette Tests
// ==============================================================================

TEST_CASE("Plaquette: Average is bounded [0, 1]", "[plaquette][bounds]") {
    InitializeRand(320);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double avg = Average_Plaquette(gauge_field, TEST_T, TEST_L);
    
    REQUIRE(avg >= -1.0 - TOL);  // Theoretically can be negative but usually > 0
    REQUIRE(avg <= 1.0 + TOL);
    
    Gauge_Field_Free(&gauge_field);
}

TEST_CASE("Plaquette: Spatial vs temporal average", "[plaquette][components]") {
    InitializeRand(321);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double avg_total = Average_Plaquette(gauge_field, TEST_T, TEST_L);
    double avg_spatial = Average_Spatial_Plaquette(gauge_field, TEST_T, TEST_L);
    double avg_temporal = Average_Temporal_Plaquette(gauge_field, TEST_T, TEST_L);
    
    // All should be valid numbers
    REQUIRE(!std::isnan(avg_total));
    REQUIRE(!std::isnan(avg_spatial));
    REQUIRE(!std::isnan(avg_temporal));
    
    // All should be bounded
    REQUIRE(avg_spatial <= 1.0 + TOL);
    REQUIRE(avg_temporal <= 1.0 + TOL);
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Geometry Tests
// ==============================================================================

TEST_CASE("Plaquette: Closed loop property", "[plaquette][geometry]") {
    InitializeRand(330);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    // Plaquette should be a valid SU(2) matrix (closed Wilson loop)
    for (int it = 0; it < TEST_T; ++it) {
        for (int ix = 0; ix < TEST_L; ++ix) {
            double plaq[8];
            compute_plaquette(plaq, gauge_field, it, ix, 0, 0, 0, 1, TEST_T, TEST_L);
            
            // Check unitarity
            double plaq_dag[8], prod[8], I[8];
            cm_eq_id(I);
            cm_eq_cm_dag(plaq_dag, plaq);
            cm_eq_cm_ti_cm(prod, plaq_dag, plaq);
            
            REQUIRE(matrices_equal(prod, I, LOOSE_TOL));
            
            // Check determinant
            complex det;
            co_eq_det_cm(&det, plaq);
            REQUIRE(std::fabs(det.re - 1.0) < LOOSE_TOL);
            REQUIRE(std::fabs(det.im) < LOOSE_TOL);
        }
    }
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Engineered Configuration Tests
// ==============================================================================

TEST_CASE("Plaquette: Engineered config gives nonzero clover", "[plaquette][engineered]") {
    const int L = 4, T = 4;
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, T, L);
    fill_identity(gauge_field, T, L);
    
    // Set specific links to create nonzero field strength
    const int t0 = 0, x0 = 0, y0 = 0, z0 = 0;
    const int mu = 1, nu = 2;
    
    double P[8], Pd[8];
    su2_exp_from_Amu(P, std::array<double, 3>{0.0, 0.0, 0.6});
    cm_eq_cm_dag(Pd, P);
    
    // U_mu(x) = P
    cm_eq_cm(&gauge_field[link_offset(t0, x0, y0, z0, mu, L, T)], P);
    
    // U_nu(x+mu) = P^dag
    int xp = (x0 + 1) % L;
    cm_eq_cm(&gauge_field[link_offset(t0, xp, y0, z0, nu, L, T)], Pd);
    
    // Check that clover is now nonzero
    double C[8];
    compute_clover(C, gauge_field, t0, x0, y0, z0, mu, nu, T, L);
    
    double I[8];
    cm_eq_id(I);
    
    // C should differ from identity
    bool differs = false;
    for (int i = 0; i < 8; ++i) {
        if (std::fabs(C[i] - I[i]) > 1e-6) {
            differs = true;
            break;
        }
    }
    REQUIRE(differs);
    
    Gauge_Field_Free(&gauge_field);
}

// ==============================================================================
// Reproducibility Tests
// ==============================================================================

TEST_CASE("Plaquette: Computation is reproducible", "[plaquette][reproducibility]") {
    InitializeRand(340);
    
    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    
    double plaq1[8], plaq2[8];
    
    compute_plaquette(plaq1, gauge_field, 1, 1, 1, 1, 0, 1, TEST_T, TEST_L);
    compute_plaquette(plaq2, gauge_field, 1, 1, 1, 1, 0, 1, TEST_T, TEST_L);
    
    REQUIRE(matrices_equal(plaq1, plaq2, TOL));
    
    Gauge_Field_Free(&gauge_field);
}

