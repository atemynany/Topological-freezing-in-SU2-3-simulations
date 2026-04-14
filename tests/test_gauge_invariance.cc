// ==============================================================================
// test_gauge_invariance.cc
// ==============================================================================
// Tests that physical observables (plaquette, topological charge, topological
// charge after APE smearing) are invariant under random SU(2) gauge
// transformations. Adapted from gauge_transformation.cc by H. Ahlborn.
// ==============================================================================

#include <catch2/catch_all.hpp>
#include <cmath>
#include <vector>

#include "fields.hh"
#include "geometry.hh"
#include "linear_algebra.hh"
#include "ranlux.hh"
#include "Plaquette.hh"
#include "topcharge_su2.hh"
#include "smearing_techniques.hh"

// ==============================================================================
// Globals required by geometry/fields
// ==============================================================================

bool open_boundary_conditions = false;

// ==============================================================================
// Gauge transformation helpers (from gauge_transformation.cc)
// ==============================================================================

// g_x: one SU(2) matrix per site, stored as vol*8 doubles
static void Generate_Random_Gauge_Transform(double *g_x, int T, int L) {
    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    int site = get_index(it, ix, iy, iz, T, L);
                    double h[4];
                    double norm;
                    while (true) {
                        h[0] = 2.0 * DRand() - 1.0;
                        h[1] = 2.0 * DRand() - 1.0;
                        h[2] = 2.0 * DRand() - 1.0;
                        h[3] = 2.0 * DRand() - 1.0;
                        norm = h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3];
                        if (norm >= 1e-10 && norm <= 1.0)
                            break;
                    }
                    norm = 1.0 / std::sqrt(norm);
                    for (int i = 0; i < 4; i++) h[i] *= norm;
                    cm_from_h(g_x + site * 8, h);
                }
            }
        }
    }
}

// U'_mu(x) = g(x) * U_mu(x) * g^dag(x+mu)
static void Apply_Gauge_Transformation(double *gauge_field, double *g_x, int T, int L) {
    double transformed_link[8], temp[8];

    for (int it = 0; it < T; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    int site = get_index(it, ix, iy, iz, T, L);
                    for (int mu = 0; mu < 4; mu++) {

                        long long link_idx = ggi(site, mu);
                        if (link_idx == -1) continue;

                        int it_new = it, ix_new = ix, iy_new = iy, iz_new = iz;
                        switch (mu) {
                            case 0: it_new++; break;
                            case 1: ix_new++; break;
                            case 2: iy_new++; break;
                            case 3: iz_new++; break;
                        }

                        int site_fwd = get_index(it_new, ix_new, iy_new, iz_new, T, L);
                        if (site_fwd == -1) continue;

                        // U'_mu(x) = g(x) * U_mu(x) * g^dag(x+mu)
                        cm_eq_cm_ti_cm(temp, g_x + site * 8, gauge_field + link_idx);
                        cm_eq_cm_ti_cm_dag(transformed_link, temp, g_x + site_fwd * 8);
                        cm_eq_cm(gauge_field + link_idx, transformed_link);
                    }
                }
            }
        }
    }
}

// ==============================================================================
// Tests
// ==============================================================================

constexpr int TEST_T = 8;
constexpr int TEST_L = 8;
constexpr double GAUGE_TOL = 1e-8;

static double *alloc_gauge_transform(int T, int L) {
    const int vol = T * L * L * L;
    double *g = static_cast<double *>(std::malloc(vol * 8 * sizeof(double)));
    return g;
}

TEST_CASE("Gauge Invariance: Average plaquette", "[gauge_invariance][plaquette]") {
    InitializeRand(42);

    double *gauge_field;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    double *g_x = alloc_gauge_transform(TEST_T, TEST_L);

    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);

    double P_before = Average_Plaquette(gauge_field, TEST_T, TEST_L);

    Generate_Random_Gauge_Transform(g_x, TEST_T, TEST_L);
    Apply_Gauge_Transformation(gauge_field, g_x, TEST_T, TEST_L);

    double P_after = Average_Plaquette(gauge_field, TEST_T, TEST_L);

    INFO("P_before = " << P_before);
    INFO("P_after  = " << P_after);
    REQUIRE(std::fabs(P_before - P_after) < GAUGE_TOL);

    Gauge_Field_Free(&gauge_field);
    std::free(g_x);
}

TEST_CASE("Gauge Invariance: Topological charge (no smearing)", "[gauge_invariance][topcharge]") {
    InitializeRand(43);

    double *gauge_field;
    double *gauge_field_transformed;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Alloc(&gauge_field_transformed, TEST_T, TEST_L);
    double *g_x = alloc_gauge_transform(TEST_T, TEST_L);

    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    Gauge_Field_Copy(gauge_field_transformed, gauge_field, TEST_T, TEST_L);

    double Q_before = compute_topological_charge(gauge_field, TEST_T, TEST_L);

    Generate_Random_Gauge_Transform(g_x, TEST_T, TEST_L);
    Apply_Gauge_Transformation(gauge_field_transformed, g_x, TEST_T, TEST_L);

    double Q_after = compute_topological_charge(gauge_field_transformed, TEST_T, TEST_L);

    INFO("Q_before = " << Q_before);
    INFO("Q_after  = " << Q_after);
    REQUIRE(std::fabs(Q_before - Q_after) < GAUGE_TOL);

    Gauge_Field_Free(&gauge_field);
    Gauge_Field_Free(&gauge_field_transformed);
    std::free(g_x);
}

TEST_CASE("Gauge Invariance: Topological charge after APE smearing", "[gauge_invariance][topcharge][smearing]") {
    InitializeRand(44);

    const int smear_steps = 10;
    const double smear_alpha = 0.5;

    double *gauge_field;
    double *gauge_field_transformed;
    double *smeared1;
    double *smeared2;
    Gauge_Field_Alloc(&gauge_field, TEST_T, TEST_L);
    Gauge_Field_Alloc(&gauge_field_transformed, TEST_T, TEST_L);
    Gauge_Field_Alloc(&smeared1, TEST_T, TEST_L);
    Gauge_Field_Alloc(&smeared2, TEST_T, TEST_L);
    double *g_x = alloc_gauge_transform(TEST_T, TEST_L);

    Gauge_Field_Random(gauge_field, TEST_T, TEST_L);
    Gauge_Field_Copy(gauge_field_transformed, gauge_field, TEST_T, TEST_L);

    // Smear original and measure Q
    Gauge_Field_Copy(smeared1, gauge_field, TEST_T, TEST_L);
    for (int i = 0; i < smear_steps; i++) {
        APE_Smearing_all(smeared1, TEST_T, TEST_L, smear_alpha);
    }
    double Q_smear_before = compute_topological_charge(smeared1, TEST_T, TEST_L);

    // Gauge transform, then smear and measure Q
    Generate_Random_Gauge_Transform(g_x, TEST_T, TEST_L);
    Apply_Gauge_Transformation(gauge_field_transformed, g_x, TEST_T, TEST_L);

    Gauge_Field_Copy(smeared2, gauge_field_transformed, TEST_T, TEST_L);
    for (int i = 0; i < smear_steps; i++) {
        APE_Smearing_all(smeared2, TEST_T, TEST_L, smear_alpha);
    }
    double Q_smear_after = compute_topological_charge(smeared2, TEST_T, TEST_L);

    INFO("Q_smear_before = " << Q_smear_before);
    INFO("Q_smear_after  = " << Q_smear_after);
    REQUIRE(std::fabs(Q_smear_before - Q_smear_after) < GAUGE_TOL);

    Gauge_Field_Free(&gauge_field);
    Gauge_Field_Free(&gauge_field_transformed);
    Gauge_Field_Free(&smeared1);
    Gauge_Field_Free(&smeared2);
    std::free(g_x);
}
