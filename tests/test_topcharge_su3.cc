// ==============================================================================
// test_topcharge_su3.cc
// ==============================================================================
// Unit tests for SU(3) topological charge measurement.
//
// Strategy: pin the SU(3) clover against three independent ground truths —
//   1. Cold configuration          → Q = 0 exactly.
//   2. SU(2) BPST instanton embedded in the upper-left 2x2 block of SU(3)
//      → Q ≈ +1  (anti ≈ -1). The remaining 3-direction is the identity,
//        so the SU(3) topcharge must equal the SU(2) topcharge for the
//        same instanton insertion.
//   3. Random SU(3) gauge transformation U_μ(x) → G(x) U_μ(x) G†(x+μ)
//      → Q unchanged (gauge invariance).
//
// Author: Alexander de Barros Noll
// Date: 2026-04-17
// ==============================================================================

#include <catch2/catch_all.hpp>
#include <array>
#include <cmath>
#include <vector>

#include "fields.hh"
#include "geometry.hh"
#include "linear_algebra.hh"
#include "ranlux.hh"
#include "topcharge_su2.hh"
#include "topcharge_su3.hh"
#include "su3_linear_algebra.hh"
#include "instanton.hh"

// ==============================================================================
// Globals required by linked translation units
// ==============================================================================
// test_topcharge.cc, fields.cc and the SU(3) header each reach for different
// externs; declare them here once so a single test binary links cleanly.
bool open_boundary_conditions = false;
int T_size = 0;
int L_size = 0;

// ==============================================================================
// Constants
// ==============================================================================

constexpr double TOL = 1e-10;
constexpr double LOOSE_TOL = 1e-6;

// ==============================================================================
// Helpers
// ==============================================================================

// Fill SU(3) gauge field with identity links (cold configuration).
static void su3_fill_identity(double *gf, int T, int L) {
    const long long vol = static_cast<long long>(T) * L * L * L;
    for (long long site = 0; site < vol; ++site)
        for (int mu = 0; mu < 4; ++mu)
            su3_eq_id(gf + ((4 * site + mu) * 18));
}

// Embed a 2x2 SU(2) link (8 doubles, row-major {Re,Im}) into the upper-left
// block of an SU(3) link (18 doubles, row-major {Re,Im}); (2,2) entry = 1.
static void su3_embed_su2_link(double *U3, const double *U2) {
    su3_eq_zero(U3);
    // Row 0: U2[0..3] → U3[0..3]
    U3[0] = U2[0]; U3[1] = U2[1];
    U3[2] = U2[2]; U3[3] = U2[3];
    // Row 1: U2[4..7] → U3[6..9]
    U3[6] = U2[4]; U3[7] = U2[5];
    U3[8] = U2[6]; U3[9] = U2[7];
    // Row 2: (0, 0, 1)
    U3[16] = 1.0;
}

// Build an SU(3) gauge field by embedding an SU(2) BPST instanton background.
static void build_embedded_instanton(std::vector<double> &gf3, double rho,
                                     double a, int L, int T, bool is_anti) {
    const long long vol = static_cast<long long>(T) * L * L * L;

    std::vector<double> gf2(static_cast<std::size_t>(vol) * 4 * 8, 0.0);
    insert_instanton(gf2.data(), rho, a, L, T, is_anti);

    gf3.assign(static_cast<std::size_t>(vol) * 4 * 18, 0.0);
    for (long long site = 0; site < vol; ++site)
        for (int mu = 0; mu < 4; ++mu)
            su3_embed_su2_link(gf3.data() + ((4 * site + mu) * 18),
                               gf2.data() + ((4 * site + mu) * 8));
}

// Deterministic "random" SU(3) matrix from seed (Gram-Schmidt projected).
static void su3_from_seed(double *U, unsigned seed) {
    for (int i = 0; i < 18; ++i) {
        double s = static_cast<double>(seed) + 0.137 * (i + 1);
        U[i] = 0.5 * std::sin(1.9 * s + 0.3 * i);
    }
    su3_proj(U);
}

// Apply U_μ(x) → G(x) U_μ(x) G†(x+μ) on the SU(3) field, periodic shifts.
static void apply_gauge_transform_su3(double *gf, int T, int L) {
    const long long vol = static_cast<long long>(T) * L * L * L;

    std::vector<double> G(static_cast<std::size_t>(vol) * 18, 0.0);
    for (long long s = 0; s < vol; ++s)
        su3_from_seed(G.data() + s * 18, static_cast<unsigned>(8191u + s));

    auto site_index = [&](int t, int x, int y, int z) {
        return ((long long)t * L + x) * L * L + (long long)y * L + z;
    };

    auto forward = [&](int t, int x, int y, int z, int mu) {
        if (mu == 0) t = (t + 1) % T;
        if (mu == 1) x = (x + 1) % L;
        if (mu == 2) y = (y + 1) % L;
        if (mu == 3) z = (z + 1) % L;
        return site_index(t, x, y, z);
    };

    alignas(32) double tmp1[18], tmp2[18], Gd[18];

    for (int t = 0; t < T; ++t)
        for (int x = 0; x < L; ++x)
            for (int y = 0; y < L; ++y)
                for (int z = 0; z < L; ++z) {
                    const long long s = site_index(t, x, y, z);
                    for (int mu = 0; mu < 4; ++mu) {
                        const long long sp = forward(t, x, y, z, mu);
                        double *U = gf + ((4 * s + mu) * 18);
                        const double *Gs  = G.data() + s  * 18;
                        const double *Gsp = G.data() + sp * 18;

                        su3_eq_su3_dag(Gd, Gsp);
                        su3_eq_su3_ti_su3(tmp1, Gs, U);
                        su3_eq_su3_ti_su3(tmp2, tmp1, Gd);
                        su3_eq_su3(U, tmp2);
                    }
                }
}

// ==============================================================================
// 1. Cold configuration
// ==============================================================================

TEST_CASE("SU(3) Topological Charge: Q=0 on cold configuration", "[topcharge_su3][cold]") {
    const int T = 4, L = 4;
    T_size = T; L_size = L;

    std::vector<double> gf(static_cast<std::size_t>(T) * L * L * L * 4 * 18, 0.0);
    su3_fill_identity(gf.data(), T, L);

    double Q = su3_topological_charge(gf.data(), T, L);
    REQUIRE(std::fabs(Q) < TOL);
}

// ==============================================================================
// 2. Embedded SU(2) instanton — the decisive normalization check
// ==============================================================================

TEST_CASE("SU(3) Topological Charge: Embedded BPST instanton gives Q ~ +1",
          "[topcharge_su3][instanton]") {
    const int T = 14, L = 14;
    const double rho = 3.5, a = 1.0;
    T_size = T; L_size = L;

    std::vector<double> gf3;
    build_embedded_instanton(gf3, rho, a, L, T, /*is_anti=*/false);

    const double Q3 = su3_topological_charge(gf3.data(), T, L);
    INFO("SU(3) embedded instanton Q = " << Q3);
    REQUIRE(std::fabs(Q3 - 1.0) < 0.1);
}

TEST_CASE("SU(3) Topological Charge: Embedded anti-instanton gives Q ~ -1",
          "[topcharge_su3][instanton]") {
    const int T = 14, L = 14;
    const double rho = 3.5, a = 1.0;
    T_size = T; L_size = L;

    std::vector<double> gf3;
    build_embedded_instanton(gf3, rho, a, L, T, /*is_anti=*/true);

    const double Q3 = su3_topological_charge(gf3.data(), T, L);
    INFO("SU(3) embedded anti-instanton Q = " << Q3);
    REQUIRE(std::fabs(Q3 + 1.0) < 0.1);
}

TEST_CASE("SU(3) Topological Charge: Embedded instanton matches SU(2) Q exactly",
          "[topcharge_su3][instanton][cross_check]") {
    const int T = 12, L = 12;
    const double rho = 3.0, a = 1.0;
    const long long vol = static_cast<long long>(T) * L * L * L;
    T_size = T; L_size = L;

    std::vector<double> gf2(static_cast<std::size_t>(vol) * 4 * 8, 0.0);
    insert_instanton(gf2.data(), rho, a, L, T, /*is_anti=*/false);
    const double Q2 = compute_topological_charge(gf2.data(), T, L);

    std::vector<double> gf3(static_cast<std::size_t>(vol) * 4 * 18, 0.0);
    for (long long s = 0; s < vol; ++s)
        for (int mu = 0; mu < 4; ++mu)
            su3_embed_su2_link(gf3.data() + ((4 * s + mu) * 18),
                               gf2.data() + ((4 * s + mu) * 8));
    const double Q3 = su3_topological_charge(gf3.data(), T, L);

    INFO("Q_SU(2) = " << Q2 << "   Q_SU(3) = " << Q3);
    REQUIRE(std::fabs(Q3 - Q2) < LOOSE_TOL);
}

// ==============================================================================
// 3. Gauge invariance
// ==============================================================================

TEST_CASE("SU(3) Topological Charge: Invariant under SU(3) gauge transform",
          "[topcharge_su3][gauge_invariance]") {
    const int T = 6, L = 6;
    const long long vol = static_cast<long long>(T) * L * L * L;
    T_size = T; L_size = L;

    // Build a mildly non-trivial background: embed a small-ρ instanton so
    // Q is measurably nonzero yet the random gauge transform exercises all
    // three color directions.
    std::vector<double> gf(static_cast<std::size_t>(vol) * 4 * 18, 0.0);
    build_embedded_instanton(gf, /*rho=*/1.8, /*a=*/1.0, L, T, /*is_anti=*/false);

    std::vector<double> gf_transformed = gf;
    apply_gauge_transform_su3(gf_transformed.data(), T, L);

    const double Q_orig = su3_topological_charge(gf.data(), T, L);
    const double Q_tr   = su3_topological_charge(gf_transformed.data(), T, L);

    INFO("Q_original = " << Q_orig << "   Q_transformed = " << Q_tr);
    REQUIRE(std::fabs(Q_orig - Q_tr) < LOOSE_TOL);
}
