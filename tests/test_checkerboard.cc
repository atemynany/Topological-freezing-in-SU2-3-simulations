#include <catch2/catch_all.hpp>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fields.hh"
#include "geometry.hh"
#include "linear_algebra.hh"
#include "ranlux.hh"

extern "C" {
#include "ranlxd.h"
}

// globals required by the utility library
int T = 4;
int L = 4;
double *gauge_field = nullptr;
bool open_boundary_conditions = false;
bool hot_start = false;

TEST_CASE("Thread-local RNG: each thread gets independent stream", "[rng][openmp]") {
#ifdef _OPENMP
    const int n_threads = 4;
    std::vector<double> first_values(n_threads, 0.0);

    #pragma omp parallel num_threads(n_threads)
    {
        int tid = omp_get_thread_num();
        rlxd_init(2, 1000 + tid);
        first_values[tid] = DRand();
    }

    // all threads should produce different first values
    for (int i = 0; i < n_threads; i++) {
        for (int j = i + 1; j < n_threads; j++) {
            REQUIRE(first_values[i] != first_values[j]);
        }
    }
#else
    // without OpenMP, just verify RNG works single-threaded
    rlxd_init(2, 1000);
    double v1 = DRand();
    rlxd_init(2, 1001);
    double v2 = DRand();
    REQUIRE(v1 != v2);
#endif
}

TEST_CASE("Checkerboard: even/odd site classification", "[checkerboard]") {
    T = 8;
    L = 8;
    open_boundary_conditions = false;
    const int volume = T * L * L * L;

    std::vector<int> ev, od;
    for (int site = 0; site < volume; site++) {
        int iz = site % L;
        int iy = (site / L) % L;
        int ix = (site / (L * L)) % L;
        int it = site / (L * L * L);
        int parity = (it + ix + iy + iz) % 2;
        if (parity == 0)
            ev.push_back(site);
        else
            od.push_back(site);
    }

    SECTION("Equal split") {
        REQUIRE(ev.size() == (size_t)(volume / 2));
        REQUIRE(od.size() == (size_t)(volume / 2));
    }

    SECTION("Every even site has only odd neighbors") {
        for (int site : ev) {
            int iz = site % L;
            int iy = (site / L) % L;
            int ix = (site / (L * L)) % L;
            int it = site / (L * L * L);
            int coords[4] = {it, ix, iy, iz};
            int sizes[4] = {T, L, L, L};
            for (int mu = 0; mu < 4; mu++) {
                // +mu neighbor
                int c[4] = {coords[0], coords[1], coords[2], coords[3]};
                c[mu] = (c[mu] + 1) % sizes[mu];
                REQUIRE((c[0] + c[1] + c[2] + c[3]) % 2 == 1);
                // -mu neighbor
                c[0] = coords[0]; c[1] = coords[1]; c[2] = coords[2]; c[3] = coords[3];
                c[mu] = (c[mu] - 1 + sizes[mu]) % sizes[mu];
                REQUIRE((c[0] + c[1] + c[2] + c[3]) % 2 == 1);
            }
        }
    }
}
