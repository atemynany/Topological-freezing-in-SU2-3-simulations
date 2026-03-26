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
