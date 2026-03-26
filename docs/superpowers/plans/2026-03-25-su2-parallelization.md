# SU(2) Parallelization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Parallelize SU(2) heatbath sweeps via checkerboard OpenMP decomposition and enable running multiple (beta, seed) jobs simultaneously on a SLURM cluster.

**Architecture:** Checkerboard even/odd site decomposition splits each sweep into two parallel half-sweeps. RANLUX RNG is made thread-safe via `_Thread_local` storage. SLURM array jobs handle inter-run parallelism with per-job error handling.

**Tech Stack:** C11/C++17, OpenMP, SLURM, Python 3 (numpy/matplotlib), Catch2

**Spec:** `docs/superpowers/specs/2026-03-25-su2-parallelization-design.md`

---

## File Structure

**Modified:**
- `_Utility/src/ranlxd.c` — Add `_Thread_local` to all static state variables (both SSE and non-SSE paths)
- `_Utility/src/ranlxs.c` — Same `_Thread_local` treatment
- `src/MC_heatbath.cc` — Checkerboard sweep with OpenMP parallel for + per-thread RNG init
- `CMakeLists.txt` — Add `test_checkerboard` target with OpenMP linking

**Created:**
- `tests/test_checkerboard.cc` — Checkerboard correctness tests
- `scripts/slurm_beta_scan.sh` — SLURM array job script
- `scripts/monitor_runs.py` — Live run monitoring with plots

---

### Task 1: Make RANLUX thread-safe

**Files:**
- Modify: `_Utility/src/ranlxd.c:374-382` (non-SSE path) and `:61-68` (SSE path)
- Modify: `_Utility/src/ranlxs.c:372-380` (non-SSE path) and `:61-68` (SSE path)

- [ ] **Step 1: Write a test that verifies thread-safety**

Create `tests/test_checkerboard.cc` with an initial test that spawns OpenMP threads, each calling `rlxd_init` + `DRand`, and verifies they get different sequences:

```cpp
#include <catch2/catch_all.hpp>
#include <omp.h>
#include <vector>
#include <cmath>

#include "fields.hh"
#include "geometry.hh"
#include "linear_algebra.hh"
#include "ranlux.hh"

// globals required by the utility library
int T = 4;
int L = 4;
double *gauge_field = nullptr;
bool open_boundary_conditions = false;
bool hot_start = false;

TEST_CASE("Thread-local RNG: each thread gets independent stream", "[rng][openmp]") {
    const int n_threads = 4;
    const int n_samples = 100;
    std::vector<double> first_values(n_threads, 0.0);

    #pragma omp parallel num_threads(n_threads)
    {
        int tid = omp_get_thread_num();
        // each thread seeds its own RANLUX stream
        rlxd_init(2, 1000 + tid);
        first_values[tid] = DRand();
    }

    // all threads should produce different first values
    for (int i = 0; i < n_threads; i++) {
        for (int j = i + 1; j < n_threads; j++) {
            REQUIRE(first_values[i] != first_values[j]);
        }
    }
}
```

- [ ] **Step 2: Add test_checkerboard to CMakeLists.txt**

After the `test_plaquette` block (around line 368), add:

```cmake
    # --- Test: Checkerboard parallelization ---
    add_executable(test_checkerboard
        tests/test_checkerboard.cc
    )

    target_include_directories(test_checkerboard
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/_Utility/include
    )

    target_link_libraries(test_checkerboard
        PRIVATE
            su2_utility
            Catch2::Catch2WithMain
            ${MATH_LIBRARY}
            $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
    )

    add_test(NAME test_checkerboard COMMAND test_checkerboard)
```

- [ ] **Step 3: Build and verify the RNG test fails (data race without _Thread_local)**

Run: `./scripts/build.sh release && cd build && ctest -R test_checkerboard --output-on-failure`

Expected: Test may fail or produce wrong results (threads clobbering shared state). This confirms the problem exists.

- [ ] **Step 4: Add `_Thread_local` to ranlxd.c static variables**

In `_Utility/src/ranlxd.c`, the non-SSE path (line 374-382 in `#else` block) has:

```c
static int init=0,pr,prm,ir,jr,is,is_old,next[96];
static double one_bit;
static vec_t carry;

static union
{
   dble_vec_t vec[12];
   int num[96];
} x;
```

Change to:

```c
// thread-local state for OpenMP checkerboard parallelization
static _Thread_local int init=0,pr,prm,ir,jr,is,is_old,next[96];
static _Thread_local double one_bit;
static _Thread_local vec_t carry;

static _Thread_local union
{
   dble_vec_t vec[12];
   int num[96];
} x;
```

The SSE path (line 61-68, inside `#if (defined SSE)`) has the same pattern with `vec_t`, `dble_vec_t`, `__attribute__((aligned))`. Apply the same `_Thread_local`:

```c
static _Thread_local int init=0,pr,prm,ir,jr,is,is_old,next[96];
static _Thread_local vec_t one,one_bit,carry;

static _Thread_local union
{
   dble_vec_t vec[12];
   float num[96];
} x __attribute__ ((aligned (16)));
```

- [ ] **Step 5: Add `_Thread_local` to ranlxs.c static variables**

Same treatment. Non-SSE path (`_Utility/src/ranlxs.c:372-380`):

```c
static _Thread_local int init=0,pr,prm,ir,jr,is,is_old,next[96];
static _Thread_local float one_bit;
static _Thread_local vec_t carry;

static _Thread_local union
{
   dble_vec_t vec[12];
   int num[96];
} x;
```

SSE path (`:61-68`):

```c
static _Thread_local int init=0,pr,prm,ir,jr,is,is_old,next[96];
static _Thread_local vec_t one,one_bit,carry;

static _Thread_local union
{
   dble_vec_t vec[12];
   float num[96];
} x __attribute__ ((aligned (16)));
```

- [ ] **Step 6: Build and verify the RNG test passes**

Run: `./scripts/build.sh release && cd build && ctest -R test_checkerboard --output-on-failure`

Expected: PASS — each thread now has its own RANLUX state.

- [ ] **Step 7: Run all existing tests to confirm no regression**

Run: `./scripts/run_tests.sh`

Expected: All 5 existing test suites pass. The single-threaded codepath is unchanged — `_Thread_local` gives each thread (including the main thread) its own copy of the static state, so sequential behavior is identical.

- [ ] **Step 8: Commit**

```bash
git add _Utility/src/ranlxd.c _Utility/src/ranlxs.c tests/test_checkerboard.cc CMakeLists.txt
git commit -m "feat: make RANLUX thread-safe with _Thread_local for OpenMP parallelization"
```

---

### Task 2: Checkerboard site lists in MC_heatbath.cc

**Files:**
- Modify: `src/MC_heatbath.cc`

- [ ] **Step 1: Add even/odd site classification test**

Append to `tests/test_checkerboard.cc`:

```cpp
TEST_CASE("Checkerboard: even/odd site classification", "[checkerboard]") {
    T = 8;
    L = 8;
    open_boundary_conditions = false;
    const int volume = T * L * L * L;

    std::vector<int> even_sites, odd_sites;
    for (int site = 0; site < volume; site++) {
        int iz = site % L;
        int iy = (site / L) % L;
        int ix = (site / (L * L)) % L;
        int it = site / (L * L * L);
        int parity = (it + ix + iy + iz) % 2;
        if (parity == 0)
            even_sites.push_back(site);
        else
            odd_sites.push_back(site);
    }

    SECTION("Equal split") {
        REQUIRE(even_sites.size() == volume / 2);
        REQUIRE(odd_sites.size() == volume / 2);
    }

    SECTION("Every even site has only odd neighbors") {
        for (int site : even_sites) {
            int iz = site % L;
            int iy = (site / L) % L;
            int ix = (site / (L * L)) % L;
            int it = site / (L * L * L);

            // check all 8 nearest neighbors
            int coords[4] = {it, ix, iy, iz};
            int sizes[4] = {T, L, L, L};
            for (int mu = 0; mu < 4; mu++) {
                int c[4] = {coords[0], coords[1], coords[2], coords[3]};
                c[mu] = (c[mu] + 1) % sizes[mu];
                int parity_plus = (c[0] + c[1] + c[2] + c[3]) % 2;
                REQUIRE(parity_plus == 1);

                c[mu] = (coords[mu] - 1 + sizes[mu]) % sizes[mu];
                c[0] = (mu == 0) ? c[mu] : coords[0];
                c[1] = (mu == 1) ? c[mu] : coords[1];
                c[2] = (mu == 2) ? c[mu] : coords[2];
                c[3] = (mu == 3) ? c[mu] : coords[3];
                int parity_minus = (c[0] + c[1] + c[2] + c[3]) % 2;
                REQUIRE(parity_minus == 1);
            }
        }
    }
}
```

- [ ] **Step 2: Run test to verify it passes**

Run: `./scripts/build.sh release && cd build && ctest -R test_checkerboard --output-on-failure`

Expected: PASS

- [ ] **Step 3: Add even/odd site arrays and initialization to MC_heatbath.cc**

After the `init_neighbor_tables` function (line 67) and before `struct SimParams` (line 69), add:

```cpp
// checkerboard site lists for parallel sweeps
std::vector<int> even_sites;
std::vector<int> odd_sites;

void init_checkerboard(int T_size, int L_size) {
    const int volume = T_size * L_size * L_size * L_size;
    even_sites.reserve(volume / 2);
    odd_sites.reserve(volume / 2);
    for (int site = 0; site < volume; site++) {
        int iz = site % L_size;
        int iy = (site / L_size) % L_size;
        int ix = (site / (L_size * L_size)) % L_size;
        int it = site / (L_size * L_size * L_size);
        if ((it + ix + iy + iz) % 2 == 0)
            even_sites.push_back(site);
        else
            odd_sites.push_back(site);
    }
}
```

In `main()`, after the call to `init_neighbor_tables(T, L)` (line 227), add:

```cpp
    init_checkerboard(T, L);
```

- [ ] **Step 4: Build and verify all tests still pass**

Run: `./scripts/build.sh release && ./scripts/run_tests.sh`

Expected: All tests pass. The checkerboard arrays are initialized but not yet used in the sweep loop.

- [ ] **Step 5: Commit**

```bash
git add src/MC_heatbath.cc tests/test_checkerboard.cc
git commit -m "feat: add checkerboard even/odd site classification"
```

---

### Task 3: Checkerboard parallel sweep

**Files:**
- Modify: `src/MC_heatbath.cc:267-371`

- [ ] **Step 1: Add test for checkerboard sweep producing valid SU(2) links**

Append to `tests/test_checkerboard.cc`:

```cpp
#include "Plaquette.hh"

TEST_CASE("Checkerboard sweep: all links remain valid SU(2)", "[checkerboard][heatbath]") {
    T = 8;
    L = 8;
    open_boundary_conditions = false;
    const int volume = T * L * L * L;

    double *test_field = nullptr;
    Gauge_Field_Alloc(&test_field, T, L);
    Gauge_Field_Unity(test_field, T, L);

    // build checkerboard
    std::vector<int> ev, od;
    for (int site = 0; site < volume; site++) {
        int iz = site % L;
        int iy = (site / L) % L;
        int ix = (site / (L * L)) % L;
        int it = site / (L * L * L);
        if ((it + ix + iy + iz) % 2 == 0) ev.push_back(site);
        else od.push_back(site);
    }

    double beta = 2.5;

    // do 10 checkerboard sweeps
    for (int sweep = 0; sweep < 10; sweep++) {
        // even half-sweep
        #pragma omp parallel
        {
            rlxd_init(2, 42 + omp_get_thread_num() + sweep * 100);
        }
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < (int)ev.size(); i++) {
            int site = ev[i];
            int iz = site % L;
            int iy = (site / L) % L;
            int ix = (site / (L * L)) % L;
            int it = site / (L * L * L);
            for (int mu = 0; mu < 4; mu++) {
                // inline heatbath update using test helper from test_heatbath.cc pattern
                double S_l[8];
                cm_eq_zero(S_l);
                double SU2_1[8], SU2_2[8];

                for (int nu = 0; nu < 4; nu++) {
                    if (nu == mu) continue;
                    int i_arr[4] = {it, ix, iy, iz};
                    int sizes[4] = {T, L, L, L};

                    // lower staple
                    i_arr[0]=it; i_arr[1]=ix; i_arr[2]=iy; i_arr[3]=iz;
                    i_arr[nu] = (i_arr[nu] - 1 + sizes[nu]) % sizes[nu];
                    int idx1 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), nu);
                    int idx2 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), mu);
                    i_arr[mu] = (i_arr[mu] + 1) % sizes[mu];
                    int idx3 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), nu);
                    if (idx1>=0 && idx2>=0 && idx3>=0) {
                        cm_eq_cm_ti_cm(SU2_1, test_field+idx2, test_field+idx3);
                        cm_eq_cm_dag_ti_cm(SU2_2, test_field+idx1, SU2_1);
                        cm_pl_eq_cm(S_l, SU2_2);
                    }

                    // upper staple
                    i_arr[0]=it; i_arr[1]=ix; i_arr[2]=iy; i_arr[3]=iz;
                    idx1 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), nu);
                    i_arr[nu] = (i_arr[nu] + 1) % sizes[nu];
                    idx2 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), mu);
                    i_arr[nu] = (i_arr[nu] - 1 + sizes[nu]) % sizes[nu];
                    i_arr[mu] = (i_arr[mu] + 1) % sizes[mu];
                    idx3 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), nu);
                    if (idx1>=0 && idx2>=0 && idx3>=0) {
                        cm_eq_cm_ti_cm_dag(SU2_1, test_field+idx2, test_field+idx3);
                        cm_eq_cm_ti_cm(SU2_2, test_field+idx1, SU2_1);
                        cm_pl_eq_cm(S_l, SU2_2);
                    }
                }

                double S_l_sum = 0;
                for (int j=0;j<8;j++) S_l_sum += std::fabs(S_l[j]);
                if (S_l_sum < 1e-15) continue;

                cm_dag_eq_cm(S_l);
                double k = std::sqrt(S_l[0]*S_l[6]-S_l[1]*S_l[7]-S_l[2]*S_l[4]+S_l[3]*S_l[5]);
                if (k < 1e-15) continue;

                double beta_k = beta * k;
                double y_min = std::exp(-beta_k);
                double y_max = std::exp(+beta_k);
                double a[4];
                while (true) {
                    double y = y_min + (y_max - y_min) * DRand();
                    a[0] = std::log(y) / beta_k;
                    if (DRand() <= std::sqrt(1.0 - a[0]*a[0])) break;
                }
                double norm;
                while (true) {
                    a[1] = 2.0*DRand()-1.0;
                    a[2] = 2.0*DRand()-1.0;
                    a[3] = 2.0*DRand()-1.0;
                    norm = a[1]*a[1]+a[2]*a[2]+a[3]*a[3];
                    if (norm>=1e-10 && norm<=1.0) break;
                }
                norm = std::sqrt((1.0-a[0]*a[0])/norm);
                a[1]*=norm; a[2]*=norm; a[3]*=norm;

                double U_0[8];
                cm_eq_cm_dag(U_0, S_l);
                cm_ti_eq_re(U_0, 1.0/k);
                double U_0l[8];
                cm_from_h(U_0l, a);
                cm_eq_cm_ti_cm(SU2_1, U_0l, U_0);
                double h[4];
                h_from_cm(h, SU2_1);
                norm = 1.0/std::sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]+h[3]*h[3]);
                h[0]*=norm; h[1]*=norm; h[2]*=norm; h[3]*=norm;
                cm_from_h(SU2_1, h);

                int link_idx = ggi(get_index(it,ix,iy,iz,T,L), mu);
                if (link_idx >= 0)
                    cm_eq_cm(test_field + link_idx, SU2_1);
            }
        }
        // odd half-sweep
        #pragma omp parallel
        {
            rlxd_init(2, 99 + omp_get_thread_num() * 997 + sweep * 31);
        }
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < (int)od.size(); i++) {
            int site = od[i];
            int iz = site % L;
            int iy = (site / L) % L;
            int ix = (site / (L * L)) % L;
            int it = site / (L * L * L);
            for (int mu = 0; mu < 4; mu++) {
                double S_l[8];
                cm_eq_zero(S_l);
                double SU2_1[8], SU2_2[8];

                for (int nu = 0; nu < 4; nu++) {
                    if (nu == mu) continue;
                    int i_arr[4] = {it, ix, iy, iz};
                    int sizes[4] = {T, L, L, L};

                    i_arr[0]=it; i_arr[1]=ix; i_arr[2]=iy; i_arr[3]=iz;
                    i_arr[nu] = (i_arr[nu] - 1 + sizes[nu]) % sizes[nu];
                    int idx1 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), nu);
                    int idx2 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), mu);
                    i_arr[mu] = (i_arr[mu] + 1) % sizes[mu];
                    int idx3 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), nu);
                    if (idx1>=0 && idx2>=0 && idx3>=0) {
                        cm_eq_cm_ti_cm(SU2_1, test_field+idx2, test_field+idx3);
                        cm_eq_cm_dag_ti_cm(SU2_2, test_field+idx1, SU2_1);
                        cm_pl_eq_cm(S_l, SU2_2);
                    }

                    i_arr[0]=it; i_arr[1]=ix; i_arr[2]=iy; i_arr[3]=iz;
                    idx1 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), nu);
                    i_arr[nu] = (i_arr[nu] + 1) % sizes[nu];
                    idx2 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), mu);
                    i_arr[nu] = (i_arr[nu] - 1 + sizes[nu]) % sizes[nu];
                    i_arr[mu] = (i_arr[mu] + 1) % sizes[mu];
                    idx3 = ggi(get_index(i_arr[0],i_arr[1],i_arr[2],i_arr[3],T,L), nu);
                    if (idx1>=0 && idx2>=0 && idx3>=0) {
                        cm_eq_cm_ti_cm_dag(SU2_1, test_field+idx2, test_field+idx3);
                        cm_eq_cm_ti_cm(SU2_2, test_field+idx1, SU2_1);
                        cm_pl_eq_cm(S_l, SU2_2);
                    }
                }

                double S_l_sum = 0;
                for (int j=0;j<8;j++) S_l_sum += std::fabs(S_l[j]);
                if (S_l_sum < 1e-15) continue;

                cm_dag_eq_cm(S_l);
                double k = std::sqrt(S_l[0]*S_l[6]-S_l[1]*S_l[7]-S_l[2]*S_l[4]+S_l[3]*S_l[5]);
                if (k < 1e-15) continue;

                double beta_k = beta * k;
                double y_min = std::exp(-beta_k);
                double y_max = std::exp(+beta_k);
                double a[4];
                while (true) {
                    double y = y_min + (y_max - y_min) * DRand();
                    a[0] = std::log(y) / beta_k;
                    if (DRand() <= std::sqrt(1.0 - a[0]*a[0])) break;
                }
                double norm;
                while (true) {
                    a[1] = 2.0*DRand()-1.0;
                    a[2] = 2.0*DRand()-1.0;
                    a[3] = 2.0*DRand()-1.0;
                    norm = a[1]*a[1]+a[2]*a[2]+a[3]*a[3];
                    if (norm>=1e-10 && norm<=1.0) break;
                }
                norm = std::sqrt((1.0-a[0]*a[0])/norm);
                a[1]*=norm; a[2]*=norm; a[3]*=norm;

                double U0[8];
                cm_eq_cm_dag(U0, S_l);
                cm_ti_eq_re(U0, 1.0/k);
                double U0l[8];
                cm_from_h(U0l, a);
                cm_eq_cm_ti_cm(SU2_1, U0l, U0);
                double h[4];
                h_from_cm(h, SU2_1);
                norm = 1.0/std::sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]+h[3]*h[3]);
                h[0]*=norm; h[1]*=norm; h[2]*=norm; h[3]*=norm;
                cm_from_h(SU2_1, h);

                int link_idx = ggi(get_index(it,ix,iy,iz,T,L), mu);
                if (link_idx >= 0)
                    cm_eq_cm(test_field + link_idx, SU2_1);
            }
        }
    }

    // verify all links valid SU(2)
    for (int site = 0; site < volume; site++) {
        for (int mu = 0; mu < 4; mu++) {
            int idx = ggi(site, mu);
            if (idx >= 0) {
                double h[4];
                h_from_cm(h, test_field + idx);
                double norm_sq = h[0]*h[0]+h[1]*h[1]+h[2]*h[2]+h[3]*h[3];
                REQUIRE(std::fabs(norm_sq - 1.0) < 1e-10);
            }
        }
    }

    // plaquette should have moved away from 1.0 (cold start)
    double P = Average_Plaquette(test_field, T, L);
    REQUIRE(P < 1.0);
    REQUIRE(P > 0.3); // reasonable range for beta=2.5

    Gauge_Field_Free(&test_field);
}
```

- [ ] **Step 2: Build and run the test**

Run: `./scripts/build.sh release && cd build && ctest -R test_checkerboard --output-on-failure`

Expected: PASS — the parallel half-sweep produces valid SU(2) links and reasonable plaquette.

- [ ] **Step 3: Replace the sequential sweep loop in MC_heatbath.cc**

Add `#include <omp.h>` at the top of `MC_heatbath.cc` (after the existing includes, around line 10).

Before the sweep loop (after `init_checkerboard(T, L)` but before `for (int sweep = 1; ...)`), seed each thread's RNG once:

```cpp
    // seed per-thread RNG once before sweep loop
    #pragma omp parallel
    {
        rlxd_init(2, params.seed + omp_get_thread_num() * 997);
    }
```

Then replace lines 267-371 (the sweep loop body, inside `for (int sweep = 1; ...)`) with:

```cpp
        // even half-sweep
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < (int)even_sites.size(); i++) {
            int site = even_sites[i];
            const int it = site / (L * L * L);

            for (int mu = 0; mu < 4; mu++) {

                alignas(32) double S_l[8];
                cm_eq_zero(S_l);
                alignas(32) double SU2_1[8], SU2_2[8];

                for (int nu = 0; nu < 4; nu++) {
                    if (nu == mu) continue;

                    const int site_minus_nu = neighbor_minus[nu][site];
                    const int idx1 = link_index[site_minus_nu * 4 + nu];
                    const int idx2 = link_index[site_minus_nu * 4 + mu];
                    const int site_minus_nu_plus_mu = neighbor_plus[mu][site_minus_nu];
                    const int idx3 = link_index[site_minus_nu_plus_mu * 4 + nu];

                    if (idx1 >= 0 && idx2 >= 0 && idx3 >= 0) {
                        cm_eq_cm_ti_cm(SU2_1, gauge_field + idx2, gauge_field + idx3);
                        cm_eq_cm_dag_ti_cm(SU2_2, gauge_field + idx1, SU2_1);
                        if (open_boundary_conditions && (it == 0 || it == T-1)) {
                            cm_ti_eq_re(SU2_2, 0.5);
                        }
                        cm_pl_eq_cm(S_l, SU2_2);
                    }

                    const int idx4 = link_index[site * 4 + nu];
                    const int site_plus_nu = neighbor_plus[nu][site];
                    const int idx5 = link_index[site_plus_nu * 4 + mu];
                    const int site_plus_mu = neighbor_plus[mu][site];
                    const int idx6 = link_index[site_plus_mu * 4 + nu];

                    if (idx4 >= 0 && idx5 >= 0 && idx6 >= 0) {
                        cm_eq_cm_ti_cm_dag(SU2_1, gauge_field + idx5, gauge_field + idx6);
                        cm_eq_cm_ti_cm(SU2_2, gauge_field + idx4, SU2_1);
                        if (open_boundary_conditions && (it == 0 || it == T-1)) {
                            cm_ti_eq_re(SU2_2, 0.5);
                        }
                        cm_pl_eq_cm(S_l, SU2_2);
                    }
                }

                double S_l_sum = 0.0;
                for (int j = 0; j < 8; j++) {
                    S_l_sum += fabs(S_l[j]);
                }
                if (S_l_sum < 1e-15) continue;

                cm_dag_eq_cm(S_l);

                const double k = sqrt(S_l[0]*S_l[6] - S_l[1]*S_l[7] - S_l[2]*S_l[4] + S_l[3]*S_l[5]);

                const double beta_k = params.beta * k;
                const double y_min = exp(-beta_k);
                const double y_max = exp(+beta_k);

                alignas(32) double a[4];

                while (true) {
                    double y = y_min + (y_max - y_min) * DRand();
                    a[0] = log(y) / beta_k;
                    if (DRand() <= sqrt(1.0 - a[0]*a[0])) break;
                }

                double norm;
                while (true) {
                    a[1] = 2.0 * DRand() - 1.0;
                    a[2] = 2.0 * DRand() - 1.0;
                    a[3] = 2.0 * DRand() - 1.0;
                    norm = a[1]*a[1] + a[2]*a[2] + a[3]*a[3];
                    if (norm >= 1e-10 && norm <= 1.0) break;
                }
                norm = sqrt((1.0 - a[0]*a[0]) / norm);
                a[1] *= norm;
                a[2] *= norm;
                a[3] *= norm;

                alignas(32) double U_0[8];
                cm_eq_cm_dag(U_0, S_l);
                cm_ti_eq_re(U_0, 1.0/k);

                alignas(32) double U_0l[8];
                cm_from_h(U_0l, a);

                cm_eq_cm_ti_cm(SU2_1, U_0l, U_0);

                alignas(32) double h[4];
                h_from_cm(h, SU2_1);
                norm = 1.0 / sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3]);
                h[0] *= norm;
                h[1] *= norm;
                h[2] *= norm;
                h[3] *= norm;
                cm_from_h(SU2_1, h);

                const int current_link_idx = link_index[site * 4 + mu];
                if (current_link_idx >= 0) {
                    cm_eq_cm(gauge_field + current_link_idx, SU2_1);
                }
            }
        }

        // odd half-sweep (identical logic, different site list)
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < (int)odd_sites.size(); i++) {
            int site = odd_sites[i];
            const int it = site / (L * L * L);

            for (int mu = 0; mu < 4; mu++) {

                alignas(32) double S_l[8];
                cm_eq_zero(S_l);
                alignas(32) double SU2_1[8], SU2_2[8];

                for (int nu = 0; nu < 4; nu++) {
                    if (nu == mu) continue;

                    const int site_minus_nu = neighbor_minus[nu][site];
                    const int idx1 = link_index[site_minus_nu * 4 + nu];
                    const int idx2 = link_index[site_minus_nu * 4 + mu];
                    const int site_minus_nu_plus_mu = neighbor_plus[mu][site_minus_nu];
                    const int idx3 = link_index[site_minus_nu_plus_mu * 4 + nu];

                    if (idx1 >= 0 && idx2 >= 0 && idx3 >= 0) {
                        cm_eq_cm_ti_cm(SU2_1, gauge_field + idx2, gauge_field + idx3);
                        cm_eq_cm_dag_ti_cm(SU2_2, gauge_field + idx1, SU2_1);
                        if (open_boundary_conditions && (it == 0 || it == T-1)) {
                            cm_ti_eq_re(SU2_2, 0.5);
                        }
                        cm_pl_eq_cm(S_l, SU2_2);
                    }

                    const int idx4 = link_index[site * 4 + nu];
                    const int site_plus_nu = neighbor_plus[nu][site];
                    const int idx5 = link_index[site_plus_nu * 4 + mu];
                    const int site_plus_mu = neighbor_plus[mu][site];
                    const int idx6 = link_index[site_plus_mu * 4 + nu];

                    if (idx4 >= 0 && idx5 >= 0 && idx6 >= 0) {
                        cm_eq_cm_ti_cm_dag(SU2_1, gauge_field + idx5, gauge_field + idx6);
                        cm_eq_cm_ti_cm(SU2_2, gauge_field + idx4, SU2_1);
                        if (open_boundary_conditions && (it == 0 || it == T-1)) {
                            cm_ti_eq_re(SU2_2, 0.5);
                        }
                        cm_pl_eq_cm(S_l, SU2_2);
                    }
                }

                double S_l_sum = 0.0;
                for (int j = 0; j < 8; j++) {
                    S_l_sum += fabs(S_l[j]);
                }
                if (S_l_sum < 1e-15) continue;

                cm_dag_eq_cm(S_l);

                const double k = sqrt(S_l[0]*S_l[6] - S_l[1]*S_l[7] - S_l[2]*S_l[4] + S_l[3]*S_l[5]);

                const double beta_k = params.beta * k;
                const double y_min = exp(-beta_k);
                const double y_max = exp(+beta_k);

                alignas(32) double a[4];

                while (true) {
                    double y = y_min + (y_max - y_min) * DRand();
                    a[0] = log(y) / beta_k;
                    if (DRand() <= sqrt(1.0 - a[0]*a[0])) break;
                }

                double norm;
                while (true) {
                    a[1] = 2.0 * DRand() - 1.0;
                    a[2] = 2.0 * DRand() - 1.0;
                    a[3] = 2.0 * DRand() - 1.0;
                    norm = a[1]*a[1] + a[2]*a[2] + a[3]*a[3];
                    if (norm >= 1e-10 && norm <= 1.0) break;
                }
                norm = sqrt((1.0 - a[0]*a[0]) / norm);
                a[1] *= norm;
                a[2] *= norm;
                a[3] *= norm;

                alignas(32) double U_0[8];
                cm_eq_cm_dag(U_0, S_l);
                cm_ti_eq_re(U_0, 1.0/k);

                alignas(32) double U_0l[8];
                cm_from_h(U_0l, a);

                cm_eq_cm_ti_cm(SU2_1, U_0l, U_0);

                alignas(32) double h[4];
                h_from_cm(h, SU2_1);
                norm = 1.0 / sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3]);
                h[0] *= norm;
                h[1] *= norm;
                h[2] *= norm;
                h[3] *= norm;
                cm_from_h(SU2_1, h);

                const int current_link_idx = link_index[site * 4 + mu];
                if (current_link_idx >= 0) {
                    cm_eq_cm(gauge_field + current_link_idx, SU2_1);
                }
            }
        }
```

Note: The `alignas(32) double SU2_1[8], SU2_2[8];` that were previously declared before the sweep loop (line 264) are now declared inside each parallel for body so each thread has its own stack copies. Remove the old declarations at line 264.

- [ ] **Step 4: Build and run all tests**

Run: `./scripts/build.sh release && ./scripts/run_tests.sh`

Expected: All tests pass including `test_checkerboard`.

- [ ] **Step 5: Commit**

```bash
git add src/MC_heatbath.cc tests/test_checkerboard.cc
git commit -m "feat: checkerboard OpenMP parallel sweep in SU(2) heatbath"
```

---

### Task 4: Statistical validation test

**Files:**
- Modify: `tests/test_checkerboard.cc`

- [ ] **Step 1: Add plaquette validation test**

Append to `tests/test_checkerboard.cc`. This runs 200 sweeps on a small lattice and checks the average plaquette is in the expected range for beta=2.5 on a 4^4 lattice:

```cpp
TEST_CASE("Checkerboard sweep: plaquette matches expected value", "[checkerboard][validation]") {
    T = 4;
    L = 4;
    open_boundary_conditions = false;
    const int volume = T * L * L * L;

    Gauge_Field_Alloc(&gauge_field, T, L);
    Gauge_Field_Unity(gauge_field, T, L);

    std::vector<int> ev, od;
    for (int site = 0; site < volume; site++) {
        int iz = site % L;
        int iy = (site / L) % L;
        int ix = (site / (L * L)) % L;
        int it = site / (L * L * L);
        if ((it + ix + iy + iz) % 2 == 0) ev.push_back(site);
        else od.push_back(site);
    }

    double beta = 2.5;
    int n_therm = 100;
    int n_meas = 100;
    double plaq_sum = 0.0;

    auto do_half_sweep = [&](const std::vector<int> &sites) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < (int)sites.size(); i++) {
            int site = sites[i];
            int iz = site % L;
            int iy = (site / L) % L;
            int ix = (site / (L * L)) % L;
            int it = site / (L * L * L);
            for (int mu = 0; mu < 4; mu++) {
                // staple + heatbath update (same as MC_heatbath.cc)
                double S_l[8], SU2_1[8], SU2_2[8];
                cm_eq_zero(S_l);
                for (int nu = 0; nu < 4; nu++) {
                    if (nu == mu) continue;
                    int c[4] = {it,ix,iy,iz};
                    int s[4] = {T,L,L,L};
                    c[nu]=(c[nu]-1+s[nu])%s[nu];
                    int i1=ggi(get_index(c[0],c[1],c[2],c[3],T,L),nu);
                    int i2=ggi(get_index(c[0],c[1],c[2],c[3],T,L),mu);
                    c[mu]=(c[mu]+1)%s[mu];
                    int i3=ggi(get_index(c[0],c[1],c[2],c[3],T,L),nu);
                    if(i1>=0&&i2>=0&&i3>=0){cm_eq_cm_ti_cm(SU2_1,gauge_field+i2,gauge_field+i3);cm_eq_cm_dag_ti_cm(SU2_2,gauge_field+i1,SU2_1);cm_pl_eq_cm(S_l,SU2_2);}
                    c[0]=it;c[1]=ix;c[2]=iy;c[3]=iz;
                    i1=ggi(get_index(c[0],c[1],c[2],c[3],T,L),nu);
                    c[nu]=(c[nu]+1)%s[nu];
                    i2=ggi(get_index(c[0],c[1],c[2],c[3],T,L),mu);
                    c[nu]=(c[nu]-1+s[nu])%s[nu];c[mu]=(c[mu]+1)%s[mu];
                    i3=ggi(get_index(c[0],c[1],c[2],c[3],T,L),nu);
                    if(i1>=0&&i2>=0&&i3>=0){cm_eq_cm_ti_cm_dag(SU2_1,gauge_field+i2,gauge_field+i3);cm_eq_cm_ti_cm(SU2_2,gauge_field+i1,SU2_1);cm_pl_eq_cm(S_l,SU2_2);}
                }
                double ss=0; for(int j=0;j<8;j++)ss+=std::fabs(S_l[j]);
                if(ss<1e-15)continue;
                cm_dag_eq_cm(S_l);
                double k=std::sqrt(S_l[0]*S_l[6]-S_l[1]*S_l[7]-S_l[2]*S_l[4]+S_l[3]*S_l[5]);
                if(k<1e-15)continue;
                double bk=beta*k,ym=std::exp(-bk),yx=std::exp(bk);
                double a[4];
                while(true){double y=ym+(yx-ym)*DRand();a[0]=std::log(y)/bk;if(DRand()<=std::sqrt(1.0-a[0]*a[0]))break;}
                double nm;
                while(true){a[1]=2*DRand()-1;a[2]=2*DRand()-1;a[3]=2*DRand()-1;nm=a[1]*a[1]+a[2]*a[2]+a[3]*a[3];if(nm>=1e-10&&nm<=1.0)break;}
                nm=std::sqrt((1.0-a[0]*a[0])/nm);a[1]*=nm;a[2]*=nm;a[3]*=nm;
                double U0[8];cm_eq_cm_dag(U0,S_l);cm_ti_eq_re(U0,1.0/k);
                double U0l[8];cm_from_h(U0l,a);cm_eq_cm_ti_cm(SU2_1,U0l,U0);
                double h[4];h_from_cm(h,SU2_1);
                nm=1.0/std::sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]+h[3]*h[3]);
                h[0]*=nm;h[1]*=nm;h[2]*=nm;h[3]*=nm;cm_from_h(SU2_1,h);
                int li=ggi(get_index(it,ix,iy,iz,T,L),mu);
                if(li>=0)cm_eq_cm(gauge_field+li,SU2_1);
            }
        }
    };

    // seed per-thread RNG once
    #pragma omp parallel
    {
        rlxd_init(2, 42 + omp_get_thread_num() * 997);
    }

    for (int sweep = 0; sweep < n_therm + n_meas; sweep++) {
        do_half_sweep(ev);
        do_half_sweep(od);
        if (sweep >= n_therm) {
            plaq_sum += Average_Plaquette(gauge_field, T, L);
        }
    }

    double avg_plaq = plaq_sum / n_meas;

    // for SU(2) at beta=2.5 on 4^4, <P> should be around 0.62-0.68
    REQUIRE(avg_plaq > 0.55);
    REQUIRE(avg_plaq < 0.75);

    Gauge_Field_Free(&gauge_field);
    gauge_field = nullptr;
}
```

- [ ] **Step 2: Build and run validation test**

Run: `./scripts/build.sh release && cd build && ctest -R test_checkerboard --output-on-failure`

Expected: PASS — average plaquette in expected range.

- [ ] **Step 3: Commit**

```bash
git add tests/test_checkerboard.cc
git commit -m "test: add statistical validation for checkerboard sweep"
```

---

### Task 5: SLURM job array script

**Files:**
- Create: `scripts/slurm_beta_scan.sh`

- [ ] **Step 1: Create the SLURM array job script**

Create `scripts/slurm_beta_scan.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=su2_beta_scan
#SBATCH --output=slurm_logs/job_%A_%a.out
#SBATCH --error=slurm_logs/job_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=8G

# adjust --cpus-per-task, --time, --mem, --partition for your cluster

set -euo pipefail

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

PARAMS_FILE="${PROJECT_DIR}/input/base_params_su2.txt"
BETA_FILE="${PROJECT_DIR}/input/beta_scan_su2.txt"
SEEDS_FILE="${PROJECT_DIR}/input/seeds_su2.txt"

if [[ ! -f "$PARAMS_FILE" ]]; then echo "Missing $PARAMS_FILE"; exit 1; fi
if [[ ! -f "$BETA_FILE" ]]; then echo "Missing $BETA_FILE"; exit 1; fi
if [[ ! -f "$SEEDS_FILE" ]]; then echo "Missing $SEEDS_FILE"; exit 1; fi

# read betas and seeds into arrays
mapfile -t BETAS < <(grep -v '^#' "$BETA_FILE" | grep -v '^$')
mapfile -t SEEDS < <(grep -v '^#' "$SEEDS_FILE" | grep -v '^$')

N_BETAS=${#BETAS[@]}
N_SEEDS=${#SEEDS[@]}

# map SLURM_ARRAY_TASK_ID to (beta_idx, seed_idx)
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
BETA_IDX=$((TASK_ID / N_SEEDS))
SEED_IDX=$((TASK_ID % N_SEEDS))

if [[ $BETA_IDX -ge $N_BETAS ]]; then
    echo "TASK_ID=$TASK_ID out of range (${N_BETAS} betas x ${N_SEEDS} seeds)"
    exit 1
fi

BETA="${BETAS[$BETA_IDX]}"
SEED="${SEEDS[$SEED_IDX]}"

# read base params
T_SIZE=$(grep '^T ' "$PARAMS_FILE" | awk '{print $2}')
L_SIZE=$(grep '^L ' "$PARAMS_FILE" | awk '{print $2}')
BOUNDARY=$(grep '^boundary' "$PARAMS_FILE" | awk '{print $2}')
BOUNDARY=${BOUNDARY:-periodic}

RUN_DIR="${PROJECT_DIR}/data/results/T${T_SIZE}_L${L_SIZE}_b${BETA}_${BOUNDARY}_seed${SEED}"
CONFIG_DIR="${RUN_DIR}/configs"
OUTPUT_DIR="${RUN_DIR}/output"

mkdir -p "$CONFIG_DIR" "$OUTPUT_DIR" "${PROJECT_DIR}/slurm_logs"

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" | tee -a "${RUN_DIR}/job_errors.log"
}

echo "=== Job $TASK_ID: beta=$BETA seed=$SEED ==="
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "Run dir: $RUN_DIR"

# generate input file
TEMP_INPUT="${RUN_DIR}/input_params.txt"
while IFS= read -r line; do
    echo "$line"
done < "$PARAMS_FILE" > "$TEMP_INPUT"
# override beta, seed, directories
sed -i.bak "s/^beta.*/beta            $BETA/" "$TEMP_INPUT"
sed -i.bak "s/^seed.*/seed            $SEED/" "$TEMP_INPUT"
echo "output_dir      ${RUN_DIR}/" >> "$TEMP_INPUT"
echo "config_dir      ${CONFIG_DIR}/" >> "$TEMP_INPUT"
rm -f "${TEMP_INPUT}.bak"

# step 1: heatbath
echo "--- Running heatbath ---"
if ! "${PROJECT_DIR}/build/bin/mc_heatbath" -i "$TEMP_INPUT"; then
    log_error "mc_heatbath failed for beta=$BETA seed=$SEED"
    exit 1
fi

# step 2: detect thermalization
echo "--- Detecting thermalization ---"
PLAQ_FILE="${RUN_DIR}/plaquette.dat"
if [[ ! -f "$PLAQ_FILE" ]]; then
    log_error "plaquette.dat not found"
    exit 1
fi

START_CONF=$(python3 "${PROJECT_DIR}/scripts/detect_thermalization.py" \
    --plaquette-file "$PLAQ_FILE" \
    --save-interval "$(grep '^save_interval' "$PARAMS_FILE" | awk '{print $2}')" \
    2>&1) || {
    log_error "detect_thermalization.py failed: $START_CONF"
    exit 1
}
echo "Thermalization detected at conf $START_CONF"

# step 3: measure topological charge
echo "--- Measuring topological charge ---"
MEAS_INPUT="${RUN_DIR}/meas_input_params.txt"
cat > "$MEAS_INPUT" << EOF
config_dir      ${CONFIG_DIR}/
output_dir      ${OUTPUT_DIR}/
T               ${T_SIZE}
L               ${L_SIZE}
beta            ${BETA}
start_conf      ${START_CONF}
end_conf        $(grep '^num_sweeps' "$PARAMS_FILE" | awk '{print $2}')
conf_step       $(grep '^save_interval' "$PARAMS_FILE" | awk '{print $2}')
smear_steps     $(grep '^smear_steps' "$PARAMS_FILE" | awk '{print $2}')
smear_interval  $(grep '^smear_interval' "$PARAMS_FILE" | awk '{print $2}')
smear_alpha     $(grep '^smear_alpha' "$PARAMS_FILE" | awk '{print $2}')
seed            ${SEED}
boundary        ${BOUNDARY}
exclude_boundary_slices $(grep '^exclude_boundary_slices' "$PARAMS_FILE" | awk '{print $2}')
EOF

if ! "${PROJECT_DIR}/build/bin/meas_topcharge" -i "$MEAS_INPUT"; then
    log_error "meas_topcharge failed for beta=$BETA seed=$SEED"
    exit 1
fi

echo "=== Job $TASK_ID complete ==="
```

- [ ] **Step 2: Create seeds input file and template**

Create `input/example_seeds.txt` (template) and copy to `input/seeds_su2.txt` (used by the SLURM script):

```
# Seeds for independent ensembles
# One seed per line (must be >= 1)
12345
67890
24680
```

- [ ] **Step 3: Make script executable**

Run: `chmod +x scripts/slurm_beta_scan.sh`

- [ ] **Step 4: Verify script syntax**

Run: `bash -n scripts/slurm_beta_scan.sh`

Expected: No output (no syntax errors).

- [ ] **Step 5: Commit**

```bash
git add scripts/slurm_beta_scan.sh input/example_seeds.txt input/seeds_su2.txt
git commit -m "feat: add SLURM array job script for parallel beta scans"
```

---

### Task 6: Monitoring script

**Files:**
- Create: `scripts/monitor_runs.py`

- [ ] **Step 1: Create the monitoring script**

Create `scripts/monitor_runs.py`:

```python
#!/usr/bin/env python3
"""Plot thermalization progress from running simulations."""

import argparse
import glob
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_dat(path):
    """Load a .dat file, skipping comment lines."""
    try:
        data = np.loadtxt(path, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except Exception:
        return None


def plot_run(run_dir, output_dir):
    """Plot plaquette and therm_topcharge for a single run."""
    plaq_path = os.path.join(run_dir, 'plaquette.dat')
    qtherm_path = os.path.join(run_dir, 'therm_topcharge.dat')

    plaq = load_dat(plaq_path)
    qtherm = load_dat(qtherm_path)

    if plaq is None or len(plaq) < 2:
        return

    run_name = os.path.basename(run_dir)
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    axes[0].plot(plaq[:, 0], plaq[:, 1], 'b-', linewidth=0.5)
    axes[0].set_ylabel('Plaquette')
    axes[0].set_title(run_name)

    if qtherm is not None and len(qtherm) > 1:
        axes[1].plot(qtherm[:, 0], qtherm[:, 1], 'r-', linewidth=0.5)
    axes[1].set_ylabel('Q (unsmeared)')
    axes[1].set_xlabel('Sweep')

    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f'{run_name}.png')
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f'Saved: {out_path}')


def main():
    parser = argparse.ArgumentParser(description='Monitor running simulations')
    parser.add_argument('--run-dir', type=str, help='Specific run directory')
    parser.add_argument('--beta', type=str, help='Filter by beta value')
    parser.add_argument('--data-dir', type=str, default='data/results',
                        help='Base data directory')
    parser.add_argument('--output-dir', type=str, default='output/monitoring',
                        help='Output directory for plots')
    args = parser.parse_args()

    if args.run_dir:
        run_dirs = [args.run_dir]
    else:
        pattern = os.path.join(args.data_dir, 'T*_L*_b*_*_seed*')
        run_dirs = sorted(glob.glob(pattern))

    if args.beta:
        run_dirs = [d for d in run_dirs if f'_b{args.beta}_' in d]

    if not run_dirs:
        print('No runs found.')
        return

    for run_dir in run_dirs:
        if os.path.exists(os.path.join(run_dir, 'plaquette.dat')):
            plot_run(run_dir, args.output_dir)


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Make script executable**

Run: `chmod +x scripts/monitor_runs.py`

- [ ] **Step 3: Commit**

```bash
git add scripts/monitor_runs.py
git commit -m "feat: add monitoring script for live run visualization"
```

---

### Task 7: Final integration test and cleanup

- [ ] **Step 1: Build everything**

Run: `./scripts/build.sh release`

Expected: Clean build with no warnings related to our changes.

- [ ] **Step 2: Run all tests**

Run: `./scripts/run_tests.sh`

Expected: All 6 test suites pass (the 5 existing + `test_checkerboard`).

- [ ] **Step 3: Quick smoke test of mc_heatbath binary**

Create a small test input and run 50 sweeps to verify the binary works end-to-end:

Run:
```bash
mkdir -p /tmp/su2_test/configs
cat > /tmp/su2_test_input.txt << 'EOF'
T               4
L               4
beta            2.5
seed            42
start_type      cold
boundary        periodic
num_sweeps      50
save_interval   10
output_dir      /tmp/su2_test/
config_dir      /tmp/su2_test/configs/
EOF
OMP_NUM_THREADS=2 ./build/bin/mc_heatbath -i /tmp/su2_test_input.txt
```

Expected: Completes without error, prints plaquette values, creates config files.

- [ ] **Step 4: Verify monitor script works on smoke test output**

Run: `python3 scripts/monitor_runs.py --run-dir /tmp/su2_test --output-dir /tmp/su2_test/monitoring`

Expected: Creates a PNG plot file.

- [ ] **Step 5: Clean up smoke test**

Run: `rm -rf /tmp/su2_test /tmp/su2_test_input.txt`

- [ ] **Step 6: Commit any final fixes**

Only if Step 1-4 revealed issues that needed fixing.
