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
#include "Plaquette.hh"

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

// heatbath update for a single link — used by the validation test
static void heatbath_link(double *gf, int site, int mu, double beta, int T_s, int L_s) {
    int iz = site % L_s;
    int iy = (site / L_s) % L_s;
    int ix = (site / (L_s * L_s)) % L_s;
    int it = site / (L_s * L_s * L_s);

    double S_l[8], SU2_1[8], SU2_2[8];
    cm_eq_zero(S_l);

    for (int nu = 0; nu < 4; nu++) {
        if (nu == mu) continue;
        int c[4] = {it, ix, iy, iz};
        int s[4] = {T_s, L_s, L_s, L_s};

        // lower staple
        c[0]=it;c[1]=ix;c[2]=iy;c[3]=iz;
        c[nu] = (c[nu]-1+s[nu])%s[nu];
        int i1=ggi(get_index(c[0],c[1],c[2],c[3],T_s,L_s),nu);
        int i2=ggi(get_index(c[0],c[1],c[2],c[3],T_s,L_s),mu);
        c[mu]=(c[mu]+1)%s[mu];
        int i3=ggi(get_index(c[0],c[1],c[2],c[3],T_s,L_s),nu);
        if(i1>=0&&i2>=0&&i3>=0) {
            cm_eq_cm_ti_cm(SU2_1,gf+i2,gf+i3);
            cm_eq_cm_dag_ti_cm(SU2_2,gf+i1,SU2_1);
            cm_pl_eq_cm(S_l,SU2_2);
        }

        // upper staple
        c[0]=it;c[1]=ix;c[2]=iy;c[3]=iz;
        i1=ggi(get_index(c[0],c[1],c[2],c[3],T_s,L_s),nu);
        c[nu]=(c[nu]+1)%s[nu];
        i2=ggi(get_index(c[0],c[1],c[2],c[3],T_s,L_s),mu);
        c[nu]=(c[nu]-1+s[nu])%s[nu]; c[mu]=(c[mu]+1)%s[mu];
        i3=ggi(get_index(c[0],c[1],c[2],c[3],T_s,L_s),nu);
        if(i1>=0&&i2>=0&&i3>=0) {
            cm_eq_cm_ti_cm_dag(SU2_1,gf+i2,gf+i3);
            cm_eq_cm_ti_cm(SU2_2,gf+i1,SU2_1);
            cm_pl_eq_cm(S_l,SU2_2);
        }
    }

    double ss=0; for(int j=0;j<8;j++) ss+=std::fabs(S_l[j]);
    if(ss<1e-15) return;
    cm_dag_eq_cm(S_l);
    double k=std::sqrt(S_l[0]*S_l[6]-S_l[1]*S_l[7]-S_l[2]*S_l[4]+S_l[3]*S_l[5]);
    if(k<1e-15) return;

    double bk=beta*k, ym=std::exp(-bk), yx=std::exp(bk);
    double a[4];
    while(true) {
        double y=ym+(yx-ym)*DRand(); a[0]=std::log(y)/bk;
        if(DRand()<=std::sqrt(1.0-a[0]*a[0])) break;
    }
    double nm;
    while(true) {
        a[1]=2*DRand()-1; a[2]=2*DRand()-1; a[3]=2*DRand()-1;
        nm=a[1]*a[1]+a[2]*a[2]+a[3]*a[3];
        if(nm>=1e-10&&nm<=1.0) break;
    }
    nm=std::sqrt((1.0-a[0]*a[0])/nm);
    a[1]*=nm; a[2]*=nm; a[3]*=nm;

    double U0[8]; cm_eq_cm_dag(U0,S_l); cm_ti_eq_re(U0,1.0/k);
    double U0l[8]; cm_from_h(U0l,a);
    cm_eq_cm_ti_cm(SU2_1,U0l,U0);
    double h[4]; h_from_cm(h,SU2_1);
    nm=1.0/std::sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]+h[3]*h[3]);
    h[0]*=nm; h[1]*=nm; h[2]*=nm; h[3]*=nm;
    cm_from_h(SU2_1,h);
    int li=ggi(get_index(it,ix,iy,iz,T_s,L_s),mu);
    if(li>=0) cm_eq_cm(gf+li,SU2_1);
}

TEST_CASE("Checkerboard sweep: plaquette in expected range", "[checkerboard][validation]") {
    T = 4;
    L = 4;
    open_boundary_conditions = false;

    Gauge_Field_Alloc(&gauge_field, T, L);
    Gauge_Field_Unity(gauge_field, T, L);

    const int volume = T * L * L * L;
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

    InitializeRand(42);

    // checkerboard sweeps
    for (int sweep = 0; sweep < n_therm + n_meas; sweep++) {
        for (int i = 0; i < (int)ev.size(); i++)
            for (int mu = 0; mu < 4; mu++)
                heatbath_link(gauge_field, ev[i], mu, beta, T, L);
        for (int i = 0; i < (int)od.size(); i++)
            for (int mu = 0; mu < 4; mu++)
                heatbath_link(gauge_field, od[i], mu, beta, T, L);
    }

    double plaq = Average_Plaquette(gauge_field, T, L);

    // SU(2) at beta=2.5 on 4^4: <P> ~ 0.62-0.68
    REQUIRE(plaq > 0.55);
    REQUIRE(plaq < 0.75);

    // verify all links are valid SU(2)
    for (int site = 0; site < volume; site++) {
        for (int mu = 0; mu < 4; mu++) {
            int idx = ggi(site, mu);
            if (idx >= 0) {
                double h[4];
                h_from_cm(h, gauge_field + idx);
                double norm_sq = h[0]*h[0]+h[1]*h[1]+h[2]*h[2]+h[3]*h[3];
                REQUIRE(std::fabs(norm_sq - 1.0) < 1e-10);
            }
        }
    }

    Gauge_Field_Free(&gauge_field);
    gauge_field = nullptr;
}
