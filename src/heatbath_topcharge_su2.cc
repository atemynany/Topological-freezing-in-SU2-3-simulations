// ==============================================================================
// heatbath_topcharge_su2.cc
// ==============================================================================
// Fused SU(2) heatbath + in-memory topological charge measurement.
//
// At each save_interval the live gauge field is copied into a scratch buffer,
// APE-smeared, and Q (global + per-timeslice) is measured at the requested
// smear checkpoints. Nothing is written to disk for the raw gauge field —
// only plaquette.dat, topcharge.dat, topcharge_timeslice.dat.
//
// Usage: heatbath_topcharge_su2 -i <input_file>
// ==============================================================================

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

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

#include "Plaquette.hh"
#include "action_density.hh"
#include "progressbar.hh"
#include "smearing_techniques.hh"
#include "topcharge_su2.hh"

int T;
int L;
double *gauge_field;
bool open_boundary_conditions = false;
bool hot_start = false;

std::vector<int> neighbor_plus[4];
std::vector<int> neighbor_minus[4];
std::vector<long long> link_index;

std::vector<int> even_sites;
std::vector<int> odd_sites;

void init_neighbor_tables(int T_size, int L_size) {
    const int volume = T_size * L_size * L_size * L_size;

    for (int mu = 0; mu < 4; mu++) {
        neighbor_plus[mu].resize(volume);
        neighbor_minus[mu].resize(volume);
    }
    link_index.resize(volume * 4);
    even_sites.clear();
    odd_sites.clear();
    even_sites.reserve(volume / 2 + 1);
    odd_sites.reserve(volume / 2 + 1);

    for (int it = 0; it < T_size; it++) {
        for (int ix = 0; ix < L_size; ix++) {
            for (int iy = 0; iy < L_size; iy++) {
                for (int iz = 0; iz < L_size; iz++) {
                    int site = get_index(it, ix, iy, iz, T_size, L_size);

                    if ((it + ix + iy + iz) % 2 == 0)
                        even_sites.push_back(site);
                    else
                        odd_sites.push_back(site);

                    int coords[4] = {it, ix, iy, iz};
                    int sizes[4] = {T_size, L_size, L_size, L_size};

                    for (int mu = 0; mu < 4; mu++) {
                        int c_plus[4] = {coords[0], coords[1], coords[2], coords[3]};
                        int c_minus[4] = {coords[0], coords[1], coords[2], coords[3]};

                        c_plus[mu] = (coords[mu] + 1) % sizes[mu];
                        c_minus[mu] = (coords[mu] - 1 + sizes[mu]) % sizes[mu];

                        neighbor_plus[mu][site] = get_index(c_plus[0], c_plus[1], c_plus[2], c_plus[3], T_size, L_size);
                        neighbor_minus[mu][site] = get_index(c_minus[0], c_minus[1], c_minus[2], c_minus[3], T_size, L_size);

                        link_index[site * 4 + mu] = ggi(site, mu);
                    }
                }
            }
        }
    }
}

struct SimParams {
    std::string output_dir;

    int T;
    int L;
    double beta;

    int seed;
    std::string start_type;
    std::string boundary;
    int num_sweeps;
    int save_interval;

    // Measurement parameters
    int smear_steps;
    int smear_interval;
    double smear_alpha;
    int exclude_boundary_slices;
};

bool read_input_file(const char *filename, SimParams &params) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open input file: " << filename << std::endl;
        return false;
    }

    params.output_dir = "output/";
    params.T = 8;
    params.L = 8;
    params.beta = 2.5;
    params.seed = 12345;
    params.start_type = "cold";
    params.boundary = "periodic";
    params.num_sweeps = 100;
    params.save_interval = 10;
    params.smear_steps = 40;
    params.smear_interval = 10;
    params.smear_alpha = 0.5;
    params.exclude_boundary_slices = 0;

    std::string line, key;

    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        iss >> key;

        if      (key == "output_dir")              iss >> params.output_dir;
        else if (key == "T")                       iss >> params.T;
        else if (key == "L")                       iss >> params.L;
        else if (key == "beta")                    iss >> params.beta;
        else if (key == "seed")                    iss >> params.seed;
        else if (key == "start_type")              iss >> params.start_type;
        else if (key == "boundary")                iss >> params.boundary;
        else if (key == "num_sweeps")              iss >> params.num_sweeps;
        else if (key == "save_interval")           iss >> params.save_interval;
        else if (key == "smear_steps")             iss >> params.smear_steps;
        else if (key == "smear_interval")          iss >> params.smear_interval;
        else if (key == "smear_alpha")             iss >> params.smear_alpha;
        else if (key == "exclude_boundary_slices") iss >> params.exclude_boundary_slices;
    }

    infile.close();
    return true;
}

bool validate_params(const SimParams &params) {
    if (params.beta <= 0.0) {
        std::cerr << "Error: beta must be > 0" << std::endl;
        return false;
    }
    if (params.T < 2 || params.L < 2) {
        std::cerr << "Error: T and L must be >= 2" << std::endl;
        return false;
    }
    if (params.seed < 1) {
        std::cerr << "Error: seed must be >= 1" << std::endl;
        return false;
    }
    if (params.start_type != "cold" && params.start_type != "hot") {
        std::cerr << "Error: start_type must be 'cold' or 'hot'" << std::endl;
        return false;
    }
    if (params.boundary != "periodic" && params.boundary != "open") {
        std::cerr << "Error: boundary must be 'periodic' or 'open'" << std::endl;
        return false;
    }
    if (params.num_sweeps < 1 || params.save_interval < 1) {
        std::cerr << "Error: num_sweeps and save_interval must be >= 1" << std::endl;
        return false;
    }
    if (params.smear_steps < 1 || params.smear_interval < 1) {
        std::cerr << "Error: smear_steps and smear_interval must be >= 1" << std::endl;
        return false;
    }
    if (params.smear_alpha < 0.0 || params.smear_alpha > 1.0) {
        std::cerr << "Warning: smear_alpha=" << params.smear_alpha << " outside typical range [0,1]" << std::endl;
    }
    return true;
}

void print_params(const SimParams &params) {
    std::cout << "========================================" << std::endl;
    std::cout << "SU(2) Heatbath + In-Memory Topcharge" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Output directory:    " << params.output_dir << std::endl;
    std::cout << "Lattice size:        " << params.T << " x " << params.L << "^3" << std::endl;
    std::cout << "Beta:                " << params.beta << std::endl;
    std::cout << "Seed:                " << params.seed << std::endl;
    std::cout << "Start type:          " << params.start_type << std::endl;
    std::cout << "Boundary conditions: " << params.boundary << std::endl;
    std::cout << "Number of sweeps:    " << params.num_sweeps << std::endl;
    std::cout << "Measure interval:    " << params.save_interval << " sweeps" << std::endl;
    std::cout << "Smear total:         " << params.smear_steps << " steps, alpha=" << params.smear_alpha << std::endl;
    std::cout << "Smear checkpoints:   every " << params.smear_interval << " steps" << std::endl;
    if (params.boundary == "open" && params.exclude_boundary_slices > 0) {
        std::cout << "Exclude boundary:    " << params.exclude_boundary_slices << " slices per end" << std::endl;
    }
#ifdef _OPENMP
    std::cout << "OpenMP threads:      " << omp_get_max_threads() << std::endl;
#endif
    std::cout << "========================================" << std::endl;
}

// ==============================================================================
// In-memory measurement: copy live field into scratch buffer, smear step-by-
// step, record Q + q(t) at every smear_interval (and the final step).
// ==============================================================================
void measure_topcharge_inmemory(
    const SimParams &params,
    double *smeared_gauge_field,
    int sweep,
    std::ofstream &tc_file,
    std::ofstream &ts_file)
{
    Gauge_Field_Copy(smeared_gauge_field, gauge_field, params.T, params.L);

    const bool open_bc = (params.boundary == "open" && params.exclude_boundary_slices > 0);
    const int t_min = open_bc ? params.exclude_boundary_slices : 0;
    const int t_max = open_bc ? params.T - params.exclude_boundary_slices : params.T;
    const int t_count = t_max - t_min;

    std::vector<double> q_per_t(t_count, 0.0);

    for (int s = 1; s <= params.smear_steps; s++) {
        APE_Smearing_all(smeared_gauge_field, params.T, params.L, params.smear_alpha);

        const bool is_checkpoint = (s % params.smear_interval == 0) || (s == params.smear_steps);
        if (!is_checkpoint) continue;

        double Q = open_bc
            ? compute_topological_charge_open(smeared_gauge_field, params.T, params.L, params.exclude_boundary_slices)
            : compute_topological_charge(smeared_gauge_field, params.T, params.L);

        double plaq = Average_Plaquette(smeared_gauge_field, params.T, params.L);

        tc_file << std::fixed << std::setprecision(6)
                << std::setw(5)  << s << "  "
                << std::setw(8)  << sweep << "  "
                << std::setw(12) << Q << "  "
                << std::setw(10) << plaq << "\n";

        // Per-timeslice density q(t) = (1/4π²) Σ_{x,y,z} q_local(t,x,y,z)
        #pragma omp parallel for schedule(static)
        for (int it_idx = 0; it_idx < t_count; it_idx++) {
            const int it = t_min + it_idx;
            double q_t = 0.0;
            for (int ix = 0; ix < params.L; ix++) {
                for (int iy = 0; iy < params.L; iy++) {
                    for (int iz = 0; iz < params.L; iz++) {
                        q_t += compute_local_topcharge_density(smeared_gauge_field, it, ix, iy, iz, params.T, params.L);
                    }
                }
            }
            q_per_t[it_idx] = q_t / (4.0 * M_PI * M_PI);
        }

        ts_file << std::fixed << std::setprecision(8);
        for (int it_idx = 0; it_idx < t_count; it_idx++) {
            ts_file << std::setw(5)  << s << "  "
                    << std::setw(8)  << sweep << "  "
                    << std::setw(4)  << (t_min + it_idx) << "  "
                    << std::setw(14) << q_per_t[it_idx] << "\n";
        }
    }
    tc_file.flush();
    ts_file.flush();
}

int main(int argc, char **argv)
{
    const char *input_file = nullptr;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            input_file = argv[i + 1];
            i++;
        }
    }

    if (input_file == nullptr) {
        std::cerr << "Usage: " << argv[0] << " -i <input_file>" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Required keys (text file, one key per line):" << std::endl;
        std::cerr << "  output_dir   T   L   beta   seed" << std::endl;
        std::cerr << "  start_type {cold|hot}   boundary {periodic|open}" << std::endl;
        std::cerr << "  num_sweeps   save_interval" << std::endl;
        std::cerr << "  smear_steps  smear_interval  smear_alpha" << std::endl;
        std::cerr << "  exclude_boundary_slices  (only used for open BC)" << std::endl;
        return EXIT_FAILURE;
    }

    SimParams params;
    if (!read_input_file(input_file, params)) return EXIT_FAILURE;
    if (!validate_params(params)) return EXIT_FAILURE;
    print_params(params);

    T = params.T;
    L = params.L;
    hot_start = (params.start_type == "hot");
    open_boundary_conditions = (params.boundary == "open");

    InitializeRand(params.seed);
#ifdef _OPENMP
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        if (tid > 0) rlxd_init(2, params.seed + tid * 997);
    }
#endif

    init_neighbor_tables(T, L);

    mkdir(params.output_dir.c_str(), 0755);

    double *smeared_gauge_field = nullptr;
    Gauge_Field_Alloc(&gauge_field, T, L);
    Gauge_Field_Alloc(&smeared_gauge_field, T, L);

    if (hot_start) {
        std::cout << "Initializing hot start..." << std::endl;
        Gauge_Field_Random(gauge_field, T, L);
    } else {
        std::cout << "Initializing cold start..." << std::endl;
        Gauge_Field_Unity(gauge_field, T, L);
    }

    double P = Average_Plaquette(gauge_field, T, L);
    double action_density = avg_action_density(params.beta, P);
    std::cout << "Initial <P> = " << P << std::endl;

    // Output files
    std::string plaq_filename = params.output_dir + "plaquette.dat";
    FILE *plaq_file = fopen(plaq_filename.c_str(), "w");
    if (plaq_file == nullptr) {
        std::cerr << "Error: Cannot open plaquette file: " << plaq_filename << std::endl;
        return EXIT_FAILURE;
    }
    fprintf(plaq_file, "# Sweep  Plaquette  ActionDensity\n");
    fprintf(plaq_file, "%5d %+.10e %+.10e\n", 0, P, action_density);

    std::string action_filename = params.output_dir + "action_density.dat";
    FILE *action_file = fopen(action_filename.c_str(), "w");
    if (action_file == nullptr) {
        std::cerr << "Error: Cannot open action density file: " << action_filename << std::endl;
        return EXIT_FAILURE;
    }
    fprintf(action_file, "# Sweep  ActionDensity\n");
    fprintf(action_file, "%5d %+.10e\n", 0, action_density);

    std::string tc_filename = params.output_dir + "topcharge.dat";
    std::ofstream tc_file(tc_filename);
    if (!tc_file.is_open()) {
        std::cerr << "Error: Cannot open topcharge file: " << tc_filename << std::endl;
        return EXIT_FAILURE;
    }
    tc_file << "# smear_steps  sweep     Q            plaquette\n";

    std::string ts_filename = params.output_dir + "topcharge_timeslice.dat";
    std::ofstream ts_file(ts_filename);
    if (!ts_file.is_open()) {
        std::cerr << "Error: Cannot open timeslice file: " << ts_filename << std::endl;
        return EXIT_FAILURE;
    }
    ts_file << "# smear_steps  sweep     t    q_t\n";

    const int L3 = L * L * L;

    for (int sweep = 1; sweep <= params.num_sweeps; sweep++) {

        // Checkerboard heatbath: even sites then odd sites
        for (int parity = 0; parity < 2; parity++) {
            const std::vector<int> &sites = (parity == 0) ? even_sites : odd_sites;
            const int n_sites = static_cast<int>(sites.size());

            #pragma omp parallel for schedule(static)
            for (int i = 0; i < n_sites; i++) {
                alignas(32) double SU2_1[8], SU2_2[8];
                const int site = sites[i];
                const int it = site / L3;

                for (int mu = 0; mu < 4; mu++) {

                    alignas(32) double S_l[8];
                    cm_eq_zero(S_l);

                    for (int nu = 0; nu < 4; nu++) {
                        if (nu == mu) continue;

                        // Lower staple: U_nu^dag(x-nu) * U_mu(x-nu) * U_nu(x-nu+mu)
                        const int site_minus_nu = neighbor_minus[nu][site];
                        const long long idx1 = link_index[site_minus_nu * 4 + nu];
                        const long long idx2 = link_index[site_minus_nu * 4 + mu];
                        const int site_minus_nu_plus_mu = neighbor_plus[mu][site_minus_nu];
                        const long long idx3 = link_index[site_minus_nu_plus_mu * 4 + nu];

                        if (idx1 >= 0 && idx2 >= 0 && idx3 >= 0) {
                            cm_eq_cm_ti_cm(SU2_1, gauge_field + idx2, gauge_field + idx3);
                            cm_eq_cm_dag_ti_cm(SU2_2, gauge_field + idx1, SU2_1);
                            if (open_boundary_conditions && (it == 0 || it == T-1)) {
                                cm_ti_eq_re(SU2_2, 0.5);
                            }
                            cm_pl_eq_cm(S_l, SU2_2);
                        }

                        // Upper staple: U_nu(x) * U_mu(x+nu) * U_nu^dag(x+mu)
                        const long long idx4 = link_index[site * 4 + nu];
                        const int site_plus_nu = neighbor_plus[nu][site];
                        const long long idx5 = link_index[site_plus_nu * 4 + mu];
                        const int site_plus_mu = neighbor_plus[mu][site];
                        const long long idx6 = link_index[site_plus_mu * 4 + nu];

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

                    const long long current_link_idx = link_index[site * 4 + mu];
                    if (current_link_idx >= 0) {
                        cm_eq_cm(gauge_field + current_link_idx, SU2_1);
                    }
                }
            }
        }

        P = Average_Plaquette(gauge_field, T, L);
        action_density = avg_action_density(params.beta, P);
        fprintf(plaq_file, "%5d %+.10e %+.10e\n", sweep, P, action_density);
        fprintf(action_file, "%5d %+.10e\n", sweep, action_density);
        fflush(plaq_file);
        fflush(action_file);

        double progress = static_cast<double>(sweep) / params.num_sweeps;
        progress_bar(progress);

        if (sweep % params.save_interval == 0) {
            measure_topcharge_inmemory(params, smeared_gauge_field, sweep, tc_file, ts_file);
        }
    }

    progress_bar_clear();
    std::cout << "Completed " << params.num_sweeps << " sweeps, final <P> = " << P << std::endl;

    fclose(plaq_file);
    fclose(action_file);
    tc_file.close();
    ts_file.close();
    Gauge_Field_Free(&gauge_field);
    Gauge_Field_Free(&smeared_gauge_field);

    std::cout << "========================================" << std::endl;
    std::cout << "Plaquette history:   " << plaq_filename << std::endl;
    std::cout << "Action density:      " << action_filename << std::endl;
    std::cout << "Topological charge:  " << tc_filename << std::endl;
    std::cout << "Per-timeslice q(t):  " << ts_filename << std::endl;
    std::cout << "========================================" << std::endl;

    return EXIT_SUCCESS;
}
