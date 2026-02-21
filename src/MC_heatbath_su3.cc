// MC_heatbath_su3.cc - SU(3) Monte Carlo heatbath simulation
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "su3_linear_algebra.hh"
#include "su3_heatbath.hh"
#include "ranlux.hh"
#include "progressbar.hh"

int T_size, L_size;
double *gauge_field_su3;
bool open_boundary_conditions = false;

std::vector<int> neighbor_plus[4];
std::vector<int> neighbor_minus[4];
std::vector<int> link_index_su3;

inline int get_site_index(int t, int x, int y, int z) {
    int tt = (t + T_size) % T_size;
    int xx = (x + L_size) % L_size;
    int yy = (y + L_size) % L_size;
    int zz = (z + L_size) % L_size;
    return ((tt * L_size + xx) * L_size + yy) * L_size + zz;
}

inline int ggi_su3(int site, int mu) {
    return (4 * site + mu) * 18;
}

void init_neighbor_tables() {
    const int volume = T_size * L_size * L_size * L_size;
    
    for (int mu = 0; mu < 4; mu++) {
        neighbor_plus[mu].resize(volume);
        neighbor_minus[mu].resize(volume);
    }
    link_index_su3.resize(volume * 4);
    
    for (int t = 0; t < T_size; t++) {
        for (int x = 0; x < L_size; x++) {
            for (int y = 0; y < L_size; y++) {
                for (int z = 0; z < L_size; z++) {
                    int site = get_site_index(t, x, y, z);
                    int coords[4] = {t, x, y, z};
                    int sizes[4] = {T_size, L_size, L_size, L_size};
                    
                    for (int mu = 0; mu < 4; mu++) {
                        int c_plus[4] = {coords[0], coords[1], coords[2], coords[3]};
                        int c_minus[4] = {coords[0], coords[1], coords[2], coords[3]};
                        c_plus[mu] = (coords[mu] + 1) % sizes[mu];
                        c_minus[mu] = (coords[mu] - 1 + sizes[mu]) % sizes[mu];
                        
                        neighbor_plus[mu][site] = get_site_index(c_plus[0], c_plus[1], c_plus[2], c_plus[3]);
                        neighbor_minus[mu][site] = get_site_index(c_minus[0], c_minus[1], c_minus[2], c_minus[3]);
                        link_index_su3[site * 4 + mu] = ggi_su3(site, mu);
                    }
                }
            }
        }
    }
}

void su3_gauge_field_alloc(double **gf, int T, int L) {
    const int volume = T * L * L * L;
    *gf = (double *)aligned_alloc(32, volume * 4 * 18 * sizeof(double));
}

void su3_gauge_field_free(double **gf) {
    free(*gf);
    *gf = nullptr;
}

void su3_gauge_field_unity(double *gf, int T, int L) {
    const int volume = T * L * L * L;
    for (int i = 0; i < volume * 4; i++) {
        su3_eq_id(gf + i * 18);
    }
}

void su3_gauge_field_random(double *gf, int T, int L) {
    const int volume = T * L * L * L;
    for (int i = 0; i < volume * 4; i++) {
        su3_random(gf + i * 18, DRand);
    }
}

double su3_average_plaquette(const double *gf, int T, int L) {
    const int volume = T * L * L * L;
    double sum = 0.0;
    int count = 0;
    
    alignas(32) double T1[18], T2[18];
    
    for (int site = 0; site < volume; site++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu + 1; nu < 4; nu++) {
                int idx_mu = link_index_su3[site * 4 + mu];
                int idx_nu = link_index_su3[site * 4 + nu];
                int site_plus_mu = neighbor_plus[mu][site];
                int site_plus_nu = neighbor_plus[nu][site];
                int idx_nu_at_mu = link_index_su3[site_plus_mu * 4 + nu];
                int idx_mu_at_nu = link_index_su3[site_plus_nu * 4 + mu];
                
                const double *U_mu = gf + idx_mu;
                const double *U_nu = gf + idx_nu;
                const double *U_nu_xmu = gf + idx_nu_at_mu;
                const double *U_mu_xnu = gf + idx_mu_at_nu;
                
                su3_eq_su3_ti_su3(T1, U_mu, U_nu_xmu);
                su3_eq_su3_ti_su3_dag(T2, T1, U_mu_xnu);
                su3_eq_su3_ti_su3_dag(T1, T2, U_nu);
                
                sum += su3_re_tr(T1) / 3.0;
                count++;
            }
        }
    }
    return sum / count;
}

// Compute staple sum for link U_mu(x)
// With open boundaries, apply 0.5 weight at temporal boundaries t=0, t=T-1
void compute_staple(double *staple, const double *gf, int site, int mu, int it, bool open_bc) {
    su3_eq_zero(staple);
    alignas(32) double T1[18], T2[18];
    
    const bool at_boundary = open_bc && (it == 0 || it == T_size - 1);
    
    for (int nu = 0; nu < 4; nu++) {
        if (nu == mu) continue;
        
        // Upper staple: U_nu(x) * U_mu(x+nu) * U_nu^dag(x+mu)
        int idx_nu = link_index_su3[site * 4 + nu];
        int site_plus_nu = neighbor_plus[nu][site];
        int idx_mu_at_nu = link_index_su3[site_plus_nu * 4 + mu];
        int site_plus_mu = neighbor_plus[mu][site];
        int idx_nu_at_mu = link_index_su3[site_plus_mu * 4 + nu];
        
        su3_eq_su3_ti_su3(T1, gf + idx_nu, gf + idx_mu_at_nu);
        su3_eq_su3_ti_su3_dag(T2, T1, gf + idx_nu_at_mu);
        if (at_boundary) su3_ti_eq_re(T2, 0.5);
        su3_pl_eq_su3(staple, T2);
        
        // Lower staple: U_nu^dag(x-nu) * U_mu(x-nu) * U_nu(x-nu+mu)
        int site_minus_nu = neighbor_minus[nu][site];
        int idx_nu_at_minus = link_index_su3[site_minus_nu * 4 + nu];
        int idx_mu_at_minus = link_index_su3[site_minus_nu * 4 + mu];
        int site_minus_nu_plus_mu = neighbor_plus[mu][site_minus_nu];
        int idx_nu_at_minus_mu = link_index_su3[site_minus_nu_plus_mu * 4 + nu];
        
        su3_eq_su3_dag_ti_su3(T1, gf + idx_nu_at_minus, gf + idx_mu_at_minus);
        su3_eq_su3_ti_su3(T2, T1, gf + idx_nu_at_minus_mu);
        if (at_boundary) su3_ti_eq_re(T2, 0.5);
        su3_pl_eq_su3(staple, T2);
    }
}

void write_su3_config(const double *gf, const char *filename, int T, int L, const char *header) {
    FILE *f = fopen(filename, "wb");
    if (!f) {
        std::cerr << "Error: Cannot write " << filename << std::endl;
        return;
    }
    fprintf(f, "# %s\n", header);
    const int volume = T * L * L * L;
    fwrite(gf, sizeof(double), volume * 4 * 18, f);
    fclose(f);
}

void read_su3_config(double *gf, const char *filename, int T, int L) {
    FILE *f = fopen(filename, "rb");
    if (!f) {
        std::cerr << "Error: Cannot read " << filename << std::endl;
        return;
    }
    char line[1024];
    if (fgets(line, sizeof(line), f) == nullptr) { fclose(f); return; }
    const int volume = T * L * L * L;
    fread(gf, sizeof(double), volume * 4 * 18, f);
    fclose(f);
}

struct SimParams {
    std::string output_dir;
    std::string config_dir;
    int T, L;
    double beta;
    int seed;
    std::string start_type;
    std::string boundary;
    int num_sweeps;
    int save_interval;
    int overrelax_steps;
};

bool read_input_file(const char *filename, SimParams &params) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open " << filename << std::endl;
        return false;
    }
    
    params.output_dir = "output/";
    params.config_dir = "output/configs_su3/";
    params.T = 8;
    params.L = 8;
    params.beta = 6.0;
    params.seed = 12345;
    params.start_type = "cold";
    params.boundary = "periodic";
    params.num_sweeps = 100;
    params.save_interval = 10;
    params.overrelax_steps = 0;
    
    std::string line, key;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        iss >> key;
        
        if (key == "output_dir") iss >> params.output_dir;
        else if (key == "config_dir") iss >> params.config_dir;
        else if (key == "T") iss >> params.T;
        else if (key == "L") iss >> params.L;
        else if (key == "beta") iss >> params.beta;
        else if (key == "seed") iss >> params.seed;
        else if (key == "start_type") iss >> params.start_type;
        else if (key == "boundary") iss >> params.boundary;
        else if (key == "num_sweeps") iss >> params.num_sweeps;
        else if (key == "save_interval") iss >> params.save_interval;
        else if (key == "overrelax_steps") iss >> params.overrelax_steps;
    }
    infile.close();
    return true;
}

void print_params(const SimParams &params) {
    std::cout << "========================================\n";
    std::cout << "SU(3) Heatbath Monte Carlo\n";
    std::cout << "========================================\n";
    std::cout << "Lattice:            " << params.T << " x " << params.L << "^3\n";
    std::cout << "Beta:               " << params.beta << "\n";
    std::cout << "Seed:               " << params.seed << "\n";
    std::cout << "Start:              " << params.start_type << "\n";
    std::cout << "Boundary:           " << params.boundary << "\n";
    std::cout << "Sweeps:             " << params.num_sweeps << "\n";
    std::cout << "Save interval:      " << params.save_interval << "\n";
    std::cout << "Overrelax steps:    " << params.overrelax_steps << "\n";
    std::cout << "========================================\n";
}

int main(int argc, char **argv) {
    const char *input_file = nullptr;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            input_file = argv[++i];
        }
    }
    
    if (!input_file) {
        std::cerr << "Usage: " << argv[0] << " -i <input_file>\n";
        return EXIT_FAILURE;
    }
    
    SimParams params;
    if (!read_input_file(input_file, params)) return EXIT_FAILURE;
    
    print_params(params);
    
    T_size = params.T;
    L_size = params.L;
    
    InitializeRand(params.seed);
    init_neighbor_tables();
    
    open_boundary_conditions = (params.boundary == "open");
    
    mkdir(params.output_dir.c_str(), 0755);
    mkdir(params.config_dir.c_str(), 0755);
    
    su3_gauge_field_alloc(&gauge_field_su3, T_size, L_size);
    
    if (params.start_type == "hot") {
        std::cout << "Hot start..." << std::endl;
        su3_gauge_field_random(gauge_field_su3, T_size, L_size);
    } else {
        std::cout << "Cold start..." << std::endl;
        su3_gauge_field_unity(gauge_field_su3, T_size, L_size);
    }
    
    double P = su3_average_plaquette(gauge_field_su3, T_size, L_size);
    std::cout << "Initial <P> = " << P << std::endl;
    
    std::string plaq_filename = params.output_dir + "plaquette_su3.dat";
    FILE *plaq_file = fopen(plaq_filename.c_str(), "w");
    fprintf(plaq_file, "# Sweep  Plaquette\n");
    fprintf(plaq_file, "%5d %+.10e\n", 0, P);
    
    const int volume = T_size * L_size * L_size * L_size;
    const int L3 = L_size * L_size * L_size;
    alignas(32) double staple[18];
    
    for (int sweep = 1; sweep <= params.num_sweeps; sweep++) {
        // Heatbath sweep
        for (int site = 0; site < volume; site++) {
            const int it = site / L3;  // time coordinate
            for (int mu = 0; mu < 4; mu++) {
                compute_staple(staple, gauge_field_su3, site, mu, it, open_boundary_conditions);
                int idx = link_index_su3[site * 4 + mu];
                su3_heatbath_link(gauge_field_su3 + idx, staple, params.beta, DRand);
            }
        }
        
        // Overrelaxation
        for (int ov = 0; ov < params.overrelax_steps; ov++) {
            for (int site = 0; site < volume; site++) {
                const int it = site / L3;
                for (int mu = 0; mu < 4; mu++) {
                    compute_staple(staple, gauge_field_su3, site, mu, it, open_boundary_conditions);
                    int idx = link_index_su3[site * 4 + mu];
                    su3_overrelax_link(gauge_field_su3 + idx, staple);
                }
            }
        }
        
        P = su3_average_plaquette(gauge_field_su3, T_size, L_size);
        fprintf(plaq_file, "%5d %+.10e\n", sweep, P);
        fflush(plaq_file);
        
        progress_bar(static_cast<double>(sweep) / params.num_sweeps);
        
        if (sweep % params.save_interval == 0) {
            char fname[1024];
            snprintf(fname, sizeof(fname), "%sconf_su3.%04d", params.config_dir.c_str(), sweep);
            char hdr[256];
            snprintf(hdr, sizeof(hdr), "SU3 T=%d L=%d beta=%.6f sweep=%d", T_size, L_size, params.beta, sweep);
            write_su3_config(gauge_field_su3, fname, T_size, L_size, hdr);
        }
    }
    
    progress_bar_clear();
    std::cout << "Final <P> = " << P << std::endl;
    
    fclose(plaq_file);
    su3_gauge_field_free(&gauge_field_su3);
    
    std::cout << "========================================\n";
    std::cout << "Plaquette: " << plaq_filename << "\n";
    std::cout << "Configs:   " << params.config_dir << "\n";
    std::cout << "========================================\n";
    
    return EXIT_SUCCESS;
}
