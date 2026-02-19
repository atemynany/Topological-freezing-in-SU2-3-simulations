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

#include "fields.hh"
#include "geometry.hh"
#include "io.hh"
#include "linear_algebra.hh"
#include "ranlux.hh"

#include "Plaquette.hh"
#include "progressbar.hh"

int T;
int L;
double *gauge_field;
bool open_boundary_conditions = false; 
bool hot_start = false;

// Precomputed neighbor index tables for faster lookup
std::vector<int> neighbor_plus[4];   
std::vector<int> neighbor_minus[4];  
std::vector<int> link_index;

void init_neighbor_tables(int T_size, int L_size) {
    const int volume = T_size * L_size * L_size * L_size;
    
    for (int mu = 0; mu < 4; mu++) {
        neighbor_plus[mu].resize(volume);
        neighbor_minus[mu].resize(volume);
    }
    link_index.resize(volume * 4);
    
    for (int it = 0; it < T_size; it++) {
        for (int ix = 0; ix < L_size; ix++) {
            for (int iy = 0; iy < L_size; iy++) {
                for (int iz = 0; iz < L_size; iz++) {
                    int site = get_index(it, ix, iy, iz, T_size, L_size);
                    
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
    std::string config_dir;
    
    int T;
    int L;
    double beta;
    
    int seed;
    std::string start_type;
    std::string boundary;
    int num_sweeps;
    int save_interval;
};

bool read_input_file(const char *filename, SimParams &params) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open input file: " << filename << std::endl;
        return false;
    }
    
    params.output_dir = "output/";
    params.config_dir = "output/configs/";
    params.T = 8;
    params.L = 8;
    params.beta = 2.5;
    params.seed = 12345;
    params.start_type = "cold";
    params.boundary = "periodic";
    params.num_sweeps = 100;
    params.save_interval = 10;
    
    std::string line, key;
    
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        iss >> key;
        
        if (key == "output_dir") {
            iss >> params.output_dir;
        } else if (key == "config_dir") {
            iss >> params.config_dir;
        } else if (key == "T") {
            iss >> params.T;
        } else if (key == "L") {
            iss >> params.L;
        } else if (key == "beta") {
            iss >> params.beta;
        } else if (key == "seed") {
            iss >> params.seed;
        } else if (key == "start_type") {
            iss >> params.start_type;
        } else if (key == "boundary") {
            iss >> params.boundary;
        } else if (key == "num_sweeps") {
            iss >> params.num_sweeps;
        } else if (key == "save_interval") {
            iss >> params.save_interval;
        }
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
    return true;
}

void print_params(const SimParams &params) {
    std::cout << "========================================" << std::endl;
    std::cout << "SU(2) Heatbath Monte Carlo" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Output directory:    " << params.output_dir << std::endl;
    std::cout << "Config directory:    " << params.config_dir << std::endl;
    std::cout << "Lattice size:        " << params.T << " x " << params.L << "^3" << std::endl;
    std::cout << "Beta:                " << params.beta << std::endl;
    std::cout << "Seed:                " << params.seed << std::endl;
    std::cout << "Start type:          " << params.start_type << std::endl;
    std::cout << "Boundary conditions: " << params.boundary << std::endl;
    std::cout << "Number of sweeps:    " << params.num_sweeps << std::endl;
    std::cout << "Save interval:       " << params.save_interval << std::endl;
    std::cout << "========================================" << std::endl;
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
        std::cerr << "Input file parameters:" << std::endl;
        std::cerr << "  output_dir      <path>       # Base output directory" << std::endl;
        std::cerr << "  config_dir      <path>       # Config output directory" << std::endl;
        std::cerr << "  T               <int>        # Temporal lattice extent" << std::endl;
        std::cerr << "  L               <int>        # Spatial lattice extent" << std::endl;
        std::cerr << "  beta            <double>     # Coupling constant" << std::endl;
        std::cerr << "  seed            <int>        # Random seed" << std::endl;
        std::cerr << "  start_type      cold|hot    # Initial configuration" << std::endl;
        std::cerr << "  boundary        periodic|open # Boundary conditions" << std::endl;
        std::cerr << "  num_sweeps      <int>        # Total MC sweeps" << std::endl;
        std::cerr << "  save_interval   <int>        # Save config every N sweeps" << std::endl;
        return EXIT_FAILURE;
    }
    
    SimParams params;
    if (!read_input_file(input_file, params)) {
        return EXIT_FAILURE;
    }
    if (!validate_params(params)) {
        return EXIT_FAILURE;
    }
    
    print_params(params);
    
    T = params.T;
    L = params.L;
    hot_start = (params.start_type == "hot");
    open_boundary_conditions = (params.boundary == "open");
    
    InitializeRand(params.seed);
    
    // Precompute neighbor tables for fast index lookup
    init_neighbor_tables(T, L);
    
    mkdir(params.output_dir.c_str(), 0755);
    mkdir(params.config_dir.c_str(), 0755);
    
    Gauge_Field_Alloc(&gauge_field, T, L);
    
    if (hot_start) {
        std::cout << "Initializing hot start..." << std::endl;
        Gauge_Field_Random(gauge_field, T, L);
    } else {
        std::cout << "Initializing cold start..." << std::endl;
        Gauge_Field_Unity(gauge_field, T, L);
    }
    
    double P = Average_Plaquette(gauge_field, T, L);
    std::cout << "Initial <P> = " << P << std::endl;
    
    std::string plaq_filename = params.output_dir + "plaquette.dat";
    FILE *plaq_file = fopen(plaq_filename.c_str(), "w");
    if (plaq_file == nullptr) {
        std::cerr << "Error: Cannot open plaquette file: " << plaq_filename << std::endl;
        return EXIT_FAILURE;
    }
    
    fprintf(plaq_file, "# Sweep  Plaquette\n");
    fprintf(plaq_file, "%5d %+.10e\n", 0, P);
    
    alignas(32) double SU2_1[8], SU2_2[8];
    const int volume = T * L * L * L;
    
    for (int sweep = 1; sweep <= params.num_sweeps; sweep++) {
        
        for (int site = 0; site < volume; site++) {
            const int it = site / (L * L * L);
            
            for (int mu = 0; mu < 4; mu++) {
                
                alignas(32) double S_l[8];
                cm_eq_zero(S_l);
                
                for (int nu = 0; nu < 4; nu++) {
                    if (nu == mu) continue;
                    
                    // Lower staple: U_nu^dag(x-nu) * U_mu(x-nu) * U_nu(x-nu+mu)
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
                    
                    // Upper staple: U_nu(x) * U_mu(x+nu) * U_nu^dag(x+mu)
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
        
        P = Average_Plaquette(gauge_field, T, L);
        fprintf(plaq_file, "%5d %+.10e\n", sweep, P);
        fflush(plaq_file);
        
        double progress = static_cast<double>(sweep) / params.num_sweeps;
        progress_bar(progress);
        
        if (sweep % params.save_interval == 0) {
            char config_filename[1024];
            snprintf(config_filename, sizeof(config_filename), 
                     "%sconf.%04d", params.config_dir.c_str(), sweep);
            
            char header[1024];
            snprintf(header, sizeof(header), 
                     "T=%d L=%d beta=%.6f sweep=%d", T, L, params.beta, sweep);
            
            write_gauge_field(gauge_field, config_filename, T, L, header);
        }
    }
    
    progress_bar_clear();
    std::cout << "Completed " << params.num_sweeps << " sweeps, final <P> = " << P << std::endl;
    
    fclose(plaq_file);
    
    Gauge_Field_Free(&gauge_field);
    
    std::cout << "========================================" << std::endl;
    std::cout << "Simulation complete." << std::endl;
    std::cout << "Plaquette history: " << plaq_filename << std::endl;
    std::cout << "Configurations:    " << params.config_dir << std::endl;
    std::cout << "========================================" << std::endl;
    
    return EXIT_SUCCESS;
}
