// meas_topcharge_su3.cc - SU(3) topological charge measurement
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "su3_linear_algebra.hh"
#include "su3_heatbath.hh"
#include "topcharge_su3.hh"
#include "smearing_techniques.hh"

int T_size, L_size;
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
                        link_index_su3[site * 4 + mu] = (4 * site + mu) * 18;
                    }
                }
            }
        }
    }
}

void read_su3_config(double *gf, const char *filename, int T, int L) {
    FILE *f = fopen(filename, "rb");
    if (!f) { std::cerr << "Error: Cannot read " << filename << std::endl; return; }
    char line[1024];
    if (fgets(line, sizeof(line), f) == nullptr) { fclose(f); return; }
    const int volume = T * L * L * L;
    size_t nread = fread(gf, sizeof(double), volume * 4 * 18, f);
    (void)nread;
    fclose(f);
}

// Simple APE smearing for SU(3)
void su3_ape_smear(double *gf_out, const double *gf_in, int T, int L, double alpha) {
    const int volume = T * L * L * L;
    alignas(32) double staple[18], smeared[18], tmp[18];
    
    for (int site = 0; site < volume; site++) {
        for (int mu = 0; mu < 4; mu++) {
            su3_eq_zero(staple);
            
            // Sum staples
            for (int nu = 0; nu < 4; nu++) {
                if (nu == mu) continue;
                
                // Upper staple
                int idx_nu = link_index_su3[site * 4 + nu];
                int site_pnu = neighbor_plus[nu][site];
                int idx_mu_pnu = link_index_su3[site_pnu * 4 + mu];
                int site_pmu = neighbor_plus[mu][site];
                int idx_nu_pmu = link_index_su3[site_pmu * 4 + nu];
                
                su3_eq_su3_ti_su3(tmp, gf_in + idx_nu, gf_in + idx_mu_pnu);
                su3_eq_su3_ti_su3_dag(smeared, tmp, gf_in + idx_nu_pmu);
                su3_pl_eq_su3(staple, smeared);
                
                // Lower staple
                int site_mnu = neighbor_minus[nu][site];
                int idx_nu_mnu = link_index_su3[site_mnu * 4 + nu];
                int idx_mu_mnu = link_index_su3[site_mnu * 4 + mu];
                int site_mnu_pmu = neighbor_plus[mu][site_mnu];
                int idx_nu_mnu_pmu = link_index_su3[site_mnu_pmu * 4 + nu];
                
                su3_eq_su3_dag_ti_su3(tmp, gf_in + idx_nu_mnu, gf_in + idx_mu_mnu);
                su3_eq_su3_ti_su3(smeared, tmp, gf_in + idx_nu_mnu_pmu);
                su3_pl_eq_su3(staple, smeared);
            }
            
            // U' = (1-alpha)*U + (alpha/6)*staple, then project
            int idx = link_index_su3[site * 4 + mu];
            for (int i = 0; i < 18; i++) {
                smeared[i] = (1.0 - alpha) * gf_in[idx + i] + (alpha / 6.0) * staple[i];
            }
            su3_proj(smeared);
            su3_eq_su3(gf_out + idx, smeared);
        }
    }
}

struct MeasParams {
    std::string config_dir;
    std::string output_file;
    int T, L;
    int start_conf, end_conf, conf_step;
    int smear_steps;
    double smear_alpha;
};

bool read_input(const char *filename, MeasParams &p) {
    std::ifstream f(filename);
    if (!f.is_open()) return false;
    
    p.config_dir = "output/configs_su3/";
    p.output_file = "output/topcharge_su3.dat";
    p.T = 8; p.L = 8;
    p.start_conf = 10; p.end_conf = 100; p.conf_step = 10;
    p.smear_steps = 20; p.smear_alpha = 0.45;
    
    std::string line, key;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        iss >> key;
        if (key == "config_dir") iss >> p.config_dir;
        else if (key == "output_file") iss >> p.output_file;
        else if (key == "T") iss >> p.T;
        else if (key == "L") iss >> p.L;
        else if (key == "start_conf") iss >> p.start_conf;
        else if (key == "end_conf") iss >> p.end_conf;
        else if (key == "conf_step") iss >> p.conf_step;
        else if (key == "smear_steps") iss >> p.smear_steps;
        else if (key == "smear_alpha") iss >> p.smear_alpha;
    }
    return true;
}

int main(int argc, char **argv) {
    const char *input_file = nullptr;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) input_file = argv[++i];
    }
    
    if (!input_file) {
        std::cerr << "Usage: " << argv[0] << " -i <input_file>\n";
        return 1;
    }
    
    MeasParams params;
    if (!read_input(input_file, params)) return 1;
    
    T_size = params.T;
    L_size = params.L;
    init_neighbor_tables();
    
    const int volume = T_size * L_size * L_size * L_size;
    double *gf = (double *)aligned_alloc(32, volume * 4 * 18 * sizeof(double));
    double *gf_smeared = (double *)aligned_alloc(32, volume * 4 * 18 * sizeof(double));
    double *gf_tmp = (double *)aligned_alloc(32, volume * 4 * 18 * sizeof(double));
    
    FILE *out = fopen(params.output_file.c_str(), "w");
    fprintf(out, "# conf smear_step Q\n");
    
    std::cout << "SU(3) Topological Charge Measurement\n";
    std::cout << "Lattice: " << T_size << "x" << L_size << "^3\n";
    std::cout << "Smearing: " << params.smear_steps << " steps, alpha=" << params.smear_alpha << "\n";
    
    for (int conf = params.start_conf; conf <= params.end_conf; conf += params.conf_step) {
        char fname[512];
        snprintf(fname, sizeof(fname), "%sconf_su3.%04d", params.config_dir.c_str(), conf);
        
        read_su3_config(gf, fname, T_size, L_size);
        
        // Copy to smeared buffer
        memcpy(gf_smeared, gf, volume * 4 * 18 * sizeof(double));
        
        // Smear and measure
        for (int step = 0; step <= params.smear_steps; step++) {
            double Q = su3_topological_charge(gf_smeared, T_size, L_size);
            
            if (step == params.smear_steps) {
                fprintf(out, "%d %d %.10e\n", conf, step, Q);
                std::cout << "conf " << conf << " smear " << step << " Q=" << Q << "\n";
            }
            
            if (step < params.smear_steps) {
                su3_ape_smear(gf_tmp, gf_smeared, T_size, L_size, params.smear_alpha);
                std::swap(gf_smeared, gf_tmp);
            }
        }
    }
    
    fclose(out);
    free(gf); free(gf_smeared); free(gf_tmp);
    
    std::cout << "Output: " << params.output_file << "\n";
    return 0;
}
