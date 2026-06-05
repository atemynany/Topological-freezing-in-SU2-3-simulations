// meas_topcharge_su3.cc - SU(3) topological charge measurement
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
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
std::vector<long long> link_index_su3;

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
                        link_index_su3[site * 4 + mu] = ((long long)4 * site + mu) * 18;
                    }
                }
            }
        }
    }
}

static uint32_t read_be32(const char *p) {
    const unsigned char *b = reinterpret_cast<const unsigned char *>(p);
    return (static_cast<uint32_t>(b[0]) << 24) |
           (static_cast<uint32_t>(b[1]) << 16) |
           (static_cast<uint32_t>(b[2]) << 8) |
            static_cast<uint32_t>(b[3]);
}

static uint64_t read_be64(const char *p) {
    const unsigned char *b = reinterpret_cast<const unsigned char *>(p);
    uint64_t value = 0;
    for (int i = 0; i < 8; i++) {
        value = (value << 8) | static_cast<uint64_t>(b[i]);
    }
    return value;
}

static double read_be_double(const char *p) {
    const uint64_t raw = read_be64(p);
    double value;
    std::memcpy(&value, &raw, sizeof(value));
    return value;
}

bool read_su3_config_raw(double *gf, const char *filename, int T, int L) {
    FILE *f = fopen(filename, "rb");
    if (!f) {
        std::cerr << "Error: Cannot read " << filename << std::endl;
        return false;
    }
    char line[1024];
    if (fgets(line, sizeof(line), f) == nullptr) {
        std::cerr << "Error: Empty config file: " << filename << std::endl;
        fclose(f);
        return false;
    }
    const int volume = T * L * L * L;
    size_t nread = fread(gf, sizeof(double), volume * 4 * 18, f);
    fclose(f);
    if (nread != static_cast<size_t>(volume) * 4 * 18) {
        std::cerr << "Error: Config file too short: " << filename << std::endl;
        return false;
    }
    return true;
}

bool read_su3_config_cl2qcd_lime(double *gf, const char *filename, int T, int L) {
    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open()) {
        std::cerr << "Error: Cannot read " << filename << std::endl;
        return false;
    }

    const size_t volume = static_cast<size_t>(T) * L * L * L;
    const size_t expected_bytes = volume * 4 * 18 * sizeof(double);

    while (true) {
        char header[144];
        f.read(header, sizeof(header));
        if (f.gcount() == 0 && f.eof()) break;
        if (f.gcount() != static_cast<std::streamsize>(sizeof(header))) {
            std::cerr << "Error: Truncated LIME header in " << filename << std::endl;
            return false;
        }

        const uint32_t magic = read_be32(header);
        if (magic != 0x456789abU) {
            std::cerr << "Error: Bad LIME magic in " << filename << std::endl;
            return false;
        }

        const uint64_t payload_bytes = read_be64(header + 8);
        std::string record_type(header + 16, strnlen(header + 16, 128));
        const uint64_t padding = (8 - (payload_bytes % 8)) % 8;

        if (record_type == "ildg-binary-data") {
            if (payload_bytes != expected_bytes) {
                std::cerr << "Error: Unexpected ildg-binary-data size in " << filename
                          << " (got " << payload_bytes << ", expected " << expected_bytes << ")" << std::endl;
                return false;
            }

            std::vector<char> payload(expected_bytes);
            f.read(payload.data(), static_cast<std::streamsize>(payload.size()));
            if (f.gcount() != static_cast<std::streamsize>(payload.size())) {
                std::cerr << "Error: Truncated ildg-binary-data payload in " << filename << std::endl;
                return false;
            }

            for (size_t site = 0; site < volume; site++) {
                for (int file_dir = 0; file_dir < 4; file_dir++) {
                    const int mu = (file_dir + 1) % 4;  // CL2QCD writes ILDG dirs as x,y,z,t.
                    const size_t src_link = site * 4 + static_cast<size_t>(file_dir);
                    const size_t dst_link = site * 4 + static_cast<size_t>(mu);
                    for (int elem = 0; elem < 18; elem++) {
                        gf[dst_link * 18 + elem] =
                            read_be_double(payload.data() + (src_link * 18 + elem) * sizeof(double));
                    }
                }
            }

            return true;
        }

        f.seekg(static_cast<std::streamoff>(payload_bytes + padding), std::ios::cur);
        if (!f.good()) {
            std::cerr << "Error: Failed while skipping LIME record " << record_type
                      << " in " << filename << std::endl;
            return false;
        }
    }

    std::cerr << "Error: No ildg-binary-data record found in " << filename << std::endl;
    return false;
}

// Simple APE smearing for SU(3)
void su3_ape_smear(double *gf_out, const double *gf_in, int T, int L, double alpha) {
    const int volume = T * L * L * L;

    // Double-buffered: reads gf_in, writes gf_out. Each (site, mu) writes a unique
    // link, so different iterations never collide — embarrassingly parallel.
    #pragma omp parallel for schedule(static)
    for (int site = 0; site < volume; site++) {
        alignas(32) double staple[18], smeared[18], tmp[18];
        for (int mu = 0; mu < 4; mu++) {
            su3_eq_zero(staple);
            
            // Sum staples
            for (int nu = 0; nu < 4; nu++) {
                if (nu == mu) continue;
                
                // Upper staple
                long long idx_nu = link_index_su3[site * 4 + nu];
                int site_pnu = neighbor_plus[nu][site];
                long long idx_mu_pnu = link_index_su3[site_pnu * 4 + mu];
                int site_pmu = neighbor_plus[mu][site];
                long long idx_nu_pmu = link_index_su3[site_pmu * 4 + nu];
                
                su3_eq_su3_ti_su3(tmp, gf_in + idx_nu, gf_in + idx_mu_pnu);
                su3_eq_su3_ti_su3_dag(smeared, tmp, gf_in + idx_nu_pmu);
                su3_pl_eq_su3(staple, smeared);
                
                // Lower staple
                int site_mnu = neighbor_minus[nu][site];
                long long idx_nu_mnu = link_index_su3[site_mnu * 4 + nu];
                long long idx_mu_mnu = link_index_su3[site_mnu * 4 + mu];
                int site_mnu_pmu = neighbor_plus[mu][site_mnu];
                long long idx_nu_mnu_pmu = link_index_su3[site_mnu_pmu * 4 + nu];
                
                su3_eq_su3_dag_ti_su3(tmp, gf_in + idx_nu_mnu, gf_in + idx_mu_mnu);
                su3_eq_su3_ti_su3(smeared, tmp, gf_in + idx_nu_mnu_pmu);
                su3_pl_eq_su3(staple, smeared);
            }
            
            // U' = (1-alpha)*U + (alpha/6)*staple, then project
            long long idx = link_index_su3[site * 4 + mu];
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
    std::string config_format;
    std::string config_prefix;
    std::string config_postfix;
    std::string output_file;
    int T, L;
    int start_conf, end_conf, conf_step;
    int conf_digits;
    int smear_steps;
    double smear_alpha;
    std::string boundary;          // "periodic" or "open"
    int exclude_boundary_slices;   // Time slices to exclude from each boundary (for open BC)
};

bool read_input(const char *filename, MeasParams &p) {
    std::ifstream f(filename);
    if (!f.is_open()) return false;
    
    p.config_dir = "output/configs_su3/";
    p.config_format = "raw";
    p.config_prefix = "conf_su3.";
    p.config_postfix = "";
    p.output_file = "output/topcharge_su3.dat";
    p.T = 8; p.L = 8;
    p.start_conf = 10; p.end_conf = 100; p.conf_step = 10;
    p.conf_digits = 4;
    p.smear_steps = 20; p.smear_alpha = 0.45;
    p.boundary = "periodic";
    p.exclude_boundary_slices = 0;
    
    std::string line, key;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        iss >> key;
        if (key == "config_dir") iss >> p.config_dir;
        else if (key == "config_format") iss >> p.config_format;
        else if (key == "config_prefix") iss >> p.config_prefix;
        else if (key == "config_postfix") iss >> p.config_postfix;
        else if (key == "output_file") iss >> p.output_file;
        else if (key == "T") iss >> p.T;
        else if (key == "L") iss >> p.L;
        else if (key == "start_conf") iss >> p.start_conf;
        else if (key == "end_conf") iss >> p.end_conf;
        else if (key == "conf_step") iss >> p.conf_step;
        else if (key == "conf_digits") iss >> p.conf_digits;
        else if (key == "smear_steps") iss >> p.smear_steps;
        else if (key == "smear_alpha") iss >> p.smear_alpha;
        else if (key == "boundary") iss >> p.boundary;
        else if (key == "exclude_boundary_slices") iss >> p.exclude_boundary_slices;
    }
    return true;
}

std::string config_filename(const MeasParams &params, int conf) {
    std::ostringstream path;
    path << params.config_dir;
    if (!params.config_dir.empty() && params.config_dir.back() != '/') {
        path << '/';
    }
    path << params.config_prefix << std::setw(params.conf_digits) << std::setfill('0')
         << conf << params.config_postfix;
    return path.str();
}

bool read_su3_config(double *gf, const std::string &filename, const MeasParams &params) {
    if (params.config_format == "raw" || params.config_format == "native") {
        return read_su3_config_raw(gf, filename.c_str(), params.T, params.L);
    }
    if (params.config_format == "cl2qcd_lime" ||
        params.config_format == "ildg" ||
        params.config_format == "lime") {
        return read_su3_config_cl2qcd_lime(gf, filename.c_str(), params.T, params.L);
    }

    std::cerr << "Error: Unknown config_format: " << params.config_format << std::endl;
    return false;
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
    double *gf = nullptr;
    double *gf_smeared = nullptr;
    double *gf_tmp = nullptr;
    posix_memalign((void **)&gf, 32, volume * 4 * 18 * sizeof(double));
    posix_memalign((void **)&gf_smeared, 32, volume * 4 * 18 * sizeof(double));
    posix_memalign((void **)&gf_tmp, 32, volume * 4 * 18 * sizeof(double));
    
    FILE *out = fopen(params.output_file.c_str(), "w");
    if (!out) {
        std::cerr << "Error: Cannot open output file: " << params.output_file << std::endl;
        free(gf); free(gf_smeared); free(gf_tmp);
        return EXIT_FAILURE;
    }
    fprintf(out, "# conf smear_step Q\n");

    // Timeslice output: same directory as output_file, fixed name
    std::string ts_filename = params.output_file;
    size_t slash = ts_filename.rfind('/');
    ts_filename = (slash != std::string::npos)
                  ? ts_filename.substr(0, slash + 1) + "topcharge_timeslice.dat"
                  : "topcharge_timeslice.dat";
    std::ofstream timeslice_file(ts_filename);
    if (!timeslice_file.is_open()) {
        std::cerr << "Error: Cannot open timeslice output file: " << ts_filename << std::endl;
        return EXIT_FAILURE;
    }
    timeslice_file << "# smear_steps  config_number  t  q_t" << std::endl;
    
    std::cout << "SU(3) Topological Charge Measurement\n";
    std::cout << "Lattice: " << T_size << "x" << L_size << "^3\n";
    std::cout << "Config format: " << params.config_format << "\n";
    std::cout << "Smearing: " << params.smear_steps << " steps, alpha=" << params.smear_alpha << "\n";
    std::cout << "Boundary: " << params.boundary;
    if (params.boundary == "open" && params.exclude_boundary_slices > 0) {
        std::cout << " (excluding " << params.exclude_boundary_slices << " slices from each end)";
    }
    std::cout << "\n";
    
    for (int conf = params.start_conf; conf <= params.end_conf; conf += params.conf_step) {
        const std::string fname = config_filename(params, conf);
        if (!read_su3_config(gf, fname, params)) {
            return EXIT_FAILURE;
        }
        
        // Copy to smeared buffer
        memcpy(gf_smeared, gf, volume * 4 * 18 * sizeof(double));
        
        // Smear and measure
        for (int step = 0; step <= params.smear_steps; step++) {
            double Q;
            if (params.boundary == "open" && params.exclude_boundary_slices > 0) {
                Q = su3_topological_charge_open(gf_smeared, T_size, L_size, params.exclude_boundary_slices);
            } else {
                Q = su3_topological_charge(gf_smeared, T_size, L_size);
            }
            
            if (step == params.smear_steps) {
                fprintf(out, "%d %d %.10e\n", conf, step, Q);
                std::cout << "conf " << conf << " smear " << step << " Q=" << Q << "\n";

                // Topological charge density per time slice: q(t) = (1/4π²) Σ_{x,y,z} q_local(t,x,y,z)
                const int t_min = (params.boundary == "open" && params.exclude_boundary_slices > 0)
                                  ? params.exclude_boundary_slices : 0;
                const int t_max = (params.boundary == "open" && params.exclude_boundary_slices > 0)
                                  ? params.T - params.exclude_boundary_slices : params.T;
                const int t_count = t_max - t_min;

                // Compute q(t) in parallel (one thread per timeslice), then write serially.
                std::vector<double> q_per_t(t_count, 0.0);
                #pragma omp parallel for schedule(static)
                for (int it_idx = 0; it_idx < t_count; it_idx++) {
                    const int it = t_min + it_idx;
                    double q_t = 0.0;
                    for (int ix = 0; ix < params.L; ix++) {
                        for (int iy = 0; iy < params.L; iy++) {
                            for (int iz = 0; iz < params.L; iz++) {
                                q_t += su3_local_topcharge_density(gf_smeared, it, ix, iy, iz);
                            }
                        }
                    }
                    q_per_t[it_idx] = q_t / (4.0 * M_PI * M_PI);
                }

                timeslice_file << std::fixed << std::setprecision(8);
                for (int it_idx = 0; it_idx < t_count; it_idx++) {
                    timeslice_file << std::setw(5)  << params.smear_steps << "  "
                                   << std::setw(6)  << conf << "  "
                                   << std::setw(4)  << (t_min + it_idx) << "  "
                                   << std::setw(14) << q_per_t[it_idx] << "\n";
                }
            }
            
            if (step < params.smear_steps) {
                su3_ape_smear(gf_tmp, gf_smeared, T_size, L_size, params.smear_alpha);
                std::swap(gf_smeared, gf_tmp);
            }
        }
    }
    
    fclose(out);
    timeslice_file.close();
    free(gf); free(gf_smeared); free(gf_tmp);
    
    std::cout << "Output: " << params.output_file << "\n";
    return 0;
}
