// ==============================================================================
// meas_action_density_su2.cc
// ==============================================================================
// Recompute SU(2) Wilson action density from saved configurations, optionally
// using only a temporal bulk window to suppress open-boundary effects.
// ==============================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "Plaquette.hh"
#include "action_density.hh"
#include "fields.hh"
#include "geometry.hh"
#include "io.hh"

int T = 16;
int L = 16;
bool open_boundary_conditions = false;

struct Params {
    std::string config_dir = "output/configs/";
    std::string output_dir = "output/";
    std::string output_file = "action_density_bulk.dat";
    std::string boundary = "periodic";
    double beta = 2.5;
    int T = 16;
    int L = 16;
    int start_conf = 100;
    int end_conf = 100;
    int conf_step = 100;
    int exclude_boundary_slices = 0;
};

static bool read_input_file(const char *filename, Params &p) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: cannot open input file: " << filename << "\n";
        return false;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "config_dir") iss >> p.config_dir;
        else if (key == "output_dir") iss >> p.output_dir;
        else if (key == "output_file") iss >> p.output_file;
        else if (key == "boundary") iss >> p.boundary;
        else if (key == "beta") iss >> p.beta;
        else if (key == "T") iss >> p.T;
        else if (key == "L") iss >> p.L;
        else if (key == "start_conf") iss >> p.start_conf;
        else if (key == "end_conf") iss >> p.end_conf;
        else if (key == "conf_step") iss >> p.conf_step;
        else if (key == "exclude_boundary_slices") iss >> p.exclude_boundary_slices;
    }
    return true;
}

static bool validate_params(const Params &p) {
    if (p.T < 2 || p.L < 2) {
        std::cerr << "Error: invalid lattice size T=" << p.T << " L=" << p.L << "\n";
        return false;
    }
    if (p.start_conf > p.end_conf || p.conf_step <= 0) {
        std::cerr << "Error: invalid configuration range\n";
        return false;
    }
    if (p.boundary != "periodic" && p.boundary != "open") {
        std::cerr << "Error: boundary must be periodic or open\n";
        return false;
    }
    if (p.exclude_boundary_slices < 0 || 2 * p.exclude_boundary_slices >= p.T) {
        std::cerr << "Error: exclude_boundary_slices must leave at least one bulk timeslice\n";
        return false;
    }
    return true;
}

static double average_bulk_plaquette(double *gauge_field, int t_min, int t_max) {
    double sum = 0.0;
    long long count = 0;

    #pragma omp parallel for reduction(+:sum,count) schedule(static)
    for (int it = t_min; it < t_max; it++) {
        for (int ix = 0; ix < L; ix++) {
            for (int iy = 0; iy < L; iy++) {
                for (int iz = 0; iz < L; iz++) {
                    for (int mu = 0; mu < 4; mu++) {
                        for (int nu = mu + 1; nu < 4; nu++) {
                            const bool temporal = (mu == 0 || nu == 0);

                            // Keep only plaquettes fully contained in the retained
                            // temporal bulk. Spatial plaquettes live on one t-slice;
                            // temporal plaquettes connect t to t+1.
                            if (temporal && it + 1 >= t_max) continue;
                            if (open_boundary_conditions && temporal && it == T - 1) continue;

                            sum += plaquette_trace(gauge_field, it, ix, iy, iz, mu, nu, T, L);
                            count++;
                        }
                    }
                }
            }
        }
    }

    return count > 0 ? sum / static_cast<double>(count) : 0.0;
}

int main(int argc, char *argv[]) {
    const char *input_file = nullptr;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            input_file = argv[++i];
        }
    }
    if (input_file == nullptr) {
        std::cerr << "Usage: " << argv[0] << " -i <input_file>\n";
        return EXIT_FAILURE;
    }

    Params params;
    if (!read_input_file(input_file, params) || !validate_params(params)) {
        return EXIT_FAILURE;
    }

    T = params.T;
    L = params.L;
    open_boundary_conditions = (params.boundary == "open");

    const int t_min = params.exclude_boundary_slices;
    const int t_max = params.T - params.exclude_boundary_slices;

    std::string output_path = params.output_dir + params.output_file;
    std::ofstream out(output_path);
    if (!out.is_open()) {
        std::cerr << "Error: cannot open output file: " << output_path << "\n";
        return EXIT_FAILURE;
    }

    out << "# config  plaquette_bulk  action_density_bulk  t_min  t_max_exclusive\n";

    double *gauge_field = nullptr;
    Gauge_Field_Alloc(&gauge_field, T, L);

    int processed = 0;
    for (int n = params.start_conf; n <= params.end_conf; n += params.conf_step) {
        char config_filename[1024];
        std::snprintf(config_filename, sizeof(config_filename), "%sconf.%04d",
                      params.config_dir.c_str(), n);

        std::cout << "Reading " << config_filename << "\n";
        read_gauge_field(gauge_field, config_filename, T, L);

        const double plaq = average_bulk_plaquette(gauge_field, t_min, t_max);
        const double action = avg_action_density(params.beta, plaq);

        out << std::setw(6) << n << "  "
            << std::scientific << std::setprecision(10) << plaq << "  "
            << std::scientific << std::setprecision(10) << action << "  "
            << t_min << "  " << t_max << "\n";
        out.flush();
        processed++;
    }

    Gauge_Field_Free(&gauge_field);
    std::cout << "Processed " << processed << " configurations\n";
    std::cout << "Output written to " << output_path << "\n";
    return EXIT_SUCCESS;
}
