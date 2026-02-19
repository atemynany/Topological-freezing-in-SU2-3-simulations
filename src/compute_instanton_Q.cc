// ********************
// compute_instanton_Q.cc
//
// Driver to compute Q for instanton and anti-instanton configurations
// ********************

#include "../include/instanton.hh"
#include "../_Utility/include/linear_algebra.hh"
#include "../_Utility/include/geometry.hh"
#include "../include/topcharge_su2.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>

// Global flag for boundary conditions (required by geometry.hh)
bool open_boundary_conditions = false;

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " [options]\n"
              << "Options:\n"
              << "  -L <int>        Spatial extent (default: 16)\n"
              << "  -T <int>        Temporal extent (default: 16)\n"
              << "  -rho <double>   Instanton size (default: 3.0)\n"
              << "  -a <double>     Lattice spacing (default: 1.0)\n"
              << "  -o <prefix>     Output file prefix (default: instanton_Q)\n"
              << "  -h              Print this help\n";
}

int main(int argc, char *argv[]) {
    // Default parameters
    int L = 16;
    int T = 16;
    double rho = 3.0;
    double a = 1.0;
    std::string output_prefix = "instanton_Q";
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-L") == 0 && i + 1 < argc) {
            L = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "-T") == 0 && i + 1 < argc) {
            T = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "-rho") == 0 && i + 1 < argc) {
            rho = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "-a") == 0 && i + 1 < argc) {
            a = std::atof(argv[++i]);
        } else if (std::strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_prefix = argv[++i];
        } else if (std::strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }
    
    std::cout << "=== Instanton Topological Charge Computation ===\n";
    std::cout << "Parameters:\n";
    std::cout << "  Lattice: " << L << "^3 x " << T << "\n";
    std::cout << "  Instanton size rho: " << rho << "\n";
    std::cout << "  Lattice spacing a: " << a << "\n";
    std::cout << "\n";
    
    const std::size_t volume = static_cast<std::size_t>(T) * L * L * L;
    
    // Allocate gauge fields for instanton and anti-instanton
    std::vector<double> gauge_inst(volume * 4 * 8, 0.0);
    std::vector<double> gauge_anti(volume * 4 * 8, 0.0);
    
    std::cout << "Building instanton configuration...\n";
    insert_instanton(gauge_inst.data(), rho, a, L, T, false);
    
    std::cout << "Building anti-instanton configuration...\n";
    insert_instanton(gauge_anti.data(), rho, a, L, T, true);
    
    // Compute topological charge
    std::cout << "Computing topological charge...\n";
    
    double Q_inst = compute_topological_charge(gauge_inst.data(), T, L);
    double Q_anti = compute_topological_charge(gauge_anti.data(), T, L);
    
    std::cout << "\nResults:\n";
    std::cout << "  Q (instanton):      " << Q_inst << " (expected: +1)\n";
    std::cout << "  Q (anti-instanton): " << Q_anti << " (expected: -1)\n";
    
    // Open output file
    std::string out_file = output_prefix + ".dat";
    std::ofstream out(out_file);
    
    if (!out) {
        std::cerr << "Error: Could not open output file\n";
        return 1;
    }
    
    out << "# Instanton topological charge computation\n";
    out << "# L = " << L << ", T = " << T << ", rho = " << rho << ", a = " << a << "\n";
    out << "# type Q\n";
    out << "instanton " << Q_inst << "\n";
    out << "anti-instanton " << Q_anti << "\n";
    
    out.close();
    
    std::cout << "\nOutput written to: " << out_file << "\n";
    std::cout << "Done.\n";
    
    return 0;
}
