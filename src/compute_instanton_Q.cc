// ********************
// compute_instanton_Q.cc
//
// Driver to compute Q vs smearing steps for instanton and anti-instanton
// ********************

#include "../include/instanton.hh"
#include "../_Utility/include/linear_algebra.hh"
#include "../_Utility/include/geometry.hh"
#include "../_Utility/include/smearing_techniques.hh"
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
              << "  -nsmear <int>   Max smearing steps (default: 100)\n"
              << "  -alpha <double> APE smearing alpha (default: 0.5)\n"
              << "  -o <prefix>     Output file prefix (default: instanton_Q)\n"
              << "  -h              Print this help\n";
}

int main(int argc, char *argv[]) {
    // Default parameters
    int L = 16;
    int T = 16;
    double rho = 3.0;
    double a = 1.0;
    int max_smear_steps = 100;
    double smearing_alpha = 0.5;
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
        } else if (std::strcmp(argv[i], "-nsmear") == 0 && i + 1 < argc) {
            max_smear_steps = std::atoi(argv[++i]);
        } else if (std::strcmp(argv[i], "-alpha") == 0 && i + 1 < argc) {
            smearing_alpha = std::atof(argv[++i]);
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
    std::cout << "  Max smearing steps: " << max_smear_steps << "\n";
    std::cout << "  Smearing alpha: " << smearing_alpha << "\n";
    std::cout << "\n";
    
    const std::size_t volume = static_cast<std::size_t>(T) * L * L * L;
    
    // Allocate gauge fields for instanton and anti-instanton
    std::vector<double> gauge_inst(volume * 4 * 8, 0.0);
    std::vector<double> gauge_anti(volume * 4 * 8, 0.0);
    
    std::cout << "Building instanton configuration...\n";
    insert_instanton(gauge_inst.data(), rho, a, L, T, false);
    
    std::cout << "Building anti-instanton configuration...\n";
    insert_instanton(gauge_anti.data(), rho, a, L, T, true);
    
    // Open output files
    std::string inst_file = output_prefix + "_instanton.dat";
    std::string anti_file = output_prefix + "_anti_instanton.dat";
    
    std::ofstream out_inst(inst_file);
    std::ofstream out_anti(anti_file);
    
    if (!out_inst || !out_anti) {
        std::cerr << "Error: Could not open output files\n";
        return 1;
    }
    
    // Write headers
    out_inst << "# Topological charge Q vs smearing steps for INSTANTON\n";
    out_inst << "# rho = " << rho << ", L = " << L << ", T = " << T 
             << ", alpha = " << smearing_alpha << "\n";
    out_inst << "# smear_step Q\n";
    
    out_anti << "# Topological charge Q vs smearing steps for ANTI-INSTANTON\n";
    out_anti << "# rho = " << rho << ", L = " << L << ", T = " << T 
             << ", alpha = " << smearing_alpha << "\n";
    out_anti << "# smear_step Q\n";
    
    std::cout << "Computing Q vs smearing steps...\n";
    
    // Measure Q at each smearing step
    for (int step = 0; step <= max_smear_steps; ++step) {
        double Q_inst = compute_topological_charge(gauge_inst.data(), T, L);
        double Q_anti = compute_topological_charge(gauge_anti.data(), T, L);
        
        out_inst << step << " " << Q_inst << "\n";
        out_anti << step << " " << Q_anti << "\n";
        
        if (step % 10 == 0) {
            std::cout << "  Step " << step << ": Q_inst = " << Q_inst 
                      << ", Q_anti = " << Q_anti << "\n";
        }
        
        if (step < max_smear_steps) {
            APE_Smearing_Step(gauge_inst.data(), T, L, smearing_alpha);
            APE_Smearing_Step(gauge_anti.data(), T, L, smearing_alpha);
        }
    }
    
    out_inst.close();
    out_anti.close();
    
    std::cout << "\nOutput written to:\n";
    std::cout << "  " << inst_file << "\n";
    std::cout << "  " << anti_file << "\n";
    std::cout << "\nDone.\n";
    
    return 0;
}
