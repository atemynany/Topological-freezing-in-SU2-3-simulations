#include "../include/instanton.hh"
#include "../_Utility/include/linear_algebra.hh"
#include "../_Utility/include/geometry.hh"
#include "../_Utility/include/smearing_techniques.hh"
#include "../include/topcharge_su2.hh"

#include <cmath>
#include <iostream>
#include <vector>
#include <array>

static inline int bar_eta(int a, int mu, int nu) {
    int val = 0;
    if (mu > 0 && nu > 0)
        val += levi3(a, mu - 1, nu - 1);
    if (mu > 0 && nu == 0)
        val += (a == mu - 1) ? -1 : 0;
    if (mu == 0 && nu > 0)
        val += (a == nu - 1) ? 1 : 0;
    return val;
}


static inline void compute_Amu_x4(std::array<double, 3> &A_mu, 
                               const std::array<double, 4> &x, 
                               int mu, double rho2, bool is_anti) {
    double x2 = 0.0;
    for (int k = 0; k < 4; ++k)
        x2 += x[k] * x[k];
    
    const double prefactor = 2.0 * rho2 / (x2 * (x2 + rho2));
    
    for (int a = 0; a < 3; ++a) {
        double sum = 0.0;
        for (int nu = 0; nu < 4; ++nu) {
            sum += (is_anti ? eta(a, mu, nu) : bar_eta(a, mu, nu)) * x[nu];
        }
        A_mu[a] = prefactor * sum;
    }
}

static inline void compute_Amu_x2(std::array<double, 3> &A_mu, 
                               const std::array<double, 4> &x, 
                               int mu, double rho2, bool is_anti) {
    double x2 = 0.0;
    for (int k = 0; k < 4; ++k)
        x2 += x[k] * x[k];
    
    const double prefactor = 2.0 / (x2 + rho2);
    
    for (int a = 0; a < 3; ++a) {
        double sum = 0.0;
        for (int nu = 0; nu < 4; ++nu) {
            sum += (is_anti ? eta(a, mu, nu) : bar_eta(a, mu, nu)) * x[nu];
        }
        A_mu[a] = prefactor * sum;
    }
}

static const int N_SEGMENTS = 8;

void insert_instanton(double *gauge_field, double rho, double a, int L, int T, bool is_anti) {
    const double rho2 = rho * rho;
    const double a_over_N = a / static_cast<double>(N_SEGMENTS);
    
    const double center_t = T / 2.0 + 0.5;
    const double center_x = L / 2.0 + 0.5;
    const double center_y = L / 2.0 + 0.5;
    const double center_z = L / 2.0 + 0.5;
    
    for (int t = 0; t < T; ++t) {
        for (int x = 0; x < L; ++x) {
            for (int y = 0; y < L; ++y) {
                for (int z = 0; z < L; ++z) {
                    std::array<double, 4> coords = {
                        static_cast<double>(t) - center_t,
                        static_cast<double>(x) - center_x,
                        static_cast<double>(y) - center_y,
                        static_cast<double>(z) - center_z
                    };
                    
                    const int site = ((t * L + x) * L + y) * L + z;
                    
                    for (int mu = 0; mu < 4; ++mu) {
                        double U[8];
                        cm_eq_id(U);
                        
                        for (int seg = 0; seg < N_SEGMENTS; ++seg) {
                            std::array<double, 4> x_seg = coords;
                            x_seg[mu] += (seg + 0.5) / static_cast<double>(N_SEGMENTS);
                            
                            std::array<double, 3> A_seg;
                            compute_Amu_x2(A_seg, x_seg, mu, rho2, is_anti);
                            
                            A_seg[0] *= a_over_N;
                            A_seg[1] *= a_over_N;
                            A_seg[2] *= a_over_N;
                            
                            double dU[8], tmp[8];
                            su2_exp_from_Amu(dU, A_seg);
                            cm_eq_cm_ti_cm(tmp, U, dU);
                            cm_eq_cm(U, tmp);
                        }
                        
                        cm_proj(U);
                        
                        const int idx = (site * 4 + mu) * 8;
                        for (int i = 0; i < 8; ++i) {
                            gauge_field[idx + i] = U[i];
                        }
                    }
                }
            }
        }
    }
}

double instanton_topological_charge(double rho, double a, int L, int T,
                                    int smearing_steps, double smearing_alpha, bool is_anti) {
    const std::size_t volume = static_cast<std::size_t>(T) * L * L * L;
    
    std::vector<double> gauge_field(volume * 4 * 8, 0.0);
    
    insert_instanton(gauge_field.data(), rho, a, L, T, is_anti);
    
    for (int step = 0; step < smearing_steps; ++step) {
        APE_Smearing_Step(gauge_field.data(), T, L, smearing_alpha);
    }
    
    return compute_topological_charge(gauge_field.data(), T, L);
}

void instanton_Q_vs_smearing(double rho, double a, int L, int T,
                             int max_smear_steps, double smearing_alpha, bool is_anti) {
    const std::size_t volume = static_cast<std::size_t>(T) * L * L * L;
    
    std::vector<double> gauge_field(volume * 4 * 8, 0.0);
    
    insert_instanton(gauge_field.data(), rho, a, L, T, is_anti);
    
    const char* type = is_anti ? "anti-instanton" : "instanton";
    
    std::cout << "# Q vs smearing steps for " << type << "\n";
    std::cout << "# rho = " << rho << ", L = " << L << ", T = " << T << "\n";
    std::cout << "# smear_step Q\n";
    
    for (int step = 0; step <= max_smear_steps; ++step) {
        double Q = compute_topological_charge(gauge_field.data(), T, L);
        std::cout << step << " " << Q << "\n";
        
        if (step < max_smear_steps) {
            APE_Smearing_Step(gauge_field.data(), T, L, smearing_alpha);
        }
    }
}
