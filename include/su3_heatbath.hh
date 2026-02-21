// su3_heatbath.hh - SU(3) heatbath with Cabibbo-Marinari
//
// Algorithm reference: CL2QCD (https://github.com/AG-Philipsen/cl2qcd)
//   - Kennedy-Pendleton: host_functionality/host_random.hpp, ocl_kernel/operations_heatbath.cl
//   - Cabibbo-Marinari decomposition with 3 SU(2) subgroups
//   - Effective coupling: beta_eff = (2*beta/Nc) * k where k = |det(W_sub)|^(1/2)
//
// References:
//   - A.D. Kennedy, B.J. Pendleton, Phys. Lett. B 156 (1985) 393
//   - N. Cabibbo, E. Marinari, Phys. Lett. B 119 (1982) 387
//
#ifndef __SU3_HEATBATH_HH__
#define __SU3_HEATBATH_HH__

#include "su3_linear_algebra.hh"
#include <cmath>

// Kennedy-Pendleton heatbath for SU(2)
// Algorithm from CL2QCD: ocl_kernel/operations_heatbath.cl::SU2Update()
// Generates a0 from P(a0) ∝ sqrt(1-a0²) * exp(2*alpha*a0)
inline void su2_kennedy_pendleton(double *a, double k, double beta, double (*drand)()) {
    double alpha = beta * k;  // effective coupling
    
    if (alpha < 1e-10) {
        // Random SU(2) for small coupling
        double x1, x2, x3, x4, norm;
        while (true) {
            x1 = 2.0 * drand() - 1.0;
            x2 = 2.0 * drand() - 1.0;
            x3 = 2.0 * drand() - 1.0;
            x4 = 2.0 * drand() - 1.0;
            norm = x1*x1 + x2*x2 + x3*x3 + x4*x4;
            if (norm > 0.01 && norm < 1.0) break;
        }
        norm = 1.0 / sqrt(norm);
        a[0] = x1 * norm;
        a[1] = x2 * norm;
        a[2] = x3 * norm;
        a[3] = x4 * norm;
        return;
    }
    
    // Kennedy-Pendleton: delta = -log(r1)/alpha * cos²(2πr2) - log(r3)/alpha
    // Accept if (1 - delta/2) > eta²
    double delta, a0;
    while (true) {
        double r1 = drand(), r2 = drand(), r3 = drand(), eta = drand();
        double c = cos(2.0 * M_PI * r2);
        delta = -log(r1) / alpha * c * c - log(r3) / alpha;
        a0 = 1.0 - delta;
        if ((1.0 - 0.5 * delta) > eta * eta) break;
    }
    
    // Uniform on sphere for (a1, a2, a3)
    double r = sqrt(1.0 - a0 * a0);
    double phi = 2.0 * M_PI * drand();
    double theta = asin(2.0 * drand() - 1.0);  // CL2QCD formula
    
    a[0] = a0;
    a[1] = r * cos(theta) * cos(phi);
    a[2] = r * cos(theta) * sin(phi);
    a[3] = r * sin(theta);
}

// Cabibbo-Marinari heatbath - updates U using staple
inline void su3_heatbath_link(double *U, const double *staple, double beta, double (*drand)()) {
    alignas(32) double W[18], R[18], T1[18];
    
    // W = U * Staple^dag (the combination we want to maximize trace of)
    su3_eq_su3_ti_su3_dag(W, U, staple);
    
    // Three SU(2) subgroups
    const int subgroups[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    
    for (int s = 0; s < 3; s++) {
        int i = subgroups[s][0];
        int j = subgroups[s][1];
        
        // Extract 2x2 block from W
        double w00_re = W[2*(i*3 + i)], w00_im = W[2*(i*3 + i) + 1];
        double w01_re = W[2*(i*3 + j)], w01_im = W[2*(i*3 + j) + 1];
        double w10_re = W[2*(j*3 + i)], w10_im = W[2*(j*3 + i) + 1];
        double w11_re = W[2*(j*3 + j)], w11_im = W[2*(j*3 + j) + 1];
        
        // Project to SU(2): V = w_block, k = sqrt(det(V))
        // SU(2) parametrization: V00 = v0 + i*v3, V01 = v2 + i*v1, V10 = -v2 + i*v1, V11 = v0 - i*v3
        double v0 = 0.5 * (w00_re + w11_re);
        double v3 = 0.5 * (w00_im - w11_im);
        double v1 = 0.5 * (w01_im + w10_im);
        double v2 = 0.5 * (w01_re - w10_re);
        
        double k = sqrt(v0*v0 + v1*v1 + v2*v2 + v3*v3);
        if (k < 1e-10) continue;
        
        // Normalize V: V_norm = V/k is in SU(2)
        double inv_k = 1.0 / k;
        v0 *= inv_k; v1 *= inv_k; v2 *= inv_k; v3 *= inv_k;
        
        // Generate new SU(2) element a[] via heatbath
        // For SU(3) Wilson action: effective coupling is (beta/3) * 2 * k = (2*beta/3) * k
        double a[4];
        su2_kennedy_pendleton(a, k, 2.0 * beta / 3.0, drand);
        
        // R = a * V^dag (in SU(2))
        // V^dag has v0, -v1, -v2, -v3
        double r0 = a[0]*v0 + a[1]*v1 + a[2]*v2 + a[3]*v3;
        double r1 = a[1]*v0 - a[0]*v1 + a[2]*v3 - a[3]*v2;
        double r2 = a[2]*v0 - a[0]*v2 + a[3]*v1 - a[1]*v3;
        double r3 = a[3]*v0 - a[0]*v3 + a[1]*v2 - a[2]*v1;
        
        // Embed R into SU(3)
        su3_eq_id(R);
        R[2*(i*3 + i)]     = r0;  R[2*(i*3 + i) + 1] = r3;
        R[2*(i*3 + j)]     = r2;  R[2*(i*3 + j) + 1] = r1;
        R[2*(j*3 + i)]     = -r2; R[2*(j*3 + i) + 1] = r1;
        R[2*(j*3 + j)]     = r0;  R[2*(j*3 + j) + 1] = -r3;
        
        // W = R * W
        su3_eq_su3_ti_su3(T1, R, W);
        su3_eq_su3(W, T1);
        
        // U = R * U
        su3_eq_su3_ti_su3(T1, R, U);
        su3_eq_su3(U, T1);
    }
    
    su3_proj(U);
}

// SU(3) overrelaxation
inline void su3_overrelax_link(double *U, const double *staple) {
    alignas(32) double W[18], T1[18], R[18];
    
    su3_eq_su3_ti_su3_dag(W, U, staple);
    
    const int subgroups[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    
    for (int s = 0; s < 3; s++) {
        int i = subgroups[s][0];
        int j = subgroups[s][1];
        
        double w00_re = W[2*(i*3 + i)], w00_im = W[2*(i*3 + i) + 1];
        double w01_re = W[2*(i*3 + j)], w01_im = W[2*(i*3 + j) + 1];
        double w10_re = W[2*(j*3 + i)], w10_im = W[2*(j*3 + i) + 1];
        double w11_re = W[2*(j*3 + j)], w11_im = W[2*(j*3 + j) + 1];
        
        double v0 = 0.5 * (w00_re + w11_re);
        double v3 = 0.5 * (w00_im - w11_im);
        double v1 = 0.5 * (w01_im + w10_im);
        double v2 = 0.5 * (w01_re - w10_re);
        
        double k = sqrt(v0*v0 + v1*v1 + v2*v2 + v3*v3);
        if (k < 1e-10) continue;
        
        double inv_k = 1.0 / k;
        v0 *= inv_k; v1 *= inv_k; v2 *= inv_k; v3 *= inv_k;
        
        // Overrelax: R = V^dag * V^dag = V^{-2}
        double r0 = v0*v0 - v1*v1 - v2*v2 - v3*v3;
        double r1 = -2.0*v0*v1;
        double r2 = -2.0*v0*v2;
        double r3 = -2.0*v0*v3;
        
        su3_eq_id(R);
        R[2*(i*3 + i)]     = r0;  R[2*(i*3 + i) + 1] = r3;
        R[2*(i*3 + j)]     = r2;  R[2*(i*3 + j) + 1] = r1;
        R[2*(j*3 + i)]     = -r2; R[2*(j*3 + i) + 1] = r1;
        R[2*(j*3 + j)]     = r0;  R[2*(j*3 + j) + 1] = -r3;
        
        su3_eq_su3_ti_su3(T1, R, W);
        su3_eq_su3(W, T1);
        
        su3_eq_su3_ti_su3(T1, R, U);
        su3_eq_su3(U, T1);
    }
    
    su3_proj(U);
}

#endif
