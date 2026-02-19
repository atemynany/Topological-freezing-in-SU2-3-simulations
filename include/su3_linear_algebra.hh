// su3_linear_algebra.hh - SU(3) matrix operations
#ifndef __SU3_LINEAR_ALGEBRA_HH__
#define __SU3_LINEAR_ALGEBRA_HH__

#include <cmath>
#include <cstddef>

// SU(3) matrix: 3x3 complex = 18 doubles
// Layout: M[2*i*3 + 2*j] = Re(M_ij), M[2*i*3 + 2*j + 1] = Im(M_ij)
// Index: M[0..5] = row 0, M[6..11] = row 1, M[12..17] = row 2

inline void su3_eq_zero(double *A) {
    for (int i = 0; i < 18; i++) A[i] = 0.0;
}

inline void su3_eq_id(double *A) {
    su3_eq_zero(A);
    A[0] = A[8] = A[16] = 1.0;
}

inline void su3_eq_su3(double * __restrict__ A, const double * __restrict__ B) {
    for (int i = 0; i < 18; i++) A[i] = B[i];
}

inline void su3_pl_eq_su3(double * __restrict__ A, const double * __restrict__ B) {
    for (int i = 0; i < 18; i++) A[i] += B[i];
}

inline void su3_ti_eq_re(double *A, double c) {
    for (int i = 0; i < 18; i++) A[i] *= c;
}

inline void su3_eq_su3_dag(double * __restrict__ A, const double * __restrict__ B) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            A[2*(i*3 + j)]     =  B[2*(j*3 + i)];
            A[2*(i*3 + j) + 1] = -B[2*(j*3 + i) + 1];
        }
    }
}

inline void su3_dag_eq_su3(double *A) {
    double tmp[18];
    su3_eq_su3_dag(tmp, A);
    su3_eq_su3(A, tmp);
}

// C = A * B
inline void su3_eq_su3_ti_su3(double * __restrict__ C, 
                               const double * __restrict__ A, 
                               const double * __restrict__ B) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double re = 0.0, im = 0.0;
            for (int k = 0; k < 3; k++) {
                double a_re = A[2*(i*3 + k)];
                double a_im = A[2*(i*3 + k) + 1];
                double b_re = B[2*(k*3 + j)];
                double b_im = B[2*(k*3 + j) + 1];
                re += a_re * b_re - a_im * b_im;
                im += a_re * b_im + a_im * b_re;
            }
            C[2*(i*3 + j)] = re;
            C[2*(i*3 + j) + 1] = im;
        }
    }
}

// C = A^dag * B
inline void su3_eq_su3_dag_ti_su3(double * __restrict__ C,
                                   const double * __restrict__ A,
                                   const double * __restrict__ B) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double re = 0.0, im = 0.0;
            for (int k = 0; k < 3; k++) {
                double a_re =  A[2*(k*3 + i)];
                double a_im = -A[2*(k*3 + i) + 1];
                double b_re = B[2*(k*3 + j)];
                double b_im = B[2*(k*3 + j) + 1];
                re += a_re * b_re - a_im * b_im;
                im += a_re * b_im + a_im * b_re;
            }
            C[2*(i*3 + j)] = re;
            C[2*(i*3 + j) + 1] = im;
        }
    }
}

// C = A * B^dag
inline void su3_eq_su3_ti_su3_dag(double * __restrict__ C,
                                   const double * __restrict__ A,
                                   const double * __restrict__ B) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double re = 0.0, im = 0.0;
            for (int k = 0; k < 3; k++) {
                double a_re = A[2*(i*3 + k)];
                double a_im = A[2*(i*3 + k) + 1];
                double b_re =  B[2*(j*3 + k)];
                double b_im = -B[2*(j*3 + k) + 1];
                re += a_re * b_re - a_im * b_im;
                im += a_re * b_im + a_im * b_re;
            }
            C[2*(i*3 + j)] = re;
            C[2*(i*3 + j) + 1] = im;
        }
    }
}

// C = A^dag * B^dag
inline void su3_eq_su3_dag_ti_su3_dag(double * __restrict__ C,
                                       const double * __restrict__ A,
                                       const double * __restrict__ B) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double re = 0.0, im = 0.0;
            for (int k = 0; k < 3; k++) {
                double a_re =  A[2*(k*3 + i)];
                double a_im = -A[2*(k*3 + i) + 1];
                double b_re =  B[2*(j*3 + k)];
                double b_im = -B[2*(j*3 + k) + 1];
                re += a_re * b_re - a_im * b_im;
                im += a_re * b_im + a_im * b_re;
            }
            C[2*(i*3 + j)] = re;
            C[2*(i*3 + j) + 1] = im;
        }
    }
}

inline double su3_re_tr(const double *A) {
    return A[0] + A[8] + A[16];
}

inline void su3_co_tr(double *re, double *im, const double *A) {
    *re = A[0] + A[8] + A[16];
    *im = A[1] + A[9] + A[17];
}

// Extract SU(2) subgroup (rows/cols i,j) from SU(3) staple
inline void su3_get_su2_subgroup(double *su2, const double *su3, int i, int j) {
    // SU(2) stored as: [a0_re, a0_im, a1_re, a1_im, a2_re, a2_im, a3_re, a3_im]
    // where U = a0*1 + i*(a1*sigma1 + a2*sigma2 + a3*sigma3)
    // = [[a0+i*a3, a2+i*a1], [-a2+i*a1, a0-i*a3]]
    double m00_re = su3[2*(i*3 + i)];
    double m00_im = su3[2*(i*3 + i) + 1];
    double m01_re = su3[2*(i*3 + j)];
    double m01_im = su3[2*(i*3 + j) + 1];
    double m10_re = su3[2*(j*3 + i)];
    double m10_im = su3[2*(j*3 + i) + 1];
    double m11_re = su3[2*(j*3 + j)];
    double m11_im = su3[2*(j*3 + j) + 1];
    
    // a0 = (m00 + m11)/2, a3 = (m00 - m11)/(2i)
    su2[0] = 0.5 * (m00_re + m11_re);  // a0_re
    su2[1] = 0.5 * (m00_im + m11_im);  // a0_im
    su2[6] = 0.5 * (m00_im - m11_im);  // a3_re  (from -i*(m00-m11)/2)
    su2[7] = 0.5 * (m11_re - m00_re);  // a3_im
    
    // a1 = (m01 + m10)/(2i), a2 = (m01 - m10)/2
    su2[2] = 0.5 * (m01_im + m10_im);  // a1_re (from -i*(m01+m10)/2)
    su2[3] = 0.5 * (-m01_re - m10_re); // a1_im
    su2[4] = 0.5 * (m01_re - m10_re);  // a2_re
    su2[5] = 0.5 * (m01_im - m10_im);  // a2_im
}

// Embed SU(2) into SU(3) at subgroup (i,j)
inline void su3_embed_su2(double *su3, const double *su2, int i, int j) {
    su3_eq_id(su3);
    // u = [[a0+i*a3, a2+i*a1], [-a2+i*a1, a0-i*a3]]
    double a0_re = su2[0], a0_im = su2[1];
    double a1_re = su2[2], a1_im = su2[3];
    double a2_re = su2[4], a2_im = su2[5];
    double a3_re = su2[6], a3_im = su2[7];
    
    su3[2*(i*3 + i)]     = a0_re + a3_im;
    su3[2*(i*3 + i) + 1] = a0_im + a3_re;
    su3[2*(i*3 + j)]     = a2_re - a1_im;
    su3[2*(i*3 + j) + 1] = a2_im + a1_re;
    su3[2*(j*3 + i)]     = -a2_re - a1_im;
    su3[2*(j*3 + i) + 1] = -a2_im + a1_re;
    su3[2*(j*3 + j)]     = a0_re - a3_im;
    su3[2*(j*3 + j) + 1] = a0_im - a3_re;
}

// Project to SU(3) via Gram-Schmidt
inline void su3_proj(double *A) {
    // Normalize row 0
    double norm = 0.0;
    for (int j = 0; j < 3; j++) {
        norm += A[2*j]*A[2*j] + A[2*j+1]*A[2*j+1];
    }
    norm = 1.0 / sqrt(norm);
    for (int j = 0; j < 3; j++) {
        A[2*j] *= norm;
        A[2*j+1] *= norm;
    }
    
    // Orthogonalize row 1 to row 0
    double dot_re = 0.0, dot_im = 0.0;
    for (int j = 0; j < 3; j++) {
        dot_re += A[6 + 2*j]*A[2*j] + A[6 + 2*j + 1]*A[2*j + 1];
        dot_im += A[6 + 2*j + 1]*A[2*j] - A[6 + 2*j]*A[2*j + 1];
    }
    for (int j = 0; j < 3; j++) {
        A[6 + 2*j]     -= dot_re*A[2*j] - dot_im*A[2*j + 1];
        A[6 + 2*j + 1] -= dot_re*A[2*j + 1] + dot_im*A[2*j];
    }
    
    // Normalize row 1
    norm = 0.0;
    for (int j = 0; j < 3; j++) {
        norm += A[6 + 2*j]*A[6 + 2*j] + A[6 + 2*j + 1]*A[6 + 2*j + 1];
    }
    norm = 1.0 / sqrt(norm);
    for (int j = 0; j < 3; j++) {
        A[6 + 2*j] *= norm;
        A[6 + 2*j + 1] *= norm;
    }
    
    // Row 2 = conj(row 0 x row 1) for SU(3)
    // (u x v)_i = eps_ijk u_j v_k  then conjugate
    // row2[0] = conj(row0[1]*row1[2] - row0[2]*row1[1])
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        double uj_re = A[2*j], uj_im = A[2*j + 1];
        double vk_re = A[6 + 2*k], vk_im = A[6 + 2*k + 1];
        double uk_re = A[2*k], uk_im = A[2*k + 1];
        double vj_re = A[6 + 2*j], vj_im = A[6 + 2*j + 1];
        
        double cross_re = (uj_re*vk_re - uj_im*vk_im) - (uk_re*vj_re - uk_im*vj_im);
        double cross_im = (uj_re*vk_im + uj_im*vk_re) - (uk_re*vj_im + uk_im*vj_re);
        
        A[12 + 2*i] = cross_re;
        A[12 + 2*i + 1] = -cross_im;  // conjugate
    }
}

// Generate random SU(3) matrix
inline void su3_random(double *A, double (*drand)()) {
    for (int i = 0; i < 18; i++) {
        A[i] = 2.0 * drand() - 1.0;
    }
    su3_proj(A);
}

#endif
