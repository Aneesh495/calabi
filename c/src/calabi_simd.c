#include "calabi_simd.h"
#include <math.h>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#define CALABI_USE_NEON 1
#elif defined(__AVX2__)
#include <immintrin.h>
#define CALABI_USE_AVX2 1
#endif

double calabi_simd_section_norm(
    const double *s_re,
    const double *s_im,
    const double *h_inv,
    size_t n
) {
    double total = 0.0;

    for (size_t a = 0; a < n; ++a) {
        double sa_r = s_re[a];
        double sa_i = s_im[a];
        if (sa_r == 0.0 && sa_i == 0.0) continue;

        const double *h_row = &h_inv[(a * n) * 2];

        for (size_t b = 0; b < n; ++b) {
            double sb_r = s_re[b];
            double sb_i = -s_im[b]; // conj

            double hr = h_row[b * 2];
            double hi = h_row[b * 2 + 1];

            double prod_r = sa_r * sb_r - sa_i * sb_i;
            double prod_i = sa_r * sb_i + sa_i * sb_r;

            total += (prod_r * hr - prod_i * hi);
        }
    }

    return total;
}

void calabi_simd_accumulate_donaldson(
    double *t_mat,
    const double *s_re_batch,
    const double *s_im_batch,
    const double *weights,
    const double *norms,
    size_t num_points,
    size_t n
) {
    for (size_t pt = 0; pt < num_points; ++pt) {
        double w = weights[pt];
        double s_norm = norms[pt];
        if (s_norm <= 1e-15) continue;
        double factor = w / s_norm;

        const double *s_r = &s_re_batch[pt * n];
        const double *s_i = &s_im_batch[pt * n];

        for (size_t a = 0; a < n; ++a) {
            double sa_r = s_r[a];
            double sa_i = s_i[a];
            if (sa_r == 0.0 && sa_i == 0.0) continue;

            double *t_row = &t_mat[(a * n) * 2];

            for (size_t b = 0; b < n; ++b) {
                double sb_r = s_r[b];
                double sb_i = -s_i[b]; // conj

                double term_r = (sa_r * sb_r - sa_i * sb_i) * factor;
                double term_i = (sa_r * sb_i + sa_i * sb_r) * factor;

                t_row[b * 2] += term_r;
                t_row[b * 2 + 1] += term_i;
            }
        }
    }
}

void calabi_simd_matrix_mult(
    double *c,
    const double *a,
    const double *b,
    size_t n
) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t k = 0; k < n; ++k) {
            size_t a_idx = (i * n + k) * 2;
            double ar = a[a_idx];
            double ai = a[a_idx + 1];
            if (ar == 0.0 && ai == 0.0) continue;

            for (size_t j = 0; j < n; ++j) {
                size_t b_idx = (k * n + j) * 2;
                double br = b[b_idx];
                double bi = b[b_idx + 1];

                size_t c_idx = (i * n + j) * 2;
                c[c_idx] += (ar * br - ai * bi);
                c[c_idx + 1] += (ar * bi + ai * br);
            }
        }
    }
}
