#ifndef CALABI_SIMD_H
#define CALABI_SIMD_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * calabi_simd_section_norm:
 * Evaluates the Hermitian quadratic form S = sum_{a, b} hInv[a, b] * s[a] * conj(s[b])
 * where s is a complex vector of length n, and hInv is an n x n complex Hermitian matrix
 * stored in interleaved real/imaginary format: [re_00, im_00, re_01, im_01, ...].
 */
double calabi_simd_section_norm(
    const double *s_re,
    const double *s_im,
    const double *h_inv,
    size_t n
);

/**
 * calabi_simd_accumulate_donaldson:
 * Batch accumulates outer products into the Hermitian matrix T:
 *   T[a, b] += (w / S) * s[a] * conj(s[b])
 * for a batch of points.
 */
void calabi_simd_accumulate_donaldson(
    double *t_mat,
    const double *s_re_batch,
    const double *s_im_batch,
    const double *weights,
    const double *norms,
    size_t num_points,
    size_t n
);

/**
 * calabi_simd_matrix_mult:
 * Complex matrix-matrix multiplication C = A * B for n x n complex matrices.
 */
void calabi_simd_matrix_mult(
    double *c,
    const double *a,
    const double *b,
    size_t n
);

#ifdef __cplusplus
}
#endif

#endif /* CALABI_SIMD_H */
