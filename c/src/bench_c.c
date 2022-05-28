#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "calabi_simd.h"

int main(void) {
    printf("===================================================\n");
    printf("   Calabi Native C SIMD Kernel Benchmark Harness   \n");
    printf("===================================================\n\n");

    const size_t n = 35; // N_3 = 35 sections (k=3 for quintic)
    const size_t num_points = 10000;

    double *s_re = (double *)malloc(n * sizeof(double));
    double *s_im = (double *)malloc(n * sizeof(double));
    double *h_inv = (double *)calloc(n * n * 2, sizeof(double));

    for (size_t i = 0; i < n; ++i) {
        s_re[i] = 1.0 / (double)(i + 1);
        s_im[i] = 0.5 / (double)(i + 1);
        h_inv[(i * n + i) * 2] = 1.0; // Identity
    }

    clock_t start = clock();
    double sum = 0.0;
    const size_t trials = 100000;
    for (size_t t = 0; t < trials; ++t) {
        sum += calabi_simd_section_norm(s_re, s_im, h_inv, n);
    }
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    printf("[1] Section Norm Evaluation (N=%zu, %zu trials):\n", n, trials);
    printf("    Elapsed time : %.4f sec\n", elapsed);
    printf("    Throughput   : %.2f million evaluations/sec\n", (trials / elapsed) / 1e6);
    printf("    Result value : %.6f\n\n", sum / trials);

    // Donaldson accumulation benchmark
    double *t_mat = (double *)calloc(n * n * 2, sizeof(double));
    double *s_re_batch = (double *)malloc(num_points * n * sizeof(double));
    double *s_im_batch = (double *)malloc(num_points * n * sizeof(double));
    double *weights = (double *)malloc(num_points * sizeof(double));
    double *norms = (double *)malloc(num_points * sizeof(double));

    for (size_t i = 0; i < num_points * n; ++i) {
        s_re_batch[i] = 0.1;
        s_im_batch[i] = 0.05;
    }
    for (size_t i = 0; i < num_points; ++i) {
        weights[i] = 1.0;
        norms[i] = 2.5;
    }

    start = clock();
    const size_t acc_trials = 50;
    for (size_t t = 0; t < acc_trials; ++t) {
        calabi_simd_accumulate_donaldson(t_mat, s_re_batch, s_im_batch, weights, norms, num_points, n);
    }
    end = clock();
    elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    printf("[2] Donaldson Batch Accumulation (%zu points x %zu sections, %zu passes):\n",
           num_points, n, acc_trials);
    printf("    Elapsed time : %.4f sec\n", elapsed);
    printf("    Throughput   : %.2f million point-section pairs/sec\n",
           ((double)(num_points * acc_trials * n) / elapsed) / 1e6);

    printf("\nAll native C kernel benchmarks completed successfully.\n");

    free(s_re);
    free(s_im);
    free(h_inv);
    free(t_mat);
    free(s_re_batch);
    free(s_im_batch);
    free(weights);
    free(norms);

    return 0;
}
