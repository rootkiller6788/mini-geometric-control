#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nonlinear_system.h"
#include "lie_derivative.h"
#include "lie_bracket.h"
#include "feedback_linearization.h"

#define BENCH_ITERATIONS 10000

double bench_vector_operations(void) {
    clock_t start = clock();
    Vector *v1 = vector_create(10);
    Vector *v2 = vector_create(10);
    Vector *v3 = vector_create(10);
    vector_set(v1, 1.0);
    vector_set(v2, 2.0);
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        vector_add(v1, v2, v3);
        vector_dot(v1, v2);
        vector_norm(v1);
    }
    vector_free(v1); vector_free(v2); vector_free(v3);
    clock_t end = clock();
    return (double)(end - start) / CLOCKS_PER_SEC;
}

double bench_matrix_operations(void) {
    clock_t start = clock();
    Matrix *A = matrix_create(4, 4);
    Matrix *B = matrix_create(4, 4);
    Matrix *C = matrix_create(4, 4);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            matrix_set(A, i, j, (double)(i+j));
            matrix_set(B, i, j, (double)(i-j));
        }
    for (int i = 0; i < BENCH_ITERATIONS; i++) {
        matrix_multiply(A, B, C);
        matrix_determinant(A);
    }
    matrix_free(A); matrix_free(B); matrix_free(C);
    clock_t end = clock();
    return (double)(end - clock()) / CLOCKS_PER_SEC;
}

int main(void) {
    printf("Feedback Linearization Benchmark
");
    printf("Iterations: %d

", BENCH_ITERATIONS);

    double t_vec = bench_vector_operations();
    printf("Vector ops (10-dim): %.4f sec
", t_vec);

    double t_mat = bench_matrix_operations();
    printf("Matrix ops (4x4):    %.4f sec
", t_mat);

    printf("
Benchmark complete.
");
    return 0;
}
