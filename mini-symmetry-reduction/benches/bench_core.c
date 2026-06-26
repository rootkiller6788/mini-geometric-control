#include "sr_core.h"
#include "sr_dynamics.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
    printf("=== Symmetry Reduction — Performance Benchmarks ===\n\n");

    /* Benchmark: Lie bracket computation (so(3)) */
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    double X[3] = {1.0, 2.0, 3.0};
    double Y[3] = {0.5, 1.5, 2.5};
    double Z[3];

    clock_t start = clock();
    int N = 10000000;
    for (int i = 0; i < N; i++)
        sr_algebra_bracket(so3, X, Y, Z);
    clock_t end = clock();

    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Lie bracket [X,Y] (so(3)):\n");
    printf("  Iterations: %d\n", N);
    printf("  Total time: %.4f sec\n", elapsed);
    printf("  Per call:   %.2f ns\n", elapsed / N * 1e9);

    /* Benchmark: Exponential map (SO(3)) */
    double omega[3] = {0.1, 0.2, 0.3};
    double R[9];
    start = clock();
    for (int i = 0; i < N; i++)
        sr_exp_map_so3(omega, R);
    end = clock();
    elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("\nexp: so(3) -> SO(3) (Rodrigues):\n");
    printf("  Total time: %.4f sec\n", elapsed);
    printf("  Per call:   %.2f ns\n", elapsed / N * 1e9);

    /* Benchmark: Euler-Poincare step */
    SREulerPoincareSystem* ep = sr_ep_create(so3, NULL, +1);
    double I[9] = {2,0,0,0,3,0,0,0,4};
    sr_ep_set_inertia(ep, I);
    ep->configuration[0] = 0.1; ep->configuration[1] = 0.05; ep->configuration[2] = 0.02;

    start = clock();
    int M = 100000;
    for (int i = 0; i < M; i++)
        sr_ep_step(ep, 0.001);
    end = clock();
    elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("\nEuler-Poincare RK4 step:\n");
    printf("  Iterations: %d\n", M);
    printf("  Total time: %.4f sec\n", elapsed);
    printf("  Per call:   %.2f us\n", elapsed / M * 1e6);

    sr_ep_free(ep);
    sr_algebra_free(so3);

    printf("\n=== Bench Complete ===\n");
    return 0;
}
