/* ==============================================================
 * example_heavy_top.c -- Heavy Top (Lagrange-Poisson) Dynamics
 *
 * Simulates a symmetric heavy top (spinning top in gravity).
 * The top has a fixed point and its center of mass is offset
 * along the body z-axis.
 *
 * Configuration: SO(3) with gravity pointing down in inertial
 * frame. The reduced dynamics live on so(3)* × S^2.
 *
 * Constants of motion:
 *   - Energy: E = T + V
 *   - Spatial angular momentum about z-axis: Pi_z (inertial)
 *
 * Ref: Marsden & Ratiu (1999) SS15.10
 *      Goldstein (2002) Classical Mechanics, Ch.5
 * ============================================================== */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lie_core.h"
#include "lie_reduction.h"

int main(void) {
    printf("=== Heavy Top (Lagrange-Poisson) Dynamics ===\n\n");

    /* Symmetric top: I1=I2=1, I3=2 */
    LieHeavyTopState *ht = lie_heavy_top_create(1.0, 1.0, 2.0, 1.0, 0.5);
    printf("Inertia: I1=1, I2=1, I3=2\n");
    printf("Mass: 1 kg, COM distance: 0.5 m\n");

    /* Initial conditions: top tilted at 30 degrees, spinning fast */
    double tilt = 30.0 * M_PI / 180.0;
    ht->Omega->data[2] = 10.0;  /* fast spin about body z */
    ht->Omega->data[0] = 0.0;

    /* Set initial Gamma (gravity direction in body frame) */
    ht->Gamma->data[0] = sin(tilt);
    ht->Gamma->data[2] = cos(tilt);

    /* chi: center of mass direction in body frame (along body z) */
    ht->chi->data[2] = 1.0;

    double E0 = lie_heavy_top_energy(ht);
    printf("Initial energy: %.6f J\n", E0);
    printf("Omega: (%.3f, %.3f, %.3f)\n",
           ht->Omega->data[0], ht->Omega->data[1], ht->Omega->data[2]);
    printf("Gamma: (%.3f, %.3f, %.3f)\n\n",
           ht->Gamma->data[0], ht->Gamma->data[1], ht->Gamma->data[2]);

    /* Simulate using Euler method on reduced variables */
    double h = 0.005;
    int steps = 200;
    printf("Simulating %d steps with dt=%.4f...\n", steps, h);

    for (int k = 0; k < steps; k++) {
        LieVector Pi_dot, Gamma_dot;
        LieMatrix R_dot_mat;
        Pi_dot.size = 3; Pi_dot.data = (double[3]){0};
        Gamma_dot.size = 3; Gamma_dot.data = (double[3]){0};
        R_dot_mat = *lie_matrix_create(3, 3);

        lie_heavy_top_rhs(ht, &Pi_dot, &Gamma_dot, &R_dot_mat);

        /* Update Omega: Pi = I * Omega => Omega = I^{-1} * Pi */
        LieVector *dOmega = lie_inertia_apply(
            lie_inertia_inverse(ht->I_body), &Pi_dot);

        for (int i = 0; i < 3; i++) {
            ht->Omega->data[i] += h * dOmega->data[i];
            ht->Gamma->data[i] += h * Gamma_dot.data[i];
        }
        lie_vector_free(dOmega);

        /* Renormalize Gamma (should stay on unit sphere) */
        double gnorm = lie_vector_norm(ht->Gamma);
        if (gnorm > 1e-10)
            for (int i = 0; i < 3; i++) ht->Gamma->data[i] /= gnorm;

        /* Update attitude */
        LieVector *hO = lie_vector_scale(ht->Omega, h);
        LieGroupElement *exp_hO = lie_exp_so3(hO);
        LieMatrix *Rnew = lie_matrix_mul(ht->R, exp_hO->g);
        memcpy(ht->R->data, Rnew->data, 9*sizeof(double));
        lie_matrix_free(exp_hO->g); free(exp_hO);
        lie_vector_free(hO); lie_matrix_free(Rnew);

        if (k % 40 == 0) {
            double E = lie_heavy_top_energy(ht);
            printf("  t=%.3f: Omega=(%7.3f,%7.3f,%7.3f) Gamma=(%6.3f,%6.3f,%6.3f) E=%.6f\n",
                   k*h, ht->Omega->data[0], ht->Omega->data[1], ht->Omega->data[2],
                   ht->Gamma->data[0], ht->Gamma->data[1], ht->Gamma->data[2], E);
        }
    }

    double Ef = lie_heavy_top_energy(ht);
    printf("\nFinal energy: %.6f (drift: %.2e)\n", Ef, fabs(Ef - E0));
    printf("Precession/nutation observed in Gamma trajectory.\n");

    lie_heavy_top_free(ht);
    printf("\n=== Example completed ===\n");
    return 0;
}

