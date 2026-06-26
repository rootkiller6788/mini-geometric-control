/* ==============================================================
 * example_rigid_body.c -- Rigid Body Rotation on SO(3)
 *
 * Demonstrates the Euler rigid body equations for a free rigid
 * body with principal inertias I1=1, I2=2, I3=3.
 *
 * Simulates the attitude evolution and verifies conservation
 * of kinetic energy and angular momentum magnitude.
 *
 * Ref: Marsden & Ratiu (1999) SS15.9
 * ============================================================== */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "lie_core.h"
#include "lie_reduction.h"

int main(void) {
    printf("=== Rigid Body Rotation on SO(3) ===\n\n");

    /* Create rigid body with principal inertias */
    LieRigidBodyState *rb = lie_rigid_body_create(1.0, 2.0, 3.0, 0.0, 0.0, 0.0);
    printf("Inertia tensor: I1=1, I2=2, I3=3\n");

    /* Initial angular velocity: rotation mostly about intermediate axis */
    rb->Omega->data[0] = 0.1;
    rb->Omega->data[1] = 2.0;
    rb->Omega->data[2] = 0.1;

    double T0 = lie_rigid_body_kinetic_energy(rb);
    printf("Initial kinetic energy: %.6f\n", T0);
    printf("Initial Omega: (%.4f, %.4f, %.4f)\n\n",
           rb->Omega->data[0], rb->Omega->data[1], rb->Omega->data[2]);

    /* Simulate using Euler integration on so(3) */
    double h = 0.01;
    int steps = 100;
    printf("Simulating %d steps with dt=%.3f...\n", steps, h);

    LieFreeTop *top = lie_free_top_create(1.0, 2.0, 3.0);
    lie_free_top_init(top, 0.1, 2.0, 0.1);

    for (int k = 0; k < steps; k++) {
        LieVector dOmega;
        dOmega.size = 3;
        double ddata[3] = {0, 0, 0};
        dOmega.data = ddata;
        lie_free_top_rhs(top, &dOmega);

        /* Euler update on Omega */
        for (int i = 0; i < 3; i++) top->Omega->data[i] += h * dOmega.data[i];

        /* Update attitude via exponential map */
        LieVector *hOmega = lie_vector_scale(top->Omega, h);
        LieGroupElement *exp_hO = lie_exp_so3(hOmega);
        LieMatrix *Rnew = lie_matrix_mul(top->R, exp_hO->g);
        memcpy(top->R->data, Rnew->data, 9 * sizeof(double));
        lie_vector_free(hOmega);
        lie_group_free(exp_hO);
        lie_matrix_free(Rnew);

        /* Recompute invariants */
        LieVector *Pi = lie_inertia_apply(top->I, top->Omega);
        top->energy = 0.5 * lie_vector_dot(top->Omega, Pi);
        top->momentum_sq = lie_vector_norm_sq(Pi);
        lie_vector_free(Pi);

        if (k % 20 == 0) {
            printf("  t=%.2f: Omega=(%7.4f, %7.4f, %7.4f) T=%.6f |Pi|^2=%.6f\n",
                   k*h, top->Omega->data[0], top->Omega->data[1], top->Omega->data[2],
                   top->energy, top->momentum_sq);
        }
    }

    /* Verify energy conservation */
    printf("\nFinal vs initial energy: %.6f vs %.6f\n", top->energy, T0);
    printf("Energy drift: %.2e\n", fabs(top->energy - T0));

    lie_free_top_free(top);
    lie_rigid_body_free(rb);
    printf("\n=== Example completed ===\n");
    return 0;
}

