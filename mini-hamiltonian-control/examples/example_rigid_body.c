/* example_rigid_body.c -- Free rigid body dynamics on so(3)* */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hamiltonian_control.h"
#include "poisson_geometry.h"

int main(void) {
    printf("=== Free Rigid Body: Euler Equations on so(3)* ===\n");

    rigid_body_params_t rb = {1.0, 2.0, 3.0}; /* I1 < I2 < I3 */
    double Pi[3] = {1.0, 0.1, 0.1};
    double dt = 0.05;
    int steps = 200;

    /* Casimir: C1 = Pi1^2 + Pi2^2 + Pi3^2, C2 = Pi1^2/I1 + Pi2^2/I2 + Pi3^2/I3 (energy) */
    double C1_0 = Pi[0]*Pi[0] + Pi[1]*Pi[1] + Pi[2]*Pi[2];
    double E_0 = rigid_body_hamiltonian(Pi, 3, &rb);

    printf("I = (%.1f, %.1f, %.1f)\n", rb.I1, rb.I2, rb.I3);
    printf("Initial: Pi=(%.3f,%.3f,%.3f), |Pi|^2=%.4f, E=%.4f\n", Pi[0],Pi[1],Pi[2],C1_0,E_0);

    for (int k = 0; k < steps; k++) {
        double dPi[3];
        rigid_body_euler_equations(&rb, Pi, dPi);
        for (int i = 0; i < 3; i++) Pi[i] += dt * dPi[i];

        if (k % 40 == 0) {
            double C1 = Pi[0]*Pi[0]+Pi[1]*Pi[1]+Pi[2]*Pi[2];
            double E = rigid_body_hamiltonian(Pi, 3, &rb);
            printf("  t=%.1f: Pi=(%.3f,%.3f,%.3f), |Pi|^2=%.4f, E=%.4f\n",
                   k*dt, Pi[0],Pi[1],Pi[2],C1,E);
        }
    }

    double C1_f = Pi[0]*Pi[0] + Pi[1]*Pi[1] + Pi[2]*Pi[2];
    double E_f = rigid_body_hamiltonian(Pi, 3, &rb);
    printf("Final: |Pi|^2=%.4f (should equal %.4f), E=%.4f (should equal %.4f)\n",
           C1_f, C1_0, E_f, E_0);
    printf("\nRotation about intermediate axis (I2) is unstable (tennis racket theorem).\n");
    return 0;
}