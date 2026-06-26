#include "sr_core.h"
#include "sr_reduction.h"
#include "sr_poisson.h"
#include "sr_control.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    printf("============================================================\n");
    printf("  Example 3: Spacecraft Attitude Control via Symmetry Reduction\n");
    printf("  Controlled Euler-Poincare equations on so(3)*\n");
    printf("============================================================\n\n");

    /* Spacecraft inertia (simulated small satellite) */
    double I_diag[3] = {12.0, 15.0, 8.0};
    double Omega[3] = {0.1, 0.0, 0.05};
    double body_torque[3] = {0.0, 0.0, 0.0};

    printf("Spacecraft inertia: I = diag(%.1f, %.1f, %.1f) kg·m^2\n",
           I_diag[0], I_diag[1], I_diag[2]);
    printf("Initial body angular velocity: Omega = (%.3f, %.3f, %.3f) rad/s\n",
           Omega[0], Omega[1], Omega[2]);

    /* Create Euler-Poincare system */
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    SREulerPoincareSystem* ep = sr_ep_create(so3, NULL, +1);
    double I_mat[9] = {I_diag[0],0,0, 0,I_diag[1],0, 0,0,I_diag[2]};
    sr_ep_set_inertia(ep, I_mat);
    memcpy(ep->configuration, Omega, 3 * sizeof(double));

    printf("\n--- Open-loop free drift (no control) ---\n");
    printf("%6s %10s %10s %10s %12s\n", "Step", "Omega1", "Omega2", "Omega3", "Energy");
    printf("------ ---------- ---------- ---------- ------------\n");

    double dt = 0.05;
    for (int s = 0; s <= 10; s++) {
        printf("%6d %10.4f %10.4f %10.4f %12.6f\n",
               s, ep->configuration[0], ep->configuration[1], ep->configuration[2],
               sr_ep_energy(ep));
        sr_ep_step(ep, dt);
    }

    /* Now apply control torque to detumble */
    printf("\n--- Controlled detumbling maneuver ---\n");
    printf("Applying body-fixed control torques to reduce angular velocity...\n");
    printf("%6s %10s %10s %10s %10s\n", "Step", "Omega1", "Omega2", "Omega3", "|Omega|");
    printf("------ ---------- ---------- ---------- ----------\n");

    /* Reset initial conditions */
    memcpy(ep->configuration, Omega, 3 * sizeof(double));

    for (int s = 0; s <= 10; s++) {
        /* Simple proportional-derivative detumbling control
         * tau = -Kp * Omega (body-fixed damping) */
        double Kp = 2.0;
        double control[3] = {-Kp * ep->configuration[0],
                             -Kp * ep->configuration[1],
                             -Kp * ep->configuration[2]};

        /* Controlled EP: I dOmega/dt = I Omega × Omega + control */
        double dOmega[3];
        sr_ctrl_ep_rhs(ep, ep->configuration, control, dOmega);

        double norm = sqrt(ep->configuration[0]*ep->configuration[0] +
                          ep->configuration[1]*ep->configuration[1] +
                          ep->configuration[2]*ep->configuration[2]);
        printf("%6d %10.4f %10.4f %10.4f %10.4f\n",
               s, ep->configuration[0], ep->configuration[1], ep->configuration[2], norm);

        /* Simple Euler step with control */
        for (int i = 0; i < 3; i++)
            ep->configuration[i] += dt * dOmega[i];
    }

    printf("\n--- Geometric Control Analysis ---\n");
    printf("SO(3) symmetry reduction: T*SO(3) -> so(3)*\n");
    printf("Reduced state space dimension: 3 (body angular velocity)\n");
    printf("Original state space: T*SO(3) dimension = 6+6 = 12 -> reduced to 3!\n");
    printf("Controllability: SO(3) is compact, body-fixed torques span so(3).\n");
    printf("Control design benefit: no attitude parametrization singularities.\n");

    printf("\n--- Application Notes (NASA/ESA) ---\n");
    printf("This reduction is used in spacecraft attitude control systems\n");
    printf("(e.g., reaction wheels, magnetorquers). The reduced Lie-Poisson\n");
    printf("structure on so(3)* preserves the geometric properties of the\n");
    printf("original Hamiltonian system.\n");

    sr_ep_free(ep);
    sr_algebra_free(so3);

    printf("\n============================================================\n");
    printf("  Example 3 Complete\n");
    printf("============================================================\n");
    return 0;
}
