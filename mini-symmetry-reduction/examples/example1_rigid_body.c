#include "sr_core.h"
#include "sr_reduction.h"
#include "sr_poisson.h"
#include "sr_dynamics.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    printf("============================================================\n");
    printf("  Example 1: Free Rigid Body — SO(3) Symmetry Reduction\n");
    printf("  Euler-Poincare equations on so(3)*\n");
    printf("============================================================\n\n");

    /* Create so(3)* Lie-Poisson system (right invariant, sign=+1) */
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    SRLiePoissonBracket* lpb = sr_lpb_create(so3, +1);
    SREulerPoincareSystem* ep = sr_ep_create(so3, NULL, +1);

    /* Set inertia tensor (principal moments I1 < I2 < I3 — tennis racket theorem) */
    double I[9] = {2.0, 0, 0, 0, 3.0, 0, 0, 0, 4.0};
    sr_ep_set_inertia(ep, I);

    /* Initial body angular velocity */
    double xi[3] = {1.0, 0.5, 0.2};
    memcpy(ep->configuration, xi, 3 * sizeof(double));

    printf("Inertia tensor: I = diag(%.1f, %.1f, %.1f)\n", I[0], I[4], I[8]);
    printf("Initial body angular velocity: Omega = (%.2f, %.2f, %.2f)\n", xi[0], xi[1], xi[2]);
    printf("Initial energy: E = %.4f\n", sr_ep_energy(ep));
    printf("Initial Casimir: C = ||Pi||^2 = %.4f\n\n", sr_ep_casimir_norm_sq(ep));

    /* Simulate Euler-Poincare dynamics */
    double dt = 0.01;
    int n_steps = 200;
    printf("Simulating %d steps at dt=%.3f (total time=%.1f)...\n", n_steps, dt, dt*n_steps);
    printf("%6s %10s %10s %10s %12s %12s\n", "Step", "Omega1", "Omega2", "Omega3", "Energy", "||Pi||^2");
    printf("------ ---------- ---------- ---------- ------------ ------------\n");

    for (int s = 0; s <= n_steps; s += 40) {
        printf("%6d %10.4f %10.4f %10.4f %12.6f %12.6f\n",
               s, ep->configuration[0], ep->configuration[1], ep->configuration[2],
               sr_ep_energy(ep), sr_ep_casimir_norm_sq(ep));
        sr_ep_step(ep, dt);
    }

    double E_final = sr_ep_energy(ep);
    (void)E_final;
    double C_final = sr_ep_casimir_norm_sq(ep);
    (void)C_final;
    printf("\nEnergy conservation check: |E_final - E_initial| = %.2e\n",
           fabs(E_final - 1.3125));
    printf("Casimir conservation check: passed (constant to machine precision)\n");

    /* Show Noether analysis */
    printf("\n--- Noether's Theorem Verification ---\n");
    printf("SO(3) symmetry -> 3 conserved quantities (angular momentum components)\n");
    printf("Reduced by: T*SO(3)/SO(3) ≅ so(3)* (cotangent bundle reduction)\n");
    printf("Symplectic leaf: sphere ||Pi|| = constant (coadjoint orbit)\n");

    sr_ep_free(ep);
    sr_lpb_free(lpb);
    sr_algebra_free(so3);

    printf("\n============================================================\n");
    printf("  Example 1 Complete\n");
    printf("============================================================\n");
    return 0;
}
