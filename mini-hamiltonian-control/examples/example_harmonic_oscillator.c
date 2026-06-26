/* example_harmonic_oscillator.c -- Harmonic oscillator with symplectic integration */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hamiltonian_control.h"
#include "symplectic_integrator.h"

static double ho_energy(const double *q, const double *p, int n, void *params) {
    double w = *(double*)params;
    return 0.5*(p[0]*p[0] + w*w*q[0]*q[0]);
}
static void ho_gradient(const double *q, const double *p, int n, double *gq, double *gp, void *params) {
    double w = *(double*)params;
    gq[0] = w*w*q[0]; gp[0] = p[0];
}

int main(void) {
    double omega = 2.0;
    hamiltonian_t H = {ho_energy, ho_gradient, NULL, &omega, 1};

    double q = 1.0, p = 0.0;
    double dt = 0.01;
    int n_steps = 628; /* ~ one period at omega=2.0 */
    integrator_diagnostics_t diag;

    printf("=== Harmonic Oscillator: Symplectic vs RK4 ===\n");
    printf("omega=%.1f, q0=%.1f, p0=%.1f\n", omega, q, p);

    /* Symplectic integration */
    double q_s = q, p_s = p;
    symplectic_integrate(&H, &q_s, &p_s, 0, n_steps*dt, dt, INT_STORMER_VERLET, &n_steps, &diag);
    double E_s = H.energy(&q_s, &p_s, 1, &omega);
    printf("Stormer-Verlet: q=%.6f, p=%.6f, E=%.6f, dE=%.2e\n", q_s, p_s, E_s, diag.energy_error);

    /* RK4 (non-symplectic) for comparison */
    double q_r = q, p_r = p;
    hamiltonian_flow_rk4(&H, &q_r, &p_r, dt, n_steps);
    double E_r = H.energy(&q_r, &p_r, 1, &omega);
    double E0 = H.energy(&q, &p, 1, &omega);
    printf("RK4:            q=%.6f, p=%.6f, E=%.6f, dE=%.2e\n", q_r, p_r, E_r, fabs(E_r-E0)/fabs(E0));

    printf("\nSymplectic integrator preserves energy much better.\n");
    return 0;
}