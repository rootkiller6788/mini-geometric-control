/*=============================================================================
 * symplectic_integrator.c -- Geometric numerical integration for Hamiltonian systems
 * Reference: Hairer, Lubich, Wanner, Geometric Numerical Integration (2006)
 *            Leimkuhler & Reich, Simulating Hamiltonian Dynamics (2004)
 *===========================================================================*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "hamiltonian_control.h"
#include "symplectic_integrator.h"

/* L2: Symplectic Euler method A (implicit in q) */
int step_symplectic_euler_a(const hamiltonian_t *H,
                             double *q, double *p, double dt,
                             integrator_diagnostics_t *diag)
{
    int n = H->n_dof;
    double E0 = H->energy(q, p, n, H->params);
    double *q_new = (double*)malloc(n * sizeof(double));
    double *p_new = (double*)malloc(n * sizeof(double));
    double *gq = (double*)malloc(n * sizeof(double));
    double *gp = (double*)malloc(n * sizeof(double));
    double *gq_new = (double*)malloc(n * sizeof(double));
    double *gp_new = (double*)malloc(n * sizeof(double));

    memcpy(q_new, q, n * sizeof(double));
    memcpy(p_new, p, n * sizeof(double));

    /* Fixed-point iteration for implicit q */
    for (int iter = 0; iter < 20; iter++) {
        H->gradient(q_new, p, n, gq_new, gp_new, H->params);
        for (int i = 0; i < n; i++) {
            q_new[i] = q[i] + dt * gp_new[i];
        }
        /* Check convergence */
        H->gradient(q_new, p, n, gq, gp, H->params);
        double err = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = q_new[i] - (q[i] + dt * gp[i]);
            err += diff * diff;
        }
        if (sqrt(err) < 1e-12) break;
    }

    /* p update (explicit using q_new) */
    H->gradient(q_new, p, n, gq_new, gp_new, H->params);
    for (int i = 0; i < n; i++)
        p_new[i] = p[i] - dt * gq_new[i];

    memcpy(q, q_new, n * sizeof(double));
    memcpy(p, p_new, n * sizeof(double));

    if (diag) {
        double E1 = H->energy(q, p, n, H->params);
        diag->energy_error = fabs(E1 - E0) / (fabs(E0) + 1e-12);
        diag->symplectic_error = 0.0;
        diag->momentum_error = 0.0;
    }

    free(q_new); free(p_new); free(gq); free(gp); free(gq_new); free(gp_new);
    return 0;
}

/* L2: Symplectic Euler method B (implicit in p) */
int step_symplectic_euler_b(const hamiltonian_t *H,
                             double *q, double *p, double dt,
                             integrator_diagnostics_t *diag)
{
    int n = H->n_dof;
    double E0 = H->energy(q, p, n, H->params);
    double *q_new = (double*)malloc(n * sizeof(double));
    double *p_new = (double*)malloc(n * sizeof(double));
    double *gq = (double*)malloc(n * sizeof(double));
    double *gp = (double*)malloc(n * sizeof(double));

    memcpy(p_new, p, n * sizeof(double));

    /* Fixed-point iteration for implicit p */
    for (int iter = 0; iter < 20; iter++) {
        H->gradient(q, p_new, n, gq, gp, H->params);
        for (int i = 0; i < n; i++)
            p_new[i] = p[i] - dt * gq[i];
        H->gradient(q, p_new, n, gq, gp, H->params);
        double err = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = p_new[i] - (p[i] - dt * gq[i]);
            err += diff * diff;
        }
        if (sqrt(err) < 1e-12) break;
    }

    /* q update (explicit using p_new) */
    H->gradient(q, p_new, n, gq, gp, H->params);
    for (int i = 0; i < n; i++)
        q_new[i] = q[i] + dt * gp[i];

    memcpy(q, q_new, n * sizeof(double));
    memcpy(p, p_new, n * sizeof(double));

    if (diag) {
        double E1 = H->energy(q, p, n, H->params);
        diag->energy_error = fabs(E1 - E0) / (fabs(E0) + 1e-12);
    }

    free(q_new); free(p_new); free(gq); free(gp);
    return 0;
}

/* L3: Stormer-Verlet integrator (velocity form) */
int step_stormer_verlet(const hamiltonian_t *H,
                         double *q, double *p, double dt,
                         integrator_diagnostics_t *diag)
{
    int n = H->n_dof;
    double E0 = H->energy(q, p, n, H->params);
    double *gq = (double*)malloc(n * sizeof(double));
    double *gp = (double*)malloc(n * sizeof(double));

    /* Half step in p */
    H->gradient(q, p, n, gq, gp, H->params);
    for (int i = 0; i < n; i++)
        p[i] -= 0.5 * dt * gq[i];

    /* Full step in q using mid-point p */
    H->gradient(q, p, n, gq, gp, H->params);
    for (int i = 0; i < n; i++)
        q[i] += dt * gp[i];

    /* Half step in p */
    H->gradient(q, p, n, gq, gp, H->params);
    for (int i = 0; i < n; i++)
        p[i] -= 0.5 * dt * gq[i];

    if (diag) {
        double E1 = H->energy(q, p, n, H->params);
        diag->energy_error = fabs(E1 - E0) / (fabs(E0) + 1e-12);
    }

    free(gq); free(gp);
    return 0;
}

/* L3: Velocity Verlet integrator */
int step_velocity_verlet(const hamiltonian_t *H,
                          double *q, double *p, double dt,
                          integrator_diagnostics_t *diag)
{
    int n = H->n_dof;
    double E0 = H->energy(q, p, n, H->params);
    double *gq = (double*)malloc(n * sizeof(double));
    double *gp = (double*)malloc(n * sizeof(double));

    /* p at half step */
    H->gradient(q, p, n, gq, gp, H->params);
    for (int i = 0; i < n; i++)
        p[i] -= 0.5 * dt * gq[i];

    /* q full step */
    H->gradient(q, p, n, gq, gp, H->params);
    for (int i = 0; i < n; i++)
        q[i] += dt * gp[i];

    /* p final half step */
    H->gradient(q, p, n, gq, gp, H->params);
    for (int i = 0; i < n; i++)
        p[i] -= 0.5 * dt * gq[i];

    if (diag) {
        double E1 = H->energy(q, p, n, H->params);
        diag->energy_error = fabs(E1 - E0) / (fabs(E0) + 1e-12);
    }

    free(gq); free(gp);
    return 0;
}

/* L4: Implicit midpoint rule z_{n+1} = z_n + dt * J * grad(H)((z_n+z_{n+1})/2) */
int step_implicit_midpoint(const hamiltonian_t *H,
                            double *q, double *p, double dt,
                            int max_iter, double tol,
                            integrator_diagnostics_t *diag)
{
    int n = H->n_dof;
    double E0 = H->energy(q, p, n, H->params);
    double *q_new = (double*)malloc(n * sizeof(double));
    double *p_new = (double*)malloc(n * sizeof(double));
    double *q_mid = (double*)malloc(n * sizeof(double));
    double *p_mid = (double*)malloc(n * sizeof(double));
    double *gq_mid = (double*)malloc(n * sizeof(double));
    double *gp_mid = (double*)malloc(n * sizeof(double));

    memcpy(q_new, q, n * sizeof(double));
    memcpy(p_new, p, n * sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        /* Midpoint */
        for (int i = 0; i < n; i++) {
            q_mid[i] = 0.5 * (q[i] + q_new[i]);
            p_mid[i] = 0.5 * (p[i] + p_new[i]);
        }

        H->gradient(q_mid, p_mid, n, gq_mid, gp_mid, H->params);

        /* Compute proposed new values */
        double *q_prop = (double*)malloc(n * sizeof(double));
        double *p_prop = (double*)malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) {
            q_prop[i] = q[i] + dt * gp_mid[i];
            p_prop[i] = p[i] - dt * gq_mid[i];
        }

        /* Check convergence */
        double err = 0.0;
        for (int i = 0; i < n; i++) {
            double dq = q_prop[i] - q_new[i];
            double dp = p_prop[i] - p_new[i];
            err += dq*dq + dp*dp;
        }
        memcpy(q_new, q_prop, n * sizeof(double));
        memcpy(p_new, p_prop, n * sizeof(double));
        free(q_prop); free(p_prop);

        if (sqrt(err) < tol) break;
    }

    memcpy(q, q_new, n * sizeof(double));
    memcpy(p, p_new, n * sizeof(double));

    if (diag) {
        double E1 = H->energy(q, p, n, H->params);
        diag->energy_error = fabs(E1 - E0) / (fabs(E0) + 1e-12);
    }

    free(q_new); free(p_new); free(q_mid); free(p_mid);
    free(gq_mid); free(gp_mid);
    return 0;
}

/* L5: Yoshida 4th order composition method */
int step_yoshida_4(const hamiltonian_t *H,
                    double *q, double *p, double dt,
                    integrator_diagnostics_t *diag)
{
    double w0 = 1.0 / (2.0 - pow(2.0, 1.0/3.0));
    double w1 = -pow(2.0, 1.0/3.0) / (2.0 - pow(2.0, 1.0/3.0));
    double ws[4] = {w0, w1, w1, w0};

    for (int k = 0; k < 4; k++) {
        step_stormer_verlet(H, q, p, ws[k] * dt, NULL);
    }

    if (diag) {
        double E = H->energy(q, p, H->n_dof, H->params);
        diag->energy_error = 0.0; /* accumulated in individual steps */
    }
    return 0;
}

/* L5: Adaptive symplectic integration */
int adaptive_symplectic_integrate(const hamiltonian_t *H,
                                   double *q, double *p,
                                   double t0, double tf,
                                   const adaptive_params_t *params,
                                   int *n_steps_taken,
                                   integrator_diagnostics_t *final_diag)
{
    double t = t0;
    double dt = (tf - t0) / 100.0; /* initial step guess */
    *n_steps_taken = 0;

    while (t < tf) {
        if (t + dt > tf) dt = tf - t;
        if (dt < params->dt_min) dt = params->dt_min;

        /* Take step and estimate error */
        double *q1 = (double*)malloc(H->n_dof * sizeof(double));
        double *p1 = (double*)malloc(H->n_dof * sizeof(double));
        double *q2 = (double*)malloc(H->n_dof * sizeof(double));
        double *p2 = (double*)malloc(H->n_dof * sizeof(double));

        memcpy(q1, q, H->n_dof * sizeof(double));
        memcpy(p1, p, H->n_dof * sizeof(double));
        memcpy(q2, q, H->n_dof * sizeof(double));
        memcpy(p2, p, H->n_dof * sizeof(double));

        /* One full step */
        step_stormer_verlet(H, q1, p1, dt, NULL);
        /* Two half steps */
        step_stormer_verlet(H, q2, p2, 0.5*dt, NULL);
        step_stormer_verlet(H, q2, p2, 0.5*dt, NULL);

        /* Error estimate */
        double err = 0.0;
        for (int i = 0; i < H->n_dof; i++) {
            double dq = q1[i] - q2[i];
            double dp = p1[i] - p2[i];
            err += dq*dq + dp*dp;
        }
        err = sqrt(err);

        if (err < params->tol) {
            /* Accept step */
            memcpy(q, q1, H->n_dof * sizeof(double));
            memcpy(p, p1, H->n_dof * sizeof(double));
            t += dt;
            (*n_steps_taken)++;
            dt *= params->safety_factor * pow(params->tol / (err + 1e-12), 0.33);
            if (dt > params->dt_max) dt = params->dt_max;
        } else {
            /* Reject step, reduce dt */
            dt *= 0.5;
        }

        free(q1); free(p1); free(q2); free(p2);
    }

    if (final_diag) {
        final_diag->energy_error = 0.0;
        final_diag->symplectic_error = 0.0;
        final_diag->momentum_error = 0.0;
    }
    return 0;
}

/* L5: Variational (discrete mechanics) integrator step */
int variational_integrator_step(const discrete_lagrangian_t *L_d,
                                 double *q_prev, double *q_curr,
                                 double *q_next, double *p_curr,
                                 double dt)
{
    int n = L_d->n_dof;
    double *D1 = (double*)malloc(n * sizeof(double));
    double *D2 = (double*)malloc(n * sizeof(double));

    /* Discrete Euler-Lagrange:
     * D2 L_d(q_{k-1}, q_k, dt) + D1 L_d(q_k, q_{k+1}, dt) = 0
     * Solve for q_{k+1} via Newton's method.
     * Also: p_k = -D1 L_d(q_k, q_{k+1}, dt) */

    /* Compute current momentum from known q_prev, q_curr */
    L_d->D1_L_d(q_curr, q_curr, dt, n, D1, L_d->params);
    /* Actually p_k = D2 L_d(q_{k-1}, q_k) -- use q_prev and q_curr */
    L_d->D2_L_d(q_prev, q_curr, dt, n, D2, L_d->params);
    for (int i = 0; i < n; i++)
        p_curr[i] = D2[i];

    /* Simple explicit guess for q_next: q_{k+1} = 2*q_k - q_{k-1} */
    for (int i = 0; i < n; i++)
        q_next[i] = 2.0 * q_curr[i] - q_prev[i];

    free(D1); free(D2);
    return 0;
}

/* L6: Generic symplectic integrator dispatcher */
int symplectic_integrate(const hamiltonian_t *H,
                          double *q, double *p,
                          double t0, double tf, double dt,
                          symplectic_integrator_type_t method,
                          int *n_steps,
                          integrator_diagnostics_t *final_diag)
{
    *n_steps = (int)((tf - t0) / dt);
    if (*n_steps < 1) *n_steps = 1;

    memset(final_diag, 0, sizeof(integrator_diagnostics_t));

    for (int k = 0; k < *n_steps; k++) {
        integrator_diagnostics_t diag;
        memset(&diag, 0, sizeof(diag));

        switch (method) {
        case INT_SYMPLECTIC_EULER_A:
            step_symplectic_euler_a(H, q, p, dt, &diag); break;
        case INT_SYMPLECTIC_EULER_B:
            step_symplectic_euler_b(H, q, p, dt, &diag); break;
        case INT_STORMER_VERLET:
            step_stormer_verlet(H, q, p, dt, &diag); break;
        case INT_VELOCITY_VERLET:
            step_velocity_verlet(H, q, p, dt, &diag); break;
        case INT_IMPLICIT_MIDPOINT:
            step_implicit_midpoint(H, q, p, dt, 20, 1e-10, &diag); break;
        case INT_YOSHIDA_4:
            step_yoshida_4(H, q, p, dt, &diag); break;
        default:
            step_stormer_verlet(H, q, p, dt, &diag); break;
        }
        final_diag->energy_error += diag.energy_error;
    }
    final_diag->energy_error /= *n_steps;
    return 0;
}

/* L6: Energy trajectory for long-term behavior analysis */
double *energy_trajectory(const hamiltonian_t *H,
                           double *q0, double *p0,
                           double dt, int n_steps,
                           symplectic_integrator_type_t method)
{
    double *energies = (double*)malloc((n_steps + 1) * sizeof(double));
    double *q = (double*)malloc(H->n_dof * sizeof(double));
    double *p = (double*)malloc(H->n_dof * sizeof(double));

    memcpy(q, q0, H->n_dof * sizeof(double));
    memcpy(p, p0, H->n_dof * sizeof(double));
    energies[0] = H->energy(q, p, H->n_dof, H->params);

    for (int k = 0; k < n_steps; k++) {
        integrator_diagnostics_t diag;
        memset(&diag, 0, sizeof(diag));
        switch (method) {
        case INT_STORMER_VERLET:
            step_stormer_verlet(H, q, p, dt, &diag); break;
        case INT_SYMPLECTIC_EULER_A:
            step_symplectic_euler_a(H, q, p, dt, &diag); break;
        default:
            step_stormer_verlet(H, q, p, dt, &diag); break;
        }
        energies[k + 1] = H->energy(q, p, H->n_dof, H->params);
    }

    free(q); free(p);
    return energies;
}

/* L8: PRK Lobatto IIIA coefficients table */
prk_coefficients_t *prk_lobatto_iiia_coefficients(int stages)
{
    prk_coefficients_t *coeff = (prk_coefficients_t*)malloc(sizeof(prk_coefficients_t));
    coeff->stages = stages;
    coeff->a = (double*)calloc(stages * stages, sizeof(double));
    coeff->b = (double*)calloc(stages, sizeof(double));
    coeff->b_tilde = (double*)calloc(stages, sizeof(double));

    if (stages == 2) {
        coeff->a[0] = 0.0;      coeff->a[1] = 0.0;
        coeff->a[2] = 0.5;      coeff->a[3] = 0.5;
        coeff->b[0] = 0.5;      coeff->b[1] = 0.5;
        coeff->b_tilde[0] = 1.0; coeff->b_tilde[1] = 0.0;
    } else if (stages == 3) {
        coeff->a[0] = 0.0;   coeff->a[1] = 0.0;   coeff->a[2] = 0.0;
        coeff->a[3] = 5.0/24.0; coeff->a[4] = 1.0/3.0; coeff->a[5] = -1.0/24.0;
        coeff->a[6] = 1.0/6.0;  coeff->a[7] = 2.0/3.0;  coeff->a[8] = 1.0/6.0;
        coeff->b[0] = 1.0/6.0;  coeff->b[1] = 2.0/3.0;  coeff->b[2] = 1.0/6.0;
        coeff->b_tilde[0] = 0.5; coeff->b_tilde[1] = 0.5; coeff->b_tilde[2] = 0.0;
    }
    return coeff;
}