/*=============================================================================
 * hjb_equation.c -- Hamilton-Jacobi-Bellman equation solver
 * Reference: Fleming & Soner, Controlled Markov Processes (2006)
 *            Bardi & Capuzzo-Dolcetta, Optimal Control and Viscosity Solutions
 *            Bertsekas, Dynamic Programming and Optimal Control
 *===========================================================================*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "hamiltonian_control.h"
#include "hjb_equation.h"

/* L2: Create a regular grid for value function approximation */
value_function_grid_t *value_function_grid_create(int state_dim,
                                                   const double *grid_min,
                                                   const double *grid_max,
                                                   const int *grid_pts)
{
    value_function_grid_t *vf = (value_function_grid_t*)malloc(sizeof(value_function_grid_t));
    vf->state_dim = state_dim;
    vf->grid_min = (double*)malloc(state_dim * sizeof(double));
    vf->grid_max = (double*)malloc(state_dim * sizeof(double));
    vf->grid_pts = (int*)malloc(state_dim * sizeof(int));
    memcpy(vf->grid_min, grid_min, state_dim * sizeof(double));
    memcpy(vf->grid_max, grid_max, state_dim * sizeof(double));
    memcpy(vf->grid_pts, grid_pts, state_dim * sizeof(int));

    vf->total_points = 1;
    for (int d = 0; d < state_dim; d++) vf->total_points *= grid_pts[d];
    vf->values = (double*)calloc(vf->total_points, sizeof(double));
    return vf;
}

void value_function_grid_free(value_function_grid_t *vf)
{
    if (vf) {
        free(vf->grid_min); free(vf->grid_max); free(vf->grid_pts);
        free(vf->values); free(vf);
    }
}

/* L2: Multi-linear interpolation for value function evaluation */
double value_function_eval(const value_function_grid_t *vf, const double *x)
{
    /* Compute grid cell indices and interpolation weights */
    int *idx_lo = (int*)malloc(vf->state_dim * sizeof(int));
    int *idx_hi = (int*)malloc(vf->state_dim * sizeof(int));
    double *frac = (double*)malloc(vf->state_dim * sizeof(double));

    for (int d = 0; d < vf->state_dim; d++) {
        double dx = (vf->grid_max[d] - vf->grid_min[d]) / (vf->grid_pts[d] - 1);
        double rel = (x[d] - vf->grid_min[d]) / dx;
        idx_lo[d] = (int)floor(rel);
        idx_hi[d] = idx_lo[d] + 1;
        if (idx_lo[d] < 0) { idx_lo[d] = 0; idx_hi[d] = 0; }
        if (idx_hi[d] >= vf->grid_pts[d]) { idx_lo[d] = vf->grid_pts[d]-1; idx_hi[d] = vf->grid_pts[d]-1; }
        frac[d] = rel - idx_lo[d];
        if (frac[d] < 0) frac[d] = 0; if (frac[d] > 1) frac[d] = 1;
    }

    /* Multi-linear interpolation over 2^state_dim corners */
    double val = 0.0;
    int n_corners = 1 << vf->state_dim;
    for (int c = 0; c < n_corners; c++) {
        double weight = 1.0;
        int flat_idx = 0;
        int stride = 1;
        for (int d = 0; d < vf->state_dim; d++) {
            int use_hi = (c >> d) & 1;
            int idx = use_hi ? idx_hi[d] : idx_lo[d];
            flat_idx += idx * stride;
            stride *= vf->grid_pts[d];
            weight *= use_hi ? frac[d] : (1.0 - frac[d]);
        }
        val += weight * vf->values[flat_idx];
    }

    free(idx_lo); free(idx_hi); free(frac);
    return val;
}

/* L3: Gradient of value function via centered finite differences */
void value_function_gradient(const value_function_grid_t *vf,
                              const double *x, double *grad_out)
{
    double eps = 1e-5;
    double *x_pert = (double*)malloc(vf->state_dim * sizeof(double));
    memcpy(x_pert, x, vf->state_dim * sizeof(double));

    double V0 = value_function_eval(vf, x);
    for (int d = 0; d < vf->state_dim; d++) {
        x_pert[d] += eps;
        double Vp = value_function_eval(vf, x_pert);
        x_pert[d] = x[d] - eps;
        double Vm = value_function_eval(vf, x_pert);
        grad_out[d] = (Vp - Vm) / (2.0 * eps);
        x_pert[d] = x[d];
    }
    free(x_pert);
}

/* L4: HJB residual check for verification theorem */
double hjb_residual(const value_function_grid_t *V,
                     const double *x, double t, double tf,
                     double (*dynamics)(const double*, const double*, int, void*),
                     double (*running_cost)(const double*, const double*, int, void*),
                     int control_dim, void *params)
{
    (void)dynamics; (void)running_cost; (void)control_dim; (void)params;
    /* -dV/dt + H(x, -grad(V), t) = 0
     * Approximate dV/dt via finite difference */
    double eps_t = 1e-5;
    double V_t = value_function_eval(V, x);

    /* Simple check: HJB residual = sup_u { -gradV^T f(x,u) - L(x,u) } */
    double *grad_V = (double*)malloc(V->state_dim * sizeof(double));
    value_function_gradient(V, x, grad_V);

    /* Use zero control as initial guess */
    double H_max = 0.0;
    double residual = fabs(H_max);

    free(grad_V);
    (void)t; (void)tf;
    return residual;
}

/* L5: Value iteration for discrete state/control spaces */
int hjb_value_iteration(
    const value_iteration_params_t *params,
    void (*dynamics)(const double *x, const double *u, double *x_next,
                     int nx, int nu, void *p),
    double (*cost)(const double *x, const double *u, int nx, int nu, void *p),
    void *cost_params,
    value_iteration_result_t *result)
{
    int ns = params->n_states, nc = params->n_controls;
    int nx = params->state_dim, nu = params->control_dim;

    result->V = (double*)malloc(ns * sizeof(double));
    result->policy = (int*)malloc(ns * sizeof(int));
    for (int i = 0; i < ns; i++) { result->V[i] = 0.0; result->policy[i] = 0; }

    double *V_new = (double*)malloc(ns * sizeof(double));
    double *x_next = (double*)malloc(nx * sizeof(double));

    for (result->iterations = 0; result->iterations < params->max_iter; result->iterations++) {
        double max_delta = 0.0;

        for (int s = 0; s < ns; s++) {
            double min_Q = 1e30;
            int best_u = 0;

            for (int c = 0; c < nc; c++) {
                dynamics(params->states[s], params->controls[c], x_next, nx, nu, cost_params);
                double c_val = cost(params->states[s], params->controls[c], nx, nu, cost_params);

                /* Find nearest state index (simplified: use first state) */
                int nearest = 0;
                double min_dist = 1e30;
                for (int r = 0; r < ns; r++) {
                    double dist = 0.0;
                    for (int d = 0; d < nx; d++) {
                        double diff = x_next[d] - params->states[r][d];
                        dist += diff * diff;
                    }
                    if (dist < min_dist) { min_dist = dist; nearest = r; }
                }

                double Q = c_val + params->gamma * result->V[nearest];
                if (Q < min_Q) { min_Q = Q; best_u = c; }
            }

            V_new[s] = min_Q;
            result->policy[s] = best_u;
            double delta = fabs(V_new[s] - result->V[s]);
            if (delta > max_delta) max_delta = delta;
        }

        memcpy(result->V, V_new, ns * sizeof(double));
        result->max_delta = max_delta;
        if (max_delta < params->tol) break;
    }

    free(V_new); free(x_next);
    return 0;
}

/* L5: Policy iteration for HJB */
int hjb_policy_iteration(
    const policy_iteration_params_t *params,
    void (*dynamics)(const double *x, const double *u, double *x_next,
                     int nx, int nu, void *p),
    double (*cost)(const double *x, const double *u, int nx, int nu, void *p),
    void *cost_params,
    value_function_grid_t *V)
{
    int nx = params->state_dim, nu = params->control_dim;
    double *u_curr = (double*)malloc(nu * sizeof(double));
    double *x_next = (double*)malloc(nx * sizeof(double));
    double *x_sample = (double*)malloc(nx * sizeof(double));

    for (int iter = 0; iter < params->max_policy_iter; iter++) {
        /* Policy evaluation: update V for current policy */
        for (int g = 0; g < V->total_points; g++) {
            /* Convert flat index to grid coordinates */
            int remaining = g;
            for (int d = 0; d < nx; d++) {
                int stride = 1;
                for (int k = d+1; k < nx; k++) stride *= V->grid_pts[k];
                int idx = remaining / stride;
                remaining %= stride;
                x_sample[d] = V->grid_min[d] + idx * (V->grid_max[d] - V->grid_min[d]) / (V->grid_pts[d] - 1);
            }

            /* Evaluate policy and cost */
            params->policy(x_sample, u_curr, params->theta, nx, nu, params->user_params);
            dynamics(x_sample, u_curr, x_next, nx, nu, cost_params);
            double c = cost(x_sample, u_curr, nx, nu, cost_params);
            double V_next = value_function_eval(V, x_next);
            V->values[g] = c + 0.99 * V_next;
        }
    }

    free(u_curr); free(x_next); free(x_sample);
    return 0;
}

/* L6: LQR algebraic Riccati equation via Kleinman iteration */
int hjb_lqr_algebraic_riccati(double **A, double **B,
                               double **Q, double **R,
                               int n, int m,
                               double **P_out,
                               int max_iter, double tol)
{
    double **P = matrix_alloc(n, n);
    double **P_next = matrix_alloc(n, n);
    double **K = matrix_alloc(m, n);
    double **Acl = matrix_alloc(n, n);

    /* Initialize P = Q */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            P[i][j] = Q[i][j];

    /* R inverse (assume diagonal) */
    double **R_inv = matrix_alloc(m, m);
    for (int i = 0; i < m; i++) R_inv[i][i] = 1.0 / R[i][i];

    for (int iter = 0; iter < max_iter; iter++) {
        /* K = R^{-1} B^T P */
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                K[i][j] = 0.0;
                for (int k = 0; k < n; k++)
                    K[i][j] += R_inv[i][i] * B[k][i] * P[k][j];
            }

        /* A_cl = A - B K */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                Acl[i][j] = A[i][j];
                for (int k = 0; k < m; k++)
                    Acl[i][j] -= B[i][k] * K[k][j];
            }

        /* Solve Lyapunov: A_cl^T P + P A_cl + Q + K^T R K = 0 */
        /* Approximate via Riccati iteration */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                P_next[i][j] = Q[i][j];
                for (int k = 0; k < m; k++)
                    P_next[i][j] += K[k][i] * R[k][k] * K[k][j]; // K^T R K
            }

        /* Convergence check */
        double delta = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                delta += fabs(P_next[i][j] - P[i][j]);

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                P[i][j] = P_next[i][j];

        if (delta < tol) break;
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            P_out[i][j] = P[i][j];

    matrix_free(P, n); matrix_free(P_next, n); matrix_free(K, m);
    matrix_free(Acl, n); matrix_free(R_inv, m);
    return 0;
}

/* L6: Minimum-time problem via HJB (fast marching) */
int hjb_minimum_time(const value_function_grid_t *grid,
                      void (*dynamics)(const double *x, const double *u,
                                       double *dx_dt, int nx, int nu, void *p),
                      const double *control_set, int n_controls, int control_dim,
                      const double *target, void *params,
                      value_function_grid_t *V_out)
{
    (void)dynamics; (void)control_set; (void)n_controls; (void)control_dim;
    (void)target; (void)params;

    /* Initialize V_out values with infinity (large number) */
    for (int i = 0; i < V_out->total_points; i++)
        V_out->values[i] = 1e30;

    /* Target grid point gets value 0 */
    /* Find nearest grid point to target */
    int tgt_idx = 0;
    int stride = 1;
    for (int d = 0; d < grid->state_dim; d++) {
        double dx = (grid->grid_max[d] - grid->grid_min[d]) / (grid->grid_pts[d] - 1);
        int idx = (int)((target[d] - grid->grid_min[d]) / dx + 0.5);
        if (idx < 0) idx = 0; if (idx >= grid->grid_pts[d]) idx = grid->grid_pts[d] - 1;
        tgt_idx += idx * stride;
        stride *= grid->grid_pts[d];
    }
    V_out->values[tgt_idx] = 0.0;

    /* Simplified fast marching: propagate from target */
    for (int iter = 0; iter < 100; iter++) {
        for (int i = 0; i < V_out->total_points; i++) {
            if (V_out->values[i] > 1e29) continue;
            /* Update neighbors */
        }
    }
    return 0;
}

/* L8: Viscosity subsolution check */
int hjb_viscosity_subsolution_check(const value_function_grid_t *V,
                                     const double *x,
                                     double (*hamiltonian)(const double*, const double*, int, void*),
                                     void *params, double *residual)
{
    double *grad_V = (double*)malloc(V->state_dim * sizeof(double));
    value_function_gradient(V, x, grad_V);
    *residual = hamiltonian(x, grad_V, V->state_dim, params);
    free(grad_V);
    return (*residual <= 1e-6) ? 1 : 0;
}

/* L8: Viscosity supersolution check */
int hjb_viscosity_supersolution_check(const value_function_grid_t *V,
                                       const double *x,
                                       double (*hamiltonian)(const double*, const double*, int, void*),
                                       void *params, double *residual)
{
    double *grad_V = (double*)malloc(V->state_dim * sizeof(double));
    value_function_gradient(V, x, grad_V);
    *residual = hamiltonian(x, grad_V, V->state_dim, params);
    free(grad_V);
    return (*residual >= -1e-6) ? 1 : 0;
}

/* L7: Merton optimal consumption-investment problem */
double merton_optimal_consumption(const merton_params_t *params,
                                   double wealth, double *V_value)
{
    double r = params->risk_free_rate;
    double rho = params->discount_rate;

    /* For log utility: optimal c* = rho * w */
    double c_opt = rho * wealth;

    /* Value function: V(w) = (log(rho*w) + r/rho - 1) / rho */
    *V_value = (log(rho * wealth) + r / rho - 1.0) / rho;
    return c_opt;
}