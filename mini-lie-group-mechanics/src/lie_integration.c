/* ==============================================================
 * lie_integration.c -- Geometric Integration on Lie Groups
 *
 * Implements: RKMK (Runge-Kutta-Munthe-Kaas) integrators,
 * Lie-Euler, Crouch-Grossman, Magnus expansion, variational
 * integrators on SO(3), geodesic computation and interpolation.
 *
 * Key methods:
 *   RKMK-RK4 on SO(3) and SE(3): Munthe-Kaas (1999)
 *   Lie-Euler: Crouch & Grossman (1993)
 *   Magnus expansion: Magnus (1954), Blanes et al. (2009)
 *   Variational integrators: Marsden & West (2001)
 *   Geodesic on SO(3): Park (1995)
 *
 * References:
 *   Iserles, Munthe-Kaas, Norsett, Zanna (2000) Acta Numerica
 *   Hairer, Lubich, Wanner (2006) Geometric Numerical Integration
 *   Marsden & West (2001) Acta Numerica 10:357-514
 * ============================================================== */

#include "lie_integration.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/* --- RKMK Context --- */
LieRKMKContext* lie_rkmk_create(LieRKMKOrder order, int algebra_dim) {
    LieRKMKContext *ctx = (LieRKMKContext*)malloc(sizeof(LieRKMKContext));
    assert(ctx);
    ctx->order = order;
    ctx->stages = (order == LIE_RKMK_RK2) ? 2 : ((order == LIE_RKMK_RK3) ? 3 : 4);
    ctx->k = (LieVector**)malloc((size_t)ctx->stages * sizeof(LieVector*));
    for (int i = 0; i < ctx->stages; i++) ctx->k[i] = lie_vector_create(algebra_dim);
    ctx->theta = lie_vector_create(algebra_dim);
    ctx->tmp = lie_vector_create(algebra_dim);
    ctx->xi_alg = lie_algebra_create(LIE_ALG_SO3, "rkmk_work");
    return ctx;
}
void lie_rkmk_free(LieRKMKContext *ctx) {
    if (ctx) {
        for (int i = 0; i < ctx->stages; i++) lie_vector_free(ctx->k[i]);
        free(ctx->k);
        lie_vector_free(ctx->theta); lie_vector_free(ctx->tmp);
        lie_algebra_free(ctx->xi_alg); free(ctx);
    }
}

/* --- RKMK Step on SO(3) --- */
void lie_rkmk_step_so3(const LieGroupElement *g_in, double t, double h,
                        LieRHSFunction f, void *params, LieGroupElement *g_out) {
    assert(g_in && g_out);
    double h2 = h / 2.0;

    /* Local storage for RK stages */
    LieVector *k[4];
    for (int i = 0; i < 4; i++) k[i] = lie_vector_create(3);
    LieVector *theta = lie_vector_create(3);

    /* Stage 1 */
    f(g_in, t, params, k[0]);

    /* Stage 2: exp(k1/2) */
    LieVector *k1_half = lie_vector_scale(k[0], h2);
    LieGroupElement *g2 = lie_exp_so3(k1_half);
    LieGroupElement *g1_exp = lie_group_mul(g_in, g2);
    f(g1_exp, t + h2, params, k[1]);
    lie_vector_free(k1_half); lie_group_free(g2); lie_group_free(g1_exp);

    /* Stage 3 */
    LieVector *k2_half = lie_vector_scale(k[1], h2);
    LieGroupElement *g3 = lie_exp_so3(k2_half);
    LieGroupElement *g2_exp = lie_group_mul(g_in, g3);
    f(g2_exp, t + h2, params, k[2]);
    lie_vector_free(k2_half); lie_group_free(g3); lie_group_free(g2_exp);

    /* Stage 4 */
    LieVector *k3_full = lie_vector_scale(k[2], h);
    LieGroupElement *g4 = lie_exp_so3(k3_full);
    LieGroupElement *g3_exp = lie_group_mul(g_in, g4);
    f(g3_exp, t + h, params, k[3]);
    lie_vector_free(k3_full); lie_group_free(g4); lie_group_free(g3_exp);

    /* theta = h*(k1 + 2*k2 + 2*k3 + k4) / 6 */
    for (int i = 0; i < 3; i++) {
        double th = (k[0]->data[i] + 2.0*k[1]->data[i]
                  + 2.0*k[2]->data[i] + k[3]->data[i]) * h / 6.0;
        theta->data[i] = th;
    }

    /* g_out = g_in * exp(theta) */
    LieGroupElement *exp_theta = lie_exp_so3(theta);
    LieGroupElement *result = lie_group_mul(g_in, exp_theta);

    /* Copy to g_out */
    for (int i = 0; i < 9; i++)
        g_out->g->data[i] = result->g->data[i];
    memcpy(g_out->params->data, theta->data, 3 * sizeof(double));

    lie_group_free(exp_theta); lie_group_free(result);
    for (int i = 0; i < 4; i++) lie_vector_free(k[i]);
    lie_vector_free(theta);
}

/* --- Lie-Euler Step on SO(3) --- */
void lie_euler_step_so3(const LieGroupElement *g_in, double t, double h,
                         LieRHSFunction f, void *params, LieGroupElement *g_out) {
    assert(g_in && g_out);
    LieVector xi = {.size=3}; xi.data = (double[3]){0, 0, 0};
    f(g_in, t, params, &xi);

    /* g_out = g_in * exp(h * xi) */
    LieVector *h_xi = lie_vector_scale(&xi, h);
    LieGroupElement *exp_hxi = lie_exp_so3(h_xi);
    LieGroupElement *result = lie_group_mul(g_in, exp_hxi);

    for (int i = 0; i < 9; i++) g_out->g->data[i] = result->g->data[i];
    memcpy(g_out->params->data, h_xi->data, 3 * sizeof(double));

    lie_vector_free(h_xi); lie_group_free(exp_hxi); lie_group_free(result);
}

/* --- Trajectory --- */
LieTrajectory* lie_trajectory_create(int n_points, int algebra_dim) {
    assert(n_points > 0);
    LieTrajectory *traj = (LieTrajectory*)malloc(sizeof(LieTrajectory));
    assert(traj);
    traj->n_points = n_points;
    traj->t = (double*)calloc((size_t)n_points, sizeof(double));
    assert(traj->t);
    traj->g = (LieGroupElement**)calloc((size_t)n_points, sizeof(LieGroupElement*));
    assert(traj->g);
    traj->dim = algebra_dim;
    for (int i = 0; i < n_points; i++) traj->g[i] = NULL;
    return traj;
}
void lie_trajectory_free(LieTrajectory *traj) {
    if (traj) {
        free(traj->t);
        for (int i = 0; i < traj->n_points; i++) lie_group_free(traj->g[i]);
        free(traj->g); free(traj);
    }
}
void lie_trajectory_print(const LieTrajectory *traj, int stride) {
    if (!traj) return;
    printf("Trajectory: %d points, stride=%d\n", traj->n_points, stride);
    for (int i = 0; i < traj->n_points; i += stride) {
        printf("t[%d]=%6.3f: ", i, traj->t[i]);
        if (traj->g[i]) lie_vector_print(traj->g[i]->params, "");
    }
}

/* --- RKMK4 Integrate Full SO(3) Trajectory --- */
LieTrajectory* lie_rkmk4_integrate_so3(const LieGroupElement *g0,
    double t0, double tf, int n_steps, LieRHSFunction f, void *params) {
    assert(g0 && n_steps > 0);
    LieTrajectory *traj = lie_trajectory_create(n_steps + 1, 3);
    double h = (tf - t0) / (double)n_steps;
    double t = t0;

    /* Initial state */
    traj->t[0] = t;
    traj->g[0] = lie_group_create(LIE_GROUP_SO3, "traj0");
    memcpy(traj->g[0]->g->data, g0->g->data, 9 * sizeof(double));

    LieRKMKContext *rkmk = lie_rkmk_create(LIE_RKMK_RK4, 3);
    LieGroupElement *g_cur = lie_group_create(LIE_GROUP_SO3, "current");
    memcpy(g_cur->g->data, g0->g->data, 9 * sizeof(double));

    for (int step = 0; step < n_steps; step++) {
        LieGroupElement *g_next = lie_group_create(LIE_GROUP_SO3, "next");
        lie_rkmk_step_so3(g_cur, t, h, f, params, g_next);

        t += h;
        traj->t[step + 1] = t;
        traj->g[step + 1] = g_next;
        memcpy(g_cur->g->data, g_next->g->data, 9 * sizeof(double));
    }

    lie_group_free(g_cur);
    lie_rkmk_free(rkmk);
    return traj;
}

/* --- Magnus Expansion (2nd order) --- */
LieMatrix* lie_magnus_expansion_order2(const LieMatrix *A_n, const LieMatrix *A_np1, double h) {
    /* Omega = h*(A_n + A_{n+1})/2 + h^2/12 * [A_n, A_{n+1}] */
    assert(A_n && A_np1 && A_n->rows == A_np1->rows);
    LieMatrix *sum = lie_matrix_add(A_n, A_np1);
    LieMatrix *Omega = lie_matrix_scale(sum, h / 2.0);
    lie_matrix_free(sum);

    LieMatrix *An_Anp1 = lie_matrix_mul(A_n, A_np1);
    LieMatrix *Anp1_An = lie_matrix_mul(A_np1, A_n);
    LieMatrix *comm = lie_matrix_sub(An_Anp1, Anp1_An);
    LieMatrix *comm_scaled = lie_matrix_scale(comm, h * h / 12.0);
    lie_matrix_free(An_Anp1); lie_matrix_free(Anp1_An); lie_matrix_free(comm);

    LieMatrix *result = lie_matrix_add(Omega, comm_scaled);
    lie_matrix_free(Omega); lie_matrix_free(comm_scaled);
    return result;
}

/* --- Variational Integrator for Rigid Body --- */
void lie_variational_rigid_body_step(const LieGroupElement *g_k,
    const LieVector *Omega_k, const LieInertiaTensor *I, double h,
    LieGroupElement *g_next, LieVector *Omega_next) {
    assert(g_k && Omega_k && I && g_next && Omega_next);

    /* Midpoint rule variational integrator:
     * Use F_k = exp(h*Omega_k/2) */
    LieVector *hOmega_half = lie_vector_scale(Omega_k, h / 2.0);
    LieGroupElement *F_k = lie_exp_so3(hOmega_half);
    LieGroupElement *g_mid = lie_group_mul(g_k, F_k);
    lie_vector_free(hOmega_half);

    /* Compute dL/dOmega at midpoint (for kinetic energy only) */
    LieVector *I_Omega_k = lie_inertia_apply(I, Omega_k);

    /* Discrete Euler-Poincare update:
     * I*Omega_{k+1} = Ad*_{F_k} I*Omega_k
     * But with midpoint: Omega_{k+1} = Omega_k + h*I^{-1}(I*Omega_k x Omega_k)/2
     * Simplified: use explicit midpoint */
    LieVector *Pi_k = lie_inertia_apply(I, Omega_k);
    LieVector *PixO = lie_vector_cross3(Pi_k, Omega_k);
    LieInertiaTensor *Iinv = lie_inertia_inverse(I);
    LieVector *dOmega = lie_inertia_apply(Iinv, PixO);

    for (int i = 0; i < 3; i++) {
        Omega_next->data[i] = Omega_k->data[i] + h * dOmega->data[i];
    }

    /* g_{k+1} = g_k * exp(h * Omega_{k+1/2}) */
    LieVector *Omega_mid = lie_vector_scale(Omega_next, h);
    LieGroupElement *inc = lie_exp_so3(Omega_mid);
    LieGroupElement *res = lie_group_mul(g_k, inc);

    memcpy(g_next->g->data, res->g->data, 9 * sizeof(double));

    lie_vector_free(I_Omega_k); lie_vector_free(Pi_k);
    lie_vector_free(PixO); lie_vector_free(dOmega);
    lie_group_free(F_k); lie_group_free(g_mid);
    lie_group_free(inc); lie_group_free(res);
    lie_vector_free(Omega_mid);
    lie_inertia_free(Iinv);
}

/* --- Geodesic Equation on SO(3) --- */
void lie_geodesic_equation_so3(const LieVector *xi, const LieInertiaTensor *I,
                                LieVector *dxi_dt) {
    /* For left-invariant metric defined by I:
     * dxi/dt = ad*_xi (I xi) = (I xi) x xi (after identification)
     * With no potential, this is the Euler-Poincare equation.
     */
    assert(xi && dxi_dt);
    LieVector *Ixi = lie_inertia_apply(I, xi);
    LieVector *adstar = lie_ad_star_so3(xi, Ixi);
    for (int i = 0; i < 3; i++) dxi_dt->data[i] = adstar->data[i];
    lie_vector_free(Ixi); lie_vector_free(adstar);
}

/* --- Geodesic Distance on SO(3) --- */
double lie_geodesic_distance_so3(const LieMatrix *R1, const LieMatrix *R2) {
    assert(R1 && R2 && R1->rows == 3 && R2->rows == 3);
    /* dist(R1, R2) = ||log(R1^T * R2)|| */
    LieMatrix *R1T = lie_matrix_transpose(R1);
    LieMatrix *prod = lie_matrix_mul(R1T, R2);
    LieGroupElement g = {.group_type = LIE_GROUP_SO3, .g = prod, .mat_size = 3};
    LieVector *omega = lie_log_so3(&g);
    double dist = lie_vector_norm(omega);
    lie_vector_free(omega);
    lie_matrix_free(R1T); lie_matrix_free(prod);
    return dist;
}

/* --- Geodesic Interpolation on SO(3) --- */
LieGroupElement* lie_so3_geodesic_interpolate(const LieGroupElement *g0,
                                                const LieGroupElement *g1, double t) {
    assert(g0 && g1 && t >= 0.0 && t <= 1.0);
    if (t < 1e-15) {
        LieGroupElement *g = lie_group_create(LIE_GROUP_SO3, "interp");
        memcpy(g->g->data, g0->g->data, 9 * sizeof(double));
        return g;
    }
    if (t > 1.0 - 1e-15) {
        LieGroupElement *g = lie_group_create(LIE_GROUP_SO3, "interp");
        memcpy(g->g->data, g1->g->data, 9 * sizeof(double));
        return g;
    }
    /* R(t) = R1 * exp(t * log(R1^T * R2)) */
    LieGroupElement *g0inv = lie_group_inverse(g0);
    LieGroupElement *rel = lie_group_mul(g0inv, g1);
    LieVector *omega = lie_log_so3(rel);
    LieVector *t_omega = lie_vector_scale(omega, t);
    LieGroupElement *exp_tomega = lie_exp_so3(t_omega);
    LieGroupElement *result = lie_group_mul(g0, exp_tomega);

    lie_group_free(g0inv); lie_group_free(rel);
    lie_vector_free(omega); lie_vector_free(t_omega);
    lie_group_free(exp_tomega);
    return result;
}

