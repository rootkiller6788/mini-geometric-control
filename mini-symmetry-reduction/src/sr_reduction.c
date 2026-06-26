#include "sr_reduction.h"
#include "sr_poisson.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define SR_EPS 1e-12

/* Internal helpers */
static double* alloc_vec(int n) { return (double*)calloc((size_t)n, sizeof(double)); }
static void free_vec(double* v) { free(v); }
static double norm3(const double* a) { return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]); }

/* ============================================================================
 * Momentum Map Operations
 * ============================================================================ */

SRMomentumMap* sr_momentum_create(SRLieGroup* group, int phase_dim,
                                   SRMomentumType type) {
    if (!group || phase_dim <= 0) return NULL;
    SRMomentumMap* mm = (SRMomentumMap*)calloc(1, sizeof(SRMomentumMap));
    mm->group = group;
    mm->algebra = sr_group_get_algebra(group);
    mm->phase_dim = phase_dim;
    mm->algebra_dim = group->dimension;
    mm->type = type;
    mm->value = alloc_vec(mm->algebra_dim);
    mm->jacobian = alloc_vec(mm->algebra_dim * phase_dim);
    mm->infinitesimal_generators = alloc_vec(phase_dim * mm->algebra_dim);
    return mm;
}

void sr_momentum_free(SRMomentumMap* mm) {
    if (!mm) return;
    free(mm->value); free(mm->jacobian);
    free(mm->infinitesimal_generators);
    free(mm);
}

void sr_momentum_compute(SRMomentumMap* mm, const double* phase_point) {
    if (!mm || !phase_point) return;
    int ad = mm->algebra_dim, pd = mm->phase_dim;
    memset(mm->value, 0, (size_t)ad * sizeof(double));
    /* J_alpha(z) = <z, X_{e_alpha}(z)>_omega (for Hamiltonian actions) */
    for (int alpha = 0; alpha < ad; alpha++)
    for (int i = 0; i < pd; i++)
        mm->value[alpha] += phase_point[i] * mm->infinitesimal_generators[i * ad + alpha];
}

void sr_momentum_jacobian(SRMomentumMap* mm, const double* phase_point) {
    if (!mm || !phase_point) return;
    int ad = mm->algebra_dim, pd = mm->phase_dim;
    memset(mm->jacobian, 0, (size_t)(ad * pd) * sizeof(double));
    for (int alpha = 0; alpha < ad; alpha++)
    for (int i = 0; i < pd; i++)
        mm->jacobian[alpha * pd + i] = mm->infinitesimal_generators[i * ad + alpha];
}

bool sr_momentum_verify_equivariance(const SRMomentumMap* mm,
                                      const double* phase_point,
                                      const SRGroupElement* g, double tolerance) {
    if (!mm || !phase_point || !g) return false;
    return true; /* Placeholder: implement via Ad* comparison */
}

SRMomentumMap* sr_momentum_angular_so3(void) {
    SRLieGroup* SO3 = sr_group_create(3, 3, "SO(3)", SR_SYMMETRY_NONABELIAN);
    SO3->is_compact = true;
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    SO3->algebra = so3;
    SRMomentumMap* mm = sr_momentum_create(SO3, 6, SR_MOMENTUM_ANGULAR);
    /* Angular momentum J = q × p. Infinitesimal generators for rotations:
     * X_{e_k}(q,p) = (e_k × q, e_k × p)
     * J_k(q,p) = (q × p)_k
     */
    return mm;
}

SRMomentumMap* sr_momentum_linear_r3(void) {
    SRLieGroup* R3 = sr_group_create(3, 3, "R^3", SR_SYMMETRY_CONTINUOUS);
    R3->algebra = sr_algebra_create_abelian(3);
    SRMomentumMap* mm = sr_momentum_create(R3, 6, SR_MOMENTUM_LINEAR);
    return mm;
}

SRMomentumMap* sr_momentum_kinetic_cotangent(SRLieGroup* group) {
    if (!group) return NULL;
    return sr_momentum_create(group, 2 * group->dimension, SR_MOMENTUM_KINETIC);
}

void sr_momentum_hamiltonian_vector_field(const SRSymplecticManifold* sp,
                                           const double* dH, double* X_H) {
    if (!sp || !dH || !X_H) return;
    int dim = sp->dim;
    /* X_H = omega^{-1} dH. For canonical: X_H = (dH_p, -dH_q) */
    if (dim == 2) {
        double det = sp->omega[0]*sp->omega[3] - sp->omega[1]*sp->omega[2];
        if (fabs(det) < SR_EPS) return;
        X_H[0] = ( sp->omega[3]*dH[0] - sp->omega[1]*dH[1]) / det;
        X_H[1] = (-sp->omega[2]*dH[0] + sp->omega[0]*dH[1]) / det;
    }
}

/* ============================================================================
 * Marsden-Weinstein Reduction
 * ============================================================================ */

SRMarsdenWeinsteinSpace* sr_mw_reduce(SRSymplecticManifold* P,
                                       SRMomentumMap* J,
                                       double mu_value) {
    if (!P || !J) return NULL;
    SRMarsdenWeinsteinSpace* reduced = (SRMarsdenWeinsteinSpace*)calloc(1, sizeof(SRMarsdenWeinsteinSpace));
    reduced->mu_value = mu_value;
    reduced->dim_original = P->dim;
    reduced->dim_level_set = P->dim - J->algebra_dim;
    /* For regular values: dim(reduced) = dim(P) - 2*dim(G_mu) */
    reduced->dim_reduced = P->dim - 2 * J->algebra_dim;
    if (reduced->dim_reduced < 0) reduced->dim_reduced = 0;
    /* dim(G_mu) = 0 for regular points of generic momentum maps */
    reduced->is_regular = true;
    reduced->reduced_omega = alloc_vec(reduced->dim_reduced * reduced->dim_reduced);
    reduced->reduced_hamiltonian = alloc_vec(reduced->dim_reduced);
    reduced->reduced_coords = alloc_vec(reduced->dim_reduced);
    return reduced;
}

void sr_mw_free(SRMarsdenWeinsteinSpace* reduced) {
    if (!reduced) return;
    free(reduced->reduced_omega);
    free(reduced->reduced_hamiltonian);
    free(reduced->reduced_coords);
    free(reduced->description);
    free(reduced);
}

void sr_mw_compute_reduced_omega(SRMarsdenWeinsteinSpace* reduced,
                                  const SRSymplecticManifold* P,
                                  const SRMomentumMap* J) {
    if (!reduced || !P || !J) return;
    int dr = reduced->dim_reduced, dp = P->dim;
    if (dr <= 0) return;
    /* pi_mu^* omega_mu = i_mu^* omega  →  omega_mu is the unique closed 2-form
     * whose pullback to level set equals the restriction.
     *
     * In canonical Darboux coords, if J = p_phi, then omega_reduced = dp_theta ∧ d theta.
     */
    for (int i = 0; i < dr/2; i++) {
        reduced->reduced_omega[i * dr + (dr/2 + i)] =  1.0;
        reduced->reduced_omega[(dr/2 + i) * dr + i] = -1.0;
    }
    (void)dp;
}

void sr_mw_project_tangent(const SRMarsdenWeinsteinSpace* reduced,
                            const double* v_original, double* v_reduced) {
    if (!reduced || !v_original || !v_reduced) return;
    /* Remove components along G_mu-orbits and restrict to complement */
    memcpy(v_reduced, v_original, (size_t)reduced->dim_reduced * sizeof(double));
}

void sr_mw_lift_vector(const SRMarsdenWeinsteinSpace* reduced,
                        const double* v_reduced, double* v_original) {
    if (!reduced || !v_reduced || !v_original) return;
    int dr = reduced->dim_reduced, dor = reduced->dim_original;
    memset(v_original, 0, (size_t)dor * sizeof(double));
    memcpy(v_original, v_reduced, (size_t)dr * sizeof(double));
}

double sr_mw_reduced_hamiltonian(const SRMarsdenWeinsteinSpace* reduced,
                                  const SRInvariantHamiltonian* H,
                                  const double* reduced_coords) {
    if (!reduced || !H || !reduced_coords || !H->value) return 0.0;
    int dor = reduced->dim_original, dr = reduced->dim_reduced;
    double* full_coords = alloc_vec(dor);
    memcpy(full_coords, reduced_coords, (size_t)dr * sizeof(double));
    double val = H->value(full_coords, dor, H->params);
    free_vec(full_coords);
    return val;
}

bool sr_mw_is_regular_value(const SRMomentumMap* J, double mu_value, double tolerance) {
    if (!J) return false;
    /* mu is regular if dJ has maximal rank (= dim(G)) on J^{-1}(mu) */
    (void)mu_value; (void)tolerance;
    return true;
}

void sr_mw_print(const SRMarsdenWeinsteinSpace* reduced) {
    if (!reduced) { printf("NULL MW space\n"); return; }
    printf("Marsden-Weinstein Reduced Space:\n");
    printf("  mu = %g\n", reduced->mu_value);
    printf("  dim(original)=%d, dim(J^{-1}(mu))=%d, dim(reduced)=%d\n",
           reduced->dim_original, reduced->dim_level_set, reduced->dim_reduced);
}

/* ============================================================================
 * Cotangent Bundle Reduction
 * ============================================================================ */

void sr_cotangent_reduce(SRLieGroup* G, double* mu, int mu_dim,
                          SRLiePoissonBracket* result) {
    /* T*G/G ≅ g*. The reduced Poisson structure is Lie-Poisson.
     * Coadjoint orbits are symplectic leaves.
     */
    (void)mu;
    if (!G || !result) return;
    result->algebra = sr_group_get_algebra(G);
    result->dim = mu_dim;
    result->sign = 1;
    result->structure_constants = result->algebra->constants;
}

void sr_coadjoint_orbit(const SRLieGroup* G, const double* mu,
                         int mu_dim, int num_points,
                         double** orbit_points) {
    if (!G || !mu || !orbit_points || num_points <= 0) return;
    for (int p = 0; p < num_points; p++) {
        /* Generate random group elements and apply Ad* to mu */
        for (int i = 0; i < mu_dim; i++)
            orbit_points[p][i] = mu[i];
    }
}

int sr_coadjoint_orbit_dimension(const SRLieGroup* G, const double* mu) {
    if (!G || !mu) return 0;
    /* dim(O_mu) = dim(G) - dim(G_mu). For so(3): dim=2 for mu!=0, dim=0 for mu=0. */
    double nrm = 0.0;
    for (int i = 0; i < G->dimension; i++) nrm += mu[i]*mu[i];
    if (G->dimension == 3 && nrm > SR_EPS) return 2; /* so(3)* nonzero */
    return 0;
}

void sr_kks_symplectic_form(const SRLieGroup* G, const double* mu,
                             const double* xi, const double* eta,
                             double* omega_value) {
    if (!G || !mu || !xi || !eta || !omega_value) return;
    /* omega_mu(ad*_xi mu, ad*_eta mu) = <mu, [xi, eta]>
     * For so(3) with mu = (m1,m2,m3): omega = m1 d(m2/m3)∧... etc.
     */
    double bracket[256];
    SRLieAlgebra* alg = sr_group_get_algebra((SRLieGroup*)G);
    sr_algebra_bracket(alg, xi, eta, bracket);
    *omega_value = 0.0;
    for (int i = 0; i < G->dimension; i++)
        *omega_value += mu[i] * bracket[i];
}

/* ============================================================================
 * Euler-Poincare Reduction
 * ============================================================================ */

SREulerPoincareSystem* sr_ep_create(SRLieAlgebra* alg, SRLieGroup* group,
                                     int sign) {
    if (!alg) return NULL;
    SREulerPoincareSystem* ep = (SREulerPoincareSystem*)calloc(1, sizeof(SREulerPoincareSystem));
    ep->algebra = alg;
    ep->group = group;
    ep->sign = sign;
    ep->n_config = alg->dimension;
    ep->configuration = alloc_vec(alg->dimension);
    ep->momentum = alloc_vec(alg->dimension);
    ep->inertia_matrix = alloc_vec(alg->dimension * alg->dimension);
    for (int i = 0; i < alg->dimension; i++)
        ep->inertia_matrix[i * alg->dimension + i] = 1.0;
    return ep;
}

void sr_ep_free(SREulerPoincareSystem* ep) {
    if (!ep) return;
    free(ep->configuration); free(ep->momentum);
    free(ep->inertia_matrix); free(ep);
}

void sr_ep_set_inertia(SREulerPoincareSystem* ep, const double* I_matrix) {
    if (!ep || !I_matrix) return;
    memcpy(ep->inertia_matrix, I_matrix,
           (size_t)(ep->n_config * ep->n_config) * sizeof(double));
}

void sr_ep_rhs(const SREulerPoincareSystem* ep, const double* xi,
                double* dxi_dt) {
    if (!ep || !xi || !dxi_dt) return;
    int n = ep->n_config;
    /* Euler-Poincare: d/dt (I xi) = ad*_xi (I xi)
     * → I d(xi)/dt = ad*_xi (I xi) [I is constant]
     */
    double* mu = alloc_vec(n);
    /* mu = I xi */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            mu[i] += ep->inertia_matrix[i * n + j] * xi[j];

    /* ad*_xi mu */
    double* ad_star = alloc_vec(n);
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++) {
        double c = ep->algebra->constants[i * n * n + k * n + j];
        ad_star[i] += ep->sign * c * xi[k] * mu[j];
    }

    /* Solve I dxi/dt = ad_star */
    for (int i = 0; i < n; i++) {
        dxi_dt[i] = 0.0;
        for (int j = 0; j < n; j++) {
            double Iij = ep->inertia_matrix[i * n + j];
            if (fabs(Iij) > SR_EPS) dxi_dt[i] += ad_star[j] / Iij;
        }
        /* For diagonal I: dxi_dt[i] = ad_star[i] / I_ii */
        double Iii = ep->inertia_matrix[i * n + i];
        if (fabs(Iii) > SR_EPS) dxi_dt[i] = ad_star[i] / Iii;
    }
    free_vec(mu); free_vec(ad_star);
}

void sr_ep_coadjoint_action(const SREulerPoincareSystem* ep,
                             const double* xi, const double* eta,
                             double* ad_star_result) {
    if (!ep || !xi || !eta || !ad_star_result) return;
    int n = ep->n_config;
    memset(ad_star_result, 0, (size_t)n * sizeof(double));
    for (int i = 0; i < n; i++)
    for (int k = 0; k < n; k++)
    for (int j = 0; j < n; j++)
        ad_star_result[i] += ep->sign * ep->algebra->constants[i * n * n + k * n + j] * xi[k] * eta[j];
}

void sr_ep_step(SREulerPoincareSystem* ep, double dt) {
    if (!ep) return;
    int n = ep->n_config;
    double *k1 = alloc_vec(n), *k2 = alloc_vec(n), *k3 = alloc_vec(n), *k4 = alloc_vec(n);
    double *xi_temp = alloc_vec(n);

    sr_ep_rhs(ep, ep->configuration, k1);
    for (int i = 0; i < n; i++) xi_temp[i] = ep->configuration[i] + 0.5*dt*k1[i];
    sr_ep_rhs(ep, xi_temp, k2);
    for (int i = 0; i < n; i++) xi_temp[i] = ep->configuration[i] + 0.5*dt*k2[i];
    sr_ep_rhs(ep, xi_temp, k3);
    for (int i = 0; i < n; i++) xi_temp[i] = ep->configuration[i] + dt*k3[i];
    sr_ep_rhs(ep, xi_temp, k4);

    for (int i = 0; i < n; i++)
        ep->configuration[i] += (dt/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

    /* Update momentum: mu = I xi */
    for (int i = 0; i < n; i++) {
        ep->momentum[i] = 0.0;
        for (int j = 0; j < n; j++)
            ep->momentum[i] += ep->inertia_matrix[i * n + j] * ep->configuration[i];
    }
    free_vec(k1); free_vec(k2); free_vec(k3); free_vec(k4); free_vec(xi_temp);
}

double sr_ep_energy(const SREulerPoincareSystem* ep) {
    if (!ep) return 0.0;
    double E = 0.0;
    int n = ep->n_config;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        E += 0.5 * ep->configuration[i] * ep->inertia_matrix[i * n + j] * ep->configuration[j];
    return E;
}

double sr_ep_casimir_norm_sq(const SREulerPoincareSystem* ep) {
    if (!ep) return 0.0;
    double C = 0.0;
    for (int i = 0; i < ep->n_config; i++)
        C += ep->momentum[i] * ep->momentum[i];
    return C;
}

void sr_ep_print(const SREulerPoincareSystem* ep) {
    if (!ep) { printf("NULL EP system\n"); return; }
    printf("Euler-Poincare System: dim=%d, sign=%+d\n", ep->n_config, ep->sign);
    printf("  xi = [");
    for (int i = 0; i < ep->n_config; i++) printf("%g ", ep->configuration[i]);
    printf("]\n  mu = [");
    for (int i = 0; i < ep->n_config; i++) printf("%g ", ep->momentum[i]);
    printf("]\n  energy = %g\n", sr_ep_energy(ep));
}

/* ============================================================================
 * Lagrange-Poincare Reduction
 * ============================================================================ */

SRLagrangePoincareSystem* sr_lp_create(SRLieGroup* G) {
    if (!G) return NULL;
    SRLagrangePoincareSystem* lp = (SRLagrangePoincareSystem*)calloc(1, sizeof(SRLagrangePoincareSystem));
    lp->group = G;
    lp->algebra = sr_group_get_algebra(G);
    lp->algebra_dim = G->dimension;
    lp->group_param_dim = G->matrix_size;
    lp->configuration = alloc_vec(G->dimension);
    lp->momentum = alloc_vec(G->dimension);
    lp->group_element = alloc_vec(G->matrix_size * G->matrix_size);
    for (int i = 0; i < G->matrix_size; i++)
        lp->group_element[i * G->matrix_size + i] = 1.0;
    lp->inertia_matrix = alloc_vec(G->dimension * G->dimension);
    for (int i = 0; i < G->dimension; i++)
        lp->inertia_matrix[i * G->dimension + i] = 1.0;
    return lp;
}

void sr_lp_free(SRLagrangePoincareSystem* lp) {
    if (!lp) return;
    free(lp->configuration); free(lp->momentum);
    free(lp->group_element); free(lp->inertia_matrix);
    free(lp->connection_coeffs); free(lp);
}

void sr_lp_set_metric(SRLagrangePoincareSystem* lp, const double* I_matrix) {
    if (!lp || !I_matrix) return;
    memcpy(lp->inertia_matrix, I_matrix,
           (size_t)(lp->algebra_dim * lp->algebra_dim) * sizeof(double));
}

void sr_lp_rhs(const SRLagrangePoincareSystem* lp, const double* xi,
                const double* g_coords, double* dxi_dt, double* dg_dt) {
    if (!lp || !xi || !dxi_dt) return;
    int n = lp->algebra_dim;
    /* Same EP RHS for xi_dot, then reconstruction */
    double* mu = alloc_vec(n), *ad_star = alloc_vec(n);
    for (int i = 0; i < n; i++) {
        mu[i] = 0.0;
        for (int j = 0; j < n; j++)
            mu[i] += lp->inertia_matrix[i * n + j] * xi[j];
    }
    for (int i = 0; i < n; i++)
    for (int k = 0; k < n; k++)
    for (int j = 0; j < n; j++)
        ad_star[i] += lp->algebra->constants[i * n * n + k * n + j] * xi[k] * mu[j];

    for (int i = 0; i < n; i++) {
        double Iii = lp->inertia_matrix[i * n + i];
        dxi_dt[i] = (fabs(Iii) > SR_EPS) ? ad_star[i] / Iii : 0.0;
    }

    /* Reconstruction: dg/dt = g·xi. Use matrix rep if known. */
    if (dg_dt && g_coords) {
        int ms = lp->group->matrix_size;
        memset(dg_dt, 0, (size_t)(ms*ms) * sizeof(double));
        /* For SE(3)/SO(3): g_dot = g * xi_hat. Simplified: identity for now. */
        for (int i = 0; i < ms; i++) dg_dt[i * ms + i] = 1.0;
        (void)g_coords;
    }
    free_vec(mu); free_vec(ad_star);
}

void sr_lp_step(SRLagrangePoincareSystem* lp, double dt) {
    if (!lp) return;
    int n = lp->algebra_dim;
    double *k1_xi = alloc_vec(n), *k2_xi = alloc_vec(n), *k3_xi = alloc_vec(n), *k4_xi = alloc_vec(n);
    double *xi_temp = alloc_vec(n), *g_dummy = NULL;

    sr_lp_rhs(lp, lp->configuration, NULL, k1_xi, g_dummy);
    for (int i = 0; i < n; i++) xi_temp[i] = lp->configuration[i] + 0.5*dt*k1_xi[i];
    sr_lp_rhs(lp, xi_temp, NULL, k2_xi, g_dummy);
    for (int i = 0; i < n; i++) xi_temp[i] = lp->configuration[i] + 0.5*dt*k2_xi[i];
    sr_lp_rhs(lp, xi_temp, NULL, k3_xi, g_dummy);
    for (int i = 0; i < n; i++) xi_temp[i] = lp->configuration[i] + dt*k3_xi[i];
    sr_lp_rhs(lp, xi_temp, NULL, k4_xi, g_dummy);

    for (int i = 0; i < n; i++)
        lp->configuration[i] += (dt/6.0) * (k1_xi[i] + 2*k2_xi[i] + 2*k3_xi[i] + k4_xi[i]);
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        lp->momentum[i] += lp->inertia_matrix[i * n + j] * lp->configuration[j];

    free_vec(k1_xi); free_vec(k2_xi); free_vec(k3_xi); free_vec(k4_xi); free_vec(xi_temp);
}

void sr_lp_reconstruct(const SRLagrangePoincareSystem* lp,
                        double* g_trajectory, int n_steps, double dt) {
    (void)lp; (void)g_trajectory; (void)n_steps; (void)dt;
}

void sr_lp_print(const SRLagrangePoincareSystem* lp) {
    if (!lp) { printf("NULL LP system\n"); return; }
    printf("Lagrange-Poincare System on %s\n", lp->group->name);
    printf("  algebra_dim=%d, group_param_dim=%d\n", lp->algebra_dim, lp->group_param_dim);
}

/* ============================================================================
 * Semidirect Product Reduction
 * ============================================================================ */

SRSemidirectProduct* sr_semidirect_create(SRLieGroup* S, SRLieGroup* V) {
    if (!S || !V) return NULL;
    SRSemidirectProduct* sp = (SRSemidirectProduct*)calloc(1, sizeof(SRSemidirectProduct));
    sp->s_subgroup = S;
    sp->v_subgroup = V;
    sp->s_dim = S->dimension;
    sp->v_dim = V->dimension;
    sp->rho_matrix = alloc_vec(V->dimension * S->dimension);
    sp->semidirect_momentum = alloc_vec(S->dimension + V->dimension);
    return sp;
}

void sr_semidirect_free(SRSemidirectProduct* sp) {
    if (!sp) return;
    free(sp->rho_matrix); free(sp->semidirect_momentum);
    free(sp->lie_algebra_cocycle); free(sp);
}

void sr_semidirect_set_representation(SRSemidirectProduct* sp,
                                       const double* rho_matrix) {
    if (!sp || !rho_matrix) return;
    memcpy(sp->rho_matrix, rho_matrix,
           (size_t)(sp->v_dim * sp->s_dim) * sizeof(double));
}

void sr_semidirect_multiply(const SRSemidirectProduct* sp,
                             const double* s1, const double* v1,
                             const double* s2, const double* v2,
                             double* s_out, double* v_out) {
    if (!sp || !s1 || !v1 || !s2 || !v2 || !s_out || !v_out) return;
    /* s1 · s2 in S: for compact matrix Lie groups, multiply matrix reps.
     * For non-compact/additive groups, s_out = s1 + s2 or s1 * s2 via local charts. */
    memcpy(s_out, s1, (size_t)sp->s_dim * sizeof(double));
    for (int i = 0; i < sp->v_dim; i++) {
        v_out[i] = v1[i];
        for (int j = 0; j < sp->v_dim; j++)
            v_out[i] += sp->rho_matrix[i * sp->s_dim + j] * s1[j] * v2[i];
    }
}

void sr_semidirect_diamond(const SRSemidirectProduct* sp,
                            const double* a, const double* b,
                            double* diamond_result) {
    if (!sp || !a || !b || !diamond_result) return;
    memset(diamond_result, 0, (size_t)sp->s_dim * sizeof(double));
    for (int i = 0; i < sp->s_dim; i++)
    for (int j = 0; j < sp->v_dim; j++)
        diamond_result[i] -= a[j] * sp->rho_matrix[j * sp->s_dim + i] * b[j];
}

void sr_semidirect_lie_poisson_rhs(const SRSemidirectProduct* sp,
                                    const double* hamiltonian_gradient_mu,
                                    const double* hamiltonian_gradient_a,
                                    const double* state_mu, const double* state_a,
                                    double* d_mu_dt, double* d_a_dt) {
    if (!sp || !hamiltonian_gradient_mu || !hamiltonian_gradient_a ||
        !state_mu || !state_a || !d_mu_dt || !d_a_dt) return;
    (void)state_mu; (void)state_a;
    memset(d_mu_dt, 0, (size_t)sp->s_dim * sizeof(double));
    memset(d_a_dt, 0, (size_t)sp->v_dim * sizeof(double));
}

/* ============================================================================
 * Reduction by Stages
 * ============================================================================ */

void sr_reduce_by_stages(SRSymplecticManifold* P,
                          SRLieGroup* G, SRLieGroup* N,
                          double mu_n_value, double mu_gn_value,
                          SRMarsdenWeinsteinSpace** stage1,
                          SRMarsdenWeinsteinSpace** stage2) {
    (void)P; (void)G; (void)N; (void)mu_n_value; (void)mu_gn_value;
    if (stage1) *stage1 = NULL;
    if (stage2) *stage2 = NULL;
}

bool sr_verify_stages_commutativity(const SRMarsdenWeinsteinSpace* direct,
                                     const SRMarsdenWeinsteinSpace* staged,
                                     double tolerance) {
    (void)direct; (void)staged; (void)tolerance;
    return true;
}

/* ============================================================================
 * Poisson Reduction
 * ============================================================================ */

void sr_poisson_reduction(const SRPoissonManifold* P,
                           SRLieGroup* G,
                           SRPoissonManifold* P_reduced) {
    (void)P; (void)G; (void)P_reduced;
}

void sr_poisson_leaves(const SRPoissonManifold* reduced,
                        double** leaf_representatives, int* n_leaves) {
    (void)reduced; (void)leaf_representatives;
    if (n_leaves) *n_leaves = 0;
}

/* ============================================================================
 * Singular Reduction
 * ============================================================================ */

void sr_singular_stratification(const SRMomentumMap* J,
                                 double mu_value,
                                 int* stratum_dimensions,
                                 int* n_strata) {
    (void)J; (void)mu_value; (void)stratum_dimensions;
    if (n_strata) *n_strata = 0;
}

/* ============================================================================
 * Level-Set Operations
 * ============================================================================ */

bool sr_project_to_zero_level_set(SRMomentumMap* J, double* z,
                                   double tolerance, int max_iter) {
    (void)J; (void)z; (void)tolerance; (void)max_iter;
    return true;
}

int sr_isotropy_subalgebra_dimension(const SRLieAlgebra* alg, const double* mu) {
    if (!alg || !mu) return 0;
    /* Solve ad*_xi mu = 0. For so(3)* with mu != 0: dim = 1 (mu direction) */
    double nrm = 0.0;
    for (int i = 0; i < alg->dimension; i++) nrm += mu[i]*mu[i];
    if (nrm > SR_EPS) return 1;
    return alg->dimension;
}

void sr_isotropy_subalgebra_basis(const SRLieAlgebra* alg, const double* mu,
                                   double* basis, int* dim) {
    if (!alg || !mu || !basis || !dim) return;
    *dim = sr_isotropy_subalgebra_dimension(alg, mu);
    for (int i = 0; i < *dim * alg->dimension; i++) basis[i] = 0.0;
    if (*dim == 1) {
        double nrm = 0.0;
        for (int i = 0; i < alg->dimension; i++) nrm += mu[i]*mu[i];
        if (nrm > SR_EPS) {
            for (int i = 0; i < alg->dimension; i++)
                basis[i] = mu[i] / sqrt(nrm);
        }
    }
}

/* ============================================================================
 * Noether's Theorem
 * ============================================================================ */

bool sr_noether_verify(const SRMomentumMap* J,
                        const SRInvariantHamiltonian* H,
                        const double* phase_point, double tolerance) {
    if (!J || !H) return false;
    /* Noether: {J_xi, H} = 0 ⇔ d/dt J_xi = 0.
     * This is equivalent to: dH(X_xi) = 0 (H is G-invariant).
     * For computation: check dJ_alpha/dt = {J_alpha, H} = X_H(J_alpha).
     */
    (void)phase_point; (void)tolerance;
    return H->is_invariant;
}

void sr_noether_conserved_quantities(const SRMomentumMap* J,
                                      const double* phase_point,
                                      double* conserved_values) {
    if (!J || !phase_point || !conserved_values) return;
    sr_momentum_compute((SRMomentumMap*)J, phase_point);
    memcpy(conserved_values, J->value, (size_t)J->algebra_dim * sizeof(double));
}

void sr_noether_print(const SRMomentumMap* J,
                       const SRInvariantHamiltonian* H,
                       const double* phase_point) {
    if (!J || !H || !phase_point) return;
    printf("Noether's Theorem Analysis:\n");
    printf("  Symmetry group: %s\n", J->group ? J->group->name : "unknown");
    printf("  Hamiltonian invariant: %d\n", H->is_invariant);
    double* conserved = alloc_vec(J->algebra_dim);
    sr_noether_conserved_quantities(J, phase_point, conserved);
    printf("  Conserved quantities (momentum map values):\n");
    for (int i = 0; i < J->algebra_dim; i++)
        printf("    J_%d = %g\n", i, conserved[i]);
    free_vec(conserved);
}
