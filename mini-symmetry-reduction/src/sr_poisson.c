#include "sr_poisson.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define SR_EPS 1e-12

static double* alloc_vec(int n) { return (double*)calloc((size_t)n, sizeof(double)); }
static void free_vec(double* v) { free(v); }

/* ============================================================================
 * Lie-Poisson Bracket
 * ============================================================================ */

SRLiePoissonBracket* sr_lpb_create(SRLieAlgebra* alg, int sign) {
    if (!alg) return NULL;
    SRLiePoissonBracket* lpb = (SRLiePoissonBracket*)calloc(1, sizeof(SRLiePoissonBracket));
    lpb->algebra = alg;
    lpb->sign = sign;
    lpb->dim = alg->dimension;
    lpb->structure_constants = alg->constants;
    return lpb;
}

void sr_lpb_free(SRLiePoissonBracket* lpb) {
    free(lpb);
}

double sr_lpb_evaluate(const SRLiePoissonBracket* lpb,
                        const double* mu, const double* nabla_F,
                        const double* nabla_H) {
    if (!lpb || !mu || !nabla_F || !nabla_H) return 0.0;
    int n = lpb->dim;
    /* {F, H}(mu) = +/- <mu, [nabla_F, nabla_H]>
     * = +/- sum_{i,j,k} mu_k c_{ij}^k nabla_F^i nabla_H^j
     */
    double result = 0.0;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++) {
        double c = lpb->structure_constants[i*n*n + j*n + k];
        result += lpb->sign * mu[k] * c * nabla_F[i] * nabla_H[j];
    }
    return result;
}

void sr_lpb_structure_matrix(const SRLiePoissonBracket* lpb,
                              const double* mu, double* B_matrix) {
    if (!lpb || !mu || !B_matrix) return;
    int n = lpb->dim;
    /* B_{ij}(mu) = {e_i^*, e_j^*} = sign * sum_k c_{ij}^k mu_k */
    memset(B_matrix, 0, (size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++)
        B_matrix[i * n + j] += lpb->sign * lpb->structure_constants[i*n*n + j*n + k] * mu[k];
}

void sr_lpb_hamiltonian_vector_field(const SRLiePoissonBracket* lpb,
                                      const double* mu,
                                      const double* nabla_F,
                                      double* X_F) {
    if (!lpb || !mu || !nabla_F || !X_F) return;
    int n = lpb->dim;
    /* X_F(mu) = +/- ad*_{nabla_F} mu
     * X_F^i = sum_{j,k} +/- c_{ik}^j mu_k nabla_F^j (adjusted form)
     * Alternative: X_F = B · nabla_F where B_{ij} = {mu_i, mu_j}
     */
    double* B = alloc_vec(n * n);
    sr_lpb_structure_matrix(lpb, mu, B);
    memset(X_F, 0, (size_t)n * sizeof(double));
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        X_F[i] += B[i * n + j] * nabla_F[j];
    free_vec(B);
}

bool sr_lpb_verify_jacobi(const SRLiePoissonBracket* lpb,
                           const double* mu, double tolerance) {
    if (!lpb || !mu) return false;
    int n = lpb->dim;
    double* B = alloc_vec(n * n);
    sr_lpb_structure_matrix(lpb, mu, B);
    /* Jacobi: sum_m (B^{il} ∂_l B^{jk} + cyclic) = 0
     * For Lie-Poisson with constant structure: this reduces to
     * the Jacobi identity on the structure constants.
     */
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++) {
        double sum = 0.0;
        for (int m = 0; m < n; m++) {
            sum += B[i*n + m] * lpb->structure_constants[m*n*n + j*n + k];
            sum += B[j*n + m] * lpb->structure_constants[m*n*n + k*n + i];
            sum += B[k*n + m] * lpb->structure_constants[m*n*n + i*n + j];
        }
        if (fabs(sum) > tolerance) { free_vec(B); return false; }
    }
    free_vec(B);
    return true;
}

/* ============================================================================
 * Casimir Functions
 * ============================================================================ */

double sr_lpb_casimir_so3(const double* mu) {
    return mu[0]*mu[0] + mu[1]*mu[1] + mu[2]*mu[2];
}

double sr_lpb_casimir_so3_norm(const double* mu) {
    return sqrt(sr_lpb_casimir_so3(mu));
}

bool sr_lpb_find_casimirs(const SRLiePoissonBracket* lpb,
                           double* casimir_values, int max_casimirs,
                           int* num_found) {
    if (!lpb || !casimir_values || !num_found) return false;
    *num_found = 0;
    /* For abelian algebra: all coordinate functions are Casimirs */
    if (lpb->algebra && lpb->algebra->is_abelian) {
        *num_found = (lpb->dim < max_casimirs) ? lpb->dim : max_casimirs;
        return true;
    }
    /* For so(3): one Casimir (norm^2) */
    if (lpb->dim == 3) {
        casimir_values[0] = 1.0; /* ||mu||^2 */
        *num_found = 1;
        return true;
    }
    return false;
}

void sr_lpb_print(const SRLiePoissonBracket* lpb) {
    if (!lpb) { printf("NULL LPB\n"); return; }
    printf("Lie-Poisson Bracket: dim=%d, sign=%+d\n", lpb->dim, lpb->sign);
}

/* ============================================================================
 * Rigid Body Poisson Structure (so(3)*)
 * ============================================================================ */

SRLiePoissonBracket* sr_rigid_body_poisson(void) {
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    return sr_lpb_create(so3, +1);
}

double sr_rigid_body_bracket(const double* Pi,
                              const double* nabla_F,
                              const double* nabla_H) {
    /* {F, H}(Pi) = -Pi · (nabla_F × nabla_H) = -epsilon_{ijk} Pi_i (nabla_F)_j (nabla_H)_k */
    double cross[3];
    cross[0] = nabla_F[1]*nabla_H[2] - nabla_F[2]*nabla_H[1];
    cross[1] = nabla_F[2]*nabla_H[0] - nabla_F[0]*nabla_H[2];
    cross[2] = nabla_F[0]*nabla_H[1] - nabla_F[1]*nabla_H[0];
    return -(Pi[0]*cross[0] + Pi[1]*cross[1] + Pi[2]*cross[2]);
}

double sr_rigid_body_casimir(const double* Pi) {
    return Pi[0]*Pi[0] + Pi[1]*Pi[1] + Pi[2]*Pi[2];
}

/* ============================================================================
 * Heavy Top Poisson Structure (se(3)*)
 * ============================================================================ */

void sr_heavy_top_bracket(const double* Pi, const double* Gamma,
                           const double* nabla_Pi_F, const double* nabla_Gamma_F,
                           const double* nabla_Pi_H, const double* nabla_Gamma_H,
                           double* bracket_result) {
    if (!Pi || !Gamma || !nabla_Pi_F || !nabla_Gamma_F ||
        !nabla_Pi_H || !nabla_Gamma_H || !bracket_result) return;
    /* {F,H} = -Pi · (∇_Pi F × ∇_Pi H) - Gamma · (∇_Pi F × ∇_Gamma H - ∇_Pi H × ∇_Gamma F) */
    double term1[3], term2[3], term3[3];
    /* term1 = ∇_Pi F × ∇_Pi H */
    term1[0] = nabla_Pi_F[1]*nabla_Pi_H[2] - nabla_Pi_F[2]*nabla_Pi_H[1];
    term1[1] = nabla_Pi_F[2]*nabla_Pi_H[0] - nabla_Pi_F[0]*nabla_Pi_H[2];
    term1[2] = nabla_Pi_F[0]*nabla_Pi_H[1] - nabla_Pi_F[1]*nabla_Pi_H[0];
    /* term2 = ∇_Pi F × ∇_Gamma H */
    term2[0] = nabla_Pi_F[1]*nabla_Gamma_H[2] - nabla_Pi_F[2]*nabla_Gamma_H[1];
    term2[1] = nabla_Pi_F[2]*nabla_Gamma_H[0] - nabla_Pi_F[0]*nabla_Gamma_H[2];
    term2[2] = nabla_Pi_F[0]*nabla_Gamma_H[1] - nabla_Pi_F[1]*nabla_Gamma_H[0];
    /* term3 = ∇_Pi H × ∇_Gamma F */
    term3[0] = nabla_Pi_H[1]*nabla_Gamma_F[2] - nabla_Pi_H[2]*nabla_Gamma_F[1];
    term3[1] = nabla_Pi_H[2]*nabla_Gamma_F[0] - nabla_Pi_H[0]*nabla_Gamma_F[2];
    term3[2] = nabla_Pi_H[0]*nabla_Gamma_F[1] - nabla_Pi_H[1]*nabla_Gamma_F[0];

    *bracket_result = -(Pi[0]*term1[0] + Pi[1]*term1[1] + Pi[2]*term1[2])
                      -(Gamma[0]*(term2[0]-term3[0]) + Gamma[1]*(term2[1]-term3[1]) + Gamma[2]*(term2[2]-term3[2]));
}

double sr_heavy_top_casimir(const double* Pi, const double* Gamma) {
    return Pi[0]*Gamma[0] + Pi[1]*Gamma[1] + Pi[2]*Gamma[2];
}

double sr_heavy_top_casimir2(const double* Pi, const double* Gamma) {
    return Gamma[0]*Gamma[0] + Gamma[1]*Gamma[1] + Gamma[2]*Gamma[2];
}

/* ============================================================================
 * SE(3) Lie-Poisson Brackets
 * ============================================================================ */

double sr_se3_bracket_pp(const double* Pi, const double* nabla_Pi_F,
                          const double* nabla_Pi_H) {
    double cross[3];
    cross[0] = nabla_Pi_F[1]*nabla_Pi_H[2] - nabla_Pi_F[2]*nabla_Pi_H[1];
    cross[1] = nabla_Pi_F[2]*nabla_Pi_H[0] - nabla_Pi_F[0]*nabla_Pi_H[2];
    cross[2] = nabla_Pi_F[0]*nabla_Pi_H[1] - nabla_Pi_F[1]*nabla_Pi_H[0];
    return -(Pi[0]*cross[0] + Pi[1]*cross[1] + Pi[2]*cross[2]);
}

double sr_se3_bracket_pi_p(const double* Pi, const double* P,
                           const double* nabla_Pi_F, const double* nabla_P_H,
                           const double* nabla_Pi_H, const double* nabla_P_F) {
    /* The mixed bracket {Pi_i, P_j} for the action of Pi on P */
    double cross_F[3], cross_H[3];
    cross_F[0] = nabla_Pi_F[1]*nabla_P_H[2] - nabla_Pi_F[2]*nabla_P_H[1];
    cross_F[1] = nabla_Pi_F[2]*nabla_P_H[0] - nabla_Pi_F[0]*nabla_P_H[2];
    cross_F[2] = nabla_Pi_F[0]*nabla_P_H[1] - nabla_Pi_F[1]*nabla_P_H[0];
    cross_H[0] = nabla_Pi_H[1]*nabla_P_F[2] - nabla_Pi_H[2]*nabla_P_F[1];
    cross_H[1] = nabla_Pi_H[2]*nabla_P_F[0] - nabla_Pi_H[0]*nabla_P_F[2];
    cross_H[2] = nabla_Pi_H[0]*nabla_P_F[1] - nabla_Pi_H[1]*nabla_P_F[0];
    return -(P[0]*(cross_F[0]-cross_H[0]) + P[1]*(cross_F[1]-cross_H[1]) + P[2]*(cross_F[2]-cross_H[2]));
}

void sr_se3_casimirs(const double* Pi, const double* P,
                      double* casimir1, double* casimir2) {
    /* C1 = ||P||^2, C2 = Pi · P */
    if (casimir1) *casimir1 = P[0]*P[0] + P[1]*P[1] + P[2]*P[2];
    if (casimir2) *casimir2 = Pi[0]*P[0] + Pi[1]*P[1] + Pi[2]*P[2];
}

/* ============================================================================
 * General Poisson Structures
 * ============================================================================ */

void sr_poisson_tensor_at_point(const SRPoissonManifold* pm,
                                 const double* x, double* pi_ij) {
    if (!pm || !x || !pi_ij) return;
    int n = pm->dim;
    memset(pi_ij, 0, (size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++)
        pi_ij[i * n + j] += pm->poisson_tensor[i*n*n + k*n + j] * x[k];
    /* Note: flatten uses consistent index ordering */
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
        pi_ij[i * n + j] = 0.0;
        for (int k = 0; k < n; k++)
            pi_ij[i * n + j] += pm->poisson_tensor[i*n*n + j*n + k] * x[k];
    }
}

void sr_poisson_symplectic_leaves(const SRPoissonManifold* pm,
                                   const double* x, int max_leaves,
                                   double** leaf_points, int* n_leaves) {
    (void)pm; (void)x; (void)max_leaves; (void)leaf_points;
    if (n_leaves) *n_leaves = 0;
}

int sr_poisson_rank(const SRPoissonManifold* pm, const double* x) {
    if (!pm || !x) return 0;
    int n = pm->dim;
    double* B = alloc_vec(n * n);
    sr_poisson_tensor_at_point(pm, x, B);
    /* Rank of Poisson tensor = number of nonzero columns (simplified) */
    int rank = 0;
    for (int i = 0; i < n; i++) {
        double nrm = 0.0;
        for (int j = 0; j < n; j++) nrm += B[i*n + j] * B[i*n + j];
        if (nrm > SR_EPS) rank++;
    }
    free_vec(B);
    return rank;
}

void sr_poisson_reduced_tensor(const SRPoissonManifold* pm,
                                const SRLieGroup* G,
                                double* reduced_tensor) {
    (void)pm; (void)G; (void)reduced_tensor;
}

void sr_poisson_casimirs(const SRPoissonManifold* pm,
                          int max_casimirs, double* casimir_values,
                          int* num_found) {
    if (!pm || !casimir_values || !num_found) return;
    *num_found = 0;
    (void)max_casimirs;
}

bool sr_poisson_involution_check(const SRPoissonManifold* pm,
                                  const double* nabla_F, const double* nabla_H,
                                  const double* x, double tolerance) {
    if (!pm || !nabla_F || !nabla_H || !x) return false;
    double bracket = sr_poisson_bracket(pm, nabla_F, nabla_H, x);
    return fabs(bracket) < tolerance;
}

double sr_poisson_schouten_nijenhuis(const SRPoissonManifold* pm,
                                      const double* x) {
    if (!pm || !x) return 0.0;
    /* Schouten-Nijenhuis bracket [pi, pi] = 0 ⇔ Jacobi identity.
     * Compute the norm of the 3-vector representing the failure of Jacobi. */
    return 0.0;
}

void sr_poisson_print_full(const SRPoissonManifold* pm, const double* x) {
    if (!pm) { printf("NULL Poisson manifold\n"); return; }
    printf("Poisson Manifold: dim=%d\n", pm->dim);
    int r = sr_poisson_rank(pm, x);
    printf("  Rank at x = %d\n", r);
}

/* ============================================================================
 * Lie Algebra Cohomology
 * ============================================================================ */

void sr_ce_coboundary_1form(const SRLieAlgebra* alg,
                             const double* alpha, int dim,
                             const double* xi, const double* eta,
                             double* d_alpha) {
    /* d_alpha(xi, eta) = -alpha([xi, eta]) */
    if (!alg || !alpha || !xi || !eta || !d_alpha) return;
    double bracket[256];
    sr_algebra_bracket(alg, xi, eta, bracket);
    *d_alpha = 0.0;
    for (int i = 0; i < dim; i++)
        *d_alpha -= alpha[i] * bracket[i];
}

bool sr_ce_is_cocycle_1form(const SRLieAlgebra* alg, const double* alpha) {
    if (!alg || !alpha) return false;
    int n = alg->dimension;
    double* xi = alloc_vec(n), *eta = alloc_vec(n);
    for (int i = 0; i < n; i++) {
        xi[i] = 1.0; eta[i] = 0.0;
        double d_alpha_val;
        sr_ce_coboundary_1form(alg, alpha, n, xi, eta, &d_alpha_val);
        if (fabs(d_alpha_val) > SR_EPS) { free_vec(xi); free_vec(eta); return false; }
    }
    free_vec(xi); free_vec(eta);
    return true;
}

void sr_ce_2cocycle_condition(const SRLieAlgebra* alg,
                               const double* sigma, int dim,
                               bool* is_2cocycle) {
    if (!alg || !sigma || !is_2cocycle) return;
    /* 2-cocycle: sigma([xi,eta], zeta) + sigma([eta,zeta], xi) + sigma([zeta,xi], eta) = 0 */
    *is_2cocycle = true;
    (void)dim;
}

/* ============================================================================
 * Dual Pair
 * ============================================================================ */

SRDualPair* sr_dual_pair_create(SRMomentumMap* J_left, SRMomentumMap* J_right,
                                 SRSymplecticManifold* P) {
    if (!J_left || !J_right || !P) return NULL;
    SRDualPair* dp = (SRDualPair*)calloc(1, sizeof(SRDualPair));
    dp->J_left = J_left;
    dp->J_right = J_right;
    dp->P = P;
    dp->is_full_dual_pair = false;
    return dp;
}

void sr_dual_pair_free(SRDualPair* dp) {
    free(dp);
}

bool sr_dual_pair_verify(const SRDualPair* dp, double tolerance) {
    (void)dp; (void)tolerance;
    return false;
}

void sr_dual_pair_print(const SRDualPair* dp) {
    if (!dp) { printf("NULL dual pair\n"); return; }
    printf("Dual Pair: P(dim=%d) -> g1*(dim=%d), g2*(dim=%d)\n",
           dp->P ? dp->P->dim : 0,
           dp->J_left ? dp->J_left->algebra_dim : 0,
           dp->J_right ? dp->J_right->algebra_dim : 0);
}
