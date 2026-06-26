#ifndef SR_POISSON_H
#define SR_POISSON_H

#include "sr_core.h"
#include <stdbool.h>

/* ============================================================================
 * Lie-Poisson and Poisson Structures
 *
 * The Lie-Poisson bracket on g* is the fundamental Poisson structure
 * arising from Lie group symmetries:
 *   {F, H}(mu) = +/- <mu, [dF/d_mu, dH/d_mu]>
 *
 * Symplectic leaves = coadjoint orbits O_mu = G/G_mu.
 * Each leaf has the KKS symplectic form: omega_mu(xi, eta) = <mu, [xi, eta]>.
 *
 * References:
 *   Kirillov (1962) "Unitary representations of nilpotent Lie groups"
 *   Kostant (1970) "Quantization and unitary representations"
 *   Souriau (1970) "Structure des systemes dynamiques"
 *   Marsden & Ratiu (1999) "Mechanics and Symmetry", Ch. 10
 * ============================================================================ */

/* === Lie-Poisson Bracket === */

/* Create a Lie-Poisson bracket on g* */
SRLiePoissonBracket* sr_lpb_create(SRLieAlgebra* alg, int sign);

void sr_lpb_free(SRLiePoissonBracket* lpb);

/* Evaluate Lie-Poisson bracket: {F, H}(mu) = +/- <mu, [nabla_F, nabla_H]>
 * where nabla_F = dF/d_mu in g (identified with g** via natural pairing).
 *
 * The bracket satisfies:
 *   (a) Bilinearity: {aF + bG, H} = a{F, H} + b{G, H}
 *   (b) Antisymmetry: {F, H} = -{H, F}
 *   (c) Jacobi identity: {{F, G}, H} + {{G, H}, F} + {{H, F}, G} = 0
 *   (d) Leibniz: {F·G, H} = F·{G, H} + {F, H}·G
 */
double sr_lpb_evaluate(const SRLiePoissonBracket* lpb,
                        const double* mu, const double* nabla_F,
                        const double* nabla_H);

/* Compute bracket values for all basis function pairs:
 *   B_{ij}(mu) = {e_i^*, e_j^*}(mu)
 * where e_i^* is the i-th coordinate function on g*.
 */
void sr_lpb_structure_matrix(const SRLiePoissonBracket* lpb,
                              const double* mu, double* B_matrix);

/* Compute the Hamiltonian vector field X_F on g*:
 * X_F[mu] = +/- ad*_{nabla_F(mu)} mu
 *
 * The flow of X_F preserves coadjoint orbits.
 */
void sr_lpb_hamiltonian_vector_field(const SRLiePoissonBracket* lpb,
                                      const double* mu,
                                      const double* nabla_F,
                                      double* X_F);

/* Verify Jacobi identity numerically at a point mu */
bool sr_lpb_verify_jacobi(const SRLiePoissonBracket* lpb,
                           const double* mu, double tolerance);

/* Compute Casimir functions: C such that {C, F} = 0 for all F.
 * Casimirs are constant on coadjoint orbits.
 * For so(3)*: C(mu) = mu_1^2 + mu_2^2 + mu_3^2
 */
double sr_lpb_casimir_so3(const double* mu);
double sr_lpb_casimir_so3_norm(const double* mu);

/* Generic Casimir: solve ad*_{dC(mu)} mu = 0 */
bool sr_lpb_find_casimirs(const SRLiePoissonBracket* lpb,
                           double* casimir_values, int max_casimirs,
                           int* num_found);

void sr_lpb_print(const SRLiePoissonBracket* lpb);

/* === Rigid Body Poisson Structure (so(3)*) === */

/* Create the so(3)* Lie-Poisson bracket (right invariant, sign = +1).
 * The bracket for the free rigid body:
 *   {Pi_i, Pi_j} = -epsilon_{ijk} Pi_k
 * where Pi is the body angular momentum.
 *
 * Euler equations: dPi/dt = {Pi, H} = Pi x omega
 * with omega = I^{-1} Pi (body angular velocity).
 */
SRLiePoissonBracket* sr_rigid_body_poisson(void);

/* Compute Poisson bracket of two functions on so(3)* */
double sr_rigid_body_bracket(const double* Pi,
                              const double* nabla_F,
                              const double* nabla_H);

/* Compute so(3) Casimir: C(Pi) = Pi_1^2 + Pi_2^2 + Pi_3^2 */
double sr_rigid_body_casimir(const double* Pi);

/* === Heavy Top Poisson Structure === */

/* The heavy top lives on se(3)* = so(3)* ⋉ R^3* (semidirect product coadjoint orbits).
 * Coordinates: (Pi, Gamma) where Pi = angular momentum, Gamma = gravity direction.
 *
 * Lie-Poisson bracket:
 *   {Pi_i, Pi_j} = -epsilon_{ijk} Pi_k
 *   {Pi_i, Gamma_j} = -epsilon_{ijk} Gamma_k
 *   {Gamma_i, Gamma_j} = 0
 *
 * Hamiltonian: H(Pi, Gamma) = 1/2 Pi·I^{-1}·Pi + mgl Gamma·chi
 * where chi is the unit vector from fixed point to center of mass.
 */
void sr_heavy_top_bracket(const double* Pi, const double* Gamma,
                           const double* nabla_Pi_F, const double* nabla_Gamma_F,
                           const double* nabla_Pi_H, const double* nabla_Gamma_H,
                           double* bracket_result);

double sr_heavy_top_casimir(const double* Pi, const double* Gamma);
double sr_heavy_top_casimir2(const double* Pi, const double* Gamma);

/* === SE(3) Lie-Poisson Brackets === */

/* SE(3)* = so(3)* ⋉ R^3*. Coordinates: (Pi, P) where Pi = angular, P = linear momentum.
 * Brackets:
 *   {Pi_i, Pi_j} = -epsilon_{ijk} Pi_k
 *   {Pi_i, P_j} = -epsilon_{ijk} P_k
 *   {P_i, P_j} = 0
 *
 * Casimirs: C1 = P·P, C2 = Pi·P
 */
double sr_se3_bracket_pp(const double* Pi, const double* nabla_Pi_F,
                          const double* nabla_Pi_H);
double sr_se3_bracket_pi_p(const double* Pi, const double* P,
                           const double* nabla_Pi_F, const double* nabla_P_H,
                           const double* nabla_Pi_H, const double* nabla_P_F);

void sr_se3_casimirs(const double* Pi, const double* P,
                      double* casimir1, double* casimir2);

/* === General Poisson Structures === */

/* Compute the Poisson tensor at a point x */
void sr_poisson_tensor_at_point(const SRPoissonManifold* pm,
                                 const double* x, double* pi_ij);

/* Compute all symplectic leaves near a point (coadjoint orbits) */
void sr_poisson_symplectic_leaves(const SRPoissonManifold* pm,
                                   const double* x, int max_leaves,
                                   double** leaf_points, int* n_leaves);

/* Rank of Poisson tensor (even, equals dimension of symplectic leaf) */
int sr_poisson_rank(const SRPoissonManifold* pm, const double* x);

/* Reduced Poisson tensor via Poisson reduction */
void sr_poisson_reduced_tensor(const SRPoissonManifold* pm,
                                const SRLieGroup* G,
                                double* reduced_tensor);

/* Compute Casimir functions for a general Poisson manifold */
void sr_poisson_casimirs(const SRPoissonManifold* pm,
                          int max_casimirs, double* casimir_values,
                          int* num_found);

/* Check if two functions are in involution: {F, H} = 0 */
bool sr_poisson_involution_check(const SRPoissonManifold* pm,
                                  const double* nabla_F, const double* nabla_H,
                                  const double* x, double tolerance);

/* Compute Nijenhuis tensor-like obstruction to integrability */
double sr_poisson_schouten_nijenhuis(const SRPoissonManifold* pm,
                                      const double* x);

void sr_poisson_print_full(const SRPoissonManifold* pm, const double* x);

/* === Lie Algebra Cohomology (for reduction context) === */

/* Compute the Chevalley-Eilenberg coboundary operator:
 * d : Alt^k(g, R) -> Alt^{k+1}(g, R)
 *
 * For k=1 (cocycles), d_alpha(xi, eta) = -alpha([xi, eta])
 * Relevant for: equivariant momentum maps, central extensions
 */
void sr_ce_coboundary_1form(const SRLieAlgebra* alg,
                             const double* alpha, int dim,
                             const double* xi, const double* eta,
                             double* d_alpha);

/* Check if alpha is a 1-cocycle: d_alpha = 0 */
bool sr_ce_is_cocycle_1form(const SRLieAlgebra* alg, const double* alpha);

/* Compute the Lie algebra cohomology H^2(g, R) for central extensions.
 * A 2-cocycle sigma satisfies:
 *   sigma([xi,eta], zeta) + sigma([eta,zeta], xi) + sigma([zeta,xi], eta) = 0
 * and is used to construct affine/extended Lie algebras.
 */
void sr_ce_2cocycle_condition(const SRLieAlgebra* alg,
                               const double* sigma, int dim,
                               bool* is_2cocycle);

/* === Dual Pair and Momentum Map Decomposition === */

/* A dual pair is a diagram:
 *       P
 *    J1/ \J2
 *    g1*  g2*
 * where the fibers of J1 and J2 are symplectic orthogonals.
 * This decomposes the symmetry.
 */
typedef struct {
    SRMomentumMap* J_left;
    SRMomentumMap* J_right;
    SRSymplecticManifold* P;
    bool is_full_dual_pair;
} SRDualPair;

SRDualPair* sr_dual_pair_create(SRMomentumMap* J_left, SRMomentumMap* J_right,
                                 SRSymplecticManifold* P);
void sr_dual_pair_free(SRDualPair* dp);
bool sr_dual_pair_verify(const SRDualPair* dp, double tolerance);
void sr_dual_pair_print(const SRDualPair* dp);

#endif /* SR_POISSON_H */
