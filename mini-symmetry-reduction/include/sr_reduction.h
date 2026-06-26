#ifndef SR_REDUCTION_H
#define SR_REDUCTION_H

#include "sr_core.h"
#include <stdbool.h>

/* ============================================================================
 * Symmetry Reduction Operations
 *
 * Implements the core reduction theorems:
 *   Marsden-Weinstein symplectic reduction
 *   Euler-Poincare reduction
 *   Lagrange-Poincare reduction
 *   Semidirect product reduction
 *   Reduction by stages
 *
 * References:
 *   Marsden & Weinstein (1974), Rep. Math. Phys. 5, 121-130
 *   Marsden, Misiolek, Ortega, Perlmutter, Ratiu (2007)
 *     "Hamiltonian Reduction by Stages", Springer LNM 1913
 * ============================================================================ */

/* === Momentum Map Operations === */

/* Create a momentum map for a given group action on phase space */
SRMomentumMap* sr_momentum_create(SRLieGroup* group, int phase_dim,
                                   SRMomentumType type);

/* Free momentum map resources */
void sr_momentum_free(SRMomentumMap* mm);

/* Compute J(z) for current phase point z */
void sr_momentum_compute(SRMomentumMap* mm, const double* phase_point);

/* Compute Jacobian of momentum map: dJ_i/dz_j */
void sr_momentum_jacobian(SRMomentumMap* mm, const double* phase_point);

/* Check Ad*-equivariance: J(g·z) = Ad*_{g^{-1}} J(z) */
bool sr_momentum_verify_equivariance(const SRMomentumMap* mm,
                                      const double* phase_point,
                                      const SRGroupElement* g, double tolerance);

/* Angular momentum map for SO(3) action on T*R^3 */
SRMomentumMap* sr_momentum_angular_so3(void);

/* Linear momentum map for R^3 translation */
SRMomentumMap* sr_momentum_linear_r3(void);

/* Kinetic momentum map: cotangent lift of left translation J_L(g,p) = T_e*L_g p */
SRMomentumMap* sr_momentum_kinetic_cotangent(SRLieGroup* group);

/* Compute the Hamiltonian vector field X_H = omega^{-1} dH */
void sr_momentum_hamiltonian_vector_field(const SRSymplecticManifold* sp,
                                           const double* dH, double* X_H);

/* === Marsden-Weinstein Reduction === */

/* Perform MW symplectic reduction at momentum value mu.
 * Theorem (Marsden-Weinstein 1974):
 *   Let (P,omega) be a symplectic manifold with free, proper Hamiltonian
 *   G-action and equivariant momentum map J. Then for each regular value
 *   mu of J, the reduced space P_mu = J^{-1}(mu)/G_mu is a symplectic
 *   manifold with a unique symplectic form omega_mu satisfying:
 *     pi_mu^* omega_mu = i_mu^* omega
 *   where pi_mu: J^{-1}(mu) -> P_mu and i_mu: J^{-1}(mu) -> P.
 */
SRMarsdenWeinsteinSpace* sr_mw_reduce(SRSymplecticManifold* P,
                                       SRMomentumMap* J,
                                       double mu_value);

void sr_mw_free(SRMarsdenWeinsteinSpace* reduced);

/* Compute the reduced symplectic form omega_mu by:
 * 1. Find a basis for T_z J^{-1}(mu) = ker(dJ(z))
 * 2. Remove gauge directions G_mu·z
 * 3. Restrict omega to the quotient
 */
void sr_mw_compute_reduced_omega(SRMarsdenWeinsteinSpace* reduced,
                                  const SRSymplecticManifold* P,
                                  const SRMomentumMap* J);

/* Project a tangent vector from original space to reduced space */
void sr_mw_project_tangent(const SRMarsdenWeinsteinSpace* reduced,
                            const double* v_original, double* v_reduced);

/* Lift a vector from reduced space back to original level set (connection-dependent) */
void sr_mw_lift_vector(const SRMarsdenWeinsteinSpace* reduced,
                        const double* v_reduced, double* v_original);

/* Compute reduced Hamiltonian: H_mu([z]) = H(z) for any z in J^{-1}(mu) */
double sr_mw_reduced_hamiltonian(const SRMarsdenWeinsteinSpace* reduced,
                                  const SRInvariantHamiltonian* H,
                                  const double* reduced_coords);

/* Check if mu is a regular value (dJ has full rank on J^{-1}(mu)) */
bool sr_mw_is_regular_value(const SRMomentumMap* J, double mu_value, double tolerance);

/* Print reduction summary */
void sr_mw_print(const SRMarsdenWeinsteinSpace* reduced);

/* === Cotangent Bundle Reduction === */

/* Theorem (Cotangent Bundle Reduction):
 *   T*G/G ≅ g*  (as Poisson manifolds, with Lie-Poisson bracket)
 *
 * This function realizes the identification:
 *   (T*G)_mu / G_mu ≅ O_mu  (coadjoint orbit of mu)
 */
void sr_cotangent_reduce(SRLieGroup* G, double* mu, int mu_dim,
                          SRLiePoissonBracket* result);

/* Compute the coadjoint orbit through mu:
 *   O_mu = {Ad*_{g^{-1}} mu : g in G} ≅ G/G_mu
 * For compact Lie groups, O_mu is a symplectic leaf of g*.
 */
void sr_coadjoint_orbit(const SRLieGroup* G, const double* mu,
                         int mu_dim, int num_points,
                         double** orbit_points);

/* Compute dimension of coadjoint orbit through mu:
 * dim(O_mu) = dim(G) - dim(G_mu) = rank of Kirillov-Kostant-Souriau symplectic form
 */
int sr_coadjoint_orbit_dimension(const SRLieGroup* G, const double* mu);

/* Compute the KKS (Kirillov-Kostant-Souriau) symplectic form on O_mu:
 *   omega_mu(ad*_xi mu, ad*_eta mu) = <mu, [xi, eta]>
 */
void sr_kks_symplectic_form(const SRLieGroup* G, const double* mu,
                             const double* xi, const double* eta,
                             double* omega_value);

/* === Euler-Poincare Reduction === */

/* Create an Euler-Poincare system for a Lie group G with Lagrangian l: g -> R.
 * Euler-Poincare equations:
 *   d/dt (delta_l/delta_xi) = ad*_xi (delta_l/delta_xi)
 * For right-invariant systems (sign=+1); sign=-1 for left-invariant.
 *
 * Source: Poincare (1901), Holm, Marsden & Ratiu (1998),
 *         "The Euler-Poincare equations and semidirect products"
 */
SREulerPoincareSystem* sr_ep_create(SRLieAlgebra* alg, SRLieGroup* group,
                                     int sign);

void sr_ep_free(SREulerPoincareSystem* ep);

/* Set the reduced Lagrangian via inertia matrix:
 * l(xi) = 1/2 <I xi, xi>  (quadratic kinetic energy)
 */
void sr_ep_set_inertia(SREulerPoincareSystem* ep, const double* I_matrix);

/* Compute the right-hand side of Euler-Poincare equations:
 * d(mu)/dt = ad*_xi mu  where mu = I xi
 */
void sr_ep_rhs(const SREulerPoincareSystem* ep, const double* xi,
                double* dxi_dt);

/* Compute coadjoint action (internal): ad*_xi eta
 * <ad*_xi eta, zeta> = <eta, [xi, zeta]>
 */
void sr_ep_coadjoint_action(const SREulerPoincareSystem* ep,
                             const double* xi, const double* eta,
                             double* ad_star_result);

/* Integrate EP equations by one timestep using RK4 */
void sr_ep_step(SREulerPoincareSystem* ep, double dt);

/* Compute conserved energy: E = 1/2 <mu, xi> */
double sr_ep_energy(const SREulerPoincareSystem* ep);

/* Compute Casimir functions (conserved for any Hamiltonian):
 * For so(3): C(mu) = ||mu||^2 (norm squared)
 * For general g*: any Ad*-invariant function
 */
double sr_ep_casimir_norm_sq(const SREulerPoincareSystem* ep);

void sr_ep_print(const SREulerPoincareSystem* ep);

/* === Lagrange-Poincare Reduction === */

/* Create Lagrange-Poincare system for TG/G reduction.
 * The LP equations consist of:
 *   (a) reduced Euler-Lagrange on g (the "dynamic" part)
 *   (b) reconstruction equation on G (the "kinematic" part)
 *
 * LP equations (without nonholonomic constraints):
 *   d/dt (dL/dxi) - ad*_xi (dL/dxi) = 0  [on g*]
 *   dg/dt = g·xi                          [on G]
 */
SRLagrangePoincareSystem* sr_lp_create(SRLieGroup* G);

void sr_lp_free(SRLagrangePoincareSystem* lp);

/* Set the reduced Lagrangian (quadratic form on g: l(xi) = 1/2 <I xi, xi>) */
void sr_lp_set_metric(SRLagrangePoincareSystem* lp, const double* I_matrix);

/* Compute LP right-hand side */
void sr_lp_rhs(const SRLagrangePoincareSystem* lp, const double* xi,
                const double* g_coords, double* dxi_dt, double* dg_dt);

/* Integrate LP system one step */
void sr_lp_step(SRLagrangePoincareSystem* lp, double dt);

/* Reconstruct configuration from reduced trajectory:
 * g(t) = g(0) * exp(∫_0^t xi(s) ds)  (for right-invariant)
 */
void sr_lp_reconstruct(const SRLagrangePoincareSystem* lp,
                        double* g_trajectory, int n_steps, double dt);

void sr_lp_print(const SRLagrangePoincareSystem* lp);

/* === Semidirect Product Reduction === */

/* Create semidirect product group G = S ⋉ V
 * Multiplication: (s1,v1)·(s2,v2) = (s1·s2, v1 + rho(s1)v2)
 *
 * Used for:
 *   SE(3) = SO(3) ⋉ R^3 (rigid body + center of mass)
 *   Heavy top: S^1 ⋉ (SO(3) ⋉ R^3)
 *   MHD: volume-preserving diffeomorphisms ⋉ magnetic fields
 *
 * Semidirect product Lie-Poisson bracket (Holm, Marsden, Ratiu 1998):
 *   {F,H}(mu,a) = <mu, [dF/d_mu, dH/d_mu]>
 *                + <a, dF/d_mu · dH/da - dH/d_mu · dF/da>
 * where mu in s*, a in V*.
 */
SRSemidirectProduct* sr_semidirect_create(SRLieGroup* S, SRLieGroup* V);

void sr_semidirect_free(SRSemidirectProduct* sp);

/* Set the representation rho: S -> Aut(V) */
void sr_semidirect_set_representation(SRSemidirectProduct* sp,
                                       const double* rho_matrix);

/* Compute semidirect product group multiplication */
void sr_semidirect_multiply(const SRSemidirectProduct* sp,
                             const double* s1, const double* v1,
                             const double* s2, const double* v2,
                             double* s_out, double* v_out);

/* Compute diamond operation (used in semidirect reduction):
 * a ◇ b ∈ s*  where <a ◇ b, xi> = <a, -xi·b>
 * xi·b = rho'(xi)·b = derivative of the representation
 */
void sr_semidirect_diamond(const SRSemidirectProduct* sp,
                            const double* a, const double* b,
                            double* diamond_result);

/* Compute semidirect Lie-Poisson bracket RHS */
void sr_semidirect_lie_poisson_rhs(const SRSemidirectProduct* sp,
                                    const double* hamiltonian_gradient_mu,
                                    const double* hamiltonian_gradient_a,
                                    const double* state_mu, const double* state_a,
                                    double* d_mu_dt, double* d_a_dt);

/* === Reduction by Stages === */

/* Reduce in two steps using a normal subgroup chain:
 *   1. Reduce by N ⊲ G  (normal subgroup)
 *   2. Reduce by G/N
 * Theorem (Marsden et al. 2007): Reduction by stages yields the same
 * reduced space as direct reduction by G.
 */
void sr_reduce_by_stages(SRSymplecticManifold* P,
                          SRLieGroup* G, SRLieGroup* N, /* N normal in G */
                          double mu_n_value, double mu_gn_value,
                          SRMarsdenWeinsteinSpace** stage1,
                          SRMarsdenWeinsteinSpace** stage2);

/* Verify commutativity of reduction by stages */
bool sr_verify_stages_commutativity(const SRMarsdenWeinsteinSpace* direct,
                                     const SRMarsdenWeinsteinSpace* staged,
                                     double tolerance);

/* === Poisson Reduction === */

/* Poisson reduction: If Poisson tensor pi is G-invariant,
 * then pi projects to a Poisson tensor on P/G (when the quotient
 * is a smooth manifold). The reduced Poisson tensor generally
 * has nonlinear structure (not Lie-Poisson unless P = T*G).
 */
void sr_poisson_reduction(const SRPoissonManifold* P,
                           SRLieGroup* G,
                           SRPoissonManifold* P_reduced);

/* Compute symplectic leaves of the reduced Poisson manifold */
void sr_poisson_leaves(const SRPoissonManifold* reduced,
                        double** leaf_representatives, int* n_leaves);

/* === Singular Reduction === */

/* For non-regular momentum values mu, J^{-1}(mu)/G_mu is not a smooth
 * manifold but a stratified symplectic space (Sjamaar-Lerman 1991).
 *
 * This function identifies the orbit type strata.
 */
void sr_singular_stratification(const SRMomentumMap* J,
                                 double mu_value,
                                 int* stratum_dimensions,
                                 int* n_strata);

/* === Momentum Map Level-Set Operations === */

/* Project a point to the zero level set J(z) = 0 using Newton's method.
 * Solve J(z + dΦ(ξ)) = 0 for ξ in g via:
 *   J_k(z) + (dJ(z) · X_{ξ})_k = 0
 * where X_ξ is the infinitesimal action of ξ on z.
 */
bool sr_project_to_zero_level_set(SRMomentumMap* J, double* z,
                                   double tolerance, int max_iter);

/* Compute the isotropy subalgebra g_mu = {xi in g : ad*_xi mu = 0} */
int sr_isotropy_subalgebra_dimension(const SRLieAlgebra* alg, const double* mu);

/* Compute a basis for g_mu */
void sr_isotropy_subalgebra_basis(const SRLieAlgebra* alg, const double* mu,
                                   double* basis, int* dim);

/* ============ Noether's Theorem ============ */

/* Noether's First Theorem (1918):
 * Every continuous symmetry of a Lagrangian system yields a conserved
 * quantity (the corresponding component of the momentum map J).
 *
 * Implementation: Verify d/dt J_ξ = {J_ξ, H} = 0 for G-invariant H.
 */
bool sr_noether_verify(const SRMomentumMap* J,
                        const SRInvariantHamiltonian* H,
                        const double* phase_point, double tolerance);

/* Compute all conserved quantities from symmetry generators */
void sr_noether_conserved_quantities(const SRMomentumMap* J,
                                      const double* phase_point,
                                      double* conserved_values);

/* Print Noether analysis */
void sr_noether_print(const SRMomentumMap* J,
                       const SRInvariantHamiltonian* H,
                       const double* phase_point);

#endif /* SR_REDUCTION_H */
