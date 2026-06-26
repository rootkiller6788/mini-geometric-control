#ifndef SR_DYNAMICS_H
#define SR_DYNAMICS_H

#include "sr_core.h"
#include "sr_reduction.h"
#include "sr_poisson.h"
#include <stdbool.h>

/* ============================================================================
 * Reduced Dynamics Engine
 *
 * Integration of reduced equations of motion, reconstruction,
 * and analysis of reduced mechanical/control systems.
 *
 * Key flows:
 *   1. Reduced Hamiltonian flow on P_mu (symplectic reduction)
 *   2. Lie-Poisson flow on g* (Euler-Poincare, rigid body)
 *   3. Reconstruction: lift reduced trajectory to full space
 *
 * References:
 *   de Leon & Rodrigues (1989) "Methods of Differential Geometry in Analytical Mechanics"
 *   Marsden & West (2001) "Discrete mechanics and variational integrators"
 *   Hairer, Lubich & Wanner (2006) "Geometric Numerical Integration"
 * ============================================================================ */

/* === Numerical Integrators for Reduced Systems === */

/* RK4 integration on coadjoint orbits (Lie-Poisson).
 * Integrates d_mu/dt = ad*_{nabla_H(mu)} mu on g*.
 * Preserves coadjoint orbit (Casimirs conserved to O(dt^4)).
 */
void sr_integrate_lie_poisson_rk4(const SRLiePoissonBracket* lpb,
                                   void (*hamiltonian_gradient)(const double*, double*),
                                   double* mu, int mu_dim,
                                   double dt, int n_steps,
                                   double** trajectory);

/* Symplectic Euler (first-order symplectic integrator) for
 * Hamiltonian system on reduced space.
 * Properties: preserves Casimirs, preserves symplectic form (up to order).
 */
void sr_integrate_symplectic_euler(const SRSymplecticManifold* sp,
                                    const SRInvariantHamiltonian* H,
                                    double* q, double* p, int n_dof,
                                    double dt, int n_steps);

/* Stormer-Verlet / leapfrog integrator for separable Hamiltonians:
 * H(q, p) = T(p) + V(q) -> symplectic, time-reversible, 2nd order.
 */
void sr_integrate_stormer_verlet(double* q, double* p, int n_dof,
                                  void (*force)(const double* q, double* f),
                                  double* mass_matrix,
                                  double dt, int n_steps);

/* Lie group variational integrator (Marsden, Pekarsky, Shkoller 1999).
 * For TG/G reduction: uses exponential map to stay on the group.
 *
 * Algorithm (one step):
 *   1. Solve discrete Euler-Lagrange eqn: D_2 L_d(g_{k-1}, g_k) + D_1 L_d(g_k, g_{k+1}) = 0
 *   2. Update using exp map: g_{k+1} = g_k exp(dt * xi_k)
 */
void sr_integrate_lie_group_var(SRLagrangePoincareSystem* lp,
                                 double dt, int n_steps,
                                 double** g_trajectory, double** xi_trajectory);

/* RATTLE integrator for constrained Hamiltonian systems (e.g., J(z)=mu).
 * Preserves constraints and symplectic structure.
 * Used for integrating on J^{-1}(mu) without explicit coordinates.
 */
void sr_integrate_rattle(const SRInvariantHamiltonian* H,
                          const SRMomentumMap* J, double mu_value,
                          double* q, double* p, int n_dof,
                          double dt, int n_steps);

/* === Reconstruction === */

/* Reconstruction of full trajectory from reduced trajectory:
 * Given xi(t) in g (reduced velocity), solve:
 *   dg/dt = g · xi(t)  (for right-invariant systems)
 * or dg/dt = xi(t) · g  (for left-invariant systems).
 *
 * This is a Lie group ODE solved by Magnus expansion or
 * piecewise exponential approximation.
 */
void sr_reconstruct_group_trajectory(SRLieGroup* G,
                                      const double** xi_history,
                                      const double* g0,
                                      int n_steps, double dt,
                                      double** g_trajectory);

/* Reconstruction in T*G (cotangent bundle):
 * Given mu(t) in g* and g0 in G, reconstruct:
 *   g(t) from dg/dt = g · xi(t) where xi(t) = I^{-1} mu(t)
 *   alpha_g(t) from alpha_g = T_e^* L_{g^{-1}} mu(t)
 */
void sr_reconstruct_cotangent_trajectory(SRLieGroup* G,
                                          const double** mu_history,
                                          const double* g0,
                                          int n_steps, double dt,
                                          double** g_trajectory,
                                          double** p_trajectory);

/* === Invariant Detection and Analysis === */

/* Test whether a Hamiltonian H is G-invariant:
 * H(g·z) = H(z) for all g in G.
 * Equivalent to: dH(X_ξ(z)) = 0 for all ξ in g (Lie derivative vanishes).
 */
bool sr_check_hamiltonian_invariance(const SRInvariantHamiltonian* H,
                                      const SRLieGroup* G,
                                      const double* z, int z_dim,
                                      double tolerance);

/* Compute symmetry group orbits of a phase point */
void sr_compute_group_orbit(const SRLieGroup* G, const double* z,
                             int z_dim, int sample_points,
                             double** orbit);

/* Check if a vector field X is G-equivariant: X(g·z) = g·X(z) */
bool sr_check_equivariance(const double* vector_field,
                            const SRLieGroup* G, const double* z,
                            int z_dim, double tolerance);

/* === Reduced Analysis === */

/* Compute the dimension of the reduced space:
 * dim(P_mu) = dim(P) - dim(G) - dim(G_mu)
 * For cotangent bundle reduction of T*G:
 *   dim(g*) = dim(G) [as vector space]
 *   But coadjoint orbits O_mu have dimension: dim(G) - dim(G_mu)
 */
int sr_reduced_dimension(int dim_P, const SRLieGroup* G, const double* mu);

/* Determine if reduction produces a relative equilibrium:
 * xi_equilibrium in g such that ad*_{xi_eq} mu = 0 (EP eqns at fixed point).
 */
bool sr_find_relative_equilibrium(const SRLiePoissonBracket* lpb,
                                   const double* mu, double* xi_eq);

/* Linearized reduced dynamics: variational equations on reduced space.
 * For small perturbations δmu about trajectory mu(t):
 *   d(δmu)/dt = D(ad*_{xi} mu)(δmu) = ad*_{δxi} mu + ad*_{xi} δmu
 */
void sr_linearize_reduced_dynamics(const SRLiePoissonBracket* lpb,
                                    const double* mu, const double* xi,
                                    double* linearized_matrix);

/* Reduced stability analysis: eigenvalues of linearized reduced dynamics */
void sr_reduced_stability(const SRLiePoissonBracket* lpb,
                           const double* mu_eq, const double* xi_eq,
                           double* eigenvalues_real, double* eigenvalues_imag);

/* === Bifurcation Analysis in Reduced Systems === */

/* Detect Hamiltonian-Hopf bifurcation in reduced system.
 * Occurrence: two purely imaginary eigenvalue pairs collide at
 * nonzero frequencies (1:1 resonance) on a symplectic leaf.
 */
bool sr_detect_hamiltonian_hopf(const SRLiePoissonBracket* lpb,
                                 const double* mu, double parameter,
                                 double* bifurcation_threshold);

/* === Integrable Reduced Systems === */

/* Test for Liouville integrability in reduced space:
 * Finding n independent involutive integrals on 2n-dimensional leaf.
 * For g*: Casimirs provide some integrals naturally.
 */
bool sr_test_liouville_integrability(const SRPoissonManifold* pm,
                                      const double* x,
                                      double* integrals, int* n_integrals);

/* Action-angle coordinates for integrable reduced system.
 * Action: I_i = 1/(2π) ∮_{γ_i} p·dq (loop integrals around invariant tori).
 */
void sr_action_angle_compute(const SRPoissonManifold* pm,
                              const double* x, int n_actions,
                              double* actions, double* angles);

/* === Energy-Momentum Analysis === */

/* Energy-momentum method for stability of relative equilibria:
 * For a relative equilibrium z_eq with momentum mu:
 *   (1) Compute the augmented Hamiltonian: H_xi(z) = H(z) - <J(z)-mu, xi>
 *   (2) Check d^2 H_xi(z_eq) on ker(dJ(z_eq)).
 * If > 0 on this subspace, the relative equilibrium is (orbitally) stable.
 *
 * Arnold's extension for Lie-Poisson systems on g*:
 * For coadjoint orbit, stability ⇔ H restricted to orbit has
 * a non-degenerate minimum at mu_eq.
 */
void sr_energy_momentum_stability(const SRInvariantHamiltonian* H,
                                   const SRMomentumMap* J,
                                   const double* z_eq, double mu_value,
                                   double* hessian_on_kernel,
                                   bool* is_stable);

/* Arnold stability on coadjoint orbits */
bool sr_arnold_stability(const SRLiePoissonBracket* lpb,
                          void (*H_func)(const double*, double*),
                          const double* mu_eq);

/* === Geometric Phase and Holonomy === */

/* Compute the geometric (Berry) phase accumulated during a closed
 * loop in the reduced space (shape space).
 * For a loop c(t) in P/G (shape space), the reconstruction:
 *   g(T) = g(0) exp(∮_c A)
 * where A is the mechanical connection 1-form on P -> P/G.
 *
 * The geometric phase = ∮_c A depends only on the loop in shape space,
 * not on the speed of traversal (holonomy of the connection).
 */
void sr_geometric_phase(const double** shape_trajectory,
                         int shape_dim, int n_steps,
                         const double* connection_form,
                         double* phase, int phase_dim);

/* Compute the dynamic phase (non-geometric contribution) */
void sr_dynamic_phase(const double** shape_trajectory,
                       const double* inertia, int n_steps,
                       double* dynamic_phase);

/* Compute the total holonomy from the curvature of the connection:
 * Phase = ∫_Σ F (Stokes theorem: F = dA + [A, A])
 * where Σ is a surface bounded by the shape loop.
 */
void sr_curvature_holonomy(const double* connection_curvature,
                            int phase_dim, int shape_dim,
                            double* holonomy);

/* === Spherical Pendulum Reduction Example === */

/* Spherical pendulum has S^1 symmetry (rotation about vertical axis).
 * Phase space: T*S^2 (dim=4), symmetry G=S^1, momentum J = p_phi.
 * Reduced space (for mu = p_phi): T*[0,pi] (dim=2).
 *
 * Reduced Hamiltonian (on reduced phase space, 2-dim):
 * Effective potential = mu^2/(2m l^2 sin^2(theta)) - mgl cos(theta)
 */
void sr_spherical_pendulum_reduced_rhs(double theta, double p_theta,
                                        double mu, double m, double l, double g,
                                        double* dtheta_dt, double* dp_dt);

double sr_spherical_pendulum_effective_potential(double theta,
                                                  double mu, double m, double l, double g);

/* === Double Spherical Pendulum with Reduction === */

/* Double spherical pendulum has S^1 symmetry (rotation about vertical).
 * Phase space: T*S^2 × S^2 (dim=8), symmetry S^1, momentum = total p_phi.
 * Reduced space: 6-dimensional.
 */
void sr_double_spherical_rhs(const double* state_6d, double mu,
                              const double* params, double* derivs);

/* ============ Mechanical Connection ============ */

/* Mechanical connection (Hannay-Berry connection for mechanical systems).
 * Defined by: A_m(v_q) = I(q)^{-1} J_L(v_q)
 * where I(q) is the locked inertia tensor and J_L is the momentum map
 * for the lifted action on TQ.
 *
 * The connection 1-form A: TQ -> g satisfies:
 *   A(ξ_Q(q)) = ξ  (for vertical vectors)
 *   A is G-equivariant
 */
void sr_mechanical_connection(const double* q, int q_dim,
                               const double* metric, SRLieGroup* G,
                               double* connection_1form);

/* Compute the locked inertia tensor I(q): g -> g* */
void sr_locked_inertia(const double* q, int q_dim,
                        const double* metric, SRLieGroup* G,
                        double* inertia_tensor);

/* === Relative Equilibria Tracking === */

/* Track a family of relative equilibria as parameters vary.
 * Solves: ad*_{xi} mu = 0, mu = I xi (for EP systems with metric I).
 */
void sr_relative_equilibria_family(const SREulerPoincareSystem* ep,
                                    double param_start, double param_end,
                                    int n_steps, double** equilibria);

/* === Efficient Coadjoint Orbit Sampling === */

/* Generate uniform samples from a coadjoint orbit O_mu.
 * For compact G, O_mu ≅ G/G_mu is a homogeneous space.
 * Sampling via Monte Carlo on G and projecting via Ad*.
 */
void sr_orbit_monte_carlo(const SRLieGroup* G, const double* mu,
                           int n_samples, int mu_dim,
                           double** samples);

#endif /* SR_DYNAMICS_H */
