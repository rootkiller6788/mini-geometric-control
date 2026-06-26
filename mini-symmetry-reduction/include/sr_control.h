#ifndef SR_CONTROL_H
#define SR_CONTROL_H

#include "sr_core.h"
#include "sr_reduction.h"
#include <stdbool.h>

/* ============================================================================
 * Symmetry Reduction in Control Systems
 *
 * Control-affine systems with symmetries:
 *   dz/dt = f(z) + Σ_i g_i(z) u_i
 * where f, g_i are G-equivariant vector fields.
 *
 * Key results:
 *   1. Controlled Euler-Poincare equations
 *   2. Reduction of optimal control problems (Bloch, Leonard, Marsden)
 *   3. Kinematic reduction for nonholonomic systems
 *   4. Controlled Lagrangians method
 *   5. Symmetry-preserving feedback design
 *
 * References:
 *   Bloch (2003) "Nonholonomic Mechanics and Control"
 *   Bullo & Lewis (2004) "Geometric Control of Mechanical Systems"
 *   Bloch, Leonard & Marsden (2000) "Controlled Lagrangians"
 * ============================================================================ */

/* === Reduced Control System === */

/* Create a reduced control system */
SRReducedControlSystem* sr_ctrl_create(int state_dim, int control_dim,
                                        int algebra_dim);

void sr_ctrl_free(SRReducedControlSystem* ctrl);

/* Set the reduced drift vector field f(x) */
void sr_ctrl_set_drift(SRReducedControlSystem* ctrl, const double* drift);

/* Set the control vector fields g_i(x) (state_dim × control_dim) */
void sr_ctrl_set_control_fields(SRReducedControlSystem* ctrl,
                                 const double* control_fields);

/* Evaluate the reduced RHS: dx/dt = f(x) + Σ g_i(x) u_i */
void sr_ctrl_rhs(const SRReducedControlSystem* ctrl,
                  const double* x, const double* u,
                  double* dx_dt);

/* === Controlled Euler-Poincare Equations === */

/* Controlled Euler-Poincare system for a mechanical system with
 * G-symmetry and body-fixed forces:
 *   d/dt (∂l/∂ξ) = ad*_ξ (∂l/∂ξ) + F_{ext}
 * where F_{ext} ∈ g* is the body-fixed external force/torque.
 *
 * For a rigid body with body-fixed torques:
 *   I dΩ/dt = I Ω × Ω + τ  (Euler equations with control torques)
 */
void sr_ctrl_ep_rhs(const SREulerPoincareSystem* ep,
                     const double* xi, const double* control_forces,
                     double* dxi_dt);

/* === Kinematic Control on Lie Groups === */

/* Driftless kinematic control system on a Lie group:
 *   dg/dt = g · (A_1 u_1 + ... + A_m u_m)
 * where A_i ∈ g are constant body-fixed control directions.
 *
 * Controllability condition (Chow's theorem / Rashevskii):
 *   Lie bracket generating condition: {A_i} and their iterated
 *   Lie brackets must span g (the Lie algebra rank condition LARC).
 */
void sr_kinematic_rhs(const SRLieGroup* G, const SRGroupElement* g,
                       const SRLieAlgebraElement* controls, int n_controls,
                       const double* u, SRGroupElement* dg);

/* Verify controllability via the Lie Algebra Rank Condition (LARC) */
bool sr_kinematic_controllable(const SRLieAlgebra* alg,
                                const SRLieAlgebraElement** A_i,
                                int n_inputs);

/* Compute the growth vector (dimensions of iterated Lie brackets):
 * dim(Delta), dim(Delta^2), dim(Delta^3), ...
 */
void sr_kinematic_growth_vector(const SRLieAlgebra* alg,
                                 const SRLieAlgebraElement** A_i,
                                 int n_inputs, int* dims, int max_iter);

/* === Optimal Control Reduction === */

/* Reduce an optimal control problem on T*G by G symmetry.
 *
 * Original problem:
 *   min_u ∫_0^T L(g, xi, u) dt  subject to  dg/dt = g·xi
 *
 * Reduced problem (via application of Pontryagin maximum principle
 * on the reduced space g*):
 *   min_u ∫_0^T l(mu, xi, u) dt  subject to mu dynamics on g*.
 *
 * Theorem (Bloch, Leonard, Marsden 2000): For G-invariant optimal
 * control problems on Lie groups, the maximum principle reduces
 * to a problem on g*.
 */
void sr_reduce_optimal_control(const SRInvariantHamiltonian* cost_integrand,
                                SRLieGroup* G,
                                SRLiePoissonBracket* reduced_pb,
                                void** reduced_data);

/* Solve LQR on a reduced space with Lie-Poisson dynamics.
 * For linearized reduced dynamics δmu_dot = A δmu + B δu,
 * find K such that u = -K δmu minimizes quadratic cost.
 */
void sr_lqr_reduced(const double* A, const double* B,
                     const double* Q, const double* R,
                     int state_dim, int control_dim,
                     double* K_gain, double* S_solution);

/* === Controlled Lagrangians === */

/* The controlled Lagrangian method for stabilization:
 *   L_controlled = L_original + L_modification
 * where L_modification adds gyroscopic forces to stabilize
 * unstable equilibria (e.g., inverted pendulum on a cart).
 *
 * Reduced version: modifies the reduced Lagrangian l(xi) on g
 * by adding a quadratic term that changes the effective metric.
 */
void sr_controlled_lagrangian_modify(SREulerPoincareSystem* ep,
                                      const double* modification_matrix,
                                      double gain);

/* Compute matching conditions: ensure the modified equations
 * correspond to a physically realizable control law.
 */
bool sr_ctrl_matching_conditions(const SREulerPoincareSystem* original,
                                  const SREulerPoincareSystem* modified,
                                  const double* original_metric,
                                  const double* modified_metric,
                                  double tolerance);

/* === Energy Shaping via Symmetry Reduction === */

/* Energy shaping using Casimir functions:
 * H_desired(mu) = H_original(mu) + Σ_i k_i C_i(mu)
 * where C_i are Casimir functions (commute with all functions
 * under the Lie-Poisson bracket).
 *
 * This preserves the Poisson structure while modifying the
 * Hamiltonian to stabilize relative equilibria.
 */
void sr_energy_shape_casimirs(const SRLiePoissonBracket* lpb,
                               const double* original_H,
                               const double* weights, int n_casimirs,
                               double* shaped_H);

/* === Nonholonomic Reduction === */

/* Nonholonomic constraints on a Lie group:
 * A(q) · dq/dt = 0  (linear velocity constraints)
 *
 * After symmetry reduction, the constraints become:
 *   xi ∈ h ⊂ g  (h is a subspace, possibly not a subalgebra)
 *
 * Lagrange-d'Alembert principle on reduced space:
 *   d/dt (∂l/∂xi) - ad*_xi (∂l/∂xi) ∈ h^0  (constraint forces)
 *   xi ∈ h  (constraint)
 */
void sr_nonholonomic_reduced_rhs(const SRLieAlgebra* alg,
                                  const double* xi, const double* mu,
                                  const double* constraint_basis, int h_dim,
                                  double* dxi_dt);

/* === Reconstruction with Feedback === */

/* Reconstruct full state from reduced state and feedback law:
 * Given reduced controller u(xi) on g and reconstruction:
 *   dxi/dt = RHS(xi, u(xi))  [reduced dynamics]
 *   dg/dt = g · xi           [reconstruction]
 *   y = h(g)                 [output in workspace]
 */
void sr_ctrl_reconstruct_feedback(const SRLagrangePoincareSystem* lp,
                                   void (*reduced_controller)(const double*, double*),
                                   double* g_state, double* xi_state,
                                   double dt, int n_steps);

/* === Inverse Optimal Control on Reduced Space === */

/* Given a stabilizing feedback law u = -K(xi) xi on g,
 * find a cost functional that this law optimizes.
 * This is the "inverse optimal" problem (Kalman 1964, Freeman & Kokotovic 1996).
 *
 * On reduced spaces, the solution involves solving a
 * Hamilton-Jacobi-Bellman equation restricted to coadjoint orbits.
 */
void sr_inverse_optimal_reduced(const SRLiePoissonBracket* lpb,
                                 const double* K, int n, double* Q, double* R);

/* === Underactuated Systems and Symmetry === */

/* For underactuated systems with symmetry (e.g., spacecraft with
 * fewer thrusters than degrees of freedom), the reduced control
 * system on g* has restricted control directions.
 *
 * Controllability on g*: Given control directions B_1,...,B_m in g*,
 * the system is controllable if the generated Lie algebra under
 * drift + controls spans the full state space.
 */
bool sr_underactuated_controllability(const SRLiePoissonBracket* lpb,
                                       const double* drift_mu,
                                       const double** control_directions,
                                       int n_controls, double tolerance);

/* === Practical Applications === */

/* Spacecraft attitude control on SO(3):
 * Reduced dynamics: I dΩ/dt = I Ω × Ω + τ (body-fixed torques)
 * State: Ω ∈ so(3) ≅ R^3 (body angular velocity)
 */
void sr_spacecraft_reduced_rhs(const double* Omega, const double* I_diag,
                                const double* tau, double* dOmega_dt);

/* Underwater vehicle on SE(3):
 * State: (v, omega) ∈ se(3) (body-fixed linear + angular velocity)
 * Reduced Kirchhoff equations with added mass and drag.
 */
void sr_underwater_vehicle_reduced_rhs(const double* state_6d,
                                        const double* mass_matrix_6x6,
                                        const double* forces_6d,
                                        double* derivs_6d);

/* Quadrotor UAV with SO(3) × R^3 symmetry:
 * Reduced attitude dynamics + translational dynamics.
 */
void sr_quadrotor_reduced_rhs(const double* Omega, const double* v,
                               const double* I_diag, double mass,
                               const double* thrust, const double* torque,
                               double* dOmega_dt, double* dv_dt);

/* Robot arm with S^1 symmetries (each joint):
 * Reduced dynamics on product of circles (toroidal geometry).
 * Symmetries arise from cyclic coordinates (joint angles that
 * don't appear in Lagrangian).
 */
void sr_robot_arm_reduced_rhs(const double* theta, const double* p_theta,
                               int n_joints, const double* inertia_matrix,
                               const double* potential_gradient,
                               double* dtheta_dt, double* dp_dt);

#endif /* SR_CONTROL_H */
