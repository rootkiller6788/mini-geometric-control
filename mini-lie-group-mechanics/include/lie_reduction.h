#ifndef LIE_REDUCTION_H
#define LIE_REDUCTION_H

#include "lie_core.h"
#include "lie_actions.h"

/* ==============================================================
 * lie_reduction.h -- Lie-Poisson and Euler-Poincare Reduction
 *
 * Symmetry reduction eliminates the symmetry directions from the
 * equations of motion, reducing dimension while preserving structure.
 *
 * Key Theorems:
 *   1. Euler-Poincare Theorem (reduction of Euler-Lagrange eqns)
 *   2. Lie-Poisson Reduction (Hamiltonian reduction on g*)
 *   3. Noether's Theorem (symmetries -> conserved quantities)
 *
 * The Lagrangian/Hamiltonian on TG (or T*G) that is G-invariant
 * induces a reduced system on g (or g*).
 *
 * References:
 *   Marsden & Ratiu (1999) Ch.13 -- Euler-Poincare and Lie-Poisson
 *   Holm, Marsden, Ratiu (1998) "Euler-Poincare Models of Ideal
 *     Fluids with Nonlinear Dispersion", Phys. Rev. Lett.
 *   Poincare (1901) "Sur une forme nouvelle des equations de
 *     la mecanique"
 * ============================================================== */

/* --- Reduced Lagrangian ---
 * l : g -> R  (Lagrangian reduced to the Lie algebra)
 * Typically: l(xi) = (1/2) <I xi, xi> - U(q)
 *
 * The Euler-Poincare equations are:
 *   d/dt (delta l / delta xi) = ad*_xi (delta l / delta xi)
 *
 * where delta l / delta xi is the functional derivative.
 */
typedef struct {
    char *name;
    LieInertiaTensor *locked_inertia;
    double potential_offset;
} LieReducedLagrangian;

LieReducedLagrangian* lie_reduced_lagrangian_create(
    const LieInertiaTensor *I, const char *name);
void lie_reduced_lagrangian_free(LieReducedLagrangian *l);

/* Evaluate reduced Lagrangian: l(xi) = (1/2) xi^T I xi */
double lie_reduced_lagrangian_eval(const LieReducedLagrangian *l,
                                    const LieVector *xi);

/* Functional derivative: delta l / delta xi = I xi = mu */
LieVector* lie_reduced_lagrangian_deriv(const LieReducedLagrangian *l,
                                         const LieVector *xi);

/* --- Euler-Poincare Equations ---
 * For SO(3) (rigid body):
 *   d/dt (I omega) = I omega x omega + tau
 *
 * In general: d/dt (delta_l/delta_xi) = ad*_xi (delta_l/delta_xi) + F
 * where F is the external force in the body frame.
 *
 * Ref: Marsden & Ratiu (1999) Theorem 13.5.2
 */
void lie_euler_poincare_rhs(const LieReducedLagrangian *l,
                             const LieVector *xi,
                             const LieVector *external_force,
                             LieVector *dxi_dt);

/* --- Lie-Poisson Bracket ---
 * On g*: {F, H}(mu) = < mu, [delta F / delta mu, delta H / delta mu] >
 *
 * For SO(3): {F, H}(Pi) = -Pi . (nabla F x nabla H)
 * This is the (-) Lie-Poisson bracket.
 *
 * The equations of motion for any Hamiltonian H : g* -> R:
 *   dF/dt = {F, H}   for all F
 *   =>  dmu/dt = ad*_{delta H / delta mu}(mu)
 *
 * Ref: Marsden & Ratiu (1999) SS13.7
 */
double lie_lie_poisson_bracket_so3(const LieVector *mu,
                                    void (*grad_F)(const LieVector*, LieVector*),
                                    void (*grad_H)(const LieVector*, LieVector*),
                                    const LieVector *mu_current);

/* --- Rigid Body on SO(3) ---
 * Reduced Hamiltonian: H(Pi) = (1/2) Pi^T I^{-1} Pi
 * Equations: dPi/dt = Pi x (I^{-1} Pi) + tau
 * (Euler equations for rigid body)
 *
 * Body angular velocity: Omega = I^{-1} Pi
 * Attitude kinematics: dR/dt = R Omega^hat
 *
 * Ref: Marsden & Ratiu (1999) SS15.9
 */
typedef struct {
    LieInertiaTensor *I_body;
    LieMatrix *R;           /* Attitude (SO(3) rotation matrix) */
    LieVector *Pi;          /* Body angular momentum */
    LieVector *Omega;       /* Body angular velocity */
    LieVector *tau;         /* External torque (body frame) */
} LieRigidBodyState;

LieRigidBodyState* lie_rigid_body_create(double Ixx, double Iyy, double Izz,
                                          double Ixy, double Ixz, double Iyz);
void lie_rigid_body_free(LieRigidBodyState *rb);
void lie_rigid_body_reset(LieRigidBodyState *rb);

/* Euler rigid body equations:
 *   I Omega_dot = I Omega x Omega + tau
 *   Pi_dot = Pi x Omega + tau
 *   R_dot = R Omega^hat
 *
 * Ref: Euler (1758), Arnold (1989) Mathematical Methods of Classical
 *      Mechanics, Appendix 2
 */
void lie_rigid_body_rhs(const LieRigidBodyState *rb,
                         LieVector *Pi_dot,
                         LieMatrix *R_dot);

/* Compute kinetic energy */
double lie_rigid_body_kinetic_energy(const LieRigidBodyState *rb);

/* --- Heavy Top (Lagrange-Poisson Equations) ---
 * System on SO(3) with a gravitational potential.
 * Reduced Hamiltonian:
 *   H(Pi, Gamma) = (1/2) Pi^T I^{-1} Pi + m g l Gamma_z
 *
 * where Gamma = R^T e_z (gravity direction in body frame).
 *
 * Equations:
 *   Pi_dot  = Pi x Omega + m g l Gamma x chi
 *   Gamma_dot = Gamma x Omega
 *
 * where chi is the center of mass direction in body frame.
 *
 * Ref: Marsden & Ratiu (1999) SS15.10
 *      Goldstein, Poole, Safko (2002) Classical Mechanics, Ch.5
 */
typedef struct {
    LieInertiaTensor *I_body;
    LieMatrix *R;
    LieVector *Pi;         /* Body angular momentum */
    LieVector *Omega;      /* Body angular velocity */
    LieVector *Gamma;      /* Gravity direction in body frame */
    LieVector *chi;        /* Center of mass direction (body frame) */
    double mass;
    double g_accel;
    double com_distance;   /* Distance from fixed point to COM */
} LieHeavyTopState;

LieHeavyTopState* lie_heavy_top_create(double Ixx, double Iyy, double Izz,
                                        double mass, double com_dist);
void lie_heavy_top_free(LieHeavyTopState *ht);
void lie_heavy_top_reset(LieHeavyTopState *ht);

/* Heavy top right-hand side:
 *   Pi_dot = Pi x Omega + m g l Gamma x chi
 *   Gamma_dot = Gamma x Omega
 *   R_dot = R Omega^hat
 */
void lie_heavy_top_rhs(const LieHeavyTopState *ht,
                        LieVector *Pi_dot,
                        LieVector *Gamma_dot,
                        LieMatrix *R_dot);

/* Compute total energy (kinetic + potential) */
double lie_heavy_top_energy(const LieHeavyTopState *ht);

/* --- Noether's Theorem ---
 * Every continuous symmetry of the Lagrangian yields a conserved quantity.
 *
 * For a left-invariant system on TG:
 *   Invariance under right translation -> conserved spatial momentum
 *   Invariance under left translation  -> conserved body momentum
 *
 * Specifically, the momentum map J : T*G -> g* is conserved
 * along the flow of a G-invariant Hamiltonian.
 *
 * Ref: Noether (1918), Marsden & Ratiu (1999) Ch.11
 */
double lie_noether_conserved_quantity(const LieReducedLagrangian *l,
                                       const LieVector *xi);

/* --- Reduction by Stages ---
 * For a semidirect product G = S semidirect V (e.g., SE(3) = SO(3) semidirect R^3),
 * reduction can be performed in stages:
 *   T*G / G  ~  g*  (Lie-Poisson on se(3)*)
 *   or: first reduce by V, then by S.
 *
 * Ref: Marsden, Misiolek, Ortega, Perlmutter, Ratiu (2007)
 *      Hamiltonian Reduction by Stages
 */
typedef struct {
    LieVector *Pi;    /* Angular momentum (so(3)* component) */
    LieVector *P;     /* Linear momentum (R^3* component)   */
} LieSE3Momentum;

LieSE3Momentum* lie_se3_momentum_create(void);
void lie_se3_momentum_free(LieSE3Momentum *mom);

/* SE(3) reduced bracket (semidirect product Lie-Poisson):
 *   dF/dt = {F, H}_{SE3}
 *
 * Ref: Marsden & Ratiu (1999) SS13.8
 */
void lie_se3_lie_poisson_rhs(const LieSE3Momentum *mu,
                              const LieVector *body_velocity,
                              const LieVector *body_angular_velocity,
                              LieVector *dmu_dt);

#endif /* LIE_REDUCTION_H */
