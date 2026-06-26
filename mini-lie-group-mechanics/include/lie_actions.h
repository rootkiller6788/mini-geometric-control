#ifndef LIE_ACTIONS_H
#define LIE_ACTIONS_H

#include "lie_core.h"

/* ==============================================================
 * lie_actions.h -- Group Actions and Symmetry
 *
 * Group Action: Phi : G x M -> M satisfying
 *   Phi(e, x) = x,  Phi(g, Phi(h, x)) = Phi(gh, x)
 *
 * References:
 *   Marsden & Ratiu (1999) Ch.4, Ch.9
 *   Abraham & Marsden (1978) Foundations of Mechanics
 * ============================================================== */

/* --- Group Action on Manifold State --- */
typedef struct {
    char *name;
    int manifold_dim;
    LieVector *position;
    LieVector *velocity;
} LieManifoldState;

LieManifoldState* lie_manifold_state_create(int dim, const char *name);
void lie_manifold_state_free(LieManifoldState *s);

/* Group action: x' = g . x  where x is in the linear representation space */
typedef struct {
    LieMatrix *rotation;
    LieVector *translation;
} LieRigidTransform;

LieRigidTransform* lie_rigid_transform_create(int dim);
void lie_rigid_transform_free(LieRigidTransform *tf);

/* Apply rigid transform to a point: p' = R p + t */
LieVector* lie_transform_point(const LieRigidTransform *tf, const LieVector *p);

/* Compose two rigid transforms: tf3 = tf1 o tf2 */
LieRigidTransform* lie_transform_compose(const LieRigidTransform *tf1,
                                          const LieRigidTransform *tf2);

/* Inverse rigid transform */
LieRigidTransform* lie_transform_inverse(const LieRigidTransform *tf);

/* Convert SE(3) element to rigid transform */
LieRigidTransform* lie_se3_to_transform(const LieGroupElement *g);

/* Convert rigid transform to SE(3) element */
LieGroupElement* lie_transform_to_se3(const LieRigidTransform *tf);

/* --- Momentum Map ---
 * J : M -> g*  (momentum map for a Hamiltonian G-action)
 * For SO(3) acting on R^3: J(q,p) = q x p  (angular momentum)
 * Satisfies: {J_xi, J_eta} = J_{[xi,eta]}
 *
 * This is the infinitesimal generator of the group action.
 *
 * Ref: Marsden & Ratiu (1999) Ch.11
 */
typedef struct {
    LieVector *linear_momentum;   /* p */
    LieVector *angular_momentum;  /* L = q x p for SO(3) */
} LieMomentum;

LieMomentum* lie_momentum_create(void);
void lie_momentum_free(LieMomentum *mom);

/* Compute angular momentum: L = r x p (SO(3) standard) */
LieVector* lie_angular_momentum_compute(const LieVector *r, const LieVector *p);

/* Momentum map for SO(3) action on T*R^3: mu = q x p */
LieVector* lie_momentum_map_so3_position(const LieVector *q, const LieVector *p);

/* --- Infinitesimal Generator ---
 * Given a Lie algebra element xi in g and a group action Phi,
 * the infinitesimal generator at x in M is:
 *   xi_M(x) = d/dt|_{t=0} Phi(exp(t xi), x)
 *
 * For SO(3) acting on R^3: xi_{R^3}(x) = omega x x
 *
 * Ref: Marsden & Ratiu (1999) SS4.1
 */
LieVector* lie_infinitesimal_generator_so3(const LieVector *omega,
                                            const LieVector *x);

/* --- Orbit and Stabilizer --- */
/* Orbit of x under G: G . x = { Phi(g, x) : g in G } */
/* Stabilizer (isotropy) subgroup: G_x = { g in G : Phi(g, x) = x } */

/* Check if two points are on the same SO(3) orbit (same norm) */
bool lie_same_so3_orbit(const LieVector *x1, const LieVector *x2, double tol);

/* --- Left/Right Translations on G ---
 * Left translation:  L_g(h) = g h
 * Right translation: R_g(h) = h g
 *
 * Left-invariant vector field: X(gh) = T_e L_g (X(e))
 * Right-invariant vector field: X(hg) = T_e R_g (X(e))
 *
 * Ref: Marsden & Ratiu (1999) SS9.1
 */
LieGroupElement* lie_left_translate(const LieGroupElement *g,
                                     const LieGroupElement *h);
LieGroupElement* lie_right_translate(const LieGroupElement *g,
                                      const LieGroupElement *h);

/* Tangent lift of left translation:
 *   T_e L_g : g -> T_g G,  xi |-> g xi  (for matrix groups)
 */
LieMatrix* lie_tangent_left_translate(const LieGroupElement *g,
                                       const LieMatrix *xi_matrix);

/* --- Maurer-Cartan Form ---
 * The canonical left-invariant g-valued 1-form on G:
 *   omega_g(v) = T_g L_{g^{-1}} (v) = g^{-1} v
 *
 * In coordinates: omega = g^{-1} dg
 * Satisfies Maurer-Cartan structure equation:
 *   d omega + (1/2)[omega, omega] = 0
 *
 * Ref: Sternberg (1964) Lectures on Differential Geometry, Ch.V
 */
LieMatrix* lie_maurer_cartan_form(const LieGroupElement *g,
                                   const LieMatrix *tangent_vector);

/* --- Coadjoint Orbit ---
 * For mu in g*, the coadjoint orbit is:
 *   O_mu = { Ad*_g(mu) : g in G }
 *
 * These are symplectic manifolds (Kirillov-Kostant-Souriau theorem).
 * For SO(3): coadjoint orbits are spheres ||mu|| = const
 *
 * Ref: Kirillov (2004) Lectures on the Orbit Method
 */
double lie_coadjoint_orbit_radius_so3(const LieVector *mu);

/* --- Locked Inertia Tensor ---
 * I(q) : g -> g*  (the locked inertia tensor for a coupled system)
 * Maps body angular velocity to body angular momentum.
 *
 * For a rigid body: I is the constant inertia tensor (in body frame).
 * For a general system: I depends on configuration q.
 *
 * Ref: Marsden & Ratiu (1999) SS13.5
 */
typedef struct {
    LieMatrix *inertia;   /* 3x3 inertia matrix */
} LieInertiaTensor;

LieInertiaTensor* lie_inertia_create(double Ixx, double Iyy, double Izz,
                                      double Ixy, double Ixz, double Iyz);
void lie_inertia_free(LieInertiaTensor *I);
LieVector* lie_inertia_apply(const LieInertiaTensor *I, const LieVector *omega);
LieInertiaTensor* lie_inertia_inverse(const LieInertiaTensor *I);
LieMatrix* lie_inertia_to_matrix(const LieInertiaTensor *I);

#endif /* LIE_ACTIONS_H */
