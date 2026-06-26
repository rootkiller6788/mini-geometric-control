#ifndef LIE_INTEGRATION_H
#define LIE_INTEGRATION_H

#include "lie_core.h"
#include "lie_actions.h"
#include "lie_reduction.h"

/* ==============================================================
 * lie_integration.h -- Geometric Integration on Lie Groups
 *
 * Standard numerical integrators (Euler, RK4) do not preserve the
 * Lie group structure: the numerical solution may leave G.
 *
 * Geometric integrators respect the manifold structure:
 *   - RKMK (Runge-Kutta-Munthe-Kaas): Uses the Lie algebra as
 *     a computational space, updating via exponential map.
 *   - Crouch-Grossman: Uses composition of exponentials.
 *   - Variational integrators: Derive from discrete variational
 *     principles; symplectic and momentum-preserving.
 *
 * For an ODE on G: dg/dt = g xi(g, t)  (right-trivialized)
 * RKMK updates: g_{n+1} = g_n exp(h * Theta)
 * where Theta is computed from Runge-Kutta stages on g.
 *
 * References:
 *   Iserles, Munthe-Kaas, Norsett, Zanna (2000) "Lie Group Methods"
 *     Acta Numerica 9:215-365
 *   Hairer, Lubich, Wanner (2006) Geometric Numerical Integration
 *   Munthe-Kaas (1999) "High Order Runge-Kutta Methods on Manifolds"
 *     Appl. Numer. Math. 29:115-127
 *   Marsden & West (2001) "Discrete Mechanics and Variational
 *     Integrators", Acta Numerica 10:357-514
 * ============================================================== */

/* --- Right-hand side function signature ---
 * f(g, t, xi_out): given g in G and time t, compute xi = g^{-1} g_dot
 * This is the right-trivialized tangent.
 */
typedef void (*LieRHSFunction)(const LieGroupElement *g, double t,
                                void *params, LieVector *xi_out);

/* --- RKMK Integrator (Runge-Kutta-Munthe-Kaas) --- */

/* RKMK Butcher tableau: uses standard RK coefficients */
typedef enum {
    LIE_RKMK_RK2 = 0,   /* Heun's method, order 2 */
    LIE_RKMK_RK3 = 1,   /* Kutta's method, order 3 */
    LIE_RKMK_RK4 = 2,   /* Classical RK4 on algebra, order 4 */
} LieRKMKOrder;

/* RKMK context for computing Theta on the algebra */
typedef struct {
    LieRKMKOrder order;
    int    stages;
    LieVector **k;       /* stages vectors on algebra */
    LieVector *theta;    /* Accumulated increment on algebra */
    LieVector *tmp;
    LieAlgebraElement *xi_alg; /* Work space */
} LieRKMKContext;

LieRKMKContext* lie_rkmk_create(LieRKMKOrder order, int algebra_dim);
void lie_rkmk_free(LieRKMKContext *ctx);

/* One RKMK step on SO(3):
 *   g_{n+1} = g_n * exp(theta)
 * where theta is computed from the RK stages of the
 * right-trivialized vector field.
 *
 * Algorithm (RK4-RKMK):
 *   k1 = h * f(g_n, t_n)
 *   k2 = h * f(g_n * exp(k1/2), t_n + h/2)
 *   k3 = h * f(g_n * exp(k2/2), t_n + h/2)
 *   k4 = h * f(g_n * exp(k3), t_n + h)
 *   theta = (k1 + 2*k2 + 2*k3 + k4) / 6
 *   g_{n+1} = g_n * exp(theta)
 *
 * Complexity: O(1) per step for SO(3).
 *
 * Ref: Munthe-Kaas (1999), Iserles et al. (2000) SS3
 */
void lie_rkmk_step_so3(const LieGroupElement *g_in, double t, double h,
                        LieRHSFunction f, void *params,
                        LieGroupElement *g_out);

/* RKMK step on SE(3) */
void lie_rkmk_step_se3(const LieGroupElement *g_in, double t, double h,
                        LieRHSFunction f, void *params,
                        LieGroupElement *g_out);

/* --- Crouch-Grossman Integrator ---
 * Uses composition of exponentials:
 *   g_{n+1} = g_n * exp(h a_1 f_1) * exp(h a_2 f_2) * ...
 *
 * For SO(3) with explicit Euler:  g_{n+1} = g_n * exp(h xi_n)
 * This is the Lie-Euler method, order 1 on the group.
 *
 * Ref: Crouch & Grossman (1993) "Numerical Integration of
 *      ODEs on Manifolds", J. Nonlinear Sci. 3:1-33
 */
void lie_euler_step_so3(const LieGroupElement *g_in, double t, double h,
                         LieRHSFunction f, void *params,
                         LieGroupElement *g_out);

void lie_euler_step_se3(const LieGroupElement *g_in, double t, double h,
                         LieRHSFunction f, void *params,
                         LieGroupElement *g_out);

/* --- Full Trajectory Integration --- */
typedef struct {
    int n_points;
    double *t;               /* Time array */
    LieGroupElement **g;     /* Trajectory on G */
    int dim;
} LieTrajectory;

LieTrajectory* lie_trajectory_create(int n_points, int algebra_dim);
void lie_trajectory_free(LieTrajectory *traj);
void lie_trajectory_print(const LieTrajectory *traj, int stride);

/* Integrate a full trajectory using RKMK-RK4 */
LieTrajectory* lie_rkmk4_integrate_so3(const LieGroupElement *g0,
                                        double t0, double tf, int n_steps,
                                        LieRHSFunction f, void *params);

LieTrajectory* lie_rkmk4_integrate_se3(const LieGroupElement *g0,
                                        double t0, double tf, int n_steps,
                                        LieRHSFunction f, void *params);

/* --- Magnus Expansion ---
 * For linear ODEs: dY/dt = A(t) Y, Y(0) = I
 * the solution is Y(t) = exp(Omega(t)) where Omega(t) is given by
 * the Magnus expansion:
 *
 *   Omega(t) = int_0^t A(t1) dt1
 *            + 1/2 int_0^t int_0^{t1} [A(t1), A(t2)] dt2 dt1
 *            + ...
 *
 * This preserves the Lie group structure automatically.
 *
 * Ref: Magnus (1954) "On the Exponential Solution of Differential
 *      Equations for a Linear Operator", CPAM 7:649-673
 *      Blanes, Casas, Oteo, Ros (2009) "The Magnus Expansion and
 *      Some of Its Applications", Phys. Rep. 470:151-238
 */
LieMatrix* lie_magnus_expansion_order2(const LieMatrix *A_n,
                                        const LieMatrix *A_np1, double h);

/* --- Variational Integrators on Lie Groups ---
 * Based on the discrete Euler-Poincare principle.
 *
 * Discrete Lagrangian: L_d : G x G -> R
 *   L_d(g_k, g_{k+1}) approximates int_{t_k}^{t_{k+1}} L(g, g_dot) dt
 *
 * Discrete Euler-Lagrange equations:
 *   D_2 L_d(g_{k-1}, g_k) + D_1 L_d(g_k, g_{k+1}) = 0
 *
 * For a left-invariant system on G, the reduced discrete equations
 * live on the Lie algebra.
 *
 * Ref: Marsden, Pekarsky, Shkoller (1999) "Discrete Euler-Poincare
 *      and Lie-Poisson Equations", Nonlinearity 12:1647-1662
 *      Bobenko & Suris (1999) "Discrete Lagrangian Reduction,
 *      Discrete Euler-Poincare Equations, and Semidirect Products",
 *      Lett. Math. Phys. 49:79-93
 */

/* Discrete Euler-Poincare step for SO(3) rigid body:
 * Uses the midpoint rule as the underlying quadrature.
 *
 * Given g_k and Omega_k (body angular velocity at step k),
 * compute g_{k+1} and Omega_{k+1}.
 *
 * The variational integrator exactly conserves a modified
 * angular momentum (the discrete Noether theorem).
 */
void lie_variational_rigid_body_step(const LieGroupElement *g_k,
                                      const LieVector *Omega_k,
                                      const LieInertiaTensor *I,
                                      double h,
                                      LieGroupElement *g_next,
                                      LieVector *Omega_next);

/* Variational integrator for heavy top */
void lie_variational_heavy_top_step(const LieGroupElement *g_k,
                                     const LieVector *Omega_k,
                                     const LieVector *Gamma_k,
                                     const LieHeavyTopState *ht,
                                     double h,
                                     LieGroupElement *g_next,
                                     LieVector *Omega_next,
                                     LieVector *Gamma_next);

/* --- Geodesic Computation on Lie Groups ---
 * The geodesic equation on G with left-invariant metric
 * defined by <.,.>_e on g is:
 *   dxi/dt = B(xi, xi)  (Euler-Poincare with zero potential)
 * where B is the bilinear form from the structure constants.
 *
 * Ref: do Carmo (1992) Riemannian Geometry, Ch.1
 */
void lie_geodesic_equation_so3(const LieVector *xi, const LieInertiaTensor *I,
                                LieVector *dxi_dt);

/* Compute geodesic distance on SO(3):
 *   dist(R1, R2) = ||log(R1^T R2)||
 *
 * Ref: Park (1995) "Distance Metrics on the Rigid-Body Motions
 *      with Applications to Mechanism Design", ASME J. Mech. Des.
 */
double lie_geodesic_distance_so3(const LieMatrix *R1, const LieMatrix *R2);

/* Interpolation on SO(3):
 *   R(t) = R1 * exp(t * log(R1^T R2))
 * for t in [0, 1] gives the geodesic between R1 and R2.
 *
 * Ref: Shoemake (1985) "Animating Rotation with Quaternion Curves"
 */
LieGroupElement* lie_so3_geodesic_interpolate(const LieGroupElement *g0,
                                                const LieGroupElement *g1,
                                                double t);

#endif /* LIE_INTEGRATION_H */
