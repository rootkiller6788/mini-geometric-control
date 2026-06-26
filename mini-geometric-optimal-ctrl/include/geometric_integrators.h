/**
 * geometric_integrators.h — Geometric Numerical Integration
 *
 * References:
 *   Hairer, Lubich, Wanner (2006) "Geometric Numerical Integration"
 *   Leimkuhler & Reich (2004) "Simulating Hamiltonian Dynamics"
 *   Marsden & West (2001) "Discrete Mechanics and Variational Integrators"
 *
 * Covers: Symplectic RK, discrete variational mechanics, energy-preserving
 * methods, Yoshida composition, DMOC for optimal control.
 *
 * Knowledge coverage: L5 Algorithms, L8 Advanced Topics
 */

#ifndef GEOMETRIC_INTEGRATORS_H
#define GEOMETRIC_INTEGRATORS_H

#include "geo_optctrl_core.h"
#include "symplectic_geometry.h"

/* -----------------------------------------------------------------------------
 * Butcher tableau for Runge-Kutta methods
 * ----------------------------------------------------------------------------- */

typedef struct {
    int     stages;
    double  a[GEO_MAX_DIM][GEO_MAX_DIM];
    double  b[GEO_MAX_DIM];
    double  c[GEO_MAX_DIM];
} GeoButcherTableau;

/** Initialize Gauss-Legendre tableau (symplectic, A-stable). s-stage = order 2s. */
void geo_gauss_legendre_tableau(GeoButcherTableau *tableau, int stages);

/**
 * Symplectic RK step: advances (q,p) using a symplectic RK method.
 * Condition: b_i*a_ij + b_j*a_ji = b_i*b_j ensures preservation of omega.
 * Reference: Sanz-Serna (1988), Lasagni (1988)
 */
void geo_symplectic_rk_step(const GeoButcherTableau *tableau,
                            GeoHamiltonian H, void *ctx,
                            double *q, double *p, double h, int n,
                            double tol, int max_iter);

/* -----------------------------------------------------------------------------
 * Discrete Lagrangian and Variational Integrators
 * ----------------------------------------------------------------------------- */

/** Discrete Lagrangian L_d: Q x Q -> R */
typedef double (*GeoDiscreteLagrangian)(const double *qk, const double *qk1,
                                         double h, int n, void *ctx);

/** Midpoint discrete Lagrangian: L_d(q_k,q_{k+1}) = h*L((q_k+q_{k+1})/2, (q_{k+1}-q_k)/h) */
double geo_discrete_lagrangian_midpoint(GeoScalarFunction L, void *ctx,
                                         const double *qk, const double *qk1,
                                         double h, int n);

/**
 * Discrete Euler-Lagrange (DEL) step.
 * DEL equation: D2 L_d(q_{k-1}, q_k) + D1 L_d(q_k, q_{k+1}) = 0
 * Properties: symplectic, momentum-preserving via discrete Noether.
 * Reference: Marsden & West (2001) Acta Numerica
 */
int geo_del_step(GeoDiscreteLagrangian L_d, void *L_d_ctx,
                 const double *q_prev, const double *q_curr,
                 double *q_next, double h, int n,
                 double tol, int max_iter);

/**
 * Discrete Legendre transform: maps (q_k, q_{k+1}) to momenta.
 * Left: p_k = -D1 L_d(q_k, q_{k+1}), Right: p_{k+1} = D2 L_d(q_k, q_{k+1})
 */
void geo_discrete_legendre(const GeoDiscreteLagrangian L_d, void *ctx,
                           const double *qk, const double *qk1,
                           double *pk, double *pk1, double h, int n);

/* -----------------------------------------------------------------------------
 * Energy-preserving methods (L8 Advanced)
 * ----------------------------------------------------------------------------- */

/**
 * Average Vector Field (AVF): exactly preserves energy H.
 * x_{n+1} = x_n + h * int_0^1 X((1-tau)x_n + tau x_{n+1}) dtau
 * Theorem (Quispel & McLaren 2008): H(x_{n+1}) = H(x_n) exactly.
 */
int geo_avf_step(GeoHamiltonian H, void *ctx, double *x,
                 double h, int n, int n_quad,
                 double tol, int max_iter);

/**
 * Discrete gradient method with Itoh-Abe gradient.
 * Exact energy conservation via H(x')-H(x) = <grad_H(x,x'), x'-x>.
 */
int geo_discrete_gradient_step(GeoHamiltonian H, void *ctx, double *x,
                               double h, int n, double tol, int max_iter);

/* -----------------------------------------------------------------------------
 * Yoshida composition: higher-order symplectic from 2nd-order base
 * ----------------------------------------------------------------------------- */

/** Yoshida weights for order 4 or 6. w0=-2^{1/(2k+1)}/(2-2^{1/(2k+1)}), w1=1/(2-2^{1/(2k+1)}) */
void geo_yoshida_weights(int order, double *w, int *n_w);

/* -----------------------------------------------------------------------------
 * Symplectic Euler (1st order) for separable H = T(p) + V(q)
 * ----------------------------------------------------------------------------- */

/** Variant A: p'=p-h*dV(q), q'=q+h*dT(p'). Variant B: q'=q+h*dT(p), p'=p-h*dV(q'). */
void geo_symplectic_euler(char variant,
                          GeoScalarFunction T, void *t_ctx,
                          GeoScalarFunction V, void *v_ctx,
                          double *q, double *p, double h, int n);

/* -----------------------------------------------------------------------------
 * Discrete Mechanics and Optimal Control (DMOC) — L8
 * ----------------------------------------------------------------------------- */

/** Discrete dynamics f_d: (q_k, u_k) -> q_{k+1} */
typedef void (*GeoDiscreteDynamics)(const double *q, const double *u,
                                     double *q_next, int n, int m, void *ctx);

/**
 * DMOC solver: structure-preserving discrete optimal control.
 * Minimize sum_k L_d(q_k,u_k) + Phi(q_N) s.t. q_{k+1}=f_d(q_k,u_k).
 * Reference: Junge, Marsden, Ober-Blobaum (2005) CDC
 */
int geo_dmoc_solve(GeoDiscreteLagrangian L_d, void *L_ctx,
                   GeoDiscreteDynamics f_d, void *f_ctx,
                   const double *q0, const double *qT, int N,
                   double *u_guess, double *q_traj, double *u_traj,
                   int n, int m, double tol, int max_iter);

#endif /* GEOMETRIC_INTEGRATORS_H */
