#ifndef LIE_DERIVATIVE_H
#define LIE_DERIVATIVE_H

#include "nonlinear_system.h"

/* ============================================================================
 * Lie Derivative Computation
 *
 * The Lie derivative of a scalar field h along a vector field f is the
 * directional derivative of h in the direction of f:
 *
 *   L_f h(x) = sum_{i=1}^n (partial h / partial x_i) * f_i(x)
 *            = (grad h)^T ? f(x)
 *            = dh ? f
 *
 * This measures the rate of change of h along trajectories of f.
 *
 * Higher-order Lie derivatives:
 *   L_f^2 h = L_f (L_f h)
 *   L_f^k h = L_f (L_f^{k-1} h)
 *
 * Mixed Lie derivatives (crucial for feedback linearization):
 *   L_g L_f h = L_g (L_f h)
 *   L_g L_f^k h = L_g (L_f^k h)
 *
 * References:
 *   Isidori (1995) Nonlinear Control Systems, Ch. 4
 *   Nijmeijer & van der Schaft (1990) Nonlinear Dynamical Control Systems
 *   Abraham, Marsden, Ratiu (1988) Manifolds, Tensor Analysis, and Applications
 * ============================================================================ */

/**
 * lie_derivative -- Compute L_f h(x).
 *
 * L_f h(x) = grad h(x) ? f(x)
 *
 * For systems with explicit gradient callback, uses grad. Otherwise falls
 * back to finite difference approximation:
 *   L_f h ? [h(x + ? f) - h(x)] / ?
 *
 * @param f  vector field
 * @param h  scalar field
 * @param x  evaluation point
 * @return   L_f h(x)
 */
double lie_derivative(const VectorField *f,
                      const ScalarField *h,
                      const Vector *x);

/**
 * lie_derivative_numeric -- Compute L_f h(x) using finite differences.
 *
 * Uses central difference for accuracy:
 *   L_f h ? [h(x + ? f(x)) - h(x - ? f(x))] / (2?)
 *
 * @param f  vector field
 * @param h  scalar field
 * @param x  evaluation point
 * @param h_step  finite difference step (suggest: 1e-6)
 * @return   L_f h(x)
 */
double lie_derivative_numeric(const VectorField *f,
                               const ScalarField *h,
                               const Vector *x,
                               double h_step);

/**
 * lie_derivative_kth -- Compute k-th order Lie derivative L_f^k h(x).
 *
 * L_f^0 h = h
 * L_f^1 h = L_f h
 * L_f^k h = L_f (L_f^{k-1} h)  for k >= 2
 *
 * Uses repeated application of the gradient of the previous Lie derivative.
 * Falls back to finite differences when gradient callback is unavailable.
 *
 * @param f  vector field
 * @param h  scalar field
 * @param x  evaluation point
 * @param k  order
 * @return   L_f^k h(x)
 */
double lie_derivative_kth(const VectorField *f,
                           const ScalarField *h,
                           const Vector *x,
                           int k);

/**
 * lie_derivative_mixed -- Compute L_g L_f^k h(x).
 *
 * This is the key computation for determining relative degree.
 * First computes L_f^k h, then takes its Lie derivative along g.
 *
 * Theorem (Isidori, 1995): A SISO system has relative degree rho iff
 *   L_g L_f^k h(x) = 0   for all k < rho-1
 *   L_g L_f^{rho-1} h(x) != 0
 *
 * @param g  vector field (control direction)
 * @param f  vector field (drift)
 * @param h  scalar field (output)
 * @param x  evaluation point
 * @param k  power of L_f before applying L_g
 * @return   L_g L_f^k h(x)
 */
double lie_derivative_mixed(const VectorField *g,
                             const VectorField *f,
                             const ScalarField *h,
                             const Vector *x,
                             int k);

/**
 * lie_derivative_vector -- Compute Lie derivative of a vector field Y along X.
 *
 * (L_X Y)(x) = limit_{t->0} [Y(x + t X) - Y(x)] / t
 *
 * This is the standard Lie derivative for vector fields, useful in
 * analyzing symmetries and invariance.
 *
 * @param X  direction vector field
 * @param Y  target vector field
 * @param x  evaluation point
 * @param out L_X Y at x
 */
void lie_derivative_vector(const VectorField *X,
                            const VectorField *Y,
                            const Vector *x,
                            Vector *out);

/**
 * lie_derivative_scalar_along_flow -- Compute L_f h(x) by following the flow.
 *
 * d/dt h(phi_t(x)) |_{t=0} = L_f h(x)
 *
 * where phi_t is the flow of f. This uses numerical integration to
 * compute the flow and then finite-differences the derivative.
 *
 * @param f  vector field
 * @param h  scalar field
 * @param x  initial point
 * @param dt  time step for flow integration
 * @return   time derivative of h along flow at t=0
 */
double lie_derivative_scalar_along_flow(const VectorField *f,
                                         const ScalarField *h,
                                         const Vector *x,
                                         double dt);

/**
 * check_relative_degree_order -- Verify L_g L_f^k h = 0 for k < order-1.
 *
 * Used in relative degree computation to check each order condition.
 * Returns true if the mixed derivative is zero within tolerance.
 *
 * @param g    control vector field
 * @param f    drift vector field
 * @param h    output map
 * @param x    evaluation point
 * @param k    power of L_f
 * @param tol  numerical tolerance
 * @return     true if derivative is zero within tolerance
 */
bool check_relative_degree_order(const VectorField *g,
                                  const VectorField *f,
                                  const ScalarField *h,
                                  const Vector *x,
                                  int k,
                                  double tol);

/**
 * compute_all_lie_derivatives -- Precompute L_f^k h(x) for k = 0,...,max_k.
 *
 * Returns array of length max_k+1 where result[k] = L_f^k h(x).
 * Caller must free the returned array.
 *
 * @param f      drift vector field
 * @param h      output scalar field
 * @param x      evaluation point
 * @param max_k  maximum order
 * @return       array of Lie derivative values
 */
double *compute_all_lie_derivatives(const VectorField *f,
                                     const ScalarField *h,
                                     const Vector *x,
                                     int max_k);

#endif /* LIE_DERIVATIVE_H */
