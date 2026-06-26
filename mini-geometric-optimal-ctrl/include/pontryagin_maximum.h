/**
 * pontryagin_maximum.h ? Pontryagin Maximum Principle for Geometric Optimal Control
 *
 * References:
 *   Pontryagin, Boltyanskii, Gamkrelidze, Mishchenko (1962) "The Mathematical
 *     Theory of Optimal Processes"
 *   Agrachev & Sachkov (2004) "Control Theory from the Geometric Viewpoint"
 *   Bressan & Piccoli (2007) "Introduction to the Mathematical Theory of Control"
 *
 * Covers:
 *   - PMP necessary conditions (state-costate dynamics)
 *   - Maximized Hamiltonian (pseudo-Hamiltonian) condition
 *   - Transversality conditions for boundary constraints
 *   - Shooting methods for PMP boundary value problems
 *   - Conjugate point analysis and second-order conditions
 *
 * Knowledge coverage: L4 Fundamental Laws, L5 Algorithms
 */

#ifndef PONTRYAGIN_MAXIMUM_H
#define PONTRYAGIN_MAXIMUM_H

#include "geo_optctrl_core.h"
#include "symplectic_geometry.h"

/* -----------------------------------------------------------------------------
 * L4: Pontryagin's Maximum Principle ? data structures
 * ----------------------------------------------------------------------------- */

/**
 * Pseudo-Hamiltonian (control Hamiltonian) for PMP:
 *
 *   H(x, p, u) = p??L(x,u) + ?p, f?(x) + ? f_a(x) u^a?
 *
 * where:
 *   p? ? {0, -1}   abnormal (p?=0) or normal (p?=-1) case
 *   p ? T*_x M   costate (adjoint) vector
 *   L(x,u)       running cost
 *   f?, f_a      control-affine vector fields
 *
 * PMP Theorem (normal case, p?=-1): If u*(t) is optimal, then ? p(t) ? 0 s.t.:
 *   1. ? = ?H/?p = f?(x) + ? f_a(x) u^a    (state equation)
 *   2. ? = -?H/?x                             (costate equation)
 *   3. H(x*(t), p(t), u*(t)) ? H(x*(t), p(t), u)
 *      for all u ? U, almost everywhere t     (maximization condition)
 *   4. H(x*(T), p(T), u*(T)) = 0 if T is free (transversality)
 */
typedef struct {
    GeoOptimalControlProblem *problem;   /**< Optimal control problem */
    double                    p0;        /**< p? = 0 (abnormal) or -1 (normal) */
    double                   *state;     /**< Current state x (length n) */
    double                   *costate;   /**< Current costate p (length n) */
    int                       n;         /**< State dimension */
    int                       m;         /**< Control dimension */
} GeoPMPState;

/**
 * Evaluate the pseudo-Hamiltonian H(x, p, u) at given state, costate, control.
 *
 * @param pmp    PMP state (x, p, p?)
 * @param u      control input (length m)
 * @return       H(x, p, u) value
 */
double geo_pmp_hamiltonian(const GeoPMPState *pmp, const double *u);

/**
 * Compute the state-costate dynamics (?, ?) given the optimal control u*.
 *
 * The optimal control u* is obtained by maximizing the Hamiltonian:
 *   u* = argmax_{u?U} H(x, p, u)
 *
 * For control-affine systems with running cost L(x,u) = ? (u^a)?/2
 * (minimum-energy cost) and no control constraints, the Hamiltonian
 * maximization condition gives:
 *   u^a = -?p, f_a(x)?  (if in interior of U)
 *
 * @param pmp       PMP state
 * @param optimal_u computed optimal control (output, length m)
 * @param dx        output state derivative ? (length n)
 * @param dp        output costate derivative ? (length n)
 */
void geo_pmp_dynamics(GeoPMPState *pmp, double *optimal_u,
                      double *dx, double *dp);

/* -----------------------------------------------------------------------------
 * L4: Transversality conditions
 * ----------------------------------------------------------------------------- */

/**
 * Transversality condition for free final state.
 *
 * If x(T) is free, then p(T) = -??(x(T)) where ? is the terminal cost.
 *
 * @param xT       final state
 * @param terminal terminal cost function
 * @param pT       output terminal costate (length n)
 * @param n        state dimension
 * @param h        finite difference step
 * @param ctx      context
 */
void geo_pmp_transversality_free(const double *xT,
                                 GeoTerminalCost terminal,
                                 void *ctx,
                                 double *pT, int n, double h);

/**
 * Transversality condition for fixed final state.
 *
 * If x(T) = x_fixed, then p(T) is free (determined by the shooting condition).
 * No constraint is placed on p(T); instead x(T) = x_fixed is enforced.
 */

/**
 * Transversality condition for free final time.
 *
 * If T is free (optimized), then the Hamiltonian vanishes at final time:
 *   H(x*(T), p(T), u*(T)) = 0
 *
 * This provides the additional equation needed to determine T.
 */
double geo_pmp_free_time_condition(const GeoPMPState *pmp, const double *u_opt);

/* -----------------------------------------------------------------------------
 * L5: Shooting methods for PMP
 * ----------------------------------------------------------------------------- */

/**
 * Single shooting for PMP boundary value problem.
 *
 * Solves: find p(0) (and possibly T) such that integrating PMP dynamics
 * from (x(0), p(0)) satisfies terminal conditions at t=T.
 *
 * The unknown variables are:
 *   p(0) ? ??            (initial costate)
 *   T (if free time)      (time horizon)
 *
 * The residuals are:
 *   x(T) - x_target       (state terminal condition)
 *   p(T) + ??(x(T))       (costate transversality, free final state)
 *   H(x(T), p(T), u*(T))  (free time condition)
 *
 * Uses Newton iteration with finite-difference Jacobian.
 *
 * @param problem    optimal control problem definition
 * @param p0_guess   initial guess for costate p(0) (length n, modified in place)
 * @param T_guess    initial guess for T (modified in place if free time)
 * @param n_steps    number of integration steps per shooting iteration
 * @param max_iter   maximum Newton iterations
 * @param tol        convergence tolerance
 * @param n          state dimension
 * @return           0 on convergence, negative on failure
 */
int geo_pmp_single_shooting(GeoOptimalControlProblem *problem,
                            double *p0_guess, double *T_guess,
                            int n_steps, int max_iter, double tol, int n);

/**
 * Multiple shooting for PMP: partitions time interval [0,T] into sub-intervals
 * to improve numerical stability for long horizons.
 *
 * Each sub-interval has its own initial state and costate guess,
 * with continuity conditions enforced at interval boundaries.
 *
 * @param problem         optimal control problem
 * @param n_subintervals  number of shooting segments
 * @param states_guess    n_subintervals ? n array (state at each node)
 * @param costates_guess  n_subintervals ? n array (costate at each node)
 * @param T               time horizon
 * @param n_steps_per     integration steps per sub-interval
 * @param max_iter        maximum Newton iterations
 * @param tol             convergence tolerance
 * @param n               state dimension
 * @return                0 on convergence
 */
int geo_pmp_multiple_shooting(GeoOptimalControlProblem *problem,
                              int n_subintervals,
                              double *states_guess,
                              double *costates_guess,
                              double T, int n_steps_per,
                              int max_iter, double tol, int n);

/* -----------------------------------------------------------------------------
 * L5: Conjugate point test (second-order optimality)
 * ----------------------------------------------------------------------------- */

/**
 * Compute the Jacobi equation along an extremal to detect conjugate points.
 *
 * The Jacobi equation is:
 *   ?? = A(t)??x + B(t)??p
 *   ?? = -C(t)??x - A(t)???p
 *
 * derived from linearizing the PMP dynamics. A conjugate point occurs
 * when a nontrivial solution satisfies ?x(0)=0 and ?x(t_c)=0.
 * This signals loss of local optimality.
 *
 * @param pmp    PMP state along the extremal
 * @param A      output A(t) (n?n matrix, linearized drift)
 * @param B      output B(t) (n?n, control effect on costate)
 * @param C      output C(t) (n?n, Hessian of Hamiltonian)
 * @param n      state dimension
 * @param h      finite difference step
 */
void geo_pmp_jacobi_matrices(const GeoPMPState *pmp,
                             const double *u_opt,
                             double *A, double *B, double *C,
                             int n, double h);

/**
 * Detect conjugate point by propagating the Jacobi equation along an extremal
 * and checking when ?x(t) crosses zero.
 *
 * @param problem       optimal control problem
 * @param p0            initial costate
 * @param T             time horizon
 * @param n_steps       integration steps
 * @param n             state dimension
 * @param conjugate_t   output: time of first conjugate point (>0), or -1 if none
 * @return              1 if conjugate point found before T, 0 if none
 */
int geo_pmp_conjugate_point_test(GeoOptimalControlProblem *problem,
                                 const double *p0, double T,
                                 int n_steps, int n,
                                 double *conjugate_t);

#endif /* PONTRYAGIN_MAXIMUM_H */
