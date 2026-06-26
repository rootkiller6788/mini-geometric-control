#ifndef PONTRYAGIN_MAXIMUM_H
#define PONTRYAGIN_MAXIMUM_H

#include "hamiltonian_control.h"

/*=============================================================================
 * Pontryagin Maximum Principle for optimal control.
 * Reference: Pontryagin et al. (1962), Liberzon §4
 *            Athans & Falb, Optimal Control
 *===========================================================================*/

/* L1: Optimal control problem data */
typedef struct {
    int     state_dim;     /* n: dimension of state x */
    int     control_dim;   /* m: dimension of control u */
    double  t0;            /* initial time */
    double  tf;            /* terminal time (fixed) */
    double *x0;            /* initial state [state_dim] */
} optimal_control_problem_t;

/* L1: Cost functional: phi(x(tf)) + integral_t0^tf L(x,u) dt */
typedef double (*terminal_cost_fn)(const double *x, int n, void *params);
typedef double (*running_cost_fn)(const double *x, const double *u,
                                   int nx, int nu, void *params);

/* L1: Dynamics: dx/dt = f(x, u) */
typedef void (*dynamics_fn)(const double *x, const double *u,
                             double *dx_dt, int nx, int nu, void *params);

/* L2: The control Hamiltonian (Pontryagin function):
 * H(x, u, lambda, lambda0) = lambda0 * L(x,u) + lambda^T f(x,u)
 * where lambda is the costate (adjoint) vector. */
typedef struct {
    dynamics_fn        dynamics;
    running_cost_fn    running_cost;
    terminal_cost_fn   terminal_cost;
    void              *params;
    int                state_dim;
    int                control_dim;
} pontryagin_hamiltonian_t;

/* L2: Compute the pre-Hamiltonian H(x, u, lambda, lambda0).
 * Reference: Liberzon §4.2.1 */
double pontryagin_prehamiltonian(const pontryagin_hamiltonian_t *ph,
                                  const double *x, const double *u,
                                  const double *lambda, double lambda0);

/* L2: Minimize pre-Hamiltonian w.r.t. u to get optimal control.
 * Returns optimal u* in u_opt. For unconstrained case:
 * dH/du = 0 -> lambda0*dL/du + lambda^T*df/du = 0
 * Uses gradient descent in control space. */
int pontryagin_minimize_control(const pontryagin_hamiltonian_t *ph,
                                 const double *x, const double *lambda,
                                 double lambda0,
                                 double *u_opt, double *H_min,
                                 int max_iter, double tol);

/* L3: Costate (adjoint) equations:
 * d(lambda)/dt = -dH/dx = -(lambda0*dL/dx + lambda^T*df/dx)
 * Computes the RHS at given (x, u, lambda). */
void pontryagin_adjoint_rhs(const pontryagin_hamiltonian_t *ph,
                             const double *x, const double *u,
                             const double *lambda, double lambda0,
                             double *dlambda_dt);

/* L3: Transversality condition: lambda(tf) = d(phi)/dx evaluated at x(tf).
 * Reference: Liberzon §4.3.1 */
void pontryagin_transversality(const pontryagin_hamiltonian_t *ph,
                                const double *x_tf,
                                double *lambda_tf, double *lambda0_out);

/* L4: Pontryagin Maximum Principle for free final time.
 * Additional condition: H(x*(tf), u*(tf), lambda*(tf), lambda0) = 0
 * Returns the terminal Hamiltonian value. */
double pontryagin_terminal_condition(const pontryagin_hamiltonian_t *ph,
                                      const double *x_tf, const double *lambda_tf,
                                      double lambda0);

/* L5: Shooting method for PMP two-point boundary value problem.
 * Guess lambda(t0), integrate forward dynamics+adjoint, iterate
 * to match terminal condition x(tf) = x_target and transversality.
 * Reference: Bryson & Ho §7.3 */
typedef struct {
    int     max_shooting_iter;
    int     integrator_steps;
    double  tol_state;       /* tolerance for state matching */
    double  tol_adjoint;     /* tolerance for adjoint mismatch */
    double  newton_damping;  /* damping factor for Newton update */
} shooting_params_t;

typedef struct {
    double *x_traj;       /* [steps][state_dim] */
    double *u_traj;       /* [steps][control_dim] */
    double *lambda_traj;  /* [steps][state_dim] */
    double  cost;
    int     converged;
    int     n_steps;
} shooting_result_t;

int pontryagin_shooting_method(const pontryagin_hamiltonian_t *ph,
                                const optimal_control_problem_t *prob,
                                const double *x_target,
                                const shooting_params_t *sparams,
                                shooting_result_t *result);

/* L5: Indirect multiple shooting.
 * Splits time interval into subintervals, matches continuity at nodes.
 * More robust than single shooting for unstable dynamics. */
int pontryagin_multiple_shooting(const pontryagin_hamiltonian_t *ph,
                                  const optimal_control_problem_t *prob,
                                  const double *x_target,
                                  int n_segments,
                                  const shooting_params_t *sparams,
                                  shooting_result_t **segment_results);

/* L6: LQR as special case of PMP.
 * A classic problem: linear dynamics, quadratic cost.
 * dx/dt = Ax + Bu, cost = x_f^T Q_f x_f + integral(x^T Q x + u^T R u)
 * PMP yields Riccati equation: dP/dt = -A^T P - PA + P B R^{-1} B^T P - Q
 * With u* = -R^{-1} B^T P x (linear state feedback).
 * Reference: Liberzon §6.1 */
typedef struct {
    int      n;
    int      m;
    double **A;     /* [n][n] */
    double **B;     /* [n][m] */
    double **Q;     /* [n][n] state weight, SPD */
    double **R;     /* [m][m] control weight, SPD */
    double **Qf;    /* [n][n] terminal state weight */
    double   tf;    /* final time */
} lqr_pmp_problem_t;

typedef struct {
    double  *P_traj;    /* [steps][n*n] Riccati matrix (flattened) */
    double  *K_traj;    /* [steps][m*n] optimal gain (flattened) */
    double  *x_traj;    /* [steps][n] */
    double  *u_traj;    /* [steps][m] */
    double   cost;
    int      n_steps;
} lqr_pmp_result_t;

int pontryagin_lqr_solve(const lqr_pmp_problem_t *prob,
                         const double *x0, int n_steps,
                         lqr_pmp_result_t *result);

/* Cleanup */
void shooting_result_free(shooting_result_t *r);

#endif /* PONTRYAGIN_MAXIMUM_H */
