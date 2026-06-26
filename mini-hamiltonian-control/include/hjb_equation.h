#ifndef HJB_EQUATION_H
#define HJB_EQUATION_H

#include "hamiltonian_control.h"

/*=============================================================================
 * Hamilton-Jacobi-Bellman equation for optimal control.
 * Reference: Fleming & Soner, Controlled Markov Processes (2006)
 *            Bardi & Capuzzo-Dolcetta, Optimal Control (1997)
 *            Evans, Partial Differential Equations §10.3
 *===========================================================================*/

/* L1: HJB equation for finite-horizon problem.
 * -dV/dt + H(x, -grad(V), t) = 0
 * V(x, T) = phi(x)  (terminal cost)
 * where H(x, p, t) = sup_{u in U} { -p^T f(x,u) - L(x,u) }
 * V(x,t) = cost-to-go from (x,t).
 * Reference: Liberzon §5.1.2 */

/* L2: Cost-to-go function approximator (value function). */
typedef struct {
    int      state_dim;
    double  *grid_min;     /* bounds for each state dimension */
    double  *grid_max;
    int     *grid_pts;     /* number of grid points per dimension */
    double  *values;       /* flat array: V(x_grid_points) */
    int      total_points; /* product of grid_pts */
} value_function_grid_t;

value_function_grid_t *value_function_grid_create(int state_dim,
                                                   const double *grid_min,
                                                   const double *grid_max,
                                                   const int *grid_pts);

void value_function_grid_free(value_function_grid_t *vf);

/* L2: Evaluate value function at arbitrary point via linear interpolation.
 * Handles points outside grid by extrapolation with warning. */
double value_function_eval(const value_function_grid_t *vf, const double *x);

/* L3: Gradient of value function via finite differences on grid.
 * grad(V)(x_i) is approximated by centered differences.
 * Returns gradient in grad_out. */
void value_function_gradient(const value_function_grid_t *vf,
                              const double *x, double *grad_out);

/* L4: HJB verification theorem.
 * If V solves HJB and is C^1, then V equals the optimal cost-to-go
 * and the optimal control is the argmin of the Hamiltonian.
 * This function checks the HJB residual at a point.
 * Reference: Liberzon §5.1.3 */
double hjb_residual(const value_function_grid_t *V,
                     const double *x, double t, double tf,
                     double (*dynamics)(const double*, const double*, int, void*),
                     double (*running_cost)(const double*, const double*, int, void*),
                     int control_dim,
                     void *params);

/* L5: Value iteration (dynamic programming) for discrete-time optimal control.
 * V_{k+1}(x) = min_u { L(x,u) + V_k(f(x,u)) }
 * Converges to optimal value function for infinite-horizon discounted problems.
 * Reference: Bertsekas, Dynamic Programming (2017) §1.3 */
typedef struct {
    int      state_dim;
    int      control_dim;
    int      n_states;      /* number of discrete states */
    double **states;        /* [n_states][state_dim] */
    double **controls;      /* [n_controls_discrete][control_dim] */
    int      n_controls;
    double   gamma;         /* discount factor */
    int      max_iter;
    double   tol;           /* convergence tolerance */
} value_iteration_params_t;

typedef struct {
    double  *V;             /* [n_states] value function */
    int     *policy;        /* [n_states] optimal control index */
    int      iterations;
    double   max_delta;     /* final Bellman error */
} value_iteration_result_t;

int hjb_value_iteration(
    const value_iteration_params_t *params,
    void (*dynamics)(const double *x, const double *u, double *x_next,
                     int nx, int nu, void *p),
    double (*cost)(const double *x, const double *u, int nx, int nu, void *p),
    void *cost_params,
    value_iteration_result_t *result);

/* L5: Policy iteration for continuous state-space HJB.
 * Alternates policy evaluation (solve linear PDE) and policy improvement.
 * Reference: Bertsekas §2.4 */
typedef struct {
    int      state_dim;
    int      control_dim;
    /* Current policy: u = pi(x, theta), parameterized by theta */
    void   (*policy)(const double *x, double *u, const double *theta,
                     int nx, int nu, void *params);
    double  *theta;          /* policy parameters */
    int      theta_dim;
    int      max_policy_iter;
    int      max_eval_iter;
    double   tol;
    void    *user_params;     /* user-defined parameters */
} policy_iteration_params_t;

int hjb_policy_iteration(
    const policy_iteration_params_t *params,
    void (*dynamics)(const double *x, const double *u, double *x_next,
                     int nx, int nu, void *p),
    double (*cost)(const double *x, const double *u, int nx, int nu, void *p),
    void *cost_params,
    value_function_grid_t *V);

/* L6: LQR solves HJB analytically.
 * For linear-quadratic problem, V(x) = x^T P x is exact solution.
 * Algebraic Riccati equation: A^T P + P A - P B R^{-1} B^T P + Q = 0.
 * Reference: Liberzon §6.1 */
int hjb_lqr_algebraic_riccati(double **A, double **B,
                               double **Q, double **R,
                               int n, int m,
                               double **P_out,
                               int max_iter, double tol);

/* L6: Minimum-time problem via HJB.
 * dynamics: dx/dt = f(x,u), cost = integral_0^T 1 dt = T.
 * HJB: sup_u { -grad(V)^T f(x,u) } = 1 (eikonal equation).
 * Solve using fast marching method on grid.
 * Reference: Bardi & Capuzzo-Dolcetta §4.2 */
int hjb_minimum_time(const value_function_grid_t *grid,
                      void (*dynamics)(const double *x, const double *u,
                                       double *dx_dt, int nx, int nu, void *p),
                      const double *control_set, int n_controls, int control_dim,
                      const double *target, void *params,
                      value_function_grid_t *V_out);

/* L8: Viscosity solutions for HJB.
 * Classical solutions may not exist (non-smooth V).
 * Viscosity solution concept: sub- and super-solution test.
 * Check if V is a viscosity subsolution/supersolution at x.
 * Reference: Crandall, Ishii, Lions (1992) */
int hjb_viscosity_subsolution_check(const value_function_grid_t *V,
                                     const double *x,
                                     double (*hamiltonian)(const double*, const double*, int, void*),
                                     void *params, double *residual);

int hjb_viscosity_supersolution_check(const value_function_grid_t *V,
                                       const double *x,
                                       double (*hamiltonian)(const double*, const double*, int, void*),
                                       void *params, double *residual);

/* L7: Application -- optimal investment/consumption (Merton problem).
 * State: wealth w, control: consumption c.
 * dw = (r w - c) dt, utility: integral e^{-rho t} U(c) dt.
 * HJB: rho V = max_c { U(c) + V'(w)(r w - c) }
 * For U(c) = log(c): optimal c* = rho w, V(w) = (log(rho w) + r/rho - 1)/rho.
 * Reference: Merton (1969) */
typedef struct {
    double risk_free_rate;
    double discount_rate;
} merton_params_t;

double merton_optimal_consumption(const merton_params_t *params,
                                   double wealth, double *V_value);

#endif /* HJB_EQUATION_H */
