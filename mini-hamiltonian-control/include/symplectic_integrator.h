#ifndef SYMPLECTIC_INTEGRATOR_H
#define SYMPLECTIC_INTEGRATOR_H

#include "hamiltonian_control.h"

/*=============================================================================
 * Symplectic (geometric) integrators for Hamiltonian systems.
 * Preserve symplectic structure, energy near-conservation, and momentum maps.
 * Reference: Hairer, Lubich, Wanner, Geometric Numerical Integration (2006)
 *            Leimkuhler & Reich, Simulating Hamiltonian Dynamics (2004)
 *            Sanz-Serna & Calvo, Numerical Hamiltonian Problems (1994)
 *===========================================================================*/

/* L1: Symplectic integrator types */
typedef enum {
    INT_SYMPLECTIC_EULER_A,   /* Explicit: q_new = q + dt*dH/dp(q_new, p) */
    INT_SYMPLECTIC_EULER_B,   /* Explicit: p_new = p - dt*dH/dq(q, p_new) */
    INT_STORMER_VERLET,       /* q_{1/2}=q+dt/2*dH/dp(p), p'=p-dt*dH/dq(q_{1/2}), ... */
    INT_VELOCITY_VERLET,      /* p_{1/2}=p-dt/2*dH/dq(q), q'=q+dt*dH/dp(p_{1/2}), ... */
    INT_IMPLICIT_MIDPOINT,    /* z_{new}=z+dt*J*grad(H)((z+z_new)/2) */
    INT_YOSHIDA_4,            /* Yoshida 4th order composition method */
    INT_RUTH_3,               /* Ruth 3rd order symplectic */
    INT_GAUSS_LEGENDRE_2      /* 2-stage Gauss-Legendre (implicit midpoint) */
} symplectic_integrator_type_t;

/* L1: Integrator step result */
typedef struct {
    double  energy_error;     /* relative energy change |H_new - H_old| / |H_old| */
    double  symplectic_error; /* ||M^T J M - J||_F for linear flow map M */
    double  momentum_error;   /* linear/angular momentum drift */
} integrator_diagnostics_t;

/* L2: Standard symplectic Euler (method A).
 * q_{n+1} = q_n + h * dH/dp(q_{n+1}, p_n)  -- implicit in q
 * p_{n+1} = p_n - h * dH/dq(q_{n+1}, p_n)
 * First-order, symplectic. Solves q via Newton iteration. */
int step_symplectic_euler_a(const hamiltonian_t *H,
                             double *q, double *p, double dt,
                             integrator_diagnostics_t *diag);

/* L2: Symplectic Euler (method B).
 * p_{n+1} = p_n - h * dH/dq(q_n, p_{n+1})  -- implicit in p
 * q_{n+1} = q_n + h * dH/dp(q_n, p_{n+1})
 * First-order, symplectic. */
int step_symplectic_euler_b(const hamiltonian_t *H,
                             double *q, double *p, double dt,
                             integrator_diagnostics_t *diag);

/* L3: Stormer-Verlet integrator.
 * p_{n+1/2} = p_n - (h/2) * dH/dq(q_n, p_{n+1/2})
 * q_{n+1}   = q_n + h * dH/dp(q_n, p_{n+1/2}) + dH/dp(q_{n+1}, p_{n+1/2}))
 * p_{n+1}   = p_{n+1/2} - (h/2) * dH/dq(q_{n+1}, p_{n+1/2})
 * Second-order, symplectic, explicit for separable H.
 * Reference: Hairer et al. §1.3 */
int step_stormer_verlet(const hamiltonian_t *H,
                         double *q, double *p, double dt,
                         integrator_diagnostics_t *diag);

/* L3: Velocity Verlet integrator.
 * For separable H = T(p) + V(q):
 * p_{n+1/2} = p_n - (h/2) * dV/dq(q_n)
 * q_{n+1}   = q_n + h * dT/dp(p_{n+1/2})
 * p_{n+1}   = p_{n+1/2} - (h/2) * dV/dq(q_{n+1})
 * Reference: Swope et al. (1982) */
int step_velocity_verlet(const hamiltonian_t *H,
                          double *q, double *p, double dt,
                          integrator_diagnostics_t *diag);

/* L4: Implicit midpoint rule.
 * z_{n+1} = z_n + h * J * grad(H)((z_n + z_{n+1})/2)
 * Second-order, symplectic for all Hamiltonians.
 * Solves implicit equation via fixed-point iteration.
 * Reference: Hairer et al. §2.4 */
int step_implicit_midpoint(const hamiltonian_t *H,
                            double *q, double *p, double dt,
                            int max_iter, double tol,
                            integrator_diagnostics_t *diag);

/* L5: Yoshida 4th order composition method.
 * Uses symmetric composition of second-order symplectic steps.
 * Coefficients: w1 = w4 = 1/(2-2^(1/3)), w2 = w3 = -2^(1/3)/(2-2^(1/3))
 * S_4(h) = S_2(w1*h) o S_2(w2*h) o S_2(w3*h) o S_2(w4*h)
 * Reference: Yoshida (1990), Hairer et al. §5.2 */
int step_yoshida_4(const hamiltonian_t *H,
                    double *q, double *p, double dt,
                    integrator_diagnostics_t *diag);

/* L5: Adaptive symplectic integration.
 * Adjusts step size based on energy error estimate.
 * Uses Stormer-Verlet with step doubling for error estimation.
 * Reference: Hairer et al. §8.2 */
typedef struct {
    double  dt_min;
    double  dt_max;
    double  tol;
    double  safety_factor;
} adaptive_params_t;

int adaptive_symplectic_integrate(const hamiltonian_t *H,
                                   double *q, double *p,
                                   double t0, double tf,
                                   const adaptive_params_t *params,
                                   int *n_steps_taken,
                                   integrator_diagnostics_t *final_diag);

/* L5: Variational integrator (discrete mechanics).
 * Derived from discrete Hamilton's principle:
 * delta sum_{k=0}^{N-1} L_d(q_k, q_{k+1}, h) = 0
 * where L_d is a discrete Lagrangian approximating the action integral.
 * Uses midpoint discrete Lagrangian: L_d(q0,q1,h) = h * L((q0+q1)/2, (q1-q0)/h)
 * Reference: Marsden & West (2001), Hairer et al. §6.2 */
typedef struct {
    /* Discrete Lagrangian L_d(q_k, q_{k+1}, h) */
    double (*L_d)(const double *qk, const double *qkp1, double h,
                  int n, void *params);
    /* D1 L_d (partial derivative w.r.t. first argument q_k) */
    void   (*D1_L_d)(const double *qk, const double *qkp1, double h,
                     int n, double *D1, void *params);
    /* D2 L_d (partial derivative w.r.t. second argument q_{k+1}) */
    void   (*D2_L_d)(const double *qk, const double *qkp1, double h,
                     int n, double *D2, void *params);
    void   *params;
    int     n_dof;
} discrete_lagrangian_t;

int variational_integrator_step(const discrete_lagrangian_t *L_d,
                                 double *q_prev, double *q_curr,
                                 double *q_next, double *p_curr,
                                 double dt);

/* L6: Symplectic integrator for separable Hamiltonian.
 * Generic dispatcher: selects appropriate integrator type. */
int symplectic_integrate(const hamiltonian_t *H,
                          double *q, double *p,
                          double t0, double tf, double dt,
                          symplectic_integrator_type_t method,
                          int *n_steps,
                          integrator_diagnostics_t *final_diag);

/* L6: Long-term energy behavior analysis.
 * For symplectic integrators, energy error remains bounded (no drift).
 * Compare with non-symplectic (RK4) which shows linear energy drift.
 * Returns array of energy values at each step. */
double *energy_trajectory(const hamiltonian_t *H,
                           double *q0, double *p0,
                           double dt, int n_steps,
                           symplectic_integrator_type_t method);

/* L8: Symplectic partitioned Runge-Kutta (PRK) coefficients.
 * For separable H = T(p) + V(q), PRK methods decouple q,p updates.
 * Table of classical coefficient sets. */
typedef struct {
    int     stages;
    double *a;    /* [stages][stages] RK matrix */
    double *b;    /* [stages] weights */
    double *b_tilde; /* [stages] weights for p-update */
} prk_coefficients_t;

prk_coefficients_t *prk_lobatto_iiia_coefficients(int stages);

#endif /* SYMPLECTIC_INTEGRATOR_H */
