#ifndef ENERGY_SHAPING_H
#define ENERGY_SHAPING_H

#include "hamiltonian_control.h"
#include "port_hamiltonian.h"

/*=============================================================================
 * Energy shaping methods for nonlinear control.
 * Reference: Ortega et al., "Interconnection and Damping Assignment
 *            Passivity-Based Control" (IDA-PBC), Automatica 2002
 *            Bloch et al., "Controlled Lagrangians" (2000)
 *            Nijmeijer & van der Schaft, Nonlinear Dynamical Control Systems
 *===========================================================================*/

/* L1: IDA-PBC problem formulation.
 * Given: dx/dt = [J(x) - R(x)] grad(H)(x) + g(x) u
 * Target: dx/dt = [J_d(x) - R_d(x)] grad(H_d)(x)
 * where J_d = -J_d^T, R_d = R_d^T >= 0 (desired structure)
 * and H_d has minimum at desired equilibrium x*.
 * Find u = beta(x) and H_d satisfying matching equations.
 * Reference: Ortega et al. (2002), §2 */
typedef struct {
    /* Desired interconnection matrix J_d(x) */
    void (*J_desired)(const double *x, int n, double **Jd_out, void *params);

    /* Desired damping matrix R_d(x) */
    void (*R_desired)(const double *x, int n, double **Rd_out, void *params);

    /* Desired Hamiltonian H_d(x) */
    hamiltonian_t  H_desired;

    /* Original system data (from port_hamiltonian.h) */
    port_hamiltonian_system_t *original_system;

    void  *params;
} ida_pbc_problem_t;

/* L2: Matching equations for IDA-PBC.
 * g_perp(x) * ([J - R] grad(H) - [J_d - R_d] grad(H_d)) = 0
 * where g_perp is a full-rank left annihilator of g(x): g_perp * g = 0.
 * Returns residual norm (should be ~0 for solvable cases). */
double ida_pbc_matching_residual(const ida_pbc_problem_t *prob,
                                  const double *x);

/* L3: Solve the matching PDE for H_d via homotopy.
 * For separable desired Hamiltonian H_d(q,p) = p^T M_d^{-1} p / 2 + V_d(q),
 * the matching conditions yield a PDE for V_d.
 * This is the Potential Energy Shaping step.
 * Reference: Ortega et al. (2002), §3.1 */
typedef struct {
    double *M_d;          /* Desired mass/inertia matrix [n][n] */
    double *M_d_inv;      /* Inverse of M_d [n][n] */
    double *K_p;          /* Potential shaping gain [n][n] */
    double *V_d_coeffs;   /* Coefficients for desired potential (polynomial) */
    int     poly_degree;   /* Degree of polynomial potential */
    int     n_dof;
} potential_shaping_t;

int potential_energy_shaping(const ida_pbc_problem_t *prob,
                              const potential_shaping_t *ps,
                              const double *x_eq,
                              double *V_d_at_eq);

/* L5: Kinetic energy shaping (mass matrix modification).
 * For mechanical systems: H = p^T M^{-1}(q) p / 2 + V(q).
 * IDA-PBC allows modifying both M(q) and V(q).
 * The matching PDEs for M_d(q) involve Christoffel symbols.
 * Reference: Blankenstein et al. (2002) */
typedef struct {
    double  **M_d;       /* Desired mass matrix */
    int       n_dof;
} kinetic_shaping_t;

int kinetic_energy_shaping_setup(const ida_pbc_problem_t *prob,
                                  const double *q,
                                  kinetic_shaping_t *ks);

/* L5: Damping assignment.
 * Given existing dissipation R(x) and desired R_d(x) = R(x) + R_a(x),
 * choose R_a = g K_v g^T for velocity feedback.
 * Reference: Ortega et al. (2002), §3.2 */
int damping_assignment(const port_hamiltonian_system_t *ph,
                        const double *x, double K_v,
                        double **R_a_out);

/* L5: Total energy shaping controller.
 * Combines potential shaping + kinetic shaping + damping injection.
 * u = u_es + u_di where:
 *   u_es: energy shaping term (from PDE solution)
 *   u_di: damping injection term (-K_d y)
 * Returns the control input u. */
int total_energy_shaping_control(const ida_pbc_problem_t *prob,
                                  const double *x,
                                  const potential_shaping_t *ps,
                                  const kinetic_shaping_t *ks,
                                  double K_d,
                                  double *u_out);

/* L6: Simple pendulum energy shaping.
 * Model: ml^2 ddot(theta) + mgl sin(theta) = u
 * Port-Hamiltonian: q = theta, p = ml^2 dot(theta), H = p^2/(2ml^2) + mgl(1-cos(q))
 * IDA-PBC: set H_d = p^2/(2ml^2 k_i) + V_d(q) with V_d minimum at desired q*.
 * k_i > 0 tunes the inertia, and V_d shapes the potential well. */
typedef struct {
    double mass;
    double length;
    double gravity;
    double k_i;          /* inertia scaling */
    double k_p;          /* potential shaping gain */
    double k_d;          /* damping injection gain */
    double q_desired;    /* desired equilibrium angle */
} pendulum_ida_params_t;

int pendulum_ida_pbc_control(const pendulum_ida_params_t *params,
                              const double *x,
                              double *u_out, double *H_d_val);

/* L6: Cart-pole (inverted pendulum) energy shaping.
 * Underactuated system: 2 DOF (cart position x_c, pendulum angle theta),
 * 1 control input (cart force F).
 * H = p^T M^{-1}(q) p / 2 + V(q) with nontrivial mass matrix.
 * IDA-PBC shapes total energy to stabilize upright position.
 * Reference: Acosta et al. (2005) */
typedef struct {
    double M_cart;       /* cart mass */
    double m_pole;       /* pole mass */
    double l_pole;       /* pole half-length */
    double gravity;
    double k_e;          /* energy shaping gain */
    double k_d;          /* damping gain */
} cartpole_ida_params_t;

int cartpole_ida_pbc_control(const cartpole_ida_params_t *params,
                              const double *q, const double *p,
                              double *u_out);

/* L8: Controlled Lagrangian method.
 * Alternative to IDA-PBC: modify the Lagrangian L_d = L + L_c
 * such that the Euler-Lagrange equations have desired equilibria.
 * L_c is designed to achieve stabilization for underactuated systems.
 * Reference: Bloch, Leonard, Marsden (2000) */
typedef struct {
    void (*L)(const double *q, const double *dq, int n,
              double *L_val, double *grad_q, double *grad_dq, void *params);
    void (*L_control)(const double *q, const double *dq, int n,
                      double *Lc_val, double *grad_q, double *grad_dq, void *params);
    int     n_dof;
    int     n_controls;
    void   *params;
} controlled_lagrangian_t;

int controlled_lagrangian_derive(const controlled_lagrangian_t *cl,
                                  const double *q, const double *dq,
                                  double *u_ctrl, double *E_total);

#endif /* ENERGY_SHAPING_H */
