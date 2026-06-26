#ifndef PORT_HAMILTONIAN_H
#define PORT_HAMILTONIAN_H

#include "hamiltonian_control.h"

/*=============================================================================
 * Port-Hamiltonian systems for modeling and control of physical systems.
 * Reference: van der Schaft & Jeltsema, Port-Hamiltonian Systems Theory (2014)
 *            Duindam et al., Modeling and Control of Complex Physical Systems
 *            Ortega et al., Passivity-based Control (2002)
 *===========================================================================*/

/* L1: Port-Hamiltonian system definition.
 * dx/dt = [J(x) - R(x)] grad(H)(x) + g(x) u
 * y     = g(x)^T grad(H)(x)
 * J = -J^T (interconnection/skew-symmetric), R = R^T >= 0 (dissipation)
 * u: input port variables, y: output port variables.
 * Power balance: dH/dt = y^T u - grad(H)^T R grad(H) <= y^T u (passivity) */
typedef struct {
    int      state_dim;     /* n: number of energy storage elements */
    int      port_dim;      /* m: number of ports */

    /* J(x): interconnection matrix [n][n], skew-symmetric */
    void   (*J_matrix)(const double *x, int n, double **J_out, void *params);

    /* R(x): dissipation matrix [n][n], symmetric positive semidefinite */
    void   (*R_matrix)(const double *x, int n, double **R_out, void *params);

    /* g(x): input matrix [n][m] */
    void   (*g_matrix)(const double *x, int n, int m,
                       double **g_out, void *params);

    /* Hamiltonian energy function */
    hamiltonian_t  H;

    void   *params;
} port_hamiltonian_system_t;

/* L1: Compute the Dirac structure underlying a port-Hamiltonian system.
 * D = {(f, e, f_R, e_R, f_P, e_P)} where:
 *   f = -dx/dt (flow of energy storage)
 *   e = grad(H) (effort of energy storage)
 *   f_R, e_R: resistive port (R)
 *   f_P, e_P: external port (u, y)
 * Power balance: e^T f + e_R^T f_R + e_P^T f_P = 0
 * Reference: van der Schaft §2.2, Courant (1990) */
typedef struct {
    int     dim_flow;       /* total dimension of flow variables */
    int     dim_storage;
    int     dim_resistive;
    int     dim_port;
    double *flow;           /* concatenated [f_storage, f_resistive, f_port] */
    double *effort;         /* concatenated [e_storage, e_resistive, e_port] */
} dirac_structure_t;

int dirac_structure_from_ph(const port_hamiltonian_system_t *ph,
                             const double *x, const double *u,
                             dirac_structure_t *ds);

/* L2: Passivity of port-Hamiltonian systems.
 * dH/dt = y^T u - grad(H)^T R grad(H) <= y^T u (passive with storage H).
 * Returns dH/dt computed from power balance. */
double port_hamiltonian_power_balance(const port_hamiltonian_system_t *ph,
                                       const double *x, const double *u,
                                       double *y_out);

/* L2: Casimir functions for port-Hamiltonian systems.
 * C(x) is a Casimir if grad(C)^T J(x) = 0 and grad(C)^T R(x) = 0.
 * Then dC/dt = 0 (conserved quantity).
 * Reference: van der Schaft §4.2 */
typedef int (*casimir_check_fn)(const double *x, int n,
                                 double *C_val, double *grad_C, void *params);

int port_hamiltonian_casimir_check(const port_hamiltonian_system_t *ph,
                                    casimir_check_fn C_fn, void *C_params,
                                    const double *x, double *violation);

/* L3: Interconnection of port-Hamiltonian systems.
 * Given two pH systems, compute their power-conserving interconnection
 * via a Dirac structure. Typical: u1 = -y2, u2 = y1 (feedback).
 * Reference: van der Schaft §6.1 */
typedef struct {
    port_hamiltonian_system_t *sys1;
    port_hamiltonian_system_t *sys2;
    /* Interconnection matrix: [u1; u2] = F * [y1; y2] + [v1; v2] */
    double **F;        /* [(m1+m2)][(m1+m2)] */
    double  *v_ext;    /* [(m1+m2)] external inputs */
} ph_interconnection_t;

int port_hamiltonian_interconnect(const ph_interconnection_t *interconn,
                                   const double *x1, const double *x2,
                                   double *dx1_dt, double *dx2_dt,
                                   double *y1, double *y2);

/* L5: Passivity-Based Control (PBC) -- basic damping injection.
 * u = -K_d y  with K_d > 0 adds dissipation.
 * dH/dt = -y^T K_d y - grad(H)^T R grad(H) <= 0 (Lyapunov stable).
 * Reference: Ortega et al. (2002) §3.1 */
int pbc_damping_injection(const port_hamiltonian_system_t *ph,
                           const double *x, double K_d,
                           double *u_out, double *y_out);

/* L5: Energy-Casimir method for stability analysis.
 * For pH systems, V(x) = H(x) + sum_i Phi_i(C_i(x)) is a Lyapunov function
 * if Phi_i are chosen such that V has a minimum at the equilibrium.
 * Find Phi_i coefficients to make grad(V)(x_eq) = 0.
 * Reference: van der Schaft §4.3, Marsden-Ratiu §10.4 */
typedef struct {
    int      n_casimirs;
    double  *phi_coeffs;  /* coefficients for Phi(C) = sum a_k C^k */
    double  *C_values;    /* Casimir values at equilibrium */
} energy_casimir_t;

int energy_casimir_stability(const port_hamiltonian_system_t *ph,
                              casimir_check_fn *casimirs,
                              void **casimir_params,
                              int n_casimirs,
                              const double *x_eq,
                              energy_casimir_t *result);

/* L5: Port-Hamiltonian simulation with symplectic integration.
 * Split into J-part (energy conserving) and R-part (dissipative).
 * Uses splitting method: flow_J(dt/2) o flow_R(dt) o flow_J(dt/2). */
int ph_symplectic_step(const port_hamiltonian_system_t *ph,
                        const double *u, double *x, double dt,
                        double *y_out);

/* L7: DC motor as port-Hamiltonian system.
 * State: x = [phi (magnetic flux), p (angular momentum)]
 * H = phi^2/(2L) + p^2/(2J), J = [0 1; -1 0], R = diag(R, b)
 * u = [V, tau_load]^T, y = [i, omega]^T */
typedef struct {
    double L;      /* inductance (H) */
    double R;      /* resistance (Ohm) */
    double J_m;    /* moment of inertia (kg.m^2) */
    double b;      /* viscous friction (N.m.s/rad) */
    double K;      /* torque/back-emf constant (N.m/A = V.s/rad) */
} dc_motor_params_t;

port_hamiltonian_system_t *dc_motor_create(const dc_motor_params_t *params);

/* L7: Mass-spring-damper as port-Hamiltonian system.
 * H = p^2/(2m) + k q^2/2, state = (q, p)
 * u = F_ext, y = v = p/m */
typedef struct {
    double mass;
    double stiffness;
    double damping;
} msd_params_t;

port_hamiltonian_system_t *mass_spring_damper_create(const msd_params_t *params);

/* Cleanup */
void port_hamiltonian_system_free(port_hamiltonian_system_t *ph);

#endif /* PORT_HAMILTONIAN_H */
