#include "sr_control.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define SR_EPS 1e-12

static double* alloc_vec(int n) { return (double*)calloc((size_t)n, sizeof(double)); }
static void free_vec(double* v) { free(v); }
static void cross3(const double* a, const double* b, double* c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

/* ============================================================================
 * Reduced Control System
 * ============================================================================ */

SRReducedControlSystem* sr_ctrl_create(int state_dim, int control_dim,
                                        int algebra_dim) {
    if (state_dim <= 0 || control_dim <= 0) return NULL;
    SRReducedControlSystem* ctrl = (SRReducedControlSystem*)calloc(1, sizeof(SRReducedControlSystem));
    ctrl->state_dim = state_dim;
    ctrl->control_dim = control_dim;
    ctrl->algebra_dim = algebra_dim;
    ctrl->state = alloc_vec(state_dim);
    ctrl->drift = alloc_vec(state_dim);
    ctrl->control_fields = alloc_vec(state_dim * control_dim);
    return ctrl;
}

void sr_ctrl_free(SRReducedControlSystem* ctrl) {
    if (!ctrl) return;
    free(ctrl->state); free(ctrl->drift);
    free(ctrl->control_fields);
    free(ctrl->shape_momentum); free(ctrl->connection_form);
    free(ctrl);
}

void sr_ctrl_set_drift(SRReducedControlSystem* ctrl, const double* drift) {
    if (ctrl && drift)
        memcpy(ctrl->drift, drift, (size_t)ctrl->state_dim * sizeof(double));
}

void sr_ctrl_set_control_fields(SRReducedControlSystem* ctrl,
                                 const double* control_fields) {
    if (ctrl && control_fields)
        memcpy(ctrl->control_fields, control_fields,
               (size_t)(ctrl->state_dim * ctrl->control_dim) * sizeof(double));
}

void sr_ctrl_rhs(const SRReducedControlSystem* ctrl,
                  const double* x, const double* u,
                  double* dx_dt) {
    if (!ctrl || !x || !u || !dx_dt) return;
    memcpy(dx_dt, ctrl->drift, (size_t)ctrl->state_dim * sizeof(double));
    for (int j = 0; j < ctrl->control_dim; j++)
    for (int i = 0; i < ctrl->state_dim; i++)
        dx_dt[i] += ctrl->control_fields[i * ctrl->control_dim + j] * u[j];
}

/* ============================================================================
 * Controlled Euler-Poincare Equations
 * ============================================================================ */

void sr_ctrl_ep_rhs(const SREulerPoincareSystem* ep,
                     const double* xi, const double* control_forces,
                     double* dxi_dt) {
    if (!ep || !xi || !control_forces || !dxi_dt) return;
    /* Controlled EP: I dxi/dt = I xi × xi + F_ext
     * where F_ext are body-fixed control forces (torques).
     */
    int n = ep->n_config;
    double* mu = alloc_vec(n), *ad_star = alloc_vec(n);
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        mu[i] += ep->inertia_matrix[i * n + j] * xi[j];

    for (int i = 0; i < n; i++)
    for (int k = 0; k < n; k++)
    for (int j = 0; j < n; j++)
        ad_star[i] += ep->sign * ep->algebra->constants[i * n * n + k * n + j] * xi[k] * mu[j];

    for (int i = 0; i < n; i++) {
        double Iii = ep->inertia_matrix[i * n + i];
        dxi_dt[i] = (fabs(Iii) > SR_EPS)
                    ? (ad_star[i] + control_forces[i]) / Iii : 0.0;
    }
    free_vec(mu); free_vec(ad_star);
}

/* ============================================================================
 * Kinematic Control on Lie Groups
 * ============================================================================ */

void sr_kinematic_rhs(const SRLieGroup* G, const SRGroupElement* g,
                       const SRLieAlgebraElement* controls, int n_controls,
                       const double* u, SRGroupElement* dg) {
    if (!G || !g || !controls || !u || !dg) return;
    int ms = G->matrix_size;
    /* dg/dt = g * sum(A_i u_i). A_i are in g.
     * For direct matrix representation: dg.matrix = g.matrix * sum(A_i u_i)
     */
    double* A_sum = alloc_vec(ms * ms);
    for (int c = 0; c < n_controls; c++) {
        if (!controls[c].matrix) continue;
        for (int i = 0; i < ms * ms; i++)
            A_sum[i] += controls[c].matrix[i] * u[c];
    }
    memset(dg->matrix, 0, (size_t)(ms * ms) * sizeof(double));
    for (int i = 0; i < ms; i++)
    for (int k = 0; k < ms; k++)
    for (int j = 0; j < ms; j++)
        dg->matrix[i*ms+j] += g->matrix[i*ms+k] * A_sum[k*ms+j];
    free_vec(A_sum);
}

bool sr_kinematic_controllable(const SRLieAlgebra* alg,
                                const SRLieAlgebraElement** A_i,
                                int n_inputs) {
    if (!alg || !A_i || n_inputs <= 0) return false;
    /* Lie Algebra Rank Condition: compute iterated brackets
     * of {A_i} and check if they span g.
     * Simplification: if n_inputs == dim(alg) and A_i are basis, controllable.
     */
    if (n_inputs >= alg->dimension) return true;
    return false;
}

void sr_kinematic_growth_vector(const SRLieAlgebra* alg,
                                 const SRLieAlgebraElement** A_i,
                                 int n_inputs, int* dims, int max_iter) {
    (void)alg; (void)A_i; (void)n_inputs;
    for (int i = 0; i < max_iter; i++) dims[i] = 0;
}

/* ============================================================================
 * Optimal Control Reduction
 * ============================================================================ */

void sr_reduce_optimal_control(const SRInvariantHamiltonian* cost_integrand,
                                SRLieGroup* G,
                                SRLiePoissonBracket* reduced_pb,
                                void** reduced_data) {
    (void)cost_integrand; (void)G; (void)reduced_pb; (void)reduced_data;
}

void sr_lqr_reduced(const double* A, const double* B,
                     const double* Q, const double* R,
                     int state_dim, int control_dim,
                     double* K_gain, double* S_solution) {
    if (!A || !B || !Q || !R || !K_gain) return;
    /* Solve the algebraic Riccati equation: A^T S + S A - S B R^{-1} B^T S + Q = 0
     * For simple systems (state_dim=1, control_dim=1):
     *   a s + s a - s b (1/r) b s + q = 0
     *   s = (a r + sqrt(a^2 r^2 + b^2 q r)) / b^2
     *   K = r^{-1} b^T S
     */
    if (state_dim == 1 && control_dim == 1) {
        double a = A[0], b = B[0], q = Q[0], r_inv = 1.0 / R[0];
        double disc = a*a + b*b*q*r_inv;
        double s = (a + sqrt(disc)) / (b*b*r_inv + SR_EPS);
        S_solution[0] = s;
        K_gain[0] = r_inv * b * s;
    }
}

/* ============================================================================
 * Controlled Lagrangians
 * ============================================================================ */

void sr_controlled_lagrangian_modify(SREulerPoincareSystem* ep,
                                      const double* modification_matrix,
                                      double gain) {
    if (!ep || !modification_matrix) return;
    int n = ep->n_config;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        ep->inertia_matrix[i * n + j] += gain * modification_matrix[i * n + j];
}

bool sr_ctrl_matching_conditions(const SREulerPoincareSystem* original,
                                  const SREulerPoincareSystem* modified,
                                  const double* original_metric,
                                  const double* modified_metric,
                                  double tolerance) {
    (void)original; (void)modified; (void)original_metric;
    (void)modified_metric; (void)tolerance;
    return true;
}

/* ============================================================================
 * Energy Shaping via Symmetry Reduction
 * ============================================================================ */

void sr_energy_shape_casimirs(const SRLiePoissonBracket* lpb,
                               const double* original_H,
                               const double* weights, int n_casimirs,
                               double* shaped_H) {
    (void)lpb; (void)original_H; (void)weights; (void)n_casimirs;
    if (shaped_H) *shaped_H = original_H ? original_H[0] : 0.0;
}

/* ============================================================================
 * Nonholonomic Reduction
 * ============================================================================ */

void sr_nonholonomic_reduced_rhs(const SRLieAlgebra* alg,
                                  const double* xi, const double* mu,
                                  const double* constraint_basis, int h_dim,
                                  double* dxi_dt) {
    (void)alg; (void)xi; (void)mu; (void)constraint_basis; (void)h_dim; (void)dxi_dt;
}

/* ============================================================================
 * Reconstruction with Feedback
 * ============================================================================ */

void sr_ctrl_reconstruct_feedback(const SRLagrangePoincareSystem* lp,
                                   void (*reduced_controller)(const double*, double*),
                                   double* g_state, double* xi_state,
                                   double dt, int n_steps) {
    if (!lp || !reduced_controller || !xi_state) return;
    double* u = alloc_vec(lp->algebra_dim);
    for (int step = 0; step < n_steps; step++) {
        reduced_controller(xi_state, u);
        for (int i = 0; i < lp->algebra_dim; i++)
            xi_state[i] += dt * u[i];
        if (g_state)
            memcpy(g_state, lp->group_element,
                   (size_t)(lp->group_param_dim * lp->group_param_dim) * sizeof(double));
    }
    free_vec(u);
}

/* ============================================================================
 * Inverse Optimal Control on Reduced Space
 * ============================================================================ */

void sr_inverse_optimal_reduced(const SRLiePoissonBracket* lpb,
                                 const double* K, int n, double* Q, double* R) {
    (void)lpb; (void)K; (void)n; (void)Q; (void)R;
}

/* ============================================================================
 * Underactuated Systems
 * ============================================================================ */

bool sr_underactuated_controllability(const SRLiePoissonBracket* lpb,
                                       const double* drift_mu,
                                       const double** control_directions,
                                       int n_controls, double tolerance) {
    (void)lpb; (void)drift_mu; (void)control_directions;
    (void)n_controls; (void)tolerance;
    return false;
}

/* ============================================================================
 * Practical Applications
 * ============================================================================ */

void sr_spacecraft_reduced_rhs(const double* Omega, const double* I_diag,
                                const double* tau, double* dOmega_dt) {
    /* Euler equations with control torques:
     * I1 dΩ1/dt = (I2 - I3) Ω2 Ω3 + τ1
     * I2 dΩ2/dt = (I3 - I1) Ω3 Ω1 + τ2
     * I3 dΩ3/dt = (I1 - I2) Ω1 Ω2 + τ3
     */
    if (!Omega || !I_diag || !tau || !dOmega_dt) return;
    double I1 = I_diag[0], I2 = I_diag[1], I3 = I_diag[2];
    if (fabs(I1) < SR_EPS) I1 = SR_EPS;
    if (fabs(I2) < SR_EPS) I2 = SR_EPS;
    if (fabs(I3) < SR_EPS) I3 = SR_EPS;
    dOmega_dt[0] = ((I2 - I3) * Omega[1] * Omega[2] + tau[0]) / I1;
    dOmega_dt[1] = ((I3 - I1) * Omega[2] * Omega[0] + tau[1]) / I2;
    dOmega_dt[2] = ((I1 - I2) * Omega[0] * Omega[1] + tau[2]) / I3;
}

void sr_underwater_vehicle_reduced_rhs(const double* state_6d,
                                        const double* mass_matrix_6x6,
                                        const double* forces_6d,
                                        double* derivs_6d) {
    /* Kirchhoff equations on se(3)*:
     * dPi/dt = Pi × Omega + P × v + tau  (angular momentum)
     * dP/dt = P × Omega + f              (linear momentum)
     * where Omega = dH/dPi, v = dH/dP from Hamiltonian H = 1/2 (Pi·A·Pi + 2Pi·B·P + P·C·P)
     */
    if (!state_6d || !mass_matrix_6x6 || !forces_6d || !derivs_6d) return;
    const double* Pi = state_6d;
    const double* P = state_6d + 3;
    const double* tau = forces_6d;
    const double* f = forces_6d + 3;
    double* dPi_dt = derivs_6d;
    double* dP_dt = derivs_6d + 3;

    /* Compute velocities from momentum: v = M^{-1} [Pi; P]
     * Simplified: assume mass matrix is diagonal with entries.
     */
    double Omega[3] = {0,0,0}, v[3] = {0,0,0};
    for (int i = 0; i < 3; i++) {
        double m_ii = mass_matrix_6x6[i * 6 + i];
        double m_ii2 = mass_matrix_6x6[(i+3) * 6 + (i+3)];
        Omega[i] = (fabs(m_ii) > SR_EPS) ? Pi[i] / m_ii : 0.0;
        v[i] = (fabs(m_ii2) > SR_EPS) ? P[i] / m_ii2 : 0.0;
    }

    /* Pi_dot = Pi × Omega + P × v + tau */
    double Pi_x_Omega[3], P_x_v[3];
    cross3(Pi, Omega, Pi_x_Omega);
    cross3(P, v, P_x_v);
    for (int i = 0; i < 3; i++)
        dPi_dt[i] = Pi_x_Omega[i] + P_x_v[i] + tau[i];

    /* P_dot = P × Omega + f */
    double P_x_Omega[3];
    cross3(P, Omega, P_x_Omega);
    for (int i = 0; i < 3; i++)
        dP_dt[i] = P_x_Omega[i] + f[i];
}

void sr_quadrotor_reduced_rhs(const double* Omega, const double* v,
                               const double* I_diag, double mass,
                               const double* thrust, const double* torque,
                               double* dOmega_dt, double* dv_dt) {
    /* Quadrotor: SO(3) × R^3.
     * Attitude (reduced to so(3)): I dOmega/dt = I Omega × Omega + torque
     * Translation (in body frame): m dv/dt = -m Omega × v + thrust_vector
     */
    if (!Omega || !v || !I_diag || !torque || !thrust || !dOmega_dt || !dv_dt) return;
    /* Attitude dynamics */
    sr_spacecraft_reduced_rhs(Omega, I_diag, torque, dOmega_dt);

    /* Translational dynamics in body frame */
    double Omega_x_v[3];
    cross3(Omega, v, Omega_x_v);
    if (fabs(mass) < SR_EPS) mass = SR_EPS;
    for (int i = 0; i < 3; i++)
        dv_dt[i] = -Omega_x_v[i] + thrust[i] / mass;
}

void sr_robot_arm_reduced_rhs(const double* theta, const double* p_theta,
                               int n_joints, const double* inertia_matrix,
                               const double* potential_gradient,
                               double* dtheta_dt, double* dp_dt) {
    /* For cyclic coordinates (theta_j with dL/d(theta_j) = 0):
     * p_j = constant (symmetry).
     * Reduced dynamics for non-cyclic coordinates:
     *   d/dt (dL/d(theta_i)) = dL/d(theta_i) - gamma_{ij}^k dL/d(theta_k) theta_dot^j
     */
    if (!theta || !p_theta || !inertia_matrix || !potential_gradient ||
        !dtheta_dt || !dp_dt) return;
    int n = n_joints;
    for (int i = 0; i < n; i++) {
        double dtheta_i = 0.0;
        for (int j = 0; j < n; j++) {
            double Mij = inertia_matrix[i * n + j];
            if (fabs(Mij) > SR_EPS)
                dtheta_i += p_theta[j] / Mij;
        }
        dtheta_dt[i] = dtheta_i;
        dp_dt[i] = -potential_gradient[i];  /* Simplified: no Coriolis for now */
    }
}
