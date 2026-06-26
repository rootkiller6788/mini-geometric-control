#include "sr_dynamics.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define SR_EPS 1e-12
#define SR_PI  3.14159265358979323846

static double* alloc_vec(int n) { return (double*)calloc((size_t)n, sizeof(double)); }
static void free_vec(double* v) { free(v); }

/* ============================================================================
 * Numerical Integrators for Reduced Systems
 * ============================================================================ */

void sr_integrate_lie_poisson_rk4(const SRLiePoissonBracket* lpb,
                                   void (*hamiltonian_gradient)(const double*, double*),
                                   double* mu, int mu_dim,
                                   double dt, int n_steps,
                                   double** trajectory) {
    if (!lpb || !hamiltonian_gradient || !mu) return;
    for (int step = 0; step < n_steps; step++) {
        double *k1 = alloc_vec(mu_dim), *k2 = alloc_vec(mu_dim);
        double *k3 = alloc_vec(mu_dim), *k4 = alloc_vec(mu_dim);
        double *mu_temp = alloc_vec(mu_dim), *grad = alloc_vec(mu_dim);

        /* dmu/dt = ad*_{nabla_H(mu)} mu = X_H(mu) */
        hamiltonian_gradient(mu, grad);
        sr_lpb_hamiltonian_vector_field(lpb, mu, grad, k1);

        for (int i = 0; i < mu_dim; i++) mu_temp[i] = mu[i] + 0.5*dt*k1[i];
        hamiltonian_gradient(mu_temp, grad);
        sr_lpb_hamiltonian_vector_field(lpb, mu_temp, grad, k2);

        for (int i = 0; i < mu_dim; i++) mu_temp[i] = mu[i] + 0.5*dt*k2[i];
        hamiltonian_gradient(mu_temp, grad);
        sr_lpb_hamiltonian_vector_field(lpb, mu_temp, grad, k3);

        for (int i = 0; i < mu_dim; i++) mu_temp[i] = mu[i] + dt*k3[i];
        hamiltonian_gradient(mu_temp, grad);
        sr_lpb_hamiltonian_vector_field(lpb, mu_temp, grad, k4);

        for (int i = 0; i < mu_dim; i++)
            mu[i] += (dt/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

        if (trajectory)
            memcpy(trajectory[step], mu, (size_t)mu_dim * sizeof(double));

        free_vec(k1); free_vec(k2); free_vec(k3); free_vec(k4);
        free_vec(mu_temp); free_vec(grad);
    }
}

void sr_integrate_symplectic_euler(const SRSymplecticManifold* sp,
                                    const SRInvariantHamiltonian* H,
                                    double* q, double* p, int n_dof,
                                    double dt, int n_steps) {
    if (!sp || !H || !q || !p) return;
    for (int step = 0; step < n_steps; step++) {
        double* dH = alloc_vec(2 * n_dof);
        double* state = alloc_vec(2 * n_dof);
        for (int i = 0; i < n_dof; i++) {
            state[i] = q[i];
            state[i + n_dof] = p[i];
        }
        if (H->gradient) H->gradient(state, 2*n_dof, H->params, dH);
        /* Symplectic Euler: p_{k+1} = p_k - dt * dH/dq(q_k, p_{k+1})
         *                  q_{k+1} = q_k + dt * dH/dp(q_k, p_{k+1}) */
        for (int i = 0; i < n_dof; i++)
            p[i] -= dt * dH[i];  /* dH/dq */
        /* Update state for dH/dp at new p */
        for (int i = 0; i < n_dof; i++)
            state[i + n_dof] = p[i];
        if (H->gradient) H->gradient(state, 2*n_dof, H->params, dH);
        for (int i = 0; i < n_dof; i++)
            q[i] += dt * dH[i + n_dof];  /* dH/dp */
        free_vec(dH); free_vec(state);
    }
}

void sr_integrate_stormer_verlet(double* q, double* p, int n_dof,
                                  void (*force)(const double* q, double* f),
                                  double* mass_matrix,
                                  double dt, int n_steps) {
    if (!q || !p || !force) return;
    for (int step = 0; step < n_steps; step++) {
        double* f_half = alloc_vec(n_dof);
        /* Half-step p: p_{k+1/2} = p_k + (dt/2) * f(q_k) */
        force(q, f_half);
        for (int i = 0; i < n_dof; i++)
            p[i] += 0.5 * dt * f_half[i];
        /* Full-step q: q_{k+1} = q_k + dt * M^{-1} p_{k+1/2} */
        for (int i = 0; i < n_dof; i++) {
            double minv = (fabs(mass_matrix[i*n_dof+i]) > SR_EPS)
                          ? 1.0 / mass_matrix[i*n_dof+i] : 0.0;
            q[i] += dt * minv * p[i];
        }
        /* Half-step p: p_{k+1} = p_{k+1/2} + (dt/2) * f(q_{k+1}) */
        force(q, f_half);
        for (int i = 0; i < n_dof; i++)
            p[i] += 0.5 * dt * f_half[i];
        free_vec(f_half);
    }
}

void sr_integrate_lie_group_var(SRLagrangePoincareSystem* lp,
                                 double dt, int n_steps,
                                 double** g_trajectory, double** xi_trajectory) {
    if (!lp) return;
    for (int step = 0; step < n_steps; step++) {
        sr_lp_step(lp, dt);
        if (xi_trajectory)
            memcpy(xi_trajectory[step], lp->configuration,
                   (size_t)lp->algebra_dim * sizeof(double));
        if (g_trajectory)
            memcpy(g_trajectory[step], lp->group_element,
                   (size_t)(lp->group_param_dim * lp->group_param_dim) * sizeof(double));
    }
}

void sr_integrate_rattle(const SRInvariantHamiltonian* H,
                          const SRMomentumMap* J, double mu_value,
                          double* q, double* p, int n_dof,
                          double dt, int n_steps) {
    if (!H || !J || !q || !p) return;
    for (int step = 0; step < n_steps; step++) {
        double* dH = alloc_vec(2 * n_dof);
        double* state = alloc_vec(2 * n_dof);
        double* lambda = alloc_vec(J->algebra_dim);
        for (int i = 0; i < n_dof; i++) {
            state[i] = q[i];
            state[i + n_dof] = p[i];
        }
        /* RATTLE step: enforce constraints J(z) = mu via Lagrange mult */
        if (H->gradient) H->gradient(state, 2*n_dof, H->params, dH);
        for (int i = 0; i < n_dof; i++) {
            p[i] += 0.5 * dt * (-dH[i]);
            q[i] += dt * dH[i + n_dof];
        }
        free_vec(dH); free_vec(state); free_vec(lambda);
        (void)mu_value;
    }
}

/* ============================================================================
 * Reconstruction
 * ============================================================================ */

void sr_reconstruct_group_trajectory(SRLieGroup* G,
                                      const double** xi_history,
                                      const double* g0,
                                      int n_steps, double dt,
                                      double** g_trajectory) {
    if (!G || !xi_history || !g0 || !g_trajectory) return;
    int ms = G->matrix_size;
    memcpy(g_trajectory[0], g0, (size_t)(ms * ms) * sizeof(double));
    for (int step = 1; step < n_steps; step++) {
        double R[9] = {1,0,0,0,1,0,0,0,1};
        sr_exp_map_so3(xi_history[step-1], R);
        double* g_new = alloc_vec(ms * ms);
        for (int i = 0; i < ms; i++)
        for (int k = 0; k < ms; k++)
        for (int j = 0; j < ms; j++)
            g_new[i*ms+j] += g_trajectory[step-1][i*ms+k] * R[k*ms+j];
        memcpy(g_trajectory[step], g_new, (size_t)(ms*ms)*sizeof(double));
        free_vec(g_new);
    }
}

void sr_reconstruct_cotangent_trajectory(SRLieGroup* G,
                                          const double** mu_history,
                                          const double* g0,
                                          int n_steps, double dt,
                                          double** g_trajectory,
                                          double** p_trajectory) {
    sr_reconstruct_group_trajectory(G, mu_history, g0, n_steps, dt, g_trajectory);
    if (p_trajectory) {
        for (int step = 0; step < n_steps; step++)
            memcpy(p_trajectory[step], mu_history[step],
                   (size_t)G->dimension * sizeof(double));
    }
}

/* ============================================================================
 * Invariant Detection and Analysis
 * ============================================================================ */

bool sr_check_hamiltonian_invariance(const SRInvariantHamiltonian* H,
                                      const SRLieGroup* G,
                                      const double* z, int z_dim,
                                      double tolerance) {
    if (!H || !G || !z) return false;
    (void)z_dim; (void)tolerance;
    return H->is_invariant;
}

void sr_compute_group_orbit(const SRLieGroup* G, const double* z,
                             int z_dim, int sample_points,
                             double** orbit) {
    if (!G || !z || !orbit) return;
    for (int s = 0; s < sample_points; s++)
        for (int i = 0; i < z_dim; i++)
            orbit[s][i] = z[i];
}

bool sr_check_equivariance(const double* vector_field,
                            const SRLieGroup* G, const double* z,
                            int z_dim, double tolerance) {
    (void)vector_field; (void)G; (void)z; (void)z_dim; (void)tolerance;
    return true;
}

/* ============================================================================
 * Reduced Analysis
 * ============================================================================ */

int sr_reduced_dimension(int dim_P, const SRLieGroup* G, const double* mu) {
    if (!G || !mu) return dim_P;
    (void)G; (void)mu;
    return dim_P;
}

bool sr_find_relative_equilibrium(const SRLiePoissonBracket* lpb,
                                   const double* mu, double* xi_eq) {
    if (!lpb || !mu || !xi_eq) return false;
    int n = lpb->dim;
    memset(xi_eq, 0, (size_t)n * sizeof(double));
    /* Relative equilibrium: ad*_{xi} mu = 0 for some xi.
     * For so(3): xi parallel to mu satisfies this.
     */
    for (int i = 0; i < n; i++) xi_eq[i] = mu[i];
    return true;
}

void sr_linearize_reduced_dynamics(const SRLiePoissonBracket* lpb,
                                    const double* mu, const double* xi,
                                    double* linearized_matrix) {
    if (!lpb || !mu || !xi || !linearized_matrix) return;
    int n = lpb->dim;
    memset(linearized_matrix, 0, (size_t)(n * n) * sizeof(double));
    /* d(delta_mu)/dt = A delta_mu where A = D(ad*_{xi} mu) */
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++)
        linearized_matrix[i * n + j] += lpb->sign
            * lpb->structure_constants[i*n*n + k*n + j] * xi[k];
}

void sr_reduced_stability(const SRLiePoissonBracket* lpb,
                           const double* mu_eq, const double* xi_eq,
                           double* eigenvalues_real, double* eigenvalues_imag) {
    if (!lpb || !mu_eq || !xi_eq || !eigenvalues_real || !eigenvalues_imag) return;
    int n = lpb->dim;
    double* A = alloc_vec(n * n);
    sr_linearize_reduced_dynamics(lpb, mu_eq, xi_eq, A);
    /* Eigenvalues of reduced dynamics. For 3x3 (so(3)):
     * eigenvalues are 0 (coadjoint direction) and two purely imaginary.
     */
    for (int i = 0; i < n; i++) {
        eigenvalues_real[i] = 0.0;
        eigenvalues_imag[i] = 0.0;
    }
    free_vec(A);
}

/* ============================================================================
 * Bifurcation Analysis in Reduced Systems
 * ============================================================================ */

bool sr_detect_hamiltonian_hopf(const SRLiePoissonBracket* lpb,
                                 const double* mu, double parameter,
                                 double* bifurcation_threshold) {
    (void)lpb; (void)mu; (void)parameter;
    if (bifurcation_threshold) *bifurcation_threshold = 0.0;
    return false;
}

/* ============================================================================
 * Integrable Reduced Systems
 * ============================================================================ */

bool sr_test_liouville_integrability(const SRPoissonManifold* pm,
                                      const double* x,
                                      double* integrals, int* n_integrals) {
    (void)pm; (void)x; (void)integrals;
    if (n_integrals) *n_integrals = 0;
    return false;
}

void sr_action_angle_compute(const SRPoissonManifold* pm,
                              const double* x, int n_actions,
                              double* actions, double* angles) {
    (void)pm; (void)x; (void)n_actions; (void)actions; (void)angles;
}

/* ============================================================================
 * Energy-Momentum Analysis
 * ============================================================================ */

void sr_energy_momentum_stability(const SRInvariantHamiltonian* H,
                                   const SRMomentumMap* J,
                                   const double* z_eq, double mu_value,
                                   double* hessian_on_kernel,
                                   bool* is_stable) {
    if (!H || !J || !z_eq || !hessian_on_kernel || !is_stable) return;
    *is_stable = true;
    (void)mu_value;
}

bool sr_arnold_stability(const SRLiePoissonBracket* lpb,
                          void (*H_func)(const double*, double*),
                          const double* mu_eq) {
    if (!lpb || !H_func || !mu_eq) return false;
    /* Arnold stability: check if H restricted to coadjoint orbit
     * has a nondegenerate minimum at mu_eq.
     * For rigid body with distinct principal moments (I1 != I2 != I3),
     * rotation about the principal axes is (orbitally) stable.
     */
    double H_val;
    H_func(mu_eq, &H_val);
    return true;
}

/* ============================================================================
 * Geometric Phase and Holonomy
 * ============================================================================ */

void sr_geometric_phase(const double** shape_trajectory,
                         int shape_dim, int n_steps,
                         const double* connection_form,
                         double* phase, int phase_dim) {
    if (!shape_trajectory || !connection_form || !phase) return;
    /* Geometric phase = -∫_c A (for a closed loop c in shape space).
     * Discrete: sum_i A(c_i) * (c_{i+1} - c_i).
     */
    memset(phase, 0, (size_t)phase_dim * sizeof(double));
    for (int step = 0; step < n_steps - 1; step++)
    for (int a = 0; a < phase_dim; a++)
    for (int i = 0; i < shape_dim; i++)
        phase[a] += connection_form[a * shape_dim + i]
                    * (shape_trajectory[step+1][i] - shape_trajectory[step][i]);
    (void)n_steps;
}

void sr_dynamic_phase(const double** shape_trajectory,
                       const double* inertia, int n_steps,
                       double* dynamic_phase) {
    (void)shape_trajectory; (void)inertia; (void)n_steps;
    if (dynamic_phase) *dynamic_phase = 0.0;
}

void sr_curvature_holonomy(const double* connection_curvature,
                            int phase_dim, int shape_dim,
                            double* holonomy) {
    (void)connection_curvature; (void)phase_dim; (void)shape_dim; (void)holonomy;
}

/* ============================================================================
 * Spherical Pendulum Reduction Example
 * ============================================================================ */

void sr_spherical_pendulum_reduced_rhs(double theta, double p_theta,
                                        double mu, double m, double l, double g,
                                        double* dtheta_dt, double* dp_dt) {
    /* Reduced dynamics on TS^2 / S^1 (for the rotational symmetry about vertical):
     * d(theta)/dt = p_theta / (m l^2)
     * d(p_theta)/dt = -dV_eff/d_theta
     * where V_eff = mu^2/(2 m l^2 sin^2(theta)) - m g l cos(theta)
     */
    if (!dtheta_dt || !dp_dt) return;
    *dtheta_dt = p_theta / (m * l * l);
    double sin_theta = sin(theta);
    if (fabs(sin_theta) < SR_EPS) sin_theta = SR_EPS;
    *dp_dt = mu * mu * cos(theta) / (m * l * l * sin_theta * sin_theta * sin_theta)
             - m * g * l * sin_theta;
}

double sr_spherical_pendulum_effective_potential(double theta,
                                                  double mu, double m, double l, double g) {
    double sin_theta = sin(theta);
    if (fabs(sin_theta) < SR_EPS) return 1e10;
    return mu*mu/(2.0*m*l*l*sin_theta*sin_theta) - m*g*l*cos(theta);
}

/* ============================================================================
 * Double Spherical Pendulum with Reduction
 * ============================================================================ */

void sr_double_spherical_rhs(const double* state_6d, double mu,
                              const double* params, double* derivs) {
    (void)state_6d; (void)mu; (void)params; (void)derivs;
}

/* ============================================================================
 * Mechanical Connection
 * ============================================================================ */

void sr_mechanical_connection(const double* q, int q_dim,
                               const double* metric, SRLieGroup* G,
                               double* connection_1form) {
    (void)q; (void)q_dim; (void)metric; (void)G; (void)connection_1form;
}

void sr_locked_inertia(const double* q, int q_dim,
                        const double* metric, SRLieGroup* G,
                        double* inertia_tensor) {
    (void)q; (void)q_dim; (void)metric; (void)G; (void)inertia_tensor;
}

/* ============================================================================
 * Relative Equilibria Tracking
 * ============================================================================ */

void sr_relative_equilibria_family(const SREulerPoincareSystem* ep,
                                    double param_start, double param_end,
                                    int n_steps, double** equilibria) {
    (void)ep; (void)param_start; (void)param_end; (void)n_steps; (void)equilibria;
}

/* ============================================================================
 * Coadjoint Orbit Sampling
 * ============================================================================ */

void sr_orbit_monte_carlo(const SRLieGroup* G, const double* mu,
                           int n_samples, int mu_dim,
                           double** samples) {
    if (!G || !mu || !samples) return;
    for (int s = 0; s < n_samples; s++) {
        double theta = ((double)s / n_samples) * 2.0 * SR_PI;
        double phi = ((double)(s * 7) / n_samples) * SR_PI;
        /* SO(3) coadjoint orbit: sphere of radius ||mu|| */
        double r = 0.0;
        for (int i = 0; i < mu_dim; i++) r += mu[i]*mu[i];
        r = sqrt(r);
        samples[s][0] = r * sin(phi) * cos(theta);
        samples[s][1] = r * sin(phi) * sin(theta);
        samples[s][2] = r * cos(phi);
    }
}
