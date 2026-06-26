/*=============================================================================
 * poisson_geometry.c -- Poisson geometry for control systems
 * Reference: Marsden & Ratiu, Introduction to Mechanics and Symmetry
 *            Weinstein, The local structure of Poisson manifolds (1983)
 *===========================================================================*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "hamiltonian_control.h"
#include "poisson_geometry.h"

/* L1: Poisson bracket on general Poisson manifold */
double poisson_bracket_general(const poisson_manifold_t *pm,
                                double (*F)(const double*, int, void*),
                                void *pF,
                                double (*G)(const double*, int, void*),
                                void *pG,
                                const double *z)
{
    int d = pm->dim;
    double **B = matrix_alloc(d, d);
    pm->B_tensor(z, d, B, pm->params);

    double eps = 1e-6;
    double *grad_F = (double*)malloc(d * sizeof(double));
    double *grad_G = (double*)malloc(d * sizeof(double));
    double *z_pert = (double*)malloc(d * sizeof(double));

    double F0 = F(z, d, pF);
    double G0 = G(z, d, pG);

    for (int i = 0; i < d; i++) {
        memcpy(z_pert, z, d * sizeof(double));
        z_pert[i] += eps;
        grad_F[i] = (F(z_pert, d, pF) - F0) / eps;
        grad_G[i] = (G(z_pert, d, pG) - G0) / eps;
    }

    double bracket = 0.0;
    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            bracket += grad_F[i] * B[i][j] * grad_G[j];

    matrix_free(B, d);
    free(grad_F); free(grad_G); free(z_pert);
    return bracket;
}

/* L1: Lie bracket [xi, eta] on Lie algebra */
void lie_bracket(const lie_algebra_t *g,
                 const double *xi, const double *eta,
                 double *bracket_out)
{
    int dim = g->dim;
    for (int k = 0; k < dim; k++) {
        bracket_out[k] = 0.0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                bracket_out[k] += g->c_ijk[i*dim*dim + j*dim + k] * xi[i] * eta[j];
    }
}

/* L1: Lie-Poisson bracket on g* */
double lie_poisson_bracket(const lie_algebra_t *g,
                            double (*F)(const double*, int, void*),
                            void *pF,
                            double (*G)(const double*, int, void*),
                            void *pG,
                            const double *mu)
{
    int dim = g->dim;
    double eps = 1e-6;
    double *grad_F = (double*)calloc(dim, sizeof(double));
    double *grad_G = (double*)calloc(dim, sizeof(double));
    double *mu_pert = (double*)malloc(dim * sizeof(double));

    double F0 = F(mu, dim, pF);
    double G0 = G(mu, dim, pG);

    for (int i = 0; i < dim; i++) {
        memcpy(mu_pert, mu, dim * sizeof(double));
        mu_pert[i] += eps;
        grad_F[i] = (F(mu_pert, dim, pF) - F0) / eps;
        grad_G[i] = (G(mu_pert, dim, pG) - G0) / eps;
    }

    /* {F,G}(mu) = <mu, [dF, dG]> */
    double *bracket_vec = (double*)malloc(dim * sizeof(double));
    lie_bracket(g, grad_F, grad_G, bracket_vec);

    double result = 0.0;
    for (int i = 0; i < dim; i++)
        result += mu[i] * bracket_vec[i];

    free(grad_F); free(grad_G); free(mu_pert); free(bracket_vec);
    return result;
}

/* L1: Check if C is a Casimir function */
int is_casimir_function(const poisson_manifold_t *pm,
                         double (*C)(const double*, int, void*),
                         void *pC,
                         int n_test_functions,
                         double (**test_fns)(const double*, int, void*),
                         void **test_params,
                         const double *z,
                         double *max_bracket_norm)
{
    *max_bracket_norm = 0.0;
    for (int t = 0; t < n_test_functions; t++) {
        double bracket = poisson_bracket_general(pm, C, pC, test_fns[t], test_params[t], z);
        if (fabs(bracket) > *max_bracket_norm) *max_bracket_norm = fabs(bracket);
    }
    return (*max_bracket_norm < 1e-8) ? 1 : 0;
}

/* L2: Dimension of symplectic leaf through z = rank(B(z)) */
int symplectic_leaf_dimension(const poisson_manifold_t *pm, const double *z)
{
    int d = pm->dim;
    double **B = matrix_alloc(d, d);
    pm->B_tensor(z, d, B, pm->params);

    /* Rank via QR decomposition with column pivoting */
    double **B_copy = matrix_alloc(d, d);
    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            B_copy[i][j] = B[i][j];

    /* Simple rank estimate: count non-zero singular values via Gaussian elimination */
    int rank = 0;
    double **U = matrix_alloc(d, d);
    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            U[i][j] = B_copy[i][j];

    for (int col = 0; col < d; col++) {
        /* Find pivot */
        int pivot = -1;
        for (int row = col; row < d; row++) {
            if (fabs(U[row][col]) > 1e-10) { pivot = row; break; }
        }
        if (pivot >= 0) {
            rank++;
            /* Swap rows */
            if (pivot != col) {
                for (int j = 0; j < d; j++) {
                    double tmp = U[col][j]; U[col][j] = U[pivot][j]; U[pivot][j] = tmp;
                }
            }
            /* Eliminate below */
            for (int row = col + 1; row < d; row++) {
                double factor = U[row][col] / U[col][col];
                for (int j = col; j < d; j++)
                    U[row][j] -= factor * U[col][j];
            }
        }
    }

    matrix_free(B, d); matrix_free(B_copy, d); matrix_free(U, d);
    return rank;
}

/* L2: Basis for symplectic leaf through z */
int symplectic_leaf_basis(const poisson_manifold_t *pm, const double *z,
                           double **leaf_basis, int *leaf_dim)
{
    int d = pm->dim;
    double **B = matrix_alloc(d, d);
    pm->B_tensor(z, d, B, pm->params);

    /* Range of B is tangent to symplectic leaf */
    /* Copy columns of B as basis, then orthogonalize */
    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            leaf_basis[i][j] = B[j][i]; /* columns */

    /* QR to extract independent columns */
    int rank = 0;
    for (int col = 0; col < d; col++) {
        /* Gram-Schmidt against previous */
        for (int k = 0; k < rank; k++) {
            double dot = 0.0;
            for (int i = 0; i < d; i++) dot += leaf_basis[i][col] * leaf_basis[i][k];
            for (int i = 0; i < d; i++) leaf_basis[i][col] -= dot * leaf_basis[i][k];
        }
        double norm = 0.0;
        for (int i = 0; i < d; i++) norm += leaf_basis[i][col] * leaf_basis[i][col];
        if (sqrt(norm) > 1e-10) {
            norm = sqrt(norm);
            for (int i = 0; i < d; i++) leaf_basis[i][col] /= norm;
            rank++;
        }
    }
    *leaf_dim = rank;
    matrix_free(B, d);
    return 0;
}

/* L3: Create so(3) Lie algebra */
lie_algebra_t *so3_lie_algebra_create(void)
{
    lie_algebra_t *g = (lie_algebra_t*)malloc(sizeof(lie_algebra_t));
    g->dim = 3;
    g->c_ijk = (double*)calloc(27, sizeof(double)); /* 3*3*3 */

    /* Levi-Civita: c_{ij}^k = epsilon_{ijk} */
    /* c_{12}^3 = 1, c_{23}^1 = 1, c_{31}^2 = 1 */
    /* c_{21}^3 = -1, c_{32}^1 = -1, c_{13}^2 = -1 */
    int idx;
    /* [e1, e2] = e3: c_{12}^3 = 1 */
    idx = 0*9 + 1*3 + 2; g->c_ijk[idx] =  1.0; /* (0,1,2) = e1,e2 -> e3 */
    idx = 1*9 + 0*3 + 2; g->c_ijk[idx] = -1.0; /* (1,0,2) = e2,e1 -> -e3 */
    /* [e2, e3] = e1: c_{23}^1 = 1 */
    idx = 1*9 + 2*3 + 0; g->c_ijk[idx] =  1.0;
    idx = 2*9 + 1*3 + 0; g->c_ijk[idx] = -1.0;
    /* [e3, e1] = e2: c_{31}^2 = 1 */
    idx = 2*9 + 0*3 + 1; g->c_ijk[idx] =  1.0;
    idx = 0*9 + 2*3 + 1; g->c_ijk[idx] = -1.0;

    return g;
}

/* L3: Rigid body Hamiltonian on so(3)* */
double rigid_body_hamiltonian(const double *Pi, int dim, void *params)
{
    rigid_body_params_t *rb = (rigid_body_params_t*)params;
    return Pi[0]*Pi[0]/(2.0*rb->I1) + Pi[1]*Pi[1]/(2.0*rb->I2) + Pi[2]*Pi[2]/(2.0*rb->I3);
}

void rigid_body_gradient(const double *Pi, int dim, double *grad, void *params)
{
    rigid_body_params_t *rb = (rigid_body_params_t*)params;
    grad[0] = Pi[0] / rb->I1;
    grad[1] = Pi[1] / rb->I2;
    grad[2] = Pi[2] / rb->I3;
}

/* L3: Euler equations: dPi/dt = Pi x omega where omega_i = Pi_i / I_i */
void rigid_body_euler_equations(const rigid_body_params_t *rb,
                                 const double *Pi,
                                 double *dPi_dt)
{
    double omega[3];
    omega[0] = Pi[0] / rb->I1;
    omega[1] = Pi[1] / rb->I2;
    omega[2] = Pi[2] / rb->I3;

    /* Pi x omega */
    dPi_dt[0] = (rb->I2 - rb->I3) * omega[1] * omega[2];
    dPi_dt[1] = (rb->I3 - rb->I1) * omega[2] * omega[0];
    dPi_dt[2] = (rb->I1 - rb->I2) * omega[0] * omega[1];
}

/* L4: Lie-Poisson simulation */
int lie_poisson_simulate(const lie_algebra_t *g,
                          double (*h)(const double*, int, void*),
                          void (*grad_h)(const double*, int, double*, void*),
                          void *params,
                          double *mu, double dt, int steps,
                          double **trajectory_out)
{
    int dim = g->dim;
    for (int k = 0; k < steps; k++) {
        /* Store */
        for (int i = 0; i < dim; i++)
            trajectory_out[k][i] = mu[i];

        /* grad_h at current mu */
        double *grad = (double*)malloc(dim * sizeof(double));
        grad_h(mu, dim, grad, params);

        /* dmu/dt = ad*_{grad(h)} mu */
        /* dmu_k/dt = sum_{ij} c_{ij}^k mu_i grad_j */
        for (int kk = 0; kk < dim; kk++) {
            double dmu = 0.0;
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    dmu += g->c_ijk[i*dim*dim + j*dim + kk] * mu[i] * grad[j];
            mu[kk] += dt * dmu;
        }
        free(grad);
    }
    return 0;
}

/* L5: Energy-Casimir Lyapunov function */
int energy_casimir_lyapunov(const lie_algebra_t *g,
                             double (*h)(const double*, int, void*),
                             void (*grad_h)(const double*, int, double*, void*),
                             double (**casimirs)(const double*, int, void*),
                             void (**grad_casimirs)(const double*, int, double*, void*),
                             int n_casimirs,
                             const double *mu_eq,
                             energy_casimir_method_t *ecm,
                             double *V_min_eigenvalue)
{
    int dim = g->dim;

    /* Construct Hessian of V = H + sum phi_k C_k at equilibrium */
    double eps = 1e-5;
    double **Hess = matrix_alloc(dim, dim);
    double *grad_V = (double*)malloc(dim * sizeof(double));
    double *mu_p = (double*)malloc(dim * sizeof(double));

    /* grad_V at equilibrium (should be zero for properly chosen phi) */
    double *grad_H = (double*)malloc(dim * sizeof(double));
    grad_h(mu_eq, dim, grad_H, NULL);
    for (int i = 0; i < dim; i++) {
        grad_V[i] = grad_H[i];
        for (int k = 0; k < n_casimirs; k++) {
            double *grad_C = (double*)malloc(dim * sizeof(double));
            grad_casimirs[k](mu_eq, dim, grad_C, NULL);
            grad_V[i] += ecm->phi_coeffs[k] * grad_C[i];
            free(grad_C);
        }
    }

    /* Hessian via finite differences */
    for (int j = 0; j < dim; j++) {
        memcpy(mu_p, mu_eq, dim * sizeof(double));
        mu_p[j] += eps;
        double *grad_p = (double*)malloc(dim * sizeof(double));
        grad_h(mu_p, dim, grad_p, NULL);
        for (int k = 0; k < n_casimirs; k++) {
            double *grad_Cp = (double*)malloc(dim * sizeof(double));
            grad_casimirs[k](mu_p, dim, grad_Cp, NULL);
            for (int i = 0; i < dim; i++)
                grad_p[i] += ecm->phi_coeffs[k] * grad_Cp[i];
            free(grad_Cp);
        }
        for (int i = 0; i < dim; i++)
            Hess[i][j] = (grad_p[i] - grad_V[i]) / eps;
        free(grad_p);
    }

    /* Minimum eigenvalue via Gershgorin circle theorem */
    double min_eig = 1e30;
    for (int i = 0; i < dim; i++) {
        double center = Hess[i][i];
        double radius = 0.0;
        for (int j = 0; j < dim; j++)
            if (j != i) radius += fabs(Hess[i][j]);
        double lower = center - radius;
        if (lower < min_eig) min_eig = lower;
    }

    *V_min_eigenvalue = min_eig;
    matrix_free(Hess, dim);
    free(grad_V); free(mu_p); free(grad_H);
    return (min_eig > 0) ? 1 : 0;
}

/* L6: Rigid body stabilization via momentum wheels */
int rigid_body_stabilize(const rigid_body_control_t *ctrl,
                          const double *Pi, const double *Pi_desired,
                          double *u_out)
{
    /* PD control: u = -K_p (Pi - Pi_desired) - K_d omega */
    double omega[3];
    omega[0] = Pi[0] / ctrl->rb.I1;
    omega[1] = Pi[1] / ctrl->rb.I2;
    omega[2] = Pi[2] / ctrl->rb.I3;

    for (int i = 0; i < 3; i++) {
        double err = Pi[i] - Pi_desired[i];
        u_out[i] = -ctrl->K_p * err - ctrl->K_d * omega[i];
    }
    return 0;
}

/* L7: Heavy top equations */
void heavy_top_equations(const heavy_top_params_t *ht,
                          const double *Pi, const double *Gamma,
                          double *dPi_dt, double *dGamma_dt)
{
    /* Euler-Poisson equations:
     * dPi/dt = Pi x omega + mgl Gamma x e3
     * dGamma/dt = Gamma x omega */
    double omega[3];
    omega[0] = Pi[0] / ht->rb.I1;
    omega[1] = Pi[1] / ht->rb.I2;
    omega[2] = Pi[2] / ht->rb.I3;

    /* Pi x omega */
    dPi_dt[0] = Pi[1]*omega[2] - Pi[2]*omega[1];
    dPi_dt[1] = Pi[2]*omega[0] - Pi[0]*omega[2];
    dPi_dt[2] = Pi[0]*omega[1] - Pi[1]*omega[0];

    /* mgl Gamma x e3 */
    double e3[3] = {0, 0, 1};
    dPi_dt[0] += ht->mass * ht->g * ht->l * (Gamma[1]*e3[2] - Gamma[2]*e3[1]);
    dPi_dt[1] += ht->mass * ht->g * ht->l * (Gamma[2]*e3[0] - Gamma[0]*e3[2]);
    dPi_dt[2] += ht->mass * ht->g * ht->l * (Gamma[0]*e3[1] - Gamma[1]*e3[0]);

    /* Gamma x omega */
    dGamma_dt[0] = Gamma[1]*omega[2] - Gamma[2]*omega[1];
    dGamma_dt[1] = Gamma[2]*omega[0] - Gamma[0]*omega[2];
    dGamma_dt[2] = Gamma[0]*omega[1] - Gamma[1]*omega[0];
}

/* L8: Poisson reduction optimal control */
int poisson_reduction_optimal_control(
    const reduced_control_system_t *rcs,
    const double *mu0, const double *mu_target,
    double tf, int n_steps,
    double **mu_traj, double **u_traj)
{
    int dim = rcs->g->dim, ctrl_dim = rcs->ctrl_dim;
    double dt = tf / n_steps;
    double *mu = (double*)malloc(dim * sizeof(double));
    double *u  = (double*)malloc(ctrl_dim * sizeof(double));
    memcpy(mu, mu0, dim * sizeof(double));

    for (int k = 0; k < n_steps; k++) {
        for (int i = 0; i < dim; i++) mu_traj[k][i] = mu[i];

        /* Simple gradient-based control to drive mu -> mu_target */
        double *grad_mu = (double*)malloc(dim * sizeof(double));
        rcs->grad_mu_h(mu, u, dim, ctrl_dim, grad_mu, rcs->params);

        /* Descent direction */
        for (int i = 0; i < dim; i++)
            mu[i] -= dt * grad_mu[i];

        for (int i = 0; i < ctrl_dim; i++) u_traj[k][i] = u[i];
        free(grad_mu);
    }

    free(mu); free(u);
    return 0;
}

void lie_algebra_free(lie_algebra_t *g)
{
    if (g) { free(g->c_ijk); free(g); }
}