/*=============================================================================
 * hamiltonian_control.c -- Core Hamiltonian mechanics implementation
 * Reference: Arnold, Mathematical Methods of Classical Mechanics
 *            Goldstein, Classical Mechanics
 *            Marsden & Ratiu, Introduction to Mechanics and Symmetry
 *===========================================================================*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "hamiltonian_control.h"

/*---------------------------------------------------------------------------*
 * L1: Hamiltonian vector field X_H = J grad(H)
 * dq_i/dt =  dH/dp_i,    dp_i/dt = -dH/dq_i
 * Reference: Arnold section 16, Goldstein section 8.1
 *---------------------------------------------------------------------------*/
void hamiltonian_vector_field(const hamiltonian_t *H,
                              const double *q, const double *p,
                              double *dq_dt, double *dp_dt)
{
    int n = H->n_dof;
    double *grad_q = (double*)malloc(n * sizeof(double));
    double *grad_p = (double*)malloc(n * sizeof(double));
    H->gradient(q, p, n, grad_q, grad_p, H->params);
    for (int i = 0; i < n; i++) {
        dq_dt[i] =  grad_p[i];
        dp_dt[i] = -grad_q[i];
    }
    free(grad_q);
    free(grad_p);
}

/*---------------------------------------------------------------------------*
 * L2: RK4 flow (non-symplectic baseline for comparison)
 * Reference: Butcher, Numerical Methods for ODEs
 *---------------------------------------------------------------------------*/
void hamiltonian_flow_rk4(const hamiltonian_t *H,
                          double *q, double *p,
                          double dt, int steps)
{
    int n = H->n_dof;
    double *k1q = (double*)malloc(n * sizeof(double));
    double *k1p = (double*)malloc(n * sizeof(double));
    double *k2q = (double*)malloc(n * sizeof(double));
    double *k2p = (double*)malloc(n * sizeof(double));
    double *k3q = (double*)malloc(n * sizeof(double));
    double *k3p = (double*)malloc(n * sizeof(double));
    double *k4q = (double*)malloc(n * sizeof(double));
    double *k4p = (double*)malloc(n * sizeof(double));
    double *qt  = (double*)malloc(n * sizeof(double));
    double *pt  = (double*)malloc(n * sizeof(double));

    for (int s = 0; s < steps; s++) {
        hamiltonian_vector_field(H, q, p, k1q, k1p);
        for (int i = 0; i < n; i++) {
            qt[i] = q[i] + 0.5 * dt * k1q[i];
            pt[i] = p[i] + 0.5 * dt * k1p[i];
        }
        hamiltonian_vector_field(H, qt, pt, k2q, k2p);
        for (int i = 0; i < n; i++) {
            qt[i] = q[i] + 0.5 * dt * k2q[i];
            pt[i] = p[i] + 0.5 * dt * k2p[i];
        }
        hamiltonian_vector_field(H, qt, pt, k3q, k3p);
        for (int i = 0; i < n; i++) {
            qt[i] = q[i] + dt * k3q[i];
            pt[i] = p[i] + dt * k3p[i];
        }
        hamiltonian_vector_field(H, qt, pt, k4q, k4p);
        for (int i = 0; i < n; i++) {
            q[i] += (dt / 6.0) * (k1q[i] + 2.0*k2q[i] + 2.0*k3q[i] + k4q[i]);
            p[i] += (dt / 6.0) * (k1p[i] + 2.0*k2p[i] + 2.0*k3p[i] + k4p[i]);
        }
    }

    free(k1q); free(k1p); free(k2q); free(k2p);
    free(k3q); free(k3p); free(k4q); free(k4p);
    free(qt); free(pt);
}

/*---------------------------------------------------------------------------*
 * L1: Energy error between two phase points
 * |H(q1,p1) - H(q0,p0)| / max(|H(q0,p0)|, 1)
 * Reference: Hairer et al., Geometric Numerical Integration
 *---------------------------------------------------------------------------*/
double hamiltonian_energy_error(const hamiltonian_t *H,
                                const double *q0, const double *p0,
                                const double *q1, const double *p1)
{
    double E0 = H->energy(q0, p0, H->n_dof, H->params);
    double E1 = H->energy(q1, p1, H->n_dof, H->params);
    double denom = fabs(E0) > 1.0 ? fabs(E0) : 1.0;
    return fabs(E1 - E0) / denom;
}

/*---------------------------------------------------------------------------*
 * L1: Poisson bracket {F, G}
 * {F, G} = sum_{i=1}^{n} (dF/dq_i * dG/dp_i - dF/dp_i * dG/dq_i)
 * Encodes the symplectic structure.
 * Reference: Arnold section 40, Marsden-Ratiu section 10.1
 *---------------------------------------------------------------------------*/
double poisson_bracket(const hamiltonian_t *H,
                       double (*F)(const double*, const double*, int, void*),
                       void *params_F,
                       double (*G)(const double*, const double*, int, void*),
                       void *params_G,
                       const double *q, const double *p, int n)
{
    (void)H;
    double eps = 1e-6;
    double *grad_F_q = (double*)malloc(n * sizeof(double));
    double *grad_F_p = (double*)malloc(n * sizeof(double));
    double *grad_G_q = (double*)malloc(n * sizeof(double));
    double *grad_G_p = (double*)malloc(n * sizeof(double));
    double *q_pert   = (double*)malloc(n * sizeof(double));
    double *p_pert   = (double*)malloc(n * sizeof(double));

    double F0 = F(q, p, n, params_F);
    for (int i = 0; i < n; i++) {
        memcpy(q_pert, q, n * sizeof(double));
        memcpy(p_pert, p, n * sizeof(double));
        q_pert[i] += eps;
        grad_F_q[i] = (F(q_pert, p, n, params_F) - F0) / eps;
        q_pert[i] = q[i];
        p_pert[i] += eps;
        grad_F_p[i] = (F(q, p_pert, n, params_F) - F0) / eps;
    }

    double G0 = G(q, p, n, params_G);
    for (int i = 0; i < n; i++) {
        memcpy(q_pert, q, n * sizeof(double));
        memcpy(p_pert, p, n * sizeof(double));
        q_pert[i] += eps;
        grad_G_q[i] = (G(q_pert, p, n, params_G) - G0) / eps;
        q_pert[i] = q[i];
        p_pert[i] += eps;
        grad_G_p[i] = (G(q, p_pert, n, params_G) - G0) / eps;
    }

    double bracket = 0.0;
    for (int i = 0; i < n; i++)
        bracket += grad_F_q[i] * grad_G_p[i] - grad_F_p[i] * grad_G_q[i];

    free(grad_F_q); free(grad_F_p); free(grad_G_q); free(grad_G_p);
    free(q_pert); free(p_pert);
    return bracket;
}

/*---------------------------------------------------------------------------*
 * L4: Jacobi identity numerical verification
 * {F,{G,K}} + {G,{K,F}} + {K,{F,G}} = 0
 * Uses numerical finite differences for nested brackets.
 * Reference: Arnold section 40.B
 *---------------------------------------------------------------------------*/
double poisson_jacobi_check(double (*F)(const double*, const double*, int, void*),
                            void *pf,
                            double (*G)(const double*, const double*, int, void*),
                            void *pg,
                            double (*K)(const double*, const double*, int, void*),
                            void *pk,
                            const double *q, const double *p, int n, double h)
{
    hamiltonian_t dummy_H;
    memset(&dummy_H, 0, sizeof(dummy_H));
    dummy_H.n_dof = n;

    double FG = poisson_bracket(&dummy_H, F, pf, G, pg, q, p, n);
    double GK = poisson_bracket(&dummy_H, G, pg, K, pk, q, p, n);
    double KF = poisson_bracket(&dummy_H, K, pk, F, pf, q, p, n);

    double *dF_dq=(double*)malloc(n*sizeof(double));
    double *dF_dp=(double*)malloc(n*sizeof(double));
    double *dG_dq=(double*)malloc(n*sizeof(double));
    double *dG_dp=(double*)malloc(n*sizeof(double));
    double *dK_dq=(double*)malloc(n*sizeof(double));
    double *dK_dp=(double*)malloc(n*sizeof(double));
    double *dGK_dq=(double*)malloc(n*sizeof(double));
    double *dGK_dp=(double*)malloc(n*sizeof(double));
    double *dKF_dq=(double*)malloc(n*sizeof(double));
    double *dKF_dp=(double*)malloc(n*sizeof(double));
    double *dFG_dq=(double*)malloc(n*sizeof(double));
    double *dFG_dp=(double*)malloc(n*sizeof(double));
    double *qt=(double*)malloc(n*sizeof(double));
    double *pt=(double*)malloc(n*sizeof(double));

    double F0=F(q,p,n,pf), G0=G(q,p,n,pg), K0=K(q,p,n,pk);

    for(int i=0;i<n;i++){
        memcpy(qt,q,n*sizeof(double)); memcpy(pt,p,n*sizeof(double));
        qt[i]+=h;
        dF_dq[i]=(F(qt,p,n,pf)-F0)/h;
        dG_dq[i]=(G(qt,p,n,pg)-G0)/h;
        dK_dq[i]=(K(qt,p,n,pk)-K0)/h;
        dGK_dq[i]=(poisson_bracket(&dummy_H,G,pg,K,pk,qt,p,n)-GK)/h;
        dKF_dq[i]=(poisson_bracket(&dummy_H,K,pk,F,pf,qt,p,n)-KF)/h;
        dFG_dq[i]=(poisson_bracket(&dummy_H,F,pf,G,pg,qt,p,n)-FG)/h;
        qt[i]=q[i]; pt[i]+=h;
        dF_dp[i]=(F(q,pt,n,pf)-F0)/h;
        dG_dp[i]=(G(q,pt,n,pg)-G0)/h;
        dK_dp[i]=(K(q,pt,n,pk)-K0)/h;
        dGK_dp[i]=(poisson_bracket(&dummy_H,G,pg,K,pk,q,pt,n)-GK)/h;
        dKF_dp[i]=(poisson_bracket(&dummy_H,K,pk,F,pf,q,pt,n)-KF)/h;
        dFG_dp[i]=(poisson_bracket(&dummy_H,F,pf,G,pg,q,pt,n)-FG)/h;
    }

    double F_GK=0.0, G_KF=0.0, K_FG=0.0;
    for(int i=0;i<n;i++){
        F_GK+=dF_dq[i]*dGK_dp[i]-dF_dp[i]*dGK_dq[i];
        G_KF+=dG_dq[i]*dKF_dp[i]-dG_dp[i]*dKF_dq[i];
        K_FG+=dK_dq[i]*dFG_dp[i]-dK_dp[i]*dFG_dq[i];
    }

    free(dF_dq);free(dF_dp);free(dG_dq);free(dG_dp);free(dK_dq);free(dK_dp);
    free(dGK_dq);free(dGK_dp);free(dKF_dq);free(dKF_dp);free(dFG_dq);free(dFG_dp);
    free(qt);free(pt);
    return fabs(F_GK+G_KF+K_FG);
}

/*---------------------------------------------------------------------------*
 * L3: Check if a transformation is canonical.
 * A transformation is canonical iff its Jacobian M satisfies M^T J M = J.
 * Returns Frobenius norm of deviation.
 * Reference: Goldstein section 9.1, Arnold section 44
 *---------------------------------------------------------------------------*/
double check_canonical_transform(const canonical_transform_t *T,
                                 const double *q, const double *p, double eps)
{
    int n = T->n_dof;
    int dim = 2 * n;
    double **M = matrix_alloc(dim, dim);
    double **J = matrix_alloc(dim, dim);

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            J[i][j] = 0.0;
    for (int i = 0; i < n; i++) {
        J[i][n + i] =  1.0;
        J[n + i][i] = -1.0;
    }

    double *qp = (double*)malloc(n*sizeof(double));
    double *pp = (double*)malloc(n*sizeof(double));
    double *Q0 = (double*)malloc(n*sizeof(double));
    double *P0 = (double*)malloc(n*sizeof(double));
    double *Qp = (double*)malloc(n*sizeof(double));
    double *Pp = (double*)malloc(n*sizeof(double));

    T->transform(q, p, n, Q0, P0, T->params);

    for (int j = 0; j < n; j++) {
        memcpy(qp, q, n*sizeof(double));
        memcpy(pp, p, n*sizeof(double));
        qp[j] += eps;
        T->transform(qp, pp, n, Qp, Pp, T->params);
        for (int i = 0; i < n; i++) {
            M[i][j]     = (Qp[i] - Q0[i]) / eps;
            M[n+i][j]   = (Pp[i] - P0[i]) / eps;
        }
        qp[j] = q[j];
        pp[j] += eps;
        T->transform(qp, pp, n, Qp, Pp, T->params);
        for (int i = 0; i < n; i++) {
            M[i][n+j]     = (Qp[i] - Q0[i]) / eps;
            M[n+i][n+j]   = (Pp[i] - P0[i]) / eps;
        }
    }

    double **MT = matrix_alloc(dim, dim);
    double **MTJ = matrix_alloc(dim, dim);
    double **MTJM = matrix_alloc(dim, dim);

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            MT[i][j] = M[j][i];

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) {
            MTJ[i][j] = 0.0;
            for (int k = 0; k < dim; k++)
                MTJ[i][j] += MT[i][k] * J[k][j];
        }

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) {
            MTJM[i][j] = 0.0;
            for (int k = 0; k < dim; k++)
                MTJM[i][j] += MTJ[i][k] * M[k][j];
        }

    double frob = 0.0;
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) {
            double diff = MTJM[i][j] - J[i][j];
            frob += diff * diff;
        }

    matrix_free(M, dim); matrix_free(J, dim);
    matrix_free(MT, dim); matrix_free(MTJ, dim); matrix_free(MTJM, dim);
    free(qp); free(pp); free(Q0); free(P0); free(Qp); free(Pp);

    return sqrt(frob);
}

/*---------------------------------------------------------------------------*
 * L3: Generating function transformation.
 * Maps old canonical variables to new via generating function derivatives.
 * Reference: Goldstein section 9.2
 *---------------------------------------------------------------------------*/
int generating_function_transform(const generating_function_t *gf,
                                  const double *old_a, const double *old_b,
                                  double *new_a, double *new_b)
{
    int n = gf->n_dof;
    double *grad_a = (double*)malloc(n*sizeof(double));
    double *grad_b = (double*)malloc(n*sizeof(double));
    gf->grad_a(old_a, old_b, n, grad_a, gf->params);
    gf->grad_b(old_a, old_b, n, grad_b, gf->params);

    switch (gf->type) {
    case GEN_F1_Qq:
        for (int i=0;i<n;i++) { new_a[i]=grad_a[i]; new_b[i]=-grad_b[i]; }
        break;
    case GEN_F2_qP:
        for (int i=0;i<n;i++) { new_a[i]=grad_a[i]; new_b[i]=grad_b[i]; }
        break;
    case GEN_F3_pQ:
        for (int i=0;i<n;i++) { new_a[i]=-grad_a[i]; new_b[i]=-grad_b[i]; }
        break;
    case GEN_F4_pP:
        for (int i=0;i<n;i++) { new_a[i]=-grad_a[i]; new_b[i]=grad_b[i]; }
        break;
    default:
        free(grad_a); free(grad_b); return -1;
    }

    free(grad_a); free(grad_b);
    return 0;
}

/*---------------------------------------------------------------------------*
 * L4: Liouville theorem -- phase volume preservation.
 * Phase volume is invariant under Hamiltonian flow.
 * Reference: Arnold section 38
 *---------------------------------------------------------------------------*/
double phase_volume_element(const phase_point_t *points, int n_points)
{
    if (n_points < 2) return 0.0;
    int n = points[0].n;
    int dim = 2 * n;

    double *centroid = (double*)calloc(dim, sizeof(double));
    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < n; j++) {
            centroid[j]     += points[i].q[j];
            centroid[n + j] += points[i].p[j];
        }
    }
    for (int j = 0; j < dim; j++)
        centroid[j] /= n_points;

    double **cov = matrix_alloc(dim, dim);
    for (int k = 0; k < n_points; k++) {
        for (int i = 0; i < n; i++) {
            double dq = points[k].q[i] - centroid[i];
            double dp = points[k].p[i] - centroid[n + i];
            for (int j = 0; j < n; j++) {
                double dqj = points[k].q[j] - centroid[j];
                double dpj = points[k].p[j] - centroid[n + j];
                cov[i][j]         += dq * dqj;
                cov[i][n + j]     += dq * dpj;
                cov[n + i][j]     += dp * dqj;
                cov[n + i][n + j] += dp * dpj;
            }
        }
    }

    double *diag = (double*)malloc(dim*sizeof(double));
    double det = 1.0;
    for (int i = 0; i < dim; i++) {
        double sum = cov[i][i];
        for (int k = 0; k < i; k++) sum -= diag[k] * cov[i][k] * cov[i][k];
        diag[i] = sqrt(fabs(sum));
        det *= diag[i];
    }

    double volume = sqrt(fabs(det));
    matrix_free(cov, dim);
    free(centroid);
    free(diag);
    return volume;
}

/*---------------------------------------------------------------------------*
 * L4: Liouville volume conservation check.
 * Evolves cloud of phase points and computes relative volume change.
 * Returns |V(t) - V(0)| / max(V(0), 1).
 * Reference: Arnold section 38
 *---------------------------------------------------------------------------*/
double liouville_volume_check(const hamiltonian_t *H,
                              const phase_point_t *cloud, int n_points,
                              double dt, int steps)
{
    double V0 = phase_volume_element(cloud, n_points);

    phase_point_t *evolved = (phase_point_t*)malloc(n_points * sizeof(phase_point_t));
    for (int i = 0; i < n_points; i++) {
        evolved[i].n = cloud[i].n;
        evolved[i].q = (double*)malloc(cloud[i].n * sizeof(double));
        evolved[i].p = (double*)malloc(cloud[i].n * sizeof(double));
        memcpy(evolved[i].q, cloud[i].q, cloud[i].n * sizeof(double));
        memcpy(evolved[i].p, cloud[i].p, cloud[i].n * sizeof(double));
        hamiltonian_flow_rk4(H, evolved[i].q, evolved[i].p, dt, steps);
    }

    double V1 = phase_volume_element(evolved, n_points);

    for (int i = 0; i < n_points; i++) {
        free(evolved[i].q);
        free(evolved[i].p);
    }
    free(evolved);

    double denom = V0 > 1.0 ? V0 : 1.0;
    return fabs(V1 - V0) / denom;
}
/*---------------------------------------------------------------------------*
 * L5: Action-angle variables for separable 1-DOF Hamiltonian.
 * Action: I = (1/pi) * integral_{q_min}^{q_max} sqrt(2m(E-V(q))) dq
 * Reference: Arnold section 50, Goldstein section 10.7
 *---------------------------------------------------------------------------*/
action_angle_t compute_action_angle_separable(const hamiltonian_t *H,
                                              double energy_level,
                                              double *q_turning_points)
{
    action_angle_t result;
    memset(&result, 0, sizeof(result));
    int n = H->n_dof;
    if (n != 1) { result.converged = 0; return result; }

    double q_min = q_turning_points[0];
    double q_max = q_turning_points[1];
    double q_mid = (q_min + q_max) / 2.0;
    double p_test = 1e-5;
    double p0_arr[1] = {0.0};
    double E_0 = H->energy(&q_mid, p0_arr, 1, H->params);
    double E_up = H->energy(&q_mid, &p_test, 1, H->params);
    double m_est = 0.5 * p_test * p_test / (E_up - E_0 + 1e-12);

    int N_pts = 200;
    double dq = (q_max - q_min) / N_pts;
    double integral = 0.0;

    for (int i = 0; i <= N_pts; i++) {
        double qi = q_min + i * dq;
        double V_qi, p_val = 0.0;
        if (i > 0 && i < N_pts) {
            p0_arr[0] = 0.0;
            V_qi = H->energy(&qi, p0_arr, 1, H->params);
            double T_val = energy_level - V_qi;
            if (T_val > 0) p_val = sqrt(2.0 * m_est * T_val);
        }
        double weight = (i == 0 || i == N_pts) ? 1.0 : (i % 2 == 0) ? 2.0 : 4.0;
        integral += weight * p_val;
    }
    integral *= dq / 3.0;
    result.action = integral / M_PI;

    /* Frequency via finite difference dE/dI */
    double E_low = energy_level * 0.99;
    action_angle_t result_low = compute_action_angle_separable(H, E_low, q_turning_points);
    if (result_low.converged && fabs(result.action - result_low.action) > 1e-12) {
        double dI = result.action - result_low.action;
        double dE = energy_level - E_low;
        result.frequency = dE / dI;
    }
    result.angle = 0.0;
    result.converged = 1;
    return result;
}

/*---------------------------------------------------------------------------*
 * L5: Normal mode analysis for Hamiltonian systems.
 * Linearizes H around equilibrium, computes eigenvalues of J*Hess(H).
 * Stability: all eigenvalues must be purely imaginary.
 * Reference: Goldstein section 6.2
 *---------------------------------------------------------------------------*/
normal_modes_t hamiltonian_normal_modes(const hamiltonian_t *H,
                                        const double *q_eq, const double *p_eq)
{
    normal_modes_t result;
    memset(&result, 0, sizeof(result));
    int n = H->n_dof;
    int dim = 2 * n;

    /* Compute Hessian via finite differences */
    double eps = 1e-5;
    double **JH = matrix_alloc(dim, dim);
    double *z0  = (double*)malloc(dim * sizeof(double));
    double *zp  = (double*)malloc(dim * sizeof(double));
    double *zm  = (double*)malloc(dim * sizeof(double));
    double *gzp_q = (double*)malloc(n * sizeof(double));
    double *gzp_p = (double*)malloc(n * sizeof(double));
    double *gzm_q = (double*)malloc(n * sizeof(double));
    double *gzm_p = (double*)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) { z0[i] = q_eq[i]; z0[n+i] = p_eq[i]; }

    /* Build J*Hess via gradient differences of Hamiltonian vector field */
    double *f0_q = (double*)malloc(n * sizeof(double));
    double *f0_p = (double*)malloc(n * sizeof(double));
    hamiltonian_vector_field(H, q_eq, p_eq, f0_q, f0_p);

    for (int j = 0; j < dim; j++) {
        memcpy(zp, z0, dim * sizeof(double));
        memcpy(zm, z0, dim * sizeof(double));
        zp[j] += eps; zm[j] -= eps;
        double *f_q_p = (double*)malloc(n*sizeof(double));
        double *f_p_p = (double*)malloc(n*sizeof(double));
        double *f_q_m = (double*)malloc(n*sizeof(double));
        double *f_p_m = (double*)malloc(n*sizeof(double));
        hamiltonian_vector_field(H, zp, zp+n, f_q_p, f_p_p);
        hamiltonian_vector_field(H, zm, zm+n, f_q_m, f_p_m);
        for (int i = 0; i < n; i++) {
            JH[i][j]     = (f_q_p[i] - f_q_m[i]) / (2.0 * eps);
            JH[n+i][j]   = (f_p_p[i] - f_p_m[i]) / (2.0 * eps);
        }
        free(f_q_p); free(f_p_p); free(f_q_m); free(f_p_m);
    }

    /* Compute eigenvalues of JH via basic QR iteration */
    double **A = matrix_alloc(dim, dim);
    double **Q = matrix_alloc(dim, dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            A[i][j] = JH[i][j];

    int max_iter = 100;
    for (int iter = 0; iter < max_iter; iter++) {
        /* Gram-Schmidt QR */
        for (int j = 0; j < dim; j++) {
            for (int i = 0; i < dim; i++) Q[i][j] = A[i][j];
            for (int k = 0; k < j; k++) {
                double dot = 0.0;
                for (int i = 0; i < dim; i++) dot += Q[i][j] * Q[i][k];
                for (int i = 0; i < dim; i++) Q[i][j] -= dot * Q[i][k];
            }
            double norm = 0.0;
            for (int i = 0; i < dim; i++) norm += Q[i][j] * Q[i][j];
            norm = sqrt(norm);
            if (norm > 1e-15)
                for (int i = 0; i < dim; i++) Q[i][j] /= norm;
        }
        /* A = Q^T * A * Q (Rayleigh quotient iteration) */
        double **temp = matrix_alloc(dim, dim);
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) {
                temp[i][j] = 0.0;
                for (int k = 0; k < dim; k++)
                    temp[i][j] += A[i][k] * Q[k][j];
            }
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) {
                A[i][j] = 0.0;
                for (int k = 0; k < dim; k++)
                    A[i][j] += Q[k][i] * temp[k][j];
            }
        matrix_free(temp, dim);
    }

    result.n_modes = n;
    result.frequencies = (double*)malloc(n * sizeof(double));
    result.eigenvectors = matrix_alloc(n, dim);

    double max_re = 0.0;
    for (int i = 0; i < n; i++) {
        result.frequencies[i] = fabs(A[i][i]);
        if (fabs(A[i][i]) > max_re) max_re = fabs(A[i][i]);
    }
    result.stability = (max_re < 1e-6) ? 1.0 : -1.0;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < dim; j++)
            result.eigenvectors[i][j] = Q[j][2*i];

    matrix_free(JH, dim); matrix_free(A, dim); matrix_free(Q, dim);
    free(z0); free(zp); free(zm);
    free(gzp_q); free(gzp_p); free(gzm_q); free(gzm_p);
    free(f0_q); free(f0_p);

    return result;
}

/*---------------------------------------------------------------------------*
 * Utility functions: phase point, symplectic matrix, linear algebra
 *---------------------------------------------------------------------------*/

phase_point_t *phase_point_alloc(int n)
{
    phase_point_t *pt = (phase_point_t*)malloc(sizeof(phase_point_t));
    pt->n = n;
    pt->q = (double*)calloc(n, sizeof(double));
    pt->p = (double*)calloc(n, sizeof(double));
    return pt;
}

void phase_point_free(phase_point_t *pt)
{
    if (pt) { free(pt->q); free(pt->p); free(pt); }
}

symplectic_matrix_t *symplectic_matrix_alloc(int n)
{
    symplectic_matrix_t *J = (symplectic_matrix_t*)malloc(sizeof(symplectic_matrix_t));
    J->dim2n = 2 * n;
    J->data = matrix_alloc(J->dim2n, J->dim2n);
    return J;
}

void symplectic_matrix_free(symplectic_matrix_t *J)
{
    if (J) { matrix_free(J->data, J->dim2n); free(J); }
}

void symplectic_matrix_standard(symplectic_matrix_t *J)
{
    int d = J->dim2n, n = d / 2;
    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
            J->data[i][j] = 0.0;
    for (int i = 0; i < n; i++) {
        J->data[i][n+i]     =  1.0;
        J->data[n+i][i]     = -1.0;
    }
}

double **matrix_alloc(int rows, int cols)
{
    double **m = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++)
        m[i] = (double*)calloc(cols, sizeof(double));
    return m;
}

void matrix_free(double **m, int rows)
{
    if (m) { for (int i = 0; i < rows; i++) free(m[i]); free(m); }
}

void matvec_mul(double **A, const double *x, double *y, int m, int n)
{
    for (int i = 0; i < m; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) y[i] += A[i][j] * x[j];
    }
}

double vec_norm2(const double *v, int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += v[i] * v[i];
    return sqrt(sum);
}

double vec_dot(const double *a, const double *b, int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += a[i] * b[i];
    return sum;
}

/*---------------------------------------------------------------------------*
 * L5: Cholesky decomposition for SPD matrices.
 * A = L * L^T. Forward/back substitution for solve.
 * Complexity: O(n^3/3).
 * Reference: Golub & Van Loan, Matrix Computations section 4.2
 *---------------------------------------------------------------------------*/
int solve_spd_cholesky(double **A, const double *b, double *x, int n)
{
    double **L = matrix_alloc(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = A[i][j];
            for (int k = 0; k < j; k++) sum -= L[i][k] * L[j][k];
            if (i == j) {
                if (sum <= 0.0) { matrix_free(L, n); return -1; }
                L[i][j] = sqrt(sum);
            } else {
                L[i][j] = sum / L[j][j];
            }
        }
    }
    double *y = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) y[i] -= L[i][j] * y[j];
        y[i] /= L[i][i];
    }
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) x[i] -= L[j][i] * x[j];
        x[i] /= L[i][i];
    }
    free(y);
    matrix_free(L, n);
    return 0;
}

/*---------------------------------------------------------------------------*
 * L5: Symmetric eigenvalue decomposition via Jacobi method.
 * Uses classical Jacobi rotation to diagonalize symmetric matrix.
 * Returns eigenvalues and eigenvectors.
 * Complexity: O(n^3) per sweep.
 * Reference: Golub & Van Loan section 8.4
 *---------------------------------------------------------------------------*/
int sym_eigen(double **A, int n, double *eigenvalues, double **eigenvectors,
              int max_iter, double tol)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            eigenvectors[i][j] = (i == j) ? 1.0 : 0.0;
    }
    double **B = matrix_alloc(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            B[i][j] = A[i][j];

    for (int iter = 0; iter < max_iter; iter++) {
        double max_off = 0.0; int p = 0, q = 1;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (fabs(B[i][j]) > max_off) {
                    max_off = fabs(B[i][j]); p = i; q = j;
                }
        if (max_off < tol) break;

        double theta;
        if (fabs(B[p][p] - B[q][q]) < 1e-15)
            theta = M_PI / 4.0;
        else
            theta = 0.5 * atan2(2.0 * B[p][q], B[p][p] - B[q][q]);

        double c = cos(theta), s = sin(theta);
        for (int i = 0; i < n; i++) {
            if (i != p && i != q) {
                double bip = B[i][p], biq = B[i][q];
                B[i][p] = c*bip - s*biq; B[p][i] = B[i][p];
                B[i][q] = s*bip + c*biq; B[q][i] = B[i][q];
            }
        }
        double bpp = B[p][p], bqq = B[q][q], bpq = B[p][q];
        B[p][p] = c*c*bpp + s*s*bqq - 2.0*c*s*bpq;
        B[q][q] = s*s*bpp + c*c*bqq + 2.0*c*s*bpq;
        B[p][q] = B[q][p] = (c*c - s*s)*bpq + c*s*(bpp - bqq);

        for (int i = 0; i < n; i++) {
            double vip = eigenvectors[i][p], viq = eigenvectors[i][q];
            eigenvectors[i][p] = c*vip - s*viq;
            eigenvectors[i][q] = s*vip + c*viq;
        }
    }
    for (int i = 0; i < n; i++) eigenvalues[i] = B[i][i];
    matrix_free(B, n);
    return 0;
}
