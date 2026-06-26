/*=============================================================================
 * energy_shaping.c -- Energy shaping control (IDA-PBC, Controlled Lagrangians)
 * Reference: Ortega et al., IDA-PBC, Automatica 2002
 *            Bloch, Leonard, Marsden, Controlled Lagrangians (2000)
 *===========================================================================*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "hamiltonian_control.h"
#include "port_hamiltonian.h"
#include "energy_shaping.h"

/* L2: IDA-PBC matching residual.
 * g_perp(x) * ([J-R]grad(H) - [J_d-R_d]grad(H_d)) = 0 */
double ida_pbc_matching_residual(const ida_pbc_problem_t *prob,
                                  const double *x)
{
    int n = prob->original_system->state_dim, ndof = n/2;
    port_hamiltonian_system_t *ph = prob->original_system;

    /* Compute grad(H) */
    double *gH_q = (double*)malloc(ndof*sizeof(double));
    double *gH_p = (double*)malloc(ndof*sizeof(double));
    ph->H.gradient(x, x+ndof, ndof, gH_q, gH_p, ph->H.params);

    /* Compute grad(H_d) */
    double *gHd_q = (double*)malloc(ndof*sizeof(double));
    double *gHd_p = (double*)malloc(ndof*sizeof(double));
    prob->H_desired.gradient(x, x+ndof, ndof, gHd_q, gHd_p, prob->H_desired.params);

    /* Build [J-R] and [J_d-R_d] */
    double **JR = matrix_alloc(n, n);
    double **JRd = matrix_alloc(n, n);
    double **J_tmp = matrix_alloc(n, n);
    double **R_tmp = matrix_alloc(n, n);

    ph->J_matrix(x, n, J_tmp, ph->params);
    ph->R_matrix(x, n, R_tmp, ph->params);
    for (int i=0;i<n;i++) for(int j=0;j<n;j++) JR[i][j]=J_tmp[i][j]-R_tmp[i][j];

    prob->J_desired(x, n, J_tmp, prob->params);
    prob->R_desired(x, n, R_tmp, prob->params);
    for (int i=0;i<n;i++) for(int j=0;j<n;j++) JRd[i][j]=J_tmp[i][j]-R_tmp[i][j];

    /* grad_H, grad_Hd as single vectors */
    double *gH = (double*)malloc(n*sizeof(double));
    double *gHd = (double*)malloc(n*sizeof(double));
    for(int i=0;i<ndof;i++){gH[i]=gH_q[i];gH[ndof+i]=gH_p[i];}
    for(int i=0;i<ndof;i++){gHd[i]=gHd_q[i];gHd[ndof+i]=gHd_p[i];}

    /* diff = [J-R]gH - [J_d-R_d]gHd */
    double *diff = (double*)malloc(n*sizeof(double));
    for(int i=0;i<n;i++){
        diff[i]=0;
        for(int j=0;j<n;j++) diff[i]+=JR[i][j]*gH[j]-JRd[i][j]*gHd[j];
    }

    /* Compute orthogonal complement of g (simplified: use identity as g_perp) */
    /* For a single-input system, g_perp = [0 1]^T if g = [1 0]^T */
    double residual = 0.0;
    for(int i=0;i<n;i++) residual += diff[i]*diff[i];

    matrix_free(JR,n);matrix_free(JRd,n);matrix_free(J_tmp,n);matrix_free(R_tmp,n);
    free(gH_q);free(gH_p);free(gHd_q);free(gHd_p);free(gH);free(gHd);free(diff);
    return sqrt(residual);
}

/* L3: Potential energy shaping for mechanical systems */
int potential_energy_shaping(const ida_pbc_problem_t *prob,
                              const potential_shaping_t *ps,
                              const double *x_eq,
                              double *V_d_at_eq)
{
    int ndof = ps->n_dof;
    /* Desired potential: V_d(q) = sum coeffs[k] * ||q - q_eq||^{2k}
     * With minimum at q_eq. */
    *V_d_at_eq = 0.0;
    for (int k = 0; k < ps->poly_degree; k++) {
        double r2 = 0.0;
        for (int i = 0; i < ndof; i++) {
            double dq = x_eq[i]; // q_eq
            r2 += dq * dq;
        }
        *V_d_at_eq += ps->V_d_coeffs[k] * pow(r2, k+1);
    }
    return 0;
}

/* L5: Kinetic energy shaping -- mass matrix modification */
int kinetic_energy_shaping_setup(const ida_pbc_problem_t *prob,
                                  const double *q,
                                  kinetic_shaping_t *ks)
{
    int ndof = ks->n_dof;
    /* For mechanical systems, desired mass M_d(q) = M(q) + M_a(q)
     * M_a is the assigned mass modification.
     * Setup M_d as identity plus shaping. */
    for (int i = 0; i < ndof; i++)
        for (int j = 0; j < ndof; j++)
            ks->M_d[i][j] = (i==j) ? 1.0 : 0.0;
    (void)prob; (void)q;
    return 0;
}

/* L5: Damping assignment R_d = R + R_a, R_a = g K_v g^T */
int damping_assignment(const port_hamiltonian_system_t *ph,
                        const double *x, double K_v,
                        double **R_a_out)
{
    int n = ph->state_dim, m = ph->port_dim;
    double **g_mat = matrix_alloc(n, m);
    ph->g_matrix(x, n, m, g_mat, ph->params);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            R_a_out[i][j] = 0.0;
            for (int k = 0; k < m; k++)
                R_a_out[i][j] += K_v * g_mat[i][k] * g_mat[j][k];
        }

    matrix_free(g_mat, n);
    return 0;
}

/* L5: Total energy shaping control */
int total_energy_shaping_control(const ida_pbc_problem_t *prob,
                                  const double *x,
                                  const potential_shaping_t *ps,
                                  const kinetic_shaping_t *ks,
                                  double K_d,
                                  double *u_out)
{
    int ndof = ps->n_dof, n = 2*ndof;
    port_hamiltonian_system_t *ph = prob->original_system;

    /* Compute grad(H_d - H) */
    double *gH_q = (double*)malloc(ndof*sizeof(double));
    double *gH_p = (double*)malloc(ndof*sizeof(double));
    double *gHd_q = (double*)malloc(ndof*sizeof(double));
    double *gHd_p = (double*)malloc(ndof*sizeof(double));
    ph->H.gradient(x, x+ndof, ndof, gH_q, gH_p, ph->H.params);
    prob->H_desired.gradient(x, x+ndof, ndof, gHd_q, gHd_p, prob->H_desired.params);

    /* Energy shaping term: u_es = g^+ * ([J_d-R_d]gHd - [J-R]gH) */
    /* Simplified for single-input: u_es = [J_d-R_d]gHd - [J-R]gH projected */
    /* For 1-DOF: u_es = -K_p(q-q*) */
    double u_es = 0.0;
    int m = ph->port_dim;
    for (int j = 0; j < m; j++) {
        u_es = 0.0;
        /* for 1-DOF: u_es = K_p*(x - x_eq) */
        for (int i = 0; i < ndof; i++) {
            u_es += (gHd_q[i] - gH_q[i]) * (x[i] - x[ndof+i]); // rough approx
        }
        u_out[j] = u_es;
    }

    /* Damping injection: u_di = -K_d * y */
    double *y_out = (double*)malloc(m * sizeof(double));
    pbc_damping_injection(ph, x, K_d, u_out, y_out);

    free(gH_q);free(gH_p);free(gHd_q);free(gHd_p);free(y_out);
    return 0;
}

/* L6: Simple pendulum IDA-PBC controller */
int pendulum_ida_pbc_control(const pendulum_ida_params_t *params,
                              const double *x,
                              double *u_out, double *H_d_val)
{
    double q = x[0], p = x[1];
    double m = params->mass, l = params->length, g = params->gravity;
    double ml2 = m * l * l;
    double q_star = params->q_desired;

    /* Original H = p^2/(2ml^2) + mgl(1-cos(q)) */
    /* Desired H_d = p^2/(2*ml^2*k_i) + k_p*mgl*(1-cos(q - q_star)) */
    double H_d_kinetic = p*p / (2.0 * ml2 * params->k_i);
    double H_d_potential = params->k_p * m * g * l * (1.0 - cos(q - q_star));
    *H_d_val = H_d_kinetic + H_d_potential;

    /* Control law: u = grad(H)(q) partial derivatives matching */
    double u_es = m*g*l*(sin(q) - params->k_p * sin(q - q_star));
    /* Damping injection: u_di = -k_d * p/(ml^2) */
    double u_di = -params->k_d * p / ml2;
    *u_out = u_es + u_di;

    return 0;
}

/* L6: Cart-pole IDA-PBC controller */
int cartpole_ida_pbc_control(const cartpole_ida_params_t *params,
                              const double *q, const double *p,
                              double *u_out)
{
    double M_c = params->M_cart, m_p = params->m_pole;
    double l = params->l_pole, g = params->gravity;
    double theta = q[1]; /* pendulum angle (0 = upright) */

    /* Total energy of the system */
    double E = 0.0;
    /* Kinetic energy from momentum */
    double Mt = M_c + m_p;
    double cos_th = cos(theta), sin_th = sin(theta);
    E += p[0]*p[0]/(2.0*Mt) + p[1]*p[1]/(2.0*m_p*l*l);

    /* Potential energy: m_p*g*l*cos(theta) */
    E += m_p * g * l * cos_th;

    /* Energy shaping: drive to homoclinic orbit (upright equilibrium) */
    double E_desired = m_p * g * l; /* energy at upright */
    double u_es = params->k_e * (E - E_desired) * p[0] / Mt;

    /* Damping injection */
    double u_di = -params->k_d * p[0] / Mt;

    *u_out = u_es + u_di;
    return 0;
}

/* L8: Controlled Lagrangian -- derive control from modified Lagrangian */
int controlled_lagrangian_derive(const controlled_lagrangian_t *cl,
                                  const double *q, const double *dq,
                                  double *u_ctrl, double *E_total)
{
    int n = cl->n_dof;
    double L_val, Lc_val;
    double *grad_q_L = (double*)malloc(n*sizeof(double));
    double *grad_dq_L = (double*)malloc(n*sizeof(double));
    double *grad_q_Lc = (double*)malloc(n*sizeof(double));
    double *grad_dq_Lc = (double*)malloc(n*sizeof(double));

    cl->L(q, dq, n, &L_val, grad_q_L, grad_dq_L, cl->params);
    cl->L_control(q, dq, n, &Lc_val, grad_q_Lc, grad_dq_Lc, cl->params);

    /* Euler-Lagrange for L_d = L + L_c:
     * d/dt(dL_d/d(dq)) - dL_d/dq = 0
     * The control appears as the difference from original EL:
     * u = d/dt(dL_c/d(dq)) - dL_c/dq */
    for (int i = 0; i < cl->n_controls; i++) {
        /* Simplified: u = grad_dq_Lc contribution */
        u_ctrl[i] = grad_dq_Lc[i];
    }

    *E_total = 0.0;
    for (int i = 0; i < n; i++)
        *E_total += dq[i] * (grad_dq_L[i] + grad_dq_Lc[i]);
    *E_total -= (L_val + Lc_val);

    free(grad_q_L);free(grad_dq_L);free(grad_q_Lc);free(grad_dq_Lc);
    return 0;
}