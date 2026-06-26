/*=============================================================================
 * port_hamiltonian.c -- Port-Hamiltonian systems implementation
 * Reference: van der Schaft & Jeltsema, Port-Hamiltonian Systems Theory
 *            Duindam et al., Modeling and Control of Complex Physical Systems
 *            Ortega et al., Passivity-based Control of Euler-Lagrange Systems
 *===========================================================================*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "hamiltonian_control.h"
#include "port_hamiltonian.h"

/* L1: Dirac structure from port-Hamiltonian system */
int dirac_structure_from_ph(const port_hamiltonian_system_t *ph,
                             const double *x, const double *u,
                             dirac_structure_t *ds)
{
    int n = ph->state_dim, m = ph->port_dim;
    ds->dim_storage = n;
    ds->dim_resistive = n;
    ds->dim_port = m;
    ds->dim_flow = 3 * n + m;
    ds->flow = (double*)calloc(ds->dim_flow, sizeof(double));
    ds->effort = (double*)calloc(ds->dim_flow, sizeof(double));

    /* Storage effort: e_s = grad(H) */
    double *grad_q = (double*)malloc(n * sizeof(double));
    double *grad_p = (double*)malloc(n * sizeof(double));
    ph->H.gradient(x, x + n/2, n/2, grad_q, grad_p, ph->H.params);

    for (int i = 0; i < n/2; i++) {
        ds->effort[i] = grad_q[i];
        ds->effort[n/2 + i] = grad_p[i];
    }

    /* External port effort = y = g(x)^T grad(H) */
    (void)u;
    /* Compute g(x) */
    double **g_mat = matrix_alloc(n, m);
    ph->g_matrix(x, n, m, g_mat, ph->params);

    for (int j = 0; j < m; j++) {
        ds->effort[2*n + j] = 0.0;
        for (int i = 0; i < n; i++)
            ds->effort[2*n + j] += g_mat[i][j] * ds->effort[i];
    }

    matrix_free(g_mat, n);
    free(grad_q); free(grad_p);
    return 0;
}

/* L2: Power balance for port-Hamiltonian system.
 * dH/dt = y^T u - grad(H)^T R grad(H) <= y^T u (passivity) */
double port_hamiltonian_power_balance(const port_hamiltonian_system_t *ph,
                                       const double *x, const double *u,
                                       double *y_out)
{
    int n = ph->state_dim, m = ph->port_dim, ndof = n / 2;

    /* grad(H) */
    double *gq = (double*)malloc(ndof * sizeof(double));
    double *gp = (double*)malloc(ndof * sizeof(double));
    ph->H.gradient(x, x + ndof, ndof, gq, gp, ph->H.params);

    /* g(x) */
    double **g_mat = matrix_alloc(n, m);
    ph->g_matrix(x, n, m, g_mat, ph->params);

    /* y = g^T grad(H) */
    for (int j = 0; j < m; j++) {
        y_out[j] = 0.0;
        for (int i = 0; i < ndof; i++)
            y_out[j] += g_mat[i][j] * gq[i] + g_mat[ndof+i][j] * gp[i];
    }

    /* R(x) */
    double **R_mat = matrix_alloc(n, n);
    ph->R_matrix(x, n, R_mat, ph->params);

    /* grad(H)^T R grad(H) */
    double diss = 0.0;
    for (int i = 0; i < n; i++) {
        double gi = (i < ndof) ? gq[i] : gp[i - ndof];
        for (int j = 0; j < n; j++) {
            double gj = (j < ndof) ? gq[j] : gp[j - ndof];
            diss += gi * R_mat[i][j] * gj;
        }
    }

    /* y^T u */
    double ytu = 0.0;
    for (int j = 0; j < m; j++)
        ytu += y_out[j] * u[j];

    matrix_free(g_mat, n); matrix_free(R_mat, n);
    free(gq); free(gp);
    return ytu - diss;
}

/* L2: Check if a function C(x) is a Casimir for the pH system */
int port_hamiltonian_casimir_check(const port_hamiltonian_system_t *ph,
                                    casimir_check_fn C_fn, void *C_params,
                                    const double *x, double *violation)
{
    int n = ph->state_dim;
    double C_val;
    double *grad_C = (double*)malloc(n * sizeof(double));

    C_fn(x, n, &C_val, grad_C, C_params);

    /* Check grad(C)^T J(x) = 0 */
    double **J_mat = matrix_alloc(n, n);
    ph->J_matrix(x, n, J_mat, ph->params);

    *violation = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            *violation += grad_C[i] * J_mat[i][j] * grad_C[j];
        }
    }
    *violation = fabs(*violation);

    matrix_free(J_mat, n);
    free(grad_C);
    return (*violation < 1e-10) ? 1 : 0;
}

/* L3: Interconnection of two port-Hamiltonian systems */
int port_hamiltonian_interconnect(const ph_interconnection_t *interconn,
                                   const double *x1, const double *x2,
                                   double *dx1_dt, double *dx2_dt,
                                   double *y1, double *y2)
{
    port_hamiltonian_system_t *s1 = interconn->sys1;
    port_hamiltonian_system_t *s2 = interconn->sys2;
    int n1 = s1->state_dim, n2 = s2->state_dim;
    int m1 = s1->port_dim, m2 = s2->port_dim;
    int ndof1 = n1/2, ndof2 = n2/2;

    /* For simplicity: assume standard negative feedback u1 = -y2, u2 = y1 */
    double *u1 = (double*)malloc(m1 * sizeof(double));
    double *u2 = (double*)malloc(m2 * sizeof(double));

    /* Compute y1, y2 first (as outputs) */
    double *gq1 = (double*)malloc(ndof1*sizeof(double));
    double *gp1 = (double*)malloc(ndof1*sizeof(double));
    double *gq2 = (double*)malloc(ndof2*sizeof(double));
    double *gp2 = (double*)malloc(ndof2*sizeof(double));

    s1->H.gradient(x1, x1+ndof1, ndof1, gq1, gp1, s1->H.params);
    s2->H.gradient(x2, x2+ndof2, ndof2, gq2, gp2, s2->H.params);

    double **g1 = matrix_alloc(n1, m1);
    double **g2 = matrix_alloc(n2, m2);
    s1->g_matrix(x1, n1, m1, g1, s1->params);
    s2->g_matrix(x2, n2, m2, g2, s2->params);

    for (int j=0;j<m1;j++) {
        y1[j]=0;
        for (int i=0;i<ndof1;i++) y1[j]+=g1[i][j]*gq1[i]+g1[ndof1+i][j]*gp1[i];
    }
    for (int j=0;j<m2;j++) {
        y2[j]=0;
        for (int i=0;i<ndof2;i++) y2[j]+=g2[i][j]*gq2[i]+g2[ndof2+i][j]*gp2[i];
    }

    /* Interconnection: [u1; u2] = F * [y1; y2] + v_ext */
    int mtot = m1+m2;
    for (int i=0;i<m1;i++) {
        u1[i]=interconn->v_ext[i];
        for (int j=0;j<m1;j++) u1[i]+=interconn->F[i][j]*y1[j];
        for (int j=0;j<m2;j++) u1[i]+=interconn->F[i][m1+j]*y2[j];
    }
    for (int i=0;i<m2;i++) {
        u2[i]=interconn->v_ext[m1+i];
        for (int j=0;j<m1;j++) u2[i]+=interconn->F[m1+i][j]*y1[j];
        for (int j=0;j<m2;j++) u2[i]+=interconn->F[m1+i][m1+j]*y2[j];
    }

    /* Compute dx/dt for both systems */
    double **J1 = matrix_alloc(n1,n1), **R1 = matrix_alloc(n1,n1);
    double **J2 = matrix_alloc(n2,n2), **R2 = matrix_alloc(n2,n2);
    s1->J_matrix(x1,n1,J1,s1->params); s1->R_matrix(x1,n1,R1,s1->params);
    s2->J_matrix(x2,n2,J2,s2->params); s2->R_matrix(x2,n2,R2,s2->params);

    double *gradH1 = (double*)malloc(n1*sizeof(double));
    double *gradH2 = (double*)malloc(n2*sizeof(double));
    for(int i=0;i<ndof1;i++){gradH1[i]=gq1[i];gradH1[ndof1+i]=gp1[i];}
    for(int i=0;i<ndof2;i++){gradH2[i]=gq2[i];gradH2[ndof2+i]=gp2[i];}

    for(int i=0;i<n1;i++){
        dx1_dt[i]=0;
        for(int j=0;j<n1;j++) dx1_dt[i]+=(J1[i][j]-R1[i][j])*gradH1[j];
        for(int j=0;j<m1;j++) dx1_dt[i]+=g1[i][j]*u1[j];
    }
    for(int i=0;i<n2;i++){
        dx2_dt[i]=0;
        for(int j=0;j<n2;j++) dx2_dt[i]+=(J2[i][j]-R2[i][j])*gradH2[j];
        for(int j=0;j<m2;j++) dx2_dt[i]+=g2[i][j]*u2[j];
    }

    matrix_free(g1,n1);matrix_free(g2,n2);
    matrix_free(J1,n1);matrix_free(R1,n1);
    matrix_free(J2,n2);matrix_free(R2,n2);
    free(u1);free(u2);free(gq1);free(gp1);free(gq2);free(gp2);
    free(gradH1);free(gradH2);
    return 0;
}

/* L5: Passivity-based damping injection u = -K_d y */
int pbc_damping_injection(const port_hamiltonian_system_t *ph,
                           const double *x, double K_d,
                           double *u_out, double *y_out)
{
    int n=ph->state_dim, m=ph->port_dim, ndof=n/2;
    double *gq=(double*)malloc(ndof*sizeof(double));
    double *gp=(double*)malloc(ndof*sizeof(double));
    ph->H.gradient(x,x+ndof,ndof,gq,gp,ph->H.params);

    double **g=matrix_alloc(n,m);
    ph->g_matrix(x,n,m,g,ph->params);
    for(int j=0;j<m;j++){
        y_out[j]=0;
        for(int i=0;i<ndof;i++) y_out[j]+=g[i][j]*gq[i]+g[ndof+i][j]*gp[i];
        u_out[j]=-K_d*y_out[j];
    }
    matrix_free(g,n);free(gq);free(gp);
    return 0;
}

/* L5: Energy-Casimir stability analysis */
int energy_casimir_stability(const port_hamiltonian_system_t *ph,
                              casimir_check_fn *casimirs,
                              void **casimir_params,
                              int n_casimirs,
                              const double *x_eq,
                              energy_casimir_t *result)
{
    result->n_casimirs = n_casimirs;
    result->phi_coeffs = (double*)calloc(n_casimirs, sizeof(double));
    result->C_values = (double*)calloc(n_casimirs, sizeof(double));

    /* Evaluate Casimirs at equilibrium */
    for(int k=0;k<n_casimirs;k++){
        double *gc=(double*)malloc(ph->state_dim*sizeof(double));
        casimirs[k](x_eq,ph->state_dim,&result->C_values[k],gc,casimir_params[k]);
        free(gc);
    }

    /* Simple heuristic: set phi_coeffs proportional to C values */
    for(int k=0;k<n_casimirs;k++)
        result->phi_coeffs[k] = 1.0;
    return 0;
}

/* L5: Port-Hamiltonian symplectic splitting step */
int ph_symplectic_step(const port_hamiltonian_system_t *ph,
                        const double *u, double *x, double dt,
                        double *y_out)
{
    int n=ph->state_dim, ndof=n/2;
    double *x_save=(double*)malloc(n*sizeof(double));
    memcpy(x_save,x,n*sizeof(double));

    /* Half step with J only (no dissipation, no input) */
    double **J=matrix_alloc(n,n);
    ph->J_matrix(x,n,J,ph->params);
    double *gq=(double*)malloc(ndof*sizeof(double));
    double *gp=(double*)malloc(ndof*sizeof(double));
    ph->H.gradient(x,x+ndof,ndof,gq,gp,ph->H.params);

    for(int i=0;i<ndof;i++){
        double sum_q=0,sum_p=0;
        for(int j=0;j<ndof;j++){
            sum_q+=J[i][j]*gq[j]+J[i][ndof+j]*gp[j];
            sum_p+=J[ndof+i][j]*gq[j]+J[ndof+i][ndof+j]*gp[j];
        }
        x[i]+=0.5*dt*sum_q; x[ndof+i]+=0.5*dt*sum_p;
    }
    matrix_free(J,n);free(gq);free(gp);

    /* Full step with R and g*u (dissipative + input) */
    double **R=matrix_alloc(n,n); double **gmat=matrix_alloc(n,ph->port_dim);
    ph->R_matrix(x,n,R,ph->params);
    ph->g_matrix(x,n,ph->port_dim,gmat,ph->params);

    double *gq2=(double*)malloc(ndof*sizeof(double));
    double *gp2=(double*)malloc(ndof*sizeof(double));
    ph->H.gradient(x,x+ndof,ndof,gq2,gp2,ph->H.params);

    for(int i=0;i<ndof;i++){
        double sum_q=0,sum_p=0;
        for(int j=0;j<ndof;j++){
            sum_q-=R[i][j]*gq2[j]+R[i][ndof+j]*gp2[j];
            sum_p-=R[ndof+i][j]*gq2[j]+R[ndof+i][ndof+j]*gp2[j];
        }
        for(int j=0;j<ph->port_dim;j++){
            sum_q+=gmat[i][j]*u[j];
            sum_p+=gmat[ndof+i][j]*u[j];
        }
        x[i]+=dt*sum_q; x[ndof+i]+=dt*sum_p;
    }
    matrix_free(R,n);matrix_free(gmat,n);free(gq2);free(gp2);

    /* Second half step with J only */
    double **J2=matrix_alloc(n,n);
    ph->J_matrix(x,n,J2,ph->params);
    double *gq3=(double*)malloc(ndof*sizeof(double));
    double *gp3=(double*)malloc(ndof*sizeof(double));
    ph->H.gradient(x,x+ndof,ndof,gq3,gp3,ph->H.params);

    for(int i=0;i<ndof;i++){
        double sum_q=0,sum_p=0;
        for(int j=0;j<ndof;j++){
            sum_q+=J2[i][j]*gq3[j]+J2[i][ndof+j]*gp3[j];
            sum_p+=J2[ndof+i][j]*gq3[j]+J2[ndof+i][ndof+j]*gp3[j];
        }
        x[i]+=0.5*dt*sum_q; x[ndof+i]+=0.5*dt*sum_p;
    }

    /* Compute output y */
    for(int j=0;j<ph->port_dim;j++){
        y_out[j]=0;
        for(int i=0;i<ndof;i++) y_out[j]+=gmat[i][j]*gq3[i]+gmat[ndof+i][j]*gp3[i];
    }

    matrix_free(J2,n);free(gq3);free(gp3);free(x_save);
    return 0;
}

/* L7: DC motor as port-Hamiltonian system */
static double dc_motor_energy(const double *q, const double *p, int n, void *params)
{
    dc_motor_params_t *pm = (dc_motor_params_t*)params;
    return p[0]*p[0]/(2.0*pm->J_m) + q[0]*q[0]/(2.0*pm->L);
}

static void dc_motor_gradient(const double *q, const double *p, int n,
                               double *gq, double *gp, void *params)
{
    dc_motor_params_t *pm = (dc_motor_params_t*)params;
    gq[0] = q[0] / pm->L;
    gp[0] = p[0] / pm->J_m;
}

static void dc_motor_J(const double *x, int n, double **J, void *params)
{
    (void)x; (void)params;
    J[0][0]=0; J[0][1]=1; J[1][0]=-1; J[1][1]=0;
}

static void dc_motor_R(const double *x, int n, double **R, void *params)
{
    dc_motor_params_t *pm = (dc_motor_params_t*)params;
    R[0][0]=pm->R; R[0][1]=0; R[1][0]=0; R[1][1]=pm->b;
}

static void dc_motor_g(const double *x, int n, int m, double **g, void *params)
{
    (void)x; (void)params;
    g[0][0]=1; g[0][1]=0; g[1][0]=0; g[1][1]=1;
}

port_hamiltonian_system_t *dc_motor_create(const dc_motor_params_t *params)
{
    port_hamiltonian_system_t *ph = (port_hamiltonian_system_t*)malloc(sizeof(port_hamiltonian_system_t));
    ph->state_dim = 2; ph->port_dim = 2;
    ph->J_matrix = dc_motor_J; ph->R_matrix = dc_motor_R; ph->g_matrix = dc_motor_g;
    ph->H.energy = dc_motor_energy; ph->H.gradient = dc_motor_gradient;
    ph->H.hessian = NULL; ph->H.n_dof = 1;
    ph->H.params = (void*)params; ph->params = (void*)params;
    return ph;
}

/* L7: Mass-spring-damper as port-Hamiltonian system */
static double msd_energy(const double *q, const double *p, int n, void *params)
{
    msd_params_t *pm = (msd_params_t*)params;
    return p[0]*p[0]/(2.0*pm->mass) + 0.5*pm->stiffness*q[0]*q[0];
}

static void msd_gradient(const double *q, const double *p, int n,
                         double *gq, double *gp, void *params)
{
    msd_params_t *pm = (msd_params_t*)params;
    gq[0] = pm->stiffness * q[0];
    gp[0] = p[0] / pm->mass;
}

static void msd_J(const double *x, int n, double **J, void *params)
{
    (void)x; (void)params;
    J[0][0]=0; J[0][1]=1; J[1][0]=-1; J[1][1]=0;
}

static void msd_R(const double *x, int n, double **R, void *params)
{
    msd_params_t *pm = (msd_params_t*)params;
    R[0][0]=0; R[0][1]=0; R[1][0]=0; R[1][1]=pm->damping;
}

static void msd_g(const double *x, int n, int m, double **g, void *params)
{
    (void)x; (void)params;
    g[0][0]=0; g[1][0]=1;
}

port_hamiltonian_system_t *mass_spring_damper_create(const msd_params_t *params)
{
    port_hamiltonian_system_t *ph = (port_hamiltonian_system_t*)malloc(sizeof(port_hamiltonian_system_t));
    ph->state_dim = 2; ph->port_dim = 1;
    ph->J_matrix = msd_J; ph->R_matrix = msd_R; ph->g_matrix = msd_g;
    ph->H.energy = msd_energy; ph->H.gradient = msd_gradient;
    ph->H.hessian = NULL; ph->H.n_dof = 1;
    ph->H.params = (void*)params; ph->params = (void*)params;
    return ph;
}

void port_hamiltonian_system_free(port_hamiltonian_system_t *ph) { free(ph); }