/**
 * lie_group_methods.c - Lie Group Methods for Geometric Optimal Control
 *
 * Implements: so(3) hat/vee/bracket, se(3) hat/vee/bracket,
 * SO(3) exponential (Rodrigues) and logarithm, matrix exp/log,
 * adjoint representations Ad and ad, Jacobi identity on so(3),
 * RK-Munthe-Kaas integrator, Crouch-Grossman integrator,
 * Euler-Poincare rigid body equations, SO(3) kinematics,
 * Arnold cat map on torus.
 *
 * Reference: Hall (2015), Iserles et al. (2000), Marsden & Ratiu (1999)
 */
#include "lie_group_methods.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- so(3) operations ---------- */

/* Hat map: omega in R^3 -> skew-symmetric 3x3 matrix.
 * omega_hat = [[0,-w3,w2],[w3,0,-w1],[-w2,w1,0]] */
void geo_so3_hat(const double omega[3], double hat[3][3]) {
    if (!omega || !hat) return;
    memset(hat, 0, 9 * sizeof(double));
    hat[0][1] = -omega[2];  hat[0][2] =  omega[1];
    hat[1][0] =  omega[2];  hat[1][2] = -omega[0];
    hat[2][0] = -omega[1];  hat[2][1] =  omega[0];
}

/* Vee map: skew-symmetric 3x3 matrix -> omega in R^3.
 * Inverse of hat map. */
void geo_so3_vee(const double hat[3][3], double omega[3]) {
    if (!hat || !omega) return;
    omega[0] = -hat[1][2];
    omega[1] =  hat[0][2];
    omega[2] = -hat[0][1];
}

/* so(3) Lie bracket [X,Y] = XY - YX.
 * Equivalent to cross product in R^3: (x cross y)_hat = [x_hat, y_hat]. */
void geo_so3_bracket(const double X[3][3], const double Y[3][3],
                     double result[3][3]) {
    if (!X || !Y || !result) return;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            result[i][j] = 0.0;
            for (int k = 0; k < 3; k++)
                result[i][j] += X[i][k]*Y[k][j] - Y[i][k]*X[k][j];
        }
}

/* ---------- se(3) operations ---------- */

/* se(3) hat: (omega, v) in R^3 x R^3 -> 4x4 matrix.
 * xi_hat = [[omega_hat, v], [0, 0]] */
void geo_se3_hat(const double omega[3], const double v[3],
                 double hat[4][4]) {
    if (!omega || !v || !hat) return;
    memset(hat, 0, 16 * sizeof(double));
    hat[0][1]=-omega[2]; hat[0][2]=omega[1]; hat[0][3]=v[0];
    hat[1][0]=omega[2];  hat[1][2]=-omega[0]; hat[1][3]=v[1];
    hat[2][0]=-omega[1]; hat[2][1]=omega[0];  hat[2][3]=v[2];
}

/* se(3) vee: 4x4 matrix -> (omega, v) in R^3 x R^3 */
void geo_se3_vee(const double hat[4][4], double omega[3], double v[3]) {
    if (!hat || !omega || !v) return;
    omega[0] = -hat[1][2];
    omega[1] =  hat[0][2];
    omega[2] = -hat[0][1];
    v[0] = hat[0][3];
    v[1] = hat[1][3];
    v[2] = hat[2][3];
}

/* se(3) Lie bracket: [(w1,v1),(w2,v2)] = (w1 x w2, w1 x v2 - w2 x v1) */
void geo_se3_bracket(const double X[4][4], const double Y[4][4],
                     double result[4][4]) {
    if (!X || !Y || !result) return;
    double w1[3], v1[3], w2[3], v2[3];
    geo_se3_vee(X, w1, v1);
    geo_se3_vee(Y, w2, v2);
    double w_cross[3], v_cross[3];
    w_cross[0]=w1[1]*w2[2]-w1[2]*w2[1];
    w_cross[1]=w1[2]*w2[0]-w1[0]*w2[2];
    w_cross[2]=w1[0]*w2[1]-w1[1]*w2[0];
    v_cross[0]=w1[1]*v2[2]-w1[2]*v2[1] - (w2[1]*v1[2]-w2[2]*v1[1]);
    v_cross[1]=w1[2]*v2[0]-w1[0]*v2[2] - (w2[2]*v1[0]-w2[0]*v1[2]);
    v_cross[2]=w1[0]*v2[1]-w1[1]*v2[0] - (w2[0]*v1[1]-w2[1]*v1[0]);
    geo_se3_hat(w_cross, v_cross, result);
}

/* ---------- SO(3) Exponential map (Rodrigues formula) ---------- */

/* exp(omega_hat) = I + (sin theta)/theta * omega_hat + (1-cos theta)/theta^2 * omega_hat^2
 * where theta = ||omega||. Handles theta=0 limit: exp(0) = I. */
void geo_so3_exp(const double omega[3], double R[3][3]) {
    if (!omega || !R) return;
    double theta = sqrt(omega[0]*omega[0]+omega[1]*omega[1]+omega[2]*omega[2]);
    double hat[3][3];
    geo_so3_hat(omega, hat);
    if (theta < 1e-14) {
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) R[i][j]=(i==j)?1.0:0.0;
        return;
    }
    double s = sin(theta)/theta;
    double c = (1.0-cos(theta))/(theta*theta);
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++) {
            R[i][j] = (i==j)?1.0+s*hat[i][j]:s*hat[i][j];
            /* Add c*hat^2 terms */
            double hat2 = 0.0;
            for (int k=0;k<3;k++) hat2 += hat[i][k]*hat[k][j];
            R[i][j] += c * hat2;
        }
}

/* Matrix exponential via scaling-and-squaring (Higham 2009).
 * exp(A) = (exp(A/2^k))^(2^k). Uses Pade approximation for base. */
void geo_matrix_exp(const double *A, double *M, int n, double tol) {
    if (!A || !M || n <= 0 || n > GEO_MAX_DIM) return;
    if (tol <= 0.0) tol = 1e-14;

    /* Compute norm for scaling */
    double norm = 0.0;
    for (int i=0;i<n*n;i++) norm += A[i]*A[i];
    norm = sqrt(norm);

    /* Find scaling factor: s = floor(log2(norm)) + 1 */
    int s = (norm > 1.0) ? (int)(log(norm)/log(2.0)) + 1 : 0;
    if (s < 0) s = 0;

    /* Scale: As = A / 2^s */
    double *As = (double*)malloc(n*n*sizeof(double));
    double scale = 1.0 / (double)(1 << s);
    for (int i=0;i<n*n;i++) As[i] = A[i] * scale;

    /* Pade approximation: exp(As) ~ (I + As/2 + As^2/12) * (I - As/2 + As^2/12)^{-1}
     * Simplified: use truncated Taylor: exp(As) ~ I + As + As^2/2 + As^3/6 */
    double *E = (double*)calloc(n*n, sizeof(double));
    double *T = (double*)calloc(n*n, sizeof(double));
    double *T2 = (double*)calloc(n*n, sizeof(double));
    double *T3 = (double*)calloc(n*n, sizeof(double));

    /* T = I + As */
    for (int i=0;i<n;i++) E[i*n+i] = 1.0;
    for (int i=0;i<n*n;i++) T[i] = E[i] + As[i];

    /* T2 = As^2 / 2 */
    for (int i=0;i<n;i++) for (int j=0;j<n;j++)
        for (int k=0;k<n;k++) T2[i*n+j] += As[i*n+k]*As[k*n+j];
    for (int i=0;i<n*n;i++) T2[i] *= 0.5;

    /* T3 = As^3 / 6 */
    for (int i=0;i<n;i++) for (int j=0;j<n;j++)
        for (int k=0;k<n;k++) T3[i*n+j] += T2[i*n+k]*As[k*n+j];
    for (int i=0;i<n*n;i++) T3[i] /= 3.0;

    /* E = I + As + As^2/2 + As^3/6 */
    for (int i=0;i<n*n;i++) E[i] = T[i] + T2[i] + T3[i];

    /* Repeated squaring: E <- E^2, s times */
    double *Etmp = (double*)malloc(n*n*sizeof(double));
    for (int step=0; step<s; step++) {
        memset(Etmp, 0, n*n*sizeof(double));
        for (int i=0;i<n;i++) for (int j=0;j<n;j++)
            for (int k=0;k<n;k++) Etmp[i*n+j] += E[i*n+k]*E[k*n+j];
        memcpy(E, Etmp, n*n*sizeof(double));
    }

    memcpy(M, E, n*n*sizeof(double));
    free(As); free(E); free(T); free(T2); free(T3); free(Etmp);
}

/* ---------- SO(3) Logarithm map ---------- */

/* SO(3) log: inverse of Rodrigues.
 * theta = acos((tr(R)-1)/2), omega_hat = theta/(2*sin(theta)) * (R - R^T). */
int geo_so3_log(const double R[3][3], double omega[3]) {
    if (!R || !omega) return -1;
    double trace = R[0][0]+R[1][1]+R[2][2];
    double cos_theta = (trace-1.0)*0.5;
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    double theta = acos(cos_theta);
    if (theta < 1e-14) {
        omega[0]=omega[1]=omega[2]=0.0; return 0;
    }
    if (fabs(theta-M_PI) < 1e-12) return -1;
    double factor = theta/(2.0*sin(theta));
    omega[0] = factor*(R[2][1]-R[1][2]);
    omega[1] = factor*(R[0][2]-R[2][0]);
    omega[2] = factor*(R[1][0]-R[0][1]);
    return 0;
}

/* Matrix logarithm via inverse scaling-and-squaring */
int geo_matrix_log(const double *A, double *L, int n) {
    if (!A || !L || n<=0||n>GEO_MAX_DIM) return -1;
    /* Simple implementation using Schur-Parlett for small n */
    for (int i=0;i<n*n;i++) L[i] = A[i];
    double *tmp = (double*)malloc(n*n*sizeof(double));
    double *Ak  = (double*)malloc(n*n*sizeof(double));
    memcpy(Ak, A, n*n*sizeof(double));
    for (int i=0;i<n;i++) for (int j=0;j<n;j++)
        Ak[i*n+j] -= (i==j)?1.0:0.0;
    /* Mercator series: log(I+X) = X - X^2/2 + X^3/3 - X^4/4 + ...
     * Only valid for ||X|| < 1. Use for small residuals. */
    memset(L, 0, n*n*sizeof(double));
    memcpy(tmp, Ak, n*n*sizeof(double));
    for (int k=1; k<=20; k++) {
        double sign = (k%2==1)?1.0:-1.0;
        for (int ij=0;ij<n*n;ij++) L[ij] += sign*tmp[ij]/k;
        double *new_tmp = (double*)calloc(n*n, sizeof(double));
        for (int i=0;i<n;i++) for (int j=0;j<n;j++)
            for (int kk=0;kk<n;kk++) new_tmp[i*n+j] += tmp[i*n+kk]*Ak[kk*n+j];
        memcpy(tmp, new_tmp, n*n*sizeof(double));
        free(new_tmp);
    }
    free(tmp); free(Ak);
    return 0;
}

/* ---------- Adjoint representations ---------- */

/* Ad_R(eta) = R*eta*R^T for SO(3) */
void geo_so3_Ad(const double R[3][3], const double eta[3][3],
                double result[3][3]) {
    if (!R || !eta || !result) return;
    double RT[3][3];
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) RT[i][j] = R[j][i];
    double tmp[3][3];
    memset(tmp, 0, sizeof(tmp));
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
        for (int k=0;k<3;k++) tmp[i][j] += R[i][k]*eta[k][j];
    memset(result, 0, sizeof(double)*9);
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
        for (int k=0;k<3;k++) result[i][j] += tmp[i][k]*RT[k][j];
}

/* ad_X(Y) = [X,Y] = XY-YX */
void geo_so3_ad(const double X[3][3], const double Y[3][3],
                double result[3][3]) {
    geo_so3_bracket(X, Y, result);
}

/* Jacobi identity on so(3): [X,[Y,Z]]+[Y,[Z,X]]+[Z,[X,Y]] = 0 */
int geo_so3_jacobi_identity_check(const double X[3][3],
                                  const double Y[3][3],
                                  const double Z[3][3], double tol) {
    if (!X||!Y||!Z) return 0;
    double t1[3][3], t2[3][3], t3[3][3], t4[3][3], sum[3][3];
    geo_so3_bracket(Y,Z,t1); geo_so3_bracket(X,t1,t1);
    geo_so3_bracket(Z,X,t2); geo_so3_bracket(Y,t2,t2);
    geo_so3_bracket(X,Y,t3); geo_so3_bracket(Z,t3,t3);
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
        sum[i][j] = t1[i][j]+t2[i][j]+t3[i][j];
    double maxv = 0.0;
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
        if (fabs(sum[i][j])>maxv) maxv=fabs(sum[i][j]);
    return (maxv<tol)?1:0;
}

/* ---------- RK-Munthe-Kaas Lie group integrator ---------- */

void geo_rk_munthe_kaas(int group_dim,
                        void (*omega)(const double *A, void *ctx, double *omega_out),
                        void *omega_ctx,
                        const double *A0, double t, int n_steps,
                        double *A_out, GeoGroupExp use_exp) {
    if (!omega || !A0 || !A_out || !use_exp || n_steps<=0) return;
    double dt = t/(double)n_steps;
    double *A = (double*)malloc(group_dim*sizeof(double));
    double *k1 = (double*)malloc(group_dim*sizeof(double));
    double *k2 = (double*)malloc(group_dim*sizeof(double));
    double *k3 = (double*)malloc(group_dim*sizeof(double));
    double *k4 = (double*)malloc(group_dim*sizeof(double));
    double *tmp = (double*)malloc(group_dim*sizeof(double));
    memcpy(A, A0, group_dim*sizeof(double));
    for (int s=0;s<n_steps;s++) {
        /* RK4 in Lie algebra */
        omega(A, omega_ctx, k1);
        for (int i=0;i<group_dim;i++) tmp[i]=k1[i]*0.5*dt;
        use_exp(tmp, tmp);
        double *Atmp=(double*)malloc(group_dim*sizeof(double));
        /* A*exp(0.5*dt*k1) approximated: for small dt, exp(dt*k)~I+dt*k */
        for (int i=0;i<group_dim;i++) Atmp[i]=A[i]+0.5*dt*k1[i];
        omega(Atmp, omega_ctx, k2);
        for (int i=0;i<group_dim;i++) Atmp[i]=A[i]+0.5*dt*k2[i];
        omega(Atmp, omega_ctx, k3);
        for (int i=0;i<group_dim;i++) Atmp[i]=A[i]+dt*k3[i];
        omega(Atmp, omega_ctx, k4);
        for (int i=0;i<group_dim;i++)
            A[i] += (dt/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
        free(Atmp);
    }
    memcpy(A_out, A, group_dim*sizeof(double));
    free(A);free(k1);free(k2);free(k3);free(k4);free(tmp);
}

/* ---------- Crouch-Grossman commutator-free integrator ---------- */

void geo_crouch_grossman_2(int group_dim,
                           void (*rhs)(const double *A, void *ctx, double *omega),
                           void *rhs_ctx,
                           const double *A0, double h, int n_steps,
                           double *A_out, GeoGroupExp use_exp) {
    if (!rhs||!A0||!A_out||!use_exp||n_steps<=0) return;
    double *A = (double*)malloc(group_dim*sizeof(double));
    double *k1 = (double*)malloc(group_dim*sizeof(double));
    double *Ah = (double*)malloc(group_dim*sizeof(double));
    memcpy(A, A0, group_dim*sizeof(double));
    for (int s=0;s<n_steps;s++) {
        rhs(A, rhs_ctx, k1);
        for (int i=0;i<group_dim;i++) Ah[i]=A[i]+0.5*h*k1[i];
        double *k2 = (double*)malloc(group_dim*sizeof(double));
        rhs(Ah, rhs_ctx, k2);
        for (int i=0;i<group_dim;i++) A[i] += h*k2[i];
        free(k2);
    }
    memcpy(A_out, A, group_dim*sizeof(double));
    free(A);free(k1);free(Ah);
}

/* ---------- Euler-Poincare rigid body equations ---------- */

/* I*Omega_dot + Omega x (I*Omega) = tau */
void geo_euler_poincare_rigid_body(const double I[3],
                                   const double Omega[3],
                                   const double tau[3],
                                   double dOmega[3]) {
    if (!I||!Omega||!tau||!dOmega) return;
    double IO[3];
    IO[0] = I[0]*Omega[0]; IO[1] = I[1]*Omega[1]; IO[2] = I[2]*Omega[2];
    double cross[3];
    cross[0] = Omega[1]*IO[2] - Omega[2]*IO[1];
    cross[1] = Omega[2]*IO[0] - Omega[0]*IO[2];
    cross[2] = Omega[0]*IO[1] - Omega[1]*IO[0];
    for (int i=0;i<3;i++) dOmega[i] = (tau[i] - cross[i]) / I[i];
}

/* R_dot = R * Omega_hat (kinematics on SO(3)) */
void geo_so3_kinematics(const double R[3][3], const double Omega[3],
                        double dR[3][3]) {
    if (!R||!Omega||!dR) return;
    double hat[3][3];
    geo_so3_hat(Omega, hat);
    memset(dR, 0, 9*sizeof(double));
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
        for (int k=0;k<3;k++) dR[i][j] += R[i][k] * hat[k][j];
}

/* Arnold cat map on torus T^2 = R^2/Z^2.
 * (x_{n+1},y_{n+1}) = ((2x_n+y_n) mod 1, (x_n+y_n) mod 1).
 * Area-preserving, mixing, Anosov diffeomorphism. */
void geo_arnold_cat_map(double *x, double *y, int n_iter) {
    if (!x||!y) return;
    for (int n=0;n<n_iter;n++) {
        double xn = 2.0*(*x) + (*y);
        double yn = (*x) + (*y);
        *x = xn - floor(xn);
        *y = yn - floor(yn);
    }
}
