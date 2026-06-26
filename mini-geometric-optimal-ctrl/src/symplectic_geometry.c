
/**
 * symplectic_geometry.c - Symplectic Geometry for Optimal Control
 * Reference: Marsden & Ratiu (1999), Arnold (1989), da Silva (2008)
 */
#include "symplectic_geometry.h"
#include <string.h>
#include <stdlib.h>

void geo_symplectic_init(GeoSymplecticSpace *s, int n) {
    if (!s || n <= 0 || 2*n > GEO_MAX_DIM) return;
    s->n = n;  s->dim_2n = 2 * n;
    memset(s->matrix, 0, sizeof(s->matrix));
    for (int i = 0; i < n; i++) {
        s->matrix[i][n + i] = 1.0;
        s->matrix[n + i][i] = -1.0;
    }
}

double geo_symplectic_form(const GeoSymplecticSpace *s,
                           const double *u, const double *v) {
    if (!s || !u || !v) return 0.0;
    int N = s->dim_2n;
    double result = 0.0;
    for (int i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < N; j++)
            row_sum += s->matrix[i][j] * v[j];
        result += u[i] * row_sum;
    }
    return result;
}

int geo_symplectic_nondegenerate(const GeoSymplecticSpace *s) {
    if (!s || s->dim_2n <= 0) return 0;
    int n = s->n;
    for (int i = 0; i < n; i++) {
        if (s->matrix[i][n+i] != 1.0) return 0;
        if (s->matrix[n+i][i] != -1.0) return 0;
    }
    return 1;
}

void geo_symplectic_flat(const GeoSymplecticSpace *s,
                         const double *covector, double *vector) {
    if (!s || !covector || !vector) return;
    int n = s->n;
    for (int i = 0; i < n; i++) {
        vector[i]     = -covector[n + i];
        vector[n + i] =  covector[i];
    }
}

void geo_symplectic_sharp(const GeoSymplecticSpace *s,
                          const double *vector, double *covector) {
    if (!s || !vector || !covector) return;
    int n = s->n;
    for (int i = 0; i < n; i++) {
        covector[i]     =  vector[n + i];
        covector[n + i] = -vector[i];
    }
}

void geo_hamiltonian_vector_field(GeoHamiltonian H,
                                  const double *state,
                                  double *result, int n,
                                  double h, void *ctx) {
    if (!H || !state || !result || n <= 0) return;
    if (h <= 0.0) h = 1e-6;
    int N = 2 * n;
    double *sp = (double *)malloc(N * sizeof(double));
    double *sm = (double *)malloc(N * sizeof(double));
    memcpy(sp, state, N * sizeof(double));
    memcpy(sm, state, N * sizeof(double));
    for (int i = 0; i < N; i++) {
        sp[i] = state[i] + h;
        double Hp = H(sp, n, ctx);
        sm[i] = state[i] - h;
        double Hm = H(sm, n, ctx);
        double dH = (Hp - Hm) / (2.0 * h);
        if (i < n) result[i] = dH;
        else       result[i] = -dH;
        sp[i] = state[i];  sm[i] = state[i];
    }
    free(sp); free(sm);
}

void geo_stormer_verlet(GeoHamiltonian H, void *ctx,
                        const double *q0, const double *p0,
                        double t, int n_steps,
                        double *q_out, double *p_out, int n) {
    if (!H || !q0 || !p0 || !q_out || !p_out || n <= 0) return;
    if (n_steps <= 0) n_steps = 100;
    double hs = t / (double)n_steps;
    double *q = (double *)malloc(n * sizeof(double));
    double *p = (double *)malloc(n * sizeof(double));
    double *state = (double *)malloc(2*n * sizeof(double));
    double *dV = (double *)malloc(n * sizeof(double));
    double *dT = (double *)malloc(n * sizeof(double));
    memcpy(q, q0, n * sizeof(double));
    memcpy(p, p0, n * sizeof(double));
    for (int s = 0; s < n_steps; s++) {
        double *qp = (double *)malloc(n * sizeof(double));
        memcpy(qp, q, n * sizeof(double));
        for (int i = 0; i < n; i++) {
            qp[i] = q[i] + 1e-6;
            for (int k = 0; k < n; k++) { state[k]=qp[k]; state[n+k]=0.0; }
            double Hp = H(state, n, ctx);
            qp[i] = q[i] - 1e-6;
            for (int k = 0; k < n; k++) { state[k]=qp[k]; state[n+k]=0.0; }
            double Hm = H(state, n, ctx);
            dV[i] = (Hp - Hm) / (2.0 * 1e-6);
            qp[i] = q[i];
        }
        free(qp);
        for (int i = 0; i < n; i++) p[i] -= 0.5 * hs * dV[i];
        double *pp = (double *)malloc(n * sizeof(double));
        memcpy(pp, p, n * sizeof(double));
        for (int i = 0; i < n; i++) {
            pp[i] = p[i] + 1e-6;
            for (int k = 0; k < n; k++) { state[k]=q[k]; state[n+k]=pp[k]; }
            double Hp = H(state, n, ctx);
            pp[i] = p[i] - 1e-6;
            for (int k = 0; k < n; k++) { state[k]=q[k]; state[n+k]=pp[k]; }
            double Hm = H(state, n, ctx);
            dT[i] = (Hp - Hm) / (2.0 * 1e-6);
            pp[i] = p[i];
        }
        free(pp);
        for (int i = 0; i < n; i++) q[i] += hs * dT[i];
        double *qp2 = (double *)malloc(n * sizeof(double));
        memcpy(qp2, q, n * sizeof(double));
        for (int i = 0; i < n; i++) {
            qp2[i] = q[i] + 1e-6;
            for (int k = 0; k < n; k++) { state[k]=qp2[k]; state[n+k]=0.0; }
            double Hp = H(state, n, ctx);
            qp2[i] = q[i] - 1e-6;
            for (int k = 0; k < n; k++) { state[k]=qp2[k]; state[n+k]=0.0; }
            double Hm = H(state, n, ctx);
            dV[i] = (Hp - Hm) / (2.0 * 1e-6);
            qp2[i] = q[i];
        }
        free(qp2);
        for (int i = 0; i < n; i++) p[i] -= 0.5 * hs * dV[i];
    }
    memcpy(q_out, q, n * sizeof(double));
    memcpy(p_out, p, n * sizeof(double));
    free(q); free(p); free(state); free(dV); free(dT);
}

double geo_poisson_bracket(GeoScalarFunction F, GeoScalarFunction G,
                           const double *state, int n,
                           double h, void *ctx) {
    if (!F || !G || !state || n <= 0) return 0.0;
    if (h <= 0.0) h = 1e-6;
    int N = 2 * n;
    double *gradF = (double *)malloc(N * sizeof(double));
    double *gradG = (double *)malloc(N * sizeof(double));
    double *sp = (double *)malloc(N * sizeof(double));
    double *sm = (double *)malloc(N * sizeof(double));
    memcpy(sp, state, N * sizeof(double));
    memcpy(sm, state, N * sizeof(double));
    for (int k = 0; k < N; k++) {
        sp[k] = state[k] + h;
        sm[k] = state[k] - h;
        gradF[k] = (F(sp, N, ctx) - F(sm, N, ctx)) / (2.0 * h);
        gradG[k] = (G(sp, N, ctx) - G(sm, N, ctx)) / (2.0 * h);
        sp[k] = state[k];  sm[k] = state[k];
    }
    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += gradF[i]*gradG[n+i] - gradF[n+i]*gradG[i];
    free(gradF); free(gradG); free(sp); free(sm);
    return result;
}

int geo_jacobi_identity_check(GeoScalarFunction F,
                              GeoScalarFunction G,
                              GeoScalarFunction H,
                              const double *state, int n,
                              double h, double tolerance, void *ctx) {
    if (!F || !G || !H || !state || n <= 0) return 0;
    if (h <= 0.0) h = 1e-6;
    int N = 2 * n;
    double *sp = (double *)malloc(N * sizeof(double));
    memcpy(sp, state, N * sizeof(double));
    double pb = 0.0;
    for (int i = 0; i < n; i++) {
        sp[i]=state[i]+h;
        double gh_qp = geo_poisson_bracket(G,H,sp,n,h,ctx);
        sp[i]=state[i]-h;
        double gh_qm = geo_poisson_bracket(G,H,sp,n,h,ctx);
        sp[i]=state[i];
        double dgh_dqi = (gh_qp-gh_qm)/(2.0*h);
        sp[n+i]=state[n+i]+h;
        double gh_pp = geo_poisson_bracket(G,H,sp,n,h,ctx);
        sp[n+i]=state[n+i]-h;
        double gh_pm = geo_poisson_bracket(G,H,sp,n,h,ctx);
        sp[n+i]=state[n+i];
        double dgh_dpi = (gh_pp-gh_pm)/(2.0*h);
        sp[i]=state[i]+h; double F_qp=F(sp,N,ctx); sp[i]=state[i]-h;
        double F_qm=F(sp,N,ctx); sp[i]=state[i];
        double dF_dqi=(F_qp-F_qm)/(2.0*h);
        sp[n+i]=state[n+i]+h; double F_pp=F(sp,N,ctx); sp[n+i]=state[n+i]-h;
        double F_pm=F(sp,N,ctx); sp[n+i]=state[n+i];
        double dF_dpi=(F_pp-F_pm)/(2.0*h);
        pb += dF_dqi*dgh_dpi - dF_dpi*dgh_dqi;
    }
    free(sp);
    return (fabs(pb) < tolerance) ? 1 : 0;
}

void geo_liouville_one_form(const double *q_coords,
                            const double *p_coords,
                            double *result, int n) {
    if (!q_coords || !p_coords || !result || n <= 0) return;
    memset(result, 0, 2 * n * sizeof(double));
    for (int i = 0; i < n; i++) result[n + i] = p_coords[i];
}

void geo_canonical_symplectic_form(int n,
                                   double matrix[GEO_MAX_DIM][GEO_MAX_DIM]) {
    if (n <= 0 || 2*n > GEO_MAX_DIM) return;
    memset(matrix, 0, GEO_MAX_DIM*GEO_MAX_DIM*sizeof(double));
    for (int i = 0; i < n; i++) {
        matrix[i][n+i] = 1.0;
        matrix[n+i][i] = -1.0;
    }
}

void geo_legendre_transform(const double *q, const double *v,
                            const double *mass,
                            GeoScalarFunction potential, void *ctx,
                            double *p_out, double *H_out, int n) {
    if (!q || !v || !p_out || !H_out || n <= 0) return;
    double m = mass ? mass[0] : 1.0;
    if (m <= 0.0) m = 1.0;
    double ke = 0.0;
    for (int i = 0; i < n; i++) {
        p_out[i] = m * v[i];
        ke += 0.5 * m * v[i] * v[i];
    }
    double pe = potential ? potential(q, n, ctx) : 0.0;
    *H_out = ke + pe;
}

void geo_symplectic_gradient(GeoScalarFunction F, void *ctx,
                             const double *state, double *result,
                             int n, double h) {
    if (!F || !state || !result || n <= 0) return;
    if (h <= 0.0) h = 1e-6;
    int N = 2 * n;
    double *sp=(double*)malloc(N*sizeof(double));
    double *sm=(double*)malloc(N*sizeof(double));
    memcpy(sp,state,N*sizeof(double)); memcpy(sm,state,N*sizeof(double));
    for (int i=0;i<N;i++) {
        sp[i]=state[i]+h; sm[i]=state[i]-h;
        double dF=(F(sp,N,ctx)-F(sm,N,ctx))/(2.0*h);
        if(i<n) result[n+i]=-dF; else result[i-n]=dF;
        sp[i]=state[i]; sm[i]=state[i];
    }
    free(sp); free(sm);
}
