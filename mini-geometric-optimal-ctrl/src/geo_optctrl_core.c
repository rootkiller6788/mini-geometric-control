/**
 * geo_optctrl_core.c - Core Implementation of Geometric Optimal Control
 *
 * Implements: Lie bracket [X,Y], iterated LB ad^k_X(Y), flow via RK4,
 * flow Jacobian DPhi^t, differential df, product manifold operations.
 * Reference: Agrachev & Sachkov (2004), Jurdjevic (1997), Lee (2012)
 */

#include "geo_optctrl_core.h"
#include <string.h>
#include <stdlib.h>

static void jacobian_vector_field(const GeoVectorField *X,
                                   const double *x, double *J,
                                   int n, double h) {
    double *xp = (double *)malloc(n * sizeof(double));
    double *xm = (double *)malloc(n * sizeof(double));
    double *fp = (double *)malloc(n * sizeof(double));
    double *fm = (double *)malloc(n * sizeof(double));
    memcpy(xp, x, n * sizeof(double));
    memcpy(xm, x, n * sizeof(double));
    for (int j = 0; j < n; j++) {
        xp[j] = x[j] + h;  xm[j] = x[j] - h;
        X->eval(xp, fp, X->ctx);
        X->eval(xm, fm, X->ctx);
        for (int i = 0; i < n; i++)
            J[i * n + j] = (fp[i] - fm[i]) / (2.0 * h);
        xp[j] = x[j];  xm[j] = x[j];
    }
    free(xp); free(xm); free(fp); free(fm);
}

void geo_lie_bracket(const GeoVectorField *X, const GeoVectorField *Y,
                     const double *at_point, double *result,
                     int dim, double h) {
    if (!X || !Y || !at_point || !result || dim <= 0 || dim > GEO_MAX_DIM)
        return;
    if (h <= 0.0) h = 1e-6;
    double *JX = (double *)calloc(dim * dim, sizeof(double));
    double *JY = (double *)calloc(dim * dim, sizeof(double));
    double *Xv = (double *)malloc(dim * sizeof(double));
    double *Yv = (double *)malloc(dim * sizeof(double));
    jacobian_vector_field(X, at_point, JX, dim, h);
    jacobian_vector_field(Y, at_point, JY, dim, h);
    X->eval(at_point, Xv, X->ctx);
    Y->eval(at_point, Yv, Y->ctx);
    for (int i = 0; i < dim; i++) {
        result[i] = 0.0;
        for (int j = 0; j < dim; j++)
            result[i] += Xv[j]*JY[i*dim+j] - Yv[j]*JX[i*dim+j];
    }
    free(JX); free(JY); free(Xv); free(Yv);
}

void geo_iterated_lie_bracket(const GeoVectorField *X,
                              const GeoVectorField *Y,
                              int k, const double *at_point,
                              double *result, int dim, double h) {
    if (k < 0 || dim <= 0 || dim > GEO_MAX_DIM) return;
    if (h <= 0.0) h = 1e-6;
    if (k == 0) { Y->eval(at_point, result, Y->ctx); return; }
    double *prev = (double *)malloc(dim * sizeof(double));
    geo_iterated_lie_bracket(X, Y, k-1, at_point, prev, dim, h);
    double *Jad = (double *)calloc(dim*dim, sizeof(double));
    double *JX  = (double *)calloc(dim*dim, sizeof(double));
    double *Xv  = (double *)malloc(dim * sizeof(double));
    double *xp  = (double *)malloc(dim * sizeof(double));
    jacobian_vector_field(X, at_point, JX, dim, h);
    X->eval(at_point, Xv, X->ctx);
    memcpy(xp, at_point, dim * sizeof(double));
    for (int j = 0; j < dim; j++) {
        xp[j] = at_point[j] + h;
        double *t1 = (double *)malloc(dim * sizeof(double));
        geo_iterated_lie_bracket(X, Y, k-1, xp, t1, dim, h);
        xp[j] = at_point[j] - h;
        double *t2 = (double *)malloc(dim * sizeof(double));
        geo_iterated_lie_bracket(X, Y, k-1, xp, t2, dim, h);
        for (int i = 0; i < dim; i++)
            Jad[i*dim+j] = (t1[i] - t2[i]) / (2.0 * h);
        xp[j] = at_point[j];
        free(t1); free(t2);
    }
    for (int i = 0; i < dim; i++) {
        result[i] = 0.0;
        for (int j = 0; j < dim; j++)
            result[i] += Xv[j]*Jad[i*dim+j] - prev[j]*JX[i*dim+j];
    }
    free(prev); free(Jad); free(JX); free(Xv); free(xp);
}

void geo_flow(const GeoVectorField *X, const double *x0, double t,
              int n_steps, double *result, int dim) {
    if (!X || !x0 || !result || dim <= 0 || dim > GEO_MAX_DIM) return;
    if (n_steps <= 0) n_steps = 100;
    double dt = t / (double)n_steps;
    double *x  = (double *)malloc(dim * sizeof(double));
    double *k1 = (double *)malloc(dim * sizeof(double));
    double *k2 = (double *)malloc(dim * sizeof(double));
    double *k3 = (double *)malloc(dim * sizeof(double));
    double *k4 = (double *)malloc(dim * sizeof(double));
    double *tp = (double *)malloc(dim * sizeof(double));
    memcpy(x, x0, dim * sizeof(double));
    for (int s = 0; s < n_steps; s++) {
        X->eval(x, k1, X->ctx);
        for (int i = 0; i < dim; i++) tp[i] = x[i] + 0.5*dt*k1[i];
        X->eval(tp, k2, X->ctx);
        for (int i = 0; i < dim; i++) tp[i] = x[i] + 0.5*dt*k2[i];
        X->eval(tp, k3, X->ctx);
        for (int i = 0; i < dim; i++) tp[i] = x[i] + dt*k3[i];
        X->eval(tp, k4, X->ctx);
        for (int i = 0; i < dim; i++)
            x[i] += (dt/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
    memcpy(result, x, dim * sizeof(double));
    free(x); free(k1); free(k2); free(k3); free(k4); free(tp);
}

void geo_flow_jacobian(const GeoVectorField *X, const double *x0,
                       double t, int n_steps, double *jac,
                       int dim, double h) {
    if (!X || !x0 || !jac || dim <= 0 || dim > GEO_MAX_DIM) return;
    if (n_steps <= 0) n_steps = 100;
    if (h <= 0.0) h = 1e-6;
    double dt = t / (double)n_steps;
    double *x    = (double *)malloc(dim * sizeof(double));
    double *Phi  = (double *)calloc(dim*dim, sizeof(double));
    double *JX   = (double *)malloc(dim*dim * sizeof(double));
    double *dPhi = (double *)calloc(dim*dim, sizeof(double));
    memcpy(x, x0, dim * sizeof(double));
    for (int i = 0; i < dim; i++) Phi[i*dim + i] = 1.0;
    for (int s = 0; s < n_steps; s++) {
        jacobian_vector_field(X, x, JX, dim, h);
        memset(dPhi, 0, dim*dim*sizeof(double));
        for (int i = 0; i < dim; i++)
            for (int k = 0; k < dim; k++)
                for (int j = 0; j < dim; j++)
                    dPhi[i*dim+j] += JX[i*dim+k] * Phi[k*dim+j];
        for (int ij = 0; ij < dim*dim; ij++) Phi[ij] += dt * dPhi[ij];
        double *kx = (double *)malloc(dim * sizeof(double));
        X->eval(x, kx, X->ctx);
        for (int i = 0; i < dim; i++) x[i] += dt * kx[i];
        free(kx);
    }
    memcpy(jac, Phi, dim*dim*sizeof(double));
    free(x); free(Phi); free(JX); free(dPhi);
}

void geo_differential(GeoScalarFunction f, void *ctx,
                      const double *at_p, double *result,
                      int dim, double h) {
    if (!f || !at_p || !result || dim <= 0 || dim > GEO_MAX_DIM) return;
    if (h <= 0.0) h = 1e-6;
    double *xp = (double *)malloc(dim * sizeof(double));
    memcpy(xp, at_p, dim * sizeof(double));
    for (int i = 0; i < dim; i++) {
        xp[i] = at_p[i] + h;
        double fp = f(xp, dim, ctx);
        xp[i] = at_p[i] - h;
        double fm = f(xp, dim, ctx);
        result[i] = (fp - fm) / (2.0 * h);
        xp[i] = at_p[i];
    }
    free(xp);
}

void geo_product_embed(const double *x1, int dim1,
                       const double *x2, int dim2,
                       double *result) {
    if (!x1 || !x2 || !result) return;
    memcpy(result, x1, dim1 * sizeof(double));
    memcpy(result + dim1, x2, dim2 * sizeof(double));
}

void geo_product_project(const double *x, int dim_total,
                         double *x1, int dim1,
                         double *x2, int dim2) {
    if (!x || !x1 || !x2) return;
    if (dim1 + dim2 != dim_total) return;
    memcpy(x1, x, dim1 * sizeof(double));
    memcpy(x2, x + dim1, dim2 * sizeof(double));
}