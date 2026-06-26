#include "lie_derivative.h"
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * Lie Derivative Implementation
 *
 * L_f h(x) = sum_{i=1}^n (partial h / partial x_i) * f_i(x)
 *          = (grad h)^T . f(x)
 *
 * For higher-order: L_f^k h = L_f (L_f^{k-1} h)
 *
 * The mixed derivative L_g L_f^k h is essential for determining
 * relative degree in feedback linearization.
 * ============================================================================ */

/* Finite difference step size */
#define DEFAULT_H 1e-6
#define DEFAULT_TOL 1e-8

/* Forward finite difference for scalar gradient */
static void finite_gradient(const ScalarField *h, const Vector *x,
                             Vector *grad, double eps) {
    int n = x->dim;
    double hx = scalar_field_eval(h, x);

    Vector *xp = vector_clone(x);
    for (int i = 0; i < n; i++) {
        double orig = xp->data[i];
        xp->data[i] = orig + eps;
        double hp = scalar_field_eval(h, xp);
        grad->data[i] = (hp - hx) / eps;
        xp->data[i] = orig;
    }
    vector_free(xp);
}

/* Central finite difference for scalar gradient (more accurate) */
static void finite_gradient_central(const ScalarField *h, const Vector *x,
                                     Vector *grad, double eps) {
    int n = x->dim;
    Vector *xp = vector_clone(x);
    Vector *xm = vector_clone(x);
    for (int i = 0; i < n; i++) {
        double orig = x->data[i];
        xp->data[i] = orig + eps;
        xm->data[i] = orig - eps;
        double hp = scalar_field_eval(h, xp);
        double hm = scalar_field_eval(h, xm);
        grad->data[i] = (hp - hm) / (2.0 * eps);
        xp->data[i] = orig;
        xm->data[i] = orig;
    }
    vector_free(xp);
    vector_free(xm);
}

/* ============================================================================
 * Public API
 * ============================================================================ */

double lie_derivative(const VectorField *f,
                       const ScalarField *h,
                       const Vector *x) {
    if (!f || !h || !x) return 0.0;
    if (f->dim != h->dim || f->dim != x->dim) return 0.0;

    Vector *grad = vector_create(x->dim);
    Vector *fx = vector_create(x->dim);

    /* Compute gradient of h */
    if (h->gradient) {
        h->gradient(x, grad);
    } else {
        finite_gradient_central(h, x, grad, DEFAULT_H);
    }

    /* Evaluate f at x */
    vector_field_eval(f, x, fx);

    /* L_f h = grad_h . f */
    double result = vector_dot(grad, fx);

    vector_free(grad);
    vector_free(fx);
    return result;
}

double lie_derivative_numeric(const VectorField *f,
                               const ScalarField *h,
                               const Vector *x,
                               double h_step) {
    if (!f || !h || !x) return 0.0;

    /* Directional derivative via central difference along f */
    Vector *fx = vector_create(x->dim);
    Vector *xp = vector_clone(x);
    Vector *xm = vector_clone(x);

    vector_field_eval(f, x, fx);

    double eps = (h_step > 0) ? h_step : DEFAULT_H;
    double nrm_f = vector_norm(fx);
    if (nrm_f < 1e-15) {
        vector_free(fx); vector_free(xp); vector_free(xm);
        return 0.0;
    }

    /* x ? eps * f_normalized */
    for (int i = 0; i < x->dim; i++) {
        double step = eps * fx->data[i] / nrm_f;
        xp->data[i] = x->data[i] + step;
        xm->data[i] = x->data[i] - step;
    }

    double hp = scalar_field_eval(h, xp);
    double hm = scalar_field_eval(h, xm);

    vector_free(fx); vector_free(xp); vector_free(xm);

    return (hp - hm) * nrm_f / (2.0 * eps);
}

double lie_derivative_kth(const VectorField *f,
                           const ScalarField *h,
                           const Vector *x,
                           int k) {
    if (k < 0) return 0.0;
    if (k == 0) return scalar_field_eval(h, x);
    if (k == 1) return lie_derivative(f, h, x);

    /* For k >= 2, we need gradient of L_f^{k-1} h */
    /* Compute via finite differences: L_f^k h(x) = limit (L_f^{k-1}h(x+eps*f) - L_f^{k-1}h(x)) / eps */

    Vector *fx = vector_create(x->dim);
    vector_field_eval(f, x, fx);

    double nrm_f = vector_norm(fx);
    if (nrm_f < 1e-15) {
        vector_free(fx);
        return 0.0;
    }

    Vector *xp = vector_clone(x);
    double eps = DEFAULT_H;
    for (int i = 0; i < x->dim; i++) {
        xp->data[i] = x->data[i] + eps * fx->data[i] / nrm_f;
    }

    double hkp = lie_derivative_kth(f, h, xp, k - 1);
    double hk = lie_derivative_kth(f, h, x, k - 1);

    vector_free(fx);
    vector_free(xp);

    return (hkp - hk) * nrm_f / eps;
}

double lie_derivative_mixed(const VectorField *g,
                             const VectorField *f,
                             const ScalarField *h,
                             const Vector *x,
                             int k) {
    if (!g || !f || !h || !x) return 0.0;

    if (k < 0) {
        /* L_g h(x) directly */
        return lie_derivative(g, h, x);
    }

    /* L_g L_f^k h(x) = L_g (L_f^k h)(x) */
    /* We treat L_f^k h as a new scalar field phi = L_f^k h, then compute L_g phi */

    /* Compute L_f^k h(x) and its gradient via finite differences */
    /* Strategy: define a helper scalar field wrapping L_f^k h */

    /* Build a local scalar field phi(x) = L_f^k h(x) */
    /* Since phi doesn't have gradient callback, lie_derivative will use finite diff */

    double Lfk_h = lie_derivative_kth(f, h, x, k);

    /* To compute L_g (L_f^k h), we need the gradient of L_f^k h */
    /* Compute via central finite differences of lie_derivative_kth */

    Vector *grad = vector_create(x->dim);
    Vector *xp = vector_clone(x);
    Vector *xm = vector_clone(x);
    double eps = DEFAULT_H;

    for (int i = 0; i < x->dim; i++) {
        double orig = x->data[i];
        xp->data[i] = orig + eps;
        xm->data[i] = orig - eps;
        double hp = lie_derivative_kth(f, h, xp, k);
        double hm = lie_derivative_kth(f, h, xm, k);
        grad->data[i] = (hp - hm) / (2.0 * eps);
        xp->data[i] = orig;
        xm->data[i] = orig;
    }

    /* Now L_g phi = grad(phi) . g(x) */
    Vector *gx = vector_create(x->dim);
    vector_field_eval(g, x, gx);
    double result = vector_dot(grad, gx);

    vector_free(grad); vector_free(xp); vector_free(xm); vector_free(gx);
    return result;
}

void lie_derivative_vector(const VectorField *X,
                            const VectorField *Y,
                            const Vector *x,
                            Vector *out) {
    if (!X || !Y || !x || !out) return;
    if (X->dim != Y->dim || X->dim != x->dim || X->dim != out->dim) return;

    /* L_X Y(x) = limit (Y(x + t X) - Y(x)) / t */
    double eps = DEFAULT_H;

    Vector *Xx = vector_create(x->dim);
    vector_field_eval(X, x, Xx);
    double nrm_X = vector_norm(Xx);

    Vector *xp = vector_clone(x);
    if (nrm_X > 1e-12) {
        for (int i = 0; i < x->dim; i++) {
            xp->data[i] = x->data[i] + eps * Xx->data[i] / nrm_X;
        }
    }

    Vector *Yxp = vector_create(x->dim);
    Vector *Yx = vector_create(x->dim);
    vector_field_eval(Y, xp, Yxp);
    vector_field_eval(Y, x, Yx);

    for (int i = 0; i < x->dim; i++) {
        out->data[i] = (Yxp->data[i] - Yx->data[i]) * nrm_X / eps;
    }

    vector_free(Xx); vector_free(xp); vector_free(Yxp); vector_free(Yx);
}

double lie_derivative_scalar_along_flow(const VectorField *f,
                                         const ScalarField *h,
                                         const Vector *x,
                                         double dt) {
    if (!f || !h || !x) return 0.0;

    /* Evolve x forward a small time step using Euler: x(t+dt) = x + dt*f(x) */
    Vector *fx = vector_create(x->dim);
    Vector *xt = vector_clone(x);
    vector_field_eval(f, x, fx);

    for (int i = 0; i < x->dim; i++) {
        xt->data[i] = x->data[i] + dt * fx->data[i];
    }

    double h0 = scalar_field_eval(h, x);
    double ht = scalar_field_eval(h, xt);

    vector_free(fx);
    vector_free(xt);

    return (ht - h0) / dt;  /* forward difference */
}

bool check_relative_degree_order(const VectorField *g,
                                  const VectorField *f,
                                  const ScalarField *h,
                                  const Vector *x,
                                  int k,
                                  double tol) {
    double val = lie_derivative_mixed(g, f, h, x, k);
    return fabs(val) < tol;
}

double *compute_all_lie_derivatives(const VectorField *f,
                                     const ScalarField *h,
                                     const Vector *x,
                                     int max_k) {
    if (!f || !h || !x || max_k < 0) return NULL;

    double *results = (double*)malloc((size_t)(max_k + 1) * sizeof(double));
    if (!results) return NULL;

    for (int k = 0; k <= max_k; k++) {
        results[k] = lie_derivative_kth(f, h, x, k);
    }
    return results;
}
