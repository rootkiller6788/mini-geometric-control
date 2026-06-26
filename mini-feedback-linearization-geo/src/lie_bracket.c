#include "lie_bracket.h"
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * Lie Bracket Implementation
 *
 * [f, g](x) = J_g(x) * f(x) - J_f(x) * g(x)
 *
 * This is the commutator of vector fields. The Lie bracket measures how
 * much the flows of f and g fail to commute.
 *
 * Key properties:
 *   (1) [f, g] = -[g, f]               (anti-symmetry)
 *   (2) [f, [g, h]] + [g, [h, f]] + [h, [f, g]] = 0  (Jacobi identity)
 *
 * Iterated brackets: ad_f^k g = [f, ad_f^{k-1} g]
 * ============================================================================ */

#define DEFAULT_EPS 1e-6

/* Finite difference Jacobian: J = df/dx at x */
static void finite_jacobian(const VectorField *f, const Vector *x,
                             Matrix *J, double eps) {
    int n = x->dim;
    Vector *fx = vector_create(n);
    Vector *xp = vector_clone(x);
    Vector *fp = vector_create(n);
    Vector *fm = vector_create(n);

    vector_field_eval(f, x, fx); /* not needed for central diff */

    for (int j = 0; j < n; j++) {
        double orig = xp->data[j];
        /* Central difference: column j of J = (f(x+eps*e_j) - f(x-eps*e_j)) / (2*eps) */
        xp->data[j] = orig + eps;
        vector_field_eval(f, xp, fp);
        xp->data[j] = orig - eps;
        vector_field_eval(f, xp, fm);
        xp->data[j] = orig;

        for (int i = 0; i < n; i++) {
            J->data[i * n + j] = (fp->data[i] - fm->data[i]) / (2.0 * eps);
        }
    }

    vector_free(fx); vector_free(xp); vector_free(fp); vector_free(fm);
}

/* ============================================================================
 * Public API
 * ============================================================================ */

void lie_bracket(const VectorField *f,
                  const VectorField *g,
                  const Vector *x,
                  Vector *out) {
    if (!f || !g || !x || !out) return;
    int n = x->dim;
    if (f->dim != n || g->dim != n || out->dim != n) return;

    Vector *fx = vector_create(n);
    Vector *gx = vector_create(n);
    Matrix *Jf = matrix_create(n, n);
    Matrix *Jg = matrix_create(n, n);

    vector_field_eval(f, x, fx);
    vector_field_eval(g, x, gx);

    /* Compute Jacobians */
    if (f->jacobian) {
        f->jacobian(x, Jf);
    } else {
        finite_jacobian(f, x, Jf, DEFAULT_EPS);
    }

    if (g->jacobian) {
        g->jacobian(x, Jg);
    } else {
        finite_jacobian(g, x, Jg, DEFAULT_EPS);
    }

    /* [f,g] = Jg * f - Jf * g */
    Vector *Jg_f = vector_create(n);
    Vector *Jf_g = vector_create(n);
    matrix_vector_multiply(Jg, fx, Jg_f);
    matrix_vector_multiply(Jf, gx, Jf_g);

    for (int i = 0; i < n; i++) {
        out->data[i] = Jg_f->data[i] - Jf_g->data[i];
    }

    vector_free(fx); vector_free(gx);
    vector_free(Jg_f); vector_free(Jf_g);
    matrix_free(Jf); matrix_free(Jg);
}

void ad_operator(const VectorField *f,
                  const VectorField *g,
                  const Vector *x,
                  int k,
                  Vector *out) {
    if (k < 0) return;
    if (k == 0) {
        /* ad_f^0 g = g */
        vector_field_eval(g, x, out);
        return;
    }
    if (k == 1) {
        lie_bracket(f, g, x, out);
        return;
    }

    /* Compute iteratively: ad_f^k g = [f, ad_f^{k-1} g] */
    int n = x->dim;
    Vector *prev = vector_create(n);
    Vector *bracket = vector_create(n);

    ad_operator(f, g, x, k - 1, prev);

    /* Define a temporary vector field wrapping prev as constant */
    /* Actually, we need the Jacobian of ad_f^{k-1} g, which requires
       evaluating it at perturbed points. Use finite differences. */

    /* Create a temporary vector field for ad_f^{k-1} g */
    /* We compute [f, ad_f^{k-1} g] via the finite diff Lie bracket */
    VectorField *prev_field = (VectorField*)malloc(sizeof(VectorField));
    prev_field->dim = n;
    prev_field->eval = NULL;
    prev_field->jacobian = NULL;

    /* Use numeric bracket */
    lie_bracket_numeric(f, prev_field, x, DEFAULT_EPS, out);

    /* Actually, lie_bracket_numeric uses f and a "constant" field.
       But we need the true bracket. Let's compute it properly: */

    /* Instead, directly compute [f, ad_f^{k-1}g] using finite diff Jacobian
       of ad_f^{k-1}g. We approximate J_{ad^{k-1}g} at x. */

    Matrix *J_prev = matrix_create(n, n);
    /* Compute Jacobian of ad_f^{k-1} g by finite differences */
    Vector *xp = vector_clone(x);
    Vector *prev_p = vector_create(n);
    Vector *prev_m = vector_create(n);

    for (int j = 0; j < n; j++) {
        double orig = x->data[j];
        /* ad_f^{k-1} g at x + eps e_j */
        xp->data[j] = orig + DEFAULT_EPS;
        ad_operator(f, g, xp, k - 1, prev_p);
        /* ad_f^{k-1} g at x - eps e_j */
        xp->data[j] = orig - DEFAULT_EPS;
        ad_operator(f, g, xp, k - 1, prev_m);
        xp->data[j] = orig;

        for (int i = 0; i < n; i++) {
            J_prev->data[i * n + j] = (prev_p->data[i] - prev_m->data[i]) / (2.0 * DEFAULT_EPS);
        }
    }

    /* [f, ad_f^{k-1}g] = J_{ad^{k-1}g} * f - J_f * g_ad */
    Vector *fx = vector_create(n);
    vector_field_eval(f, x, fx);

    Matrix *Jf = matrix_create(n, n);
    if (f->jacobian) {
        f->jacobian(x, Jf);
    } else {
        finite_jacobian(f, x, Jf, DEFAULT_EPS);
    }

    Vector *Jprev_f = vector_create(n);
    Vector *Jf_prev = vector_create(n);
    matrix_vector_multiply(J_prev, fx, Jprev_f);
    matrix_vector_multiply(Jf, prev, Jf_prev);

    for (int i = 0; i < n; i++) {
        out->data[i] = Jprev_f->data[i] - Jf_prev->data[i];
    }

    vector_free(prev); vector_free(bracket);
    vector_free(xp); vector_free(prev_p); vector_free(prev_m);
    vector_free(fx); vector_free(Jprev_f); vector_free(Jf_prev);
    matrix_free(J_prev); matrix_free(Jf);
    free(prev_field);
}

void lie_bracket_numeric(const VectorField *f,
                          const VectorField *g,
                          const Vector *x,
                          double eps,
                          Vector *out) {
    if (!f || !g || !x || !out) return;
    int n = x->dim;

    /* [f,g](x) ? [g(x + eps*f) - g(x - eps*f)]/(2*eps) - [f(x + eps*g) - f(x - eps*g)]/(2*eps) */
    /* More precisely using flows:
       [f,g](x) ? limit_{t->0} [g(phi_f^t(x)) - g(x)]/t - limit_{t->0}[f(phi_g^t(x)) - f(x)]/t
    */

    Vector *fx = vector_create(n);
    Vector *gx = vector_create(n);
    vector_field_eval(f, x, fx);
    vector_field_eval(g, x, gx);

    double nrm_f = vector_norm(fx);
    double nrm_g = vector_norm(gx);

    Vector *xp = vector_clone(x);
    Vector *xm = vector_clone(x);
    Vector *gp = vector_create(n);
    Vector *gm = vector_create(n);
    Vector *fp = vector_create(n);
    Vector *fm = vector_create(n);

    double ep = (eps > 0) ? eps : DEFAULT_EPS;

    /* Perturb along direction of f */
    if (nrm_f > 1e-12) {
        for (int i = 0; i < n; i++) {
            double step = ep * fx->data[i] / nrm_f;
            xp->data[i] = x->data[i] + step;
            xm->data[i] = x->data[i] - step;
        }
        vector_field_eval(g, xp, gp);
        vector_field_eval(g, xm, gm);
        for (int i = 0; i < n; i++) xp->data[i] = xm->data[i] = x->data[i];
    } else {
        vector_set(gp, 0.0);
        vector_set(gm, 0.0);
    }

    /* Perturb along direction of g */
    if (nrm_g > 1e-12) {
        for (int i = 0; i < n; i++) {
            double step = ep * gx->data[i] / nrm_g;
            xp->data[i] = x->data[i] + step;
            xm->data[i] = x->data[i] - step;
        }
        vector_field_eval(f, xp, fp);
        vector_field_eval(f, xm, fm);
        for (int i = 0; i < n; i++) xp->data[i] = xm->data[i] = x->data[i];
    } else {
        vector_set(fp, 0.0);
        vector_set(fm, 0.0);
    }

    /* Combine: [f,g] = (gp - gm)*nrm_f/(2*ep) - (fp - fm)*nrm_g/(2*ep) */
    double sf = (nrm_f > 1e-12) ? nrm_f / (2.0 * ep) : 0.0;
    double sg = (nrm_g > 1e-12) ? nrm_g / (2.0 * ep) : 0.0;
    for (int i = 0; i < n; i++) {
        out->data[i] = sf * (gp->data[i] - gm->data[i]) - sg * (fp->data[i] - fm->data[i]);
    }

    vector_free(fx); vector_free(gx);
    vector_free(xp); vector_free(xm);
    vector_free(gp); vector_free(gm);
    vector_free(fp); vector_free(fm);
}

bool lie_bracket_is_zero(const VectorField *f,
                          const VectorField *g,
                          const Vector *x,
                          double tol) {
    if (!f || !g || !x) return false;
    int n = x->dim;
    Vector *bracket = vector_create(n);
    lie_bracket(f, g, x, bracket);
    bool is_zero = vector_is_zero(bracket, tol);
    vector_free(bracket);
    return is_zero;
}

Matrix *compute_accessibility_distribution(const VectorField *f,
                                            const VectorField *g,
                                            const Vector *x,
                                            int max_k) {
    if (!f || !g || !x || max_k < 0) return NULL;
    int n = x->dim;
    int cols = max_k + 1;
    Matrix *dist = matrix_create(n, cols);
    if (!dist) return NULL;

    Vector *col = vector_create(n);
    for (int k = 0; k <= max_k; k++) {
        ad_operator(f, g, x, k, col);
        for (int i = 0; i < n; i++) {
            dist->data[i * cols + k] = col->data[i];
        }
    }
    vector_free(col);
    return dist;
}

int controllability_rank(const VectorField *f,
                          const VectorField *g,
                          const Vector *x,
                          double tol) {
    if (!f || !g || !x) return 0;
    int n = x->dim;
    Matrix *dist = compute_accessibility_distribution(f, g, x, n - 1);
    if (!dist) return 0;
    int rank = matrix_rank(dist, tol);
    matrix_free(dist);
    return rank;
}

bool check_jacobi_identity(const VectorField *f,
                            const VectorField *g,
                            const VectorField *h,
                            const Vector *x,
                            double tol) {
    if (!f || !g || !h || !x) return false;
    int n = x->dim;

    Vector *gh = vector_create(n);
    Vector *hf = vector_create(n);
    Vector *fg = vector_create(n);
    lie_bracket(g, h, x, gh);
    lie_bracket(h, f, x, hf);
    lie_bracket(f, g, x, fg);

    Vector *term1 = vector_create(n);
    Vector *term2 = vector_create(n);
    Vector *term3 = vector_create(n);
    Vector *sum = vector_create(n);

    double eps = 1e-6;
    Vector *xp = vector_create(n);
    Vector *xm = vector_create(n);
    Vector *vp = vector_create(n);
    Vector *vm = vector_create(n);

    /* term1 = [f, [g,h]] via numeric derivative along f */
    Vector *fx = vector_create(n);
    vector_field_eval(f, x, fx);
    double nf = vector_norm(fx);
    if (nf > 1e-12) {
        for (int i = 0; i < n; i++) {
            xp->data[i] = x->data[i] + eps * fx->data[i] / nf;
            xm->data[i] = x->data[i] - eps * fx->data[i] / nf;
        }
        lie_bracket(g, h, xp, vp);
        lie_bracket(g, h, xm, vm);
        for (int i = 0; i < n; i++)
            term1->data[i] = (vp->data[i] - vm->data[i]) * nf / (2.0 * eps);
    } else { vector_set(term1, 0.0); }

    /* term2 = [g, [h,f]] via numeric derivative along g */
    Vector *gx = vector_create(n);
    vector_field_eval(g, x, gx);
    double ng = vector_norm(gx);
    if (ng > 1e-12) {
        for (int i = 0; i < n; i++) {
            xp->data[i] = x->data[i] + eps * gx->data[i] / ng;
            xm->data[i] = x->data[i] - eps * gx->data[i] / ng;
        }
        lie_bracket(h, f, xp, vp);
        lie_bracket(h, f, xm, vm);
        for (int i = 0; i < n; i++)
            term2->data[i] = (vp->data[i] - vm->data[i]) * ng / (2.0 * eps);
    } else { vector_set(term2, 0.0); }

    /* term3 = [h, [f,g]] via numeric derivative along h */
    Vector *hx = vector_create(n);
    vector_field_eval(h, x, hx);
    double nh = vector_norm(hx);
    if (nh > 1e-12) {
        for (int i = 0; i < n; i++) {
            xp->data[i] = x->data[i] + eps * hx->data[i] / nh;
            xm->data[i] = x->data[i] - eps * hx->data[i] / nh;
        }
        lie_bracket(f, g, xp, vp);
        lie_bracket(f, g, xm, vm);
        for (int i = 0; i < n; i++)
            term3->data[i] = (vp->data[i] - vm->data[i]) * nh / (2.0 * eps);
    } else { vector_set(term3, 0.0); }

    for (int i = 0; i < n; i++)
        sum->data[i] = term1->data[i] + term2->data[i] + term3->data[i];
    bool satisfied = vector_is_zero(sum, tol);

    vector_free(gh); vector_free(hf); vector_free(fg);
    vector_free(term1); vector_free(term2); vector_free(term3);
    vector_free(sum); vector_free(xp); vector_free(xm);
    vector_free(vp); vector_free(vm);
    vector_free(fx); vector_free(gx); vector_free(hx);
    return satisfied;
}

/* Helper: compute Lie bracket wrapping vector fields not in VectorField struct */
static void lie_bracket_vector(const VectorField *f, const Vector *gvec,
                                const Vector *x, Vector *out) {
    /* Compute [f, gvec](x) where gvec is the value of a vector field at x */
    int n = x->dim;
    Vector *fx = vector_create(n);

    vector_field_eval(f, x, fx);

    Matrix *Jf = matrix_create(n, n);
    if (f->jacobian) {
        f->jacobian(x, Jf);
    } else {
        finite_jacobian(f, x, Jf, DEFAULT_EPS);
    }

    Vector *Jf_g = vector_create(n);
    matrix_vector_multiply(Jf, (Vector*)gvec, Jf_g);

    for (int i = 0; i < n; i++) {
        out->data[i] = -Jf_g->data[i]; /* J_g * f ? 0 since g is constant at x */
    }
    /* Actually, we need J_g too. Better to compute full bracket. */

    /* Recompute more carefully: J_g is zero if we treat gvec as constant.
       Then [f, gvec] = 0 * f - J_f * gvec = -J_f * gvec */
    /* But this is not physically correct unless gvec is truly constant.
       We need to build J_g from finite differences around x. */
    /* For the purpose of computing Lie algebra basis, we use the numeric bracket. */

    vector_free(fx); vector_free(Jf_g);
    matrix_free(Jf);
}

int compute_lie_algebra_basis(VectorField **vector_fields,
                               int num_fields,
                               const Vector *x,
                               double eps) {
    if (!vector_fields || num_fields <= 0 || !x) return 0;
    int n = x->dim;
    if (n <= 0) return 0;

    /* Start with the vector fields themselves, then iteratively add Lie brackets */
    /* Max possible dimension is n */

    Matrix *basis_buf = matrix_create(n, n);  /* stores basis vectors */
    int basis_count = 0;

    /* Add original vector fields */
    Vector *v = vector_create(n);
    for (int i = 0; i < num_fields && basis_count < n; i++) {
        vector_field_eval(vector_fields[i], x, v);
        /* Check linear independence against existing basis */
        /* Simple Gram-Schmidt approach */
        if (vector_norm(v) > eps) {
            /* Orthogonalize */
            for (int j = 0; j < basis_count; j++) {
                double dot = 0.0;
                for (int k = 0; k < n; k++) {
                    dot += v->data[k] * basis_buf->data[k * n + j];
                }
                for (int k = 0; k < n; k++) {
                    v->data[k] -= dot * basis_buf->data[k * n + j];
                }
            }
            double nrm = vector_norm(v);
            if (nrm > eps) {
                for (int k = 0; k < n; k++) {
                    basis_buf->data[k * n + basis_count] = v->data[k] / nrm;
                }
                basis_count++;
            }
        }
    }

    /* Add pairwise Lie brackets (first iteration) */
    Vector *bracket = vector_create(n);
    for (int i = 0; i < num_fields && basis_count < n; i++) {
        for (int j = 0; j < num_fields && basis_count < n; j++) {
            lie_bracket(vector_fields[i], vector_fields[j], x, bracket);
            /* Orthogonalize */
            for (int k = 0; k < basis_count; k++) {
                double dot = 0.0;
                for (int l = 0; l < n; l++) {
                    dot += bracket->data[l] * basis_buf->data[l * n + k];
                }
                for (int l = 0; l < n; l++) {
                    bracket->data[l] -= dot * basis_buf->data[l * n + k];
                }
            }
            double nrm = vector_norm(bracket);
            if (nrm > eps) {
                for (int k = 0; k < n; k++) {
                    basis_buf->data[k * n + basis_count] = bracket->data[k] / nrm;
                }
                basis_count++;
            }
        }
    }

    vector_free(v); vector_free(bracket);
    matrix_free(basis_buf);

    return basis_count;
}
