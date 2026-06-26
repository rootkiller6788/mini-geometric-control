#include "nonlinear_system.h"
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * Vector Operations Implementation
 * ============================================================================ */

Vector *vector_create(int dim) {
    if (dim <= 0) return NULL;
    Vector *v = (Vector*)malloc(sizeof(Vector));
    if (!v) return NULL;
    v->dim = dim;
    v->data = (double*)calloc((size_t)dim, sizeof(double));
    if (!v->data) {
        free(v);
        return NULL;
    }
    return v;
}

Vector *vector_clone(const Vector *src) {
    if (!src) return NULL;
    Vector *v = vector_create(src->dim);
    if (!v) return NULL;
    memcpy(v->data, src->data, (size_t)src->dim * sizeof(double));
    return v;
}

void vector_set(Vector *v, double value) {
    if (!v) return;
    for (int i = 0; i < v->dim; i++) {
        v->data[i] = value;
    }
}

void vector_set_element(Vector *v, int i, double value) {
    if (!v || i < 0 || i >= v->dim) return;
    v->data[i] = value;
}

double vector_get(const Vector *v, int i) {
    if (!v || i < 0 || i >= v->dim) return 0.0;
    return v->data[i];
}

void vector_add(const Vector *v1, const Vector *v2, Vector *v_out) {
    if (!v1 || !v2 || !v_out) return;
    if (v1->dim != v2->dim || v1->dim != v_out->dim) return;
    for (int i = 0; i < v1->dim; i++) {
        v_out->data[i] = v1->data[i] + v2->data[i];
    }
}

void vector_sub(const Vector *v1, const Vector *v2, Vector *v_out) {
    if (!v1 || !v2 || !v_out) return;
    if (v1->dim != v2->dim || v1->dim != v_out->dim) return;
    for (int i = 0; i < v1->dim; i++) {
        v_out->data[i] = v1->data[i] - v2->data[i];
    }
}

void vector_scale(const Vector *v, double scalar, Vector *v_out) {
    if (!v || !v_out || v->dim != v_out->dim) return;
    for (int i = 0; i < v->dim; i++) {
        v_out->data[i] = v->data[i] * scalar;
    }
}

double vector_dot(const Vector *v1, const Vector *v2) {
    if (!v1 || !v2 || v1->dim != v2->dim) return 0.0;
    double d = 0.0;
    for (int i = 0; i < v1->dim; i++) {
        d += v1->data[i] * v2->data[i];
    }
    return d;
}

double vector_norm(const Vector *v) {
    if (!v) return 0.0;
    double s = 0.0;
    for (int i = 0; i < v->dim; i++) {
        s += v->data[i] * v->data[i];
    }
    return sqrt(s);
}

double vector_max_norm(const Vector *v) {
    if (!v || v->dim == 0) return 0.0;
    double m = fabs(v->data[0]);
    for (int i = 1; i < v->dim; i++) {
        double a = fabs(v->data[i]);
        if (a > m) m = a;
    }
    return m;
}

double vector_normalize(Vector *v) {
    if (!v) return 0.0;
    double nrm = vector_norm(v);
    if (nrm < 1e-15) return 0.0;
    for (int i = 0; i < v->dim; i++) {
        v->data[i] /= nrm;
    }
    return nrm;
}

double vector_distance(const Vector *v1, const Vector *v2) {
    if (!v1 || !v2 || v1->dim != v2->dim) return INFINITY;
    double s = 0.0;
    for (int i = 0; i < v1->dim; i++) {
        double d = v1->data[i] - v2->data[i];
        s += d * d;
    }
    return sqrt(s);
}

void vector_lincomb(double a, const Vector *v1, double b, const Vector *v2, Vector *v_out) {
    if (!v1 || !v2 || !v_out) return;
    if (v1->dim != v2->dim || v1->dim != v_out->dim) return;
    for (int i = 0; i < v1->dim; i++) {
        v_out->data[i] = a * v1->data[i] + b * v2->data[i];
    }
}

bool vector_equals(const Vector *v1, const Vector *v2, double tol) {
    if (!v1 || !v2 || v1->dim != v2->dim) return false;
    for (int i = 0; i < v1->dim; i++) {
        if (fabs(v1->data[i] - v2->data[i]) > tol) return false;
    }
    return true;
}

bool vector_is_zero(const Vector *v, double tol) {
    if (!v) return true;
    for (int i = 0; i < v->dim; i++) {
        if (fabs(v->data[i]) > tol) return false;
    }
    return true;
}

void vector_print(const Vector *v, const char *label) {
    if (!v) return;
    printf("%s [%d]: ", label ? label : "vector", v->dim);
    for (int i = 0; i < v->dim; i++) {
        printf("%12.6f ", v->data[i]);
    }
    printf("\n");
}

void vector_free(Vector *v) {
    if (!v) return;
    free(v->data);
    free(v);
}

/* ============================================================================
 * Matrix Operations Implementation
 * ============================================================================ */

Matrix *matrix_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return NULL;
    Matrix *m = (Matrix*)malloc(sizeof(Matrix));
    if (!m) return NULL;
    m->rows = rows;
    m->cols = cols;
    m->data = (double*)calloc((size_t)(rows * cols), sizeof(double));
    if (!m->data) {
        free(m);
        return NULL;
    }
    return m;
}

Matrix *matrix_create_identity(int n) {
    Matrix *m = matrix_create(n, n);
    if (!m) return NULL;
    for (int i = 0; i < n; i++) {
        m->data[i * n + i] = 1.0;
    }
    return m;
}

Matrix *matrix_clone(const Matrix *src) {
    if (!src) return NULL;
    Matrix *m = matrix_create(src->rows, src->cols);
    if (!m) return NULL;
    memcpy(m->data, src->data, (size_t)(src->rows * src->cols) * sizeof(double));
    return m;
}

double matrix_get(const Matrix *m, int i, int j) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return 0.0;
    return m->data[i * m->cols + j];
}

void matrix_set(Matrix *m, int i, int j, double value) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return;
    m->data[i * m->cols + j] = value;
}

void matrix_add(const Matrix *A, const Matrix *B, Matrix *C) {
    if (!A || !B || !C) return;
    if (A->rows != B->rows || A->cols != B->cols) return;
    if (A->rows != C->rows || A->cols != C->cols) return;
    int n = A->rows * A->cols;
    for (int i = 0; i < n; i++) {
        C->data[i] = A->data[i] + B->data[i];
    }
}

void matrix_sub(const Matrix *A, const Matrix *B, Matrix *C) {
    if (!A || !B || !C) return;
    if (A->rows != B->rows || A->cols != B->cols) return;
    if (A->rows != C->rows || A->cols != C->cols) return;
    int n = A->rows * A->cols;
    for (int i = 0; i < n; i++) {
        C->data[i] = A->data[i] - B->data[i];
    }
}

void matrix_scale(const Matrix *A, double scalar, Matrix *B) {
    if (!A || !B) return;
    if (A->rows != B->rows || A->cols != B->cols) return;
    int n = A->rows * A->cols;
    for (int i = 0; i < n; i++) {
        B->data[i] = A->data[i] * scalar;
    }
}

void matrix_multiply(const Matrix *A, const Matrix *B, Matrix *C) {
    if (!A || !B || !C) return;
    if (A->cols != B->rows) return;
    if (C->rows != A->rows || C->cols != B->cols) return;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->cols; j++) {
            double s = 0.0;
            for (int k = 0; k < A->cols; k++) {
                s += A->data[i * A->cols + k] * B->data[k * B->cols + j];
            }
            C->data[i * C->cols + j] = s;
        }
    }
}

void matrix_vector_multiply(const Matrix *A, const Vector *x, Vector *y) {
    if (!A || !x || !y) return;
    if (A->cols != x->dim || A->rows != y->dim) return;
    for (int i = 0; i < A->rows; i++) {
        double s = 0.0;
        for (int j = 0; j < A->cols; j++) {
            s += A->data[i * A->cols + j] * x->data[j];
        }
        y->data[i] = s;
    }
}

void matrix_transpose(const Matrix *A, Matrix *B) {
    if (!A || !B) return;
    if (A->rows != B->cols || A->cols != B->rows) return;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            B->data[j * B->cols + i] = A->data[i * A->cols + j];
        }
    }
}

double matrix_determinant(const Matrix *A) {
    if (!A || A->rows != A->cols) return 0.0;
    int n = A->rows;
    /* Use LU decomposition with partial pivoting */
    Matrix *LU = matrix_clone(A);
    if (!LU) return 0.0;

    int sign = 1;
    /* Perform LU decomposition in-place */
    for (int k = 0; k < n; k++) {
        /* Find pivot */
        double max_val = fabs(LU->data[k * n + k]);
        int max_row = k;
        for (int i = k + 1; i < n; i++) {
            double v = fabs(LU->data[i * n + k]);
            if (v > max_val) {
                max_val = v;
                max_row = i;
            }
        }
        /* Swap rows if needed */
        if (max_row != k) {
            sign = -sign;
            for (int j = 0; j < n; j++) {
                double tmp = LU->data[k * n + j];
                LU->data[k * n + j] = LU->data[max_row * n + j];
                LU->data[max_row * n + j] = tmp;
            }
        }
        if (fabs(LU->data[k * n + k]) < 1e-15) {
            matrix_free(LU);
            return 0.0;
        }
        /* Eliminate below */
        for (int i = k + 1; i < n; i++) {
            double factor = LU->data[i * n + k] / LU->data[k * n + k];
            LU->data[i * n + k] = factor;
            for (int j = k + 1; j < n; j++) {
                LU->data[i * n + j] -= factor * LU->data[k * n + j];
            }
        }
    }
    double det = sign;
    for (int i = 0; i < n; i++) {
        det *= LU->data[i * n + i];
    }
    matrix_free(LU);
    return det;
}

bool matrix_inverse(const Matrix *A, Matrix *A_inv) {
    if (!A || !A_inv || A->rows != A->cols) return false;
    if (A->rows != A_inv->rows || A->cols != A_inv->cols) return false;
    int n = A->rows;

    /* Augment [A | I] */
    Matrix *aug = matrix_create(n, 2 * n);
    if (!aug) return false;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug->data[i * 2 * n + j] = A->data[i * n + j];
        }
        aug->data[i * 2 * n + n + i] = 1.0;
    }

    /* Gauss-Jordan elimination */
    for (int k = 0; k < n; k++) {
        /* Partial pivot */
        double max_val = fabs(aug->data[k * 2 * n + k]);
        int max_row = k;
        for (int i = k + 1; i < n; i++) {
            double v = fabs(aug->data[i * 2 * n + k]);
            if (v > max_val) { max_val = v; max_row = i; }
        }
        if (max_val < 1e-15) { matrix_free(aug); return false; }
        if (max_row != k) {
            for (int j = 0; j < 2 * n; j++) {
                double tmp = aug->data[k * 2 * n + j];
                aug->data[k * 2 * n + j] = aug->data[max_row * 2 * n + j];
                aug->data[max_row * 2 * n + j] = tmp;
            }
        }
        /* Normalize pivot row */
        double pivot = aug->data[k * 2 * n + k];
        for (int j = 0; j < 2 * n; j++) {
            aug->data[k * 2 * n + j] /= pivot;
        }
        /* Eliminate other rows */
        for (int i = 0; i < n; i++) {
            if (i == k) continue;
            double factor = aug->data[i * 2 * n + k];
            for (int j = 0; j < 2 * n; j++) {
                aug->data[i * 2 * n + j] -= factor * aug->data[k * 2 * n + j];
            }
        }
    }
    /* Extract inverse */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_inv->data[i * n + j] = aug->data[i * 2 * n + n + j];
        }
    }
    matrix_free(aug);
    return true;
}

double matrix_trace(const Matrix *A) {
    if (!A || A->rows != A->cols) return 0.0;
    double tr = 0.0;
    int n = (A->rows < A->cols) ? A->rows : A->cols;
    for (int i = 0; i < n; i++) tr += A->data[i * A->cols + i];
    return tr;
}

double matrix_frobenius_norm(const Matrix *A) {
    if (!A) return 0.0;
    double s = 0.0;
    int n = A->rows * A->cols;
    for (int i = 0; i < n; i++) s += A->data[i] * A->data[i];
    return sqrt(s);
}

int matrix_rank(const Matrix *A, double tol) {
    if (!A) return 0;
    /* Use QR decomposition via modified Gram-Schmidt */
    Matrix *R = matrix_clone(A);
    if (!R) return 0;
    int m = A->rows, n = A->cols;

    int rank = 0;
    /* Modified Gram-Schmidt on columns */
    double *col_norms = (double*)malloc((size_t)n * sizeof(double));
    if (!col_norms) { matrix_free(R); return 0; }

    for (int j = 0; j < n; j++) {
        double nsq = 0.0;
        for (int i = 0; i < m; i++) {
            double v = R->data[i * n + j];
            nsq += v * v;
        }
        col_norms[j] = sqrt(nsq);
    }

    for (int j = 0; j < n; j++) {
        if (col_norms[j] > tol) {
            rank++;
            /* Orthogonalize remaining columns against column j */
            for (int k = j + 1; k < n; k++) {
                double dot = 0.0;
                for (int i = 0; i < m; i++) {
                    dot += R->data[i * n + j] * R->data[i * n + k];
                }
                double factor = dot / (col_norms[j] * col_norms[j]);
                for (int i = 0; i < m; i++) {
                    R->data[i * n + k] -= factor * R->data[i * n + j];
                }
                /* Update norm */
                double nsq = 0.0;
                for (int i = 0; i < m; i++) {
                    double v = R->data[i * n + k];
                    nsq += v * v;
                }
                col_norms[k] = sqrt(nsq);
            }
        }
    }
    free(col_norms);
    matrix_free(R);
    return rank;
}

int matrix_nullspace(const Matrix *A, Matrix **basis, double tol) {
    if (!A || !basis) return -1;
    int m = A->rows, n = A->cols;

    /* Compute SVD via the Jacobi method for small matrices, or */
    /* Use augmented QR on A^T. Compute nullspace dimension = n - rank */
    int r = matrix_rank(A, tol);
    int nullity = n - r;
    if (nullity <= 0) { *basis = NULL; return 0; }

    /* Build [A; eps*I] and compute null vectors via elimination */
    /* Simplified: compute nullspace vectors via random perturbation + orthogonalization */
    Matrix *U = matrix_clone(A);
    *basis = matrix_create(n, nullity);

    /* Use row-reduced echelon form to find nullspace basis */
    /* Augment with identity and do Gaussian elimination */
    Matrix *aug = matrix_create(m, n + n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            aug->data[i * (n + n) + j] = A->data[i * n + j];
        }
    }

    /* Track pivot columns */
    bool *pivot_col = (bool*)calloc((size_t)n, sizeof(bool));
    bool *free_col = (bool*)calloc((size_t)n, sizeof(bool));

    Matrix *work = matrix_clone(A);
    int piv_row = 0;
    for (int col = 0; col < n && piv_row < m; col++) {
        /* Find pivot */
        double max_v = 0.0;
        int max_r = piv_row;
        for (int r = piv_row; r < m; r++) {
            double v = fabs(work->data[r * n + col]);
            if (v > max_v) { max_v = v; max_r = r; }
        }
        if (max_v < tol) continue; /* linearly dependent column */

        pivot_col[col] = true;
        /* Swap rows */
        if (max_r != piv_row) {
            for (int j = 0; j < n; j++) {
                double t = work->data[piv_row * n + j];
                work->data[piv_row * n + j] = work->data[max_r * n + j];
                work->data[max_r * n + j] = t;
            }
        }
        /* Normalize */
        double pv = work->data[piv_row * n + col];
        for (int j = col; j < n; j++) work->data[piv_row * n + j] /= pv;

        /* Eliminate other rows */
        for (int r = 0; r < m; r++) {
            if (r == piv_row) continue;
            double factor = work->data[r * n + col];
            if (fabs(factor) < tol) continue;
            for (int j = col; j < n; j++) {
                work->data[r * n + j] -= factor * work->data[piv_row * n + j];
            }
        }
        piv_row++;
    }

    /* Free columns are those without pivot */
    for (int j = 0; j < n; j++) {
        if (!pivot_col[j]) free_col[j] = true;
    }

    /* Build nullspace basis */
    int basis_idx = 0;
    for (int col = 0; col < n && basis_idx < nullity; col++) {
        if (!free_col[col]) continue;
        for (int i = 0; i < n; i++) {
            double val = (i == col) ? 1.0 : 0.0;
            if (pivot_col[i]) {
                /* Find the row where column i is pivot */
                val = -work->data[i * n + col]; /* approximate */
            }
            (*basis)->data[i * nullity + basis_idx] = val;
        }
        basis_idx++;
    }

    /* Orthonormalize basis using Gram-Schmidt */
    for (int k = 0; k < nullity; k++) {
        /* Normalize column k */
        double nrm = 0.0;
        for (int i = 0; i < n; i++) {
            double v = (*basis)->data[i * nullity + k];
            nrm += v * v;
        }
        nrm = sqrt(nrm);
        if (nrm > tol) {
            for (int i = 0; i < n; i++) {
                (*basis)->data[i * nullity + k] /= nrm;
            }
        }
        /* Orthogonalize against previous */
        for (int prev = 0; prev < k; prev++) {
            double dot = 0.0;
            for (int i = 0; i < n; i++) {
                dot += (*basis)->data[i * nullity + k] * (*basis)->data[i * nullity + prev];
            }
            for (int i = 0; i < n; i++) {
                (*basis)->data[i * nullity + k] -= dot * (*basis)->data[i * nullity + prev];
            }
        }
    }

    free(pivot_col);
    free(free_col);
    matrix_free(work);
    matrix_free(U);
    matrix_free(aug);
    return nullity;
}

int matrix_range(const Matrix *A, Matrix **basis, double tol) {
    if (!A || !basis) return -1;
    int r = matrix_rank(A, tol);
    if (r <= 0) { *basis = NULL; return 0; }

    int m = A->rows, n = A->cols;
    *basis = matrix_create(m, r);

    /* Use column pivoted QR to select independent columns */
    Matrix *R = matrix_clone(A);
    bool *selected = (bool*)calloc((size_t)n, sizeof(bool));
    int basis_cols = 0;

    for (int j = 0; j < n && basis_cols < r; j++) {
        double max_nrm = 0.0;
        int best_j = -1;
        for (int c = 0; c < n; c++) {
            if (selected[c]) continue;
            double nrm = 0.0;
            for (int i = 0; i < m; i++) {
                double v = R->data[i * n + c];
                nrm += v * v;
            }
            if (nrm > max_nrm) { max_nrm = nrm; best_j = c; }
        }
        if (best_j < 0 || max_nrm < tol * tol) break;

        selected[best_j] = true;
        /* Copy to basis */
        for (int i = 0; i < m; i++) {
            (*basis)->data[i * r + basis_cols] = A->data[i * n + best_j];
        }
        basis_cols++;

        /* Orthogonalize remaining */
        for (int c = 0; c < n; c++) {
            if (selected[c]) continue;
            double dot = 0.0;
            for (int i = 0; i < m; i++) {
                dot += R->data[i * n + best_j] * R->data[i * n + c];
            }
            double factor = dot / max_nrm;
            for (int i = 0; i < m; i++) {
                R->data[i * n + c] -= factor * R->data[i * n + best_j];
            }
        }
    }

    free(selected);
    matrix_free(R);
    return r;
}

bool matrix_equals(const Matrix *A, const Matrix *B, double tol) {
    if (!A || !B) return false;
    if (A->rows != B->rows || A->cols != B->cols) return false;
    int n = A->rows * A->cols;
    for (int i = 0; i < n; i++) {
        if (fabs(A->data[i] - B->data[i]) > tol) return false;
    }
    return true;
}

void matrix_print(const Matrix *m, const char *label) {
    if (!m) return;
    printf("%s [%d x %d]:\n", label ? label : "matrix", m->rows, m->cols);
    for (int i = 0; i < m->rows; i++) {
        printf("  ");
        for (int j = 0; j < m->cols; j++) {
            printf("%10.4f ", m->data[i * m->cols + j]);
        }
        printf("\n");
    }
}

void matrix_free(Matrix *m) {
    if (!m) return;
    free(m->data);
    free(m);
}

/* ============================================================================
 * Vector Field and Scalar Field Operations
 * ============================================================================ */

VectorField *vector_field_create(int dim,
                                  void (*eval)(const Vector*, Vector*),
                                  void (*jacobian)(const Vector*, Matrix*)) {
    if (dim <= 0) return NULL;
    VectorField *vf = (VectorField*)malloc(sizeof(VectorField));
    if (!vf) return NULL;
    vf->dim = dim;
    vf->eval = eval;
    vf->jacobian = jacobian;
    return vf;
}

ScalarField *scalar_field_create(int dim,
                                  double (*eval)(const Vector*),
                                  void (*gradient)(const Vector*, Vector*)) {
    if (dim <= 0) return NULL;
    ScalarField *sf = (ScalarField*)malloc(sizeof(ScalarField));
    if (!sf) return NULL;
    sf->dim = dim;
    sf->eval = eval;
    sf->gradient = gradient;
    return sf;
}

void vector_field_eval(const VectorField *vf, const Vector *x, Vector *out) {
    if (!vf || !x || !out) return;
    if (vf->eval) {
        vf->eval(x, out);
    } else {
        vector_set(out, 0.0);
    }
}

double scalar_field_eval(const ScalarField *sf, const Vector *x) {
    if (!sf || !x) return 0.0;
    if (sf->eval) return sf->eval(x);
    return 0.0;
}

void vector_field_free(VectorField *vf) {
    /* Note: does not free callbacks, which are typically static functions */
    free(vf);
}

void scalar_field_free(ScalarField *sf) {
    free(sf);
}

/* ============================================================================
 * Nonlinear System Operations
 * ============================================================================ */

NonlinearSystem *nonlinear_system_create(int n, int m, int p, const char *name) {
    NonlinearSystem *sys = (NonlinearSystem*)malloc(sizeof(NonlinearSystem));
    if (!sys) return NULL;
    sys->n = n;
    sys->m = m;
    sys->p = p;
    sys->f = NULL;
    sys->g = (m > 0) ? (VectorField**)calloc((size_t)m, sizeof(VectorField*)) : NULL;
    sys->h = (p > 0) ? (ScalarField**)calloc((size_t)p, sizeof(ScalarField*)) : NULL;
    sys->name = name ? strdup(name) : strdup("unnamed");
    sys->equilibrium = NULL;
    return sys;
}

void nonlinear_system_set_drift(NonlinearSystem *sys, VectorField *f) {
    if (!sys) return;
    if (sys->f) vector_field_free(sys->f);
    sys->f = f;
}

void nonlinear_system_set_control(NonlinearSystem *sys, int i, VectorField *g_i) {
    if (!sys || i < 0 || i >= sys->m) return;
    sys->g[i] = g_i;
}

void nonlinear_system_set_output(NonlinearSystem *sys, int i, ScalarField *h_i) {
    if (!sys || i < 0 || i >= sys->p) return;
    sys->h[i] = h_i;
}

void nonlinear_system_set_equilibrium(NonlinearSystem *sys, const Vector *x_eq) {
    if (!sys) return;
    if (sys->equilibrium) vector_free(sys->equilibrium);
    sys->equilibrium = vector_clone(x_eq);
}

void nonlinear_system_eval_dynamics(const NonlinearSystem *sys,
                                     const Vector *x,
                                     const Vector *u,
                                     Vector *x_dot) {
    if (!sys || !x || !x_dot) return;
    if (x->dim != sys->n || x_dot->dim != sys->n) return;

    /* Evaluate drift */
    if (sys->f) {
        vector_field_eval(sys->f, x, x_dot);
    } else {
        vector_set(x_dot, 0.0);
    }

    /* Add control contributions */
    for (int i = 0; i < sys->m; i++) {
        if (!sys->g[i]) continue;
        Vector *gx = vector_create(sys->n);
        vector_field_eval(sys->g[i], x, gx);
        double ui = (u && i < u->dim) ? u->data[i] : 0.0;
        for (int j = 0; j < sys->n; j++) {
            x_dot->data[j] += ui * gx->data[j];
        }
        vector_free(gx);
    }
}

void nonlinear_system_eval_output(const NonlinearSystem *sys,
                                   const Vector *x,
                                   Vector *y) {
    if (!sys || !x || !y) return;
    for (int i = 0; i < sys->p; i++) {
        if (sys->h[i]) {
            y->data[i] = scalar_field_eval(sys->h[i], x);
        } else {
            y->data[i] = x->data[i]; /* default: first p states */
        }
    }
}

void nonlinear_system_eval_without_control(const NonlinearSystem *sys,
                                            const Vector *x,
                                            Vector *f_out) {
    if (!sys || !x || !f_out) return;
    if (sys->f) {
        vector_field_eval(sys->f, x, f_out);
    } else {
        vector_set(f_out, 0.0);
    }
}

void nonlinear_system_linearize_at(const NonlinearSystem *sys,
                                    const Vector *x_eq,
                                    Matrix *A,
                                    Matrix *B,
                                    Matrix *C) {
    if (!sys || !x_eq) return;

    /* A = df/dx at equilibrium */
    if (A && sys->f && sys->f->jacobian) {
        sys->f->jacobian(x_eq, A);
    }

    /* B = [g_1(x_eq) ... g_m(x_eq)] */
    if (B) {
        for (int j = 0; j < sys->m; j++) {
            if (sys->g[j]) {
                Vector *gx = vector_create(sys->n);
                vector_field_eval(sys->g[j], x_eq, gx);
                for (int i = 0; i < sys->n; i++) {
                    B->data[i * B->cols + j] = gx->data[i];
                }
                vector_free(gx);
            }
        }
    }

    /* C = dh/dx at equilibrium */
    if (C) {
        for (int i = 0; i < sys->p; i++) {
            if (sys->h[i] && sys->h[i]->gradient) {
                Vector *grad = vector_create(sys->n);
                sys->h[i]->gradient(x_eq, grad);
                for (int j = 0; j < sys->n; j++) {
                    C->data[i * C->cols + j] = grad->data[j];
                }
                vector_free(grad);
            }
        }
    }
}

bool nonlinear_system_is_equilibrium(const NonlinearSystem *sys,
                                      const Vector *x_eq,
                                      double tol) {
    if (!sys || !x_eq) return false;
    Vector *fx = vector_create(sys->n);
    nonlinear_system_eval_without_control(sys, x_eq, fx);
    bool is_eq = vector_is_zero(fx, tol);
    vector_free(fx);
    return is_eq;
}

void nonlinear_system_print_info(const NonlinearSystem *sys) {
    if (!sys) return;
    printf("NonlinearSystem: %s\n", sys->name);
    printf("  State dim (n): %d\n", sys->n);
    printf("  Input  dim (m): %d\n", sys->m);
    printf("  Output dim (p): %d\n", sys->p);
    printf("  Drift f: %s\n", sys->f ? "defined" : "NULL (autonomous)");
    for (int i = 0; i < sys->m; i++) {
        printf("  Control g_%d: %s\n", i, sys->g[i] ? "defined" : "NULL");
    }
    for (int i = 0; i < sys->p; i++) {
        printf("  Output h_%d: %s\n", i, sys->h[i] ? "defined" : "NULL");
    }
    if (sys->equilibrium) {
        vector_print(sys->equilibrium, "  Equilibrium x_eq");
    }
}

void nonlinear_system_free(NonlinearSystem *sys) {
    if (!sys) return;
    if (sys->f) vector_field_free(sys->f);
    if (sys->g) {
        for (int i = 0; i < sys->m; i++) {
            if (sys->g[i]) vector_field_free(sys->g[i]);
        }
        free(sys->g);
    }
    if (sys->h) {
        for (int i = 0; i < sys->p; i++) {
            if (sys->h[i]) scalar_field_free(sys->h[i]);
        }
        free(sys->h);
    }
    free(sys->name);
    if (sys->equilibrium) vector_free(sys->equilibrium);
    free(sys);
}
