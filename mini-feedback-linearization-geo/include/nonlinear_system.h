#ifndef NONLINEAR_SYSTEM_H
#define NONLINEAR_SYSTEM_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Nonlinear Control-Affine System Representation
 *
 * Standard form: x' = f(x) + sum_{i=1}^m g_i(x) u_i
 *                y  = h(x)
 *
 * where:
 *   x in R^n        state vector
 *   u in R^m        input vector
 *   y in R^p        output vector
 *   f: R^n -> R^n   drift vector field
 *   g_i: R^n -> R^n control vector fields (i = 1,...,m)
 *   h: R^n -> R^p   output map
 *
 * References:
 *   Isidori (1995) Nonlinear Control Systems, 3rd ed.
 *   Khalil (2002) Nonlinear Systems, 3rd ed.
 *   Sastry (1999) Nonlinear Systems: Analysis, Stability, and Control
 * ============================================================================ */

/* --- Vector type (dense, fixed dimension) --- */

typedef struct {
    double *data;
    int dim;
} Vector;

/* --- Matrix type (dense, row-major) --- */

typedef struct {
    double *data;
    int rows;
    int cols;
} Matrix;

/* --- Vector field: a smooth map R^n -> R^n --- */

typedef struct {
    /* Function pointer evaluating f at x: out = f(x) */
    void (*eval)(const Vector *x, Vector *out);
    /* Jacobian of f at x: J = df/dx */
    void (*jacobian)(const Vector *x, Matrix *J);
    int dim;   /* dimension of domain/codomain */
} VectorField;

/* --- Scalar field: a smooth map R^n -> R --- */

typedef struct {
    /* Function pointer evaluating h at x */
    double (*eval)(const Vector *x);
    /* Gradient of h at x: grad_h = dh/dx */
    void (*gradient)(const Vector *x, Vector *grad);
    int dim;   /* domain dimension (n) */
} ScalarField;

/* --- Nonlinear control-affine system --- */

typedef struct {
    int n;                     /* State dimension */
    int m;                     /* Number of inputs */
    int p;                     /* Number of outputs */

    VectorField *f;            /* Drift vector field f(x) */
    VectorField **g;           /* Control vector fields g_1(x), ..., g_m(x) */
    ScalarField **h;           /* Output maps h_1(x), ..., h_p(x) */

    char *name;                /* System identifier */
    Vector *equilibrium;       /* Equilibrium point (f(x_eq) = 0 typically) */
} NonlinearSystem;

/* ============================================================================
 * Vector Operations
 * ============================================================================ */

/** vector_create -- Allocate zero vector of dimension dim. */
Vector *vector_create(int dim);

/** vector_clone -- Deep copy a vector. */
Vector *vector_clone(const Vector *src);

/** vector_set -- Set all elements to value. */
void vector_set(Vector *v, double value);

/** vector_set_element -- Set a single element. */
void vector_set_element(Vector *v, int i, double value);

/** vector_get -- Get element value (bounds-checked). */
double vector_get(const Vector *v, int i);

/** vector_add -- v3 = v1 + v2 (in-place: v_out can equal v1). */
void vector_add(const Vector *v1, const Vector *v2, Vector *v_out);

/** vector_sub -- v3 = v1 - v2. */
void vector_sub(const Vector *v1, const Vector *v2, Vector *v_out);

/** vector_scale -- v_out = scalar * v. */
void vector_scale(const Vector *v, double scalar, Vector *v_out);

/** vector_dot -- Compute dot product v1 ? v2. */
double vector_dot(const Vector *v1, const Vector *v2);

/** vector_norm -- Compute Euclidean norm ||v||_2. */
double vector_norm(const Vector *v);

/** vector_max_norm -- Compute infinity norm ||v||_inf. */
double vector_max_norm(const Vector *v);

/** vector_normalize -- Normalize v to unit length. Returns original norm. */
double vector_normalize(Vector *v);

/** vector_distance -- Euclidean distance ||v1 - v2||_2. */
double vector_distance(const Vector *v1, const Vector *v2);

/** vector_lincomb -- v_out = a * v1 + b * v2. */
void vector_lincomb(double a, const Vector *v1, double b, const Vector *v2, Vector *v_out);

/** vector_equals -- Check element-wise equality within tolerance. */
bool vector_equals(const Vector *v1, const Vector *v2, double tol);

/** vector_is_zero -- Check if all elements are within tol of zero. */
bool vector_is_zero(const Vector *v, double tol);

/** vector_print -- Print vector to stdout. */
void vector_print(const Vector *v, const char *label);

/** vector_free -- Free memory. */
void vector_free(Vector *v);

/* ============================================================================
 * Matrix Operations
 * ============================================================================ */

/** matrix_create -- Allocate zero matrix of size rows x cols. */
Matrix *matrix_create(int rows, int cols);

/** matrix_create_identity -- Allocate n x n identity matrix. */
Matrix *matrix_create_identity(int n);

/** matrix_clone -- Deep copy a matrix. */
Matrix *matrix_clone(const Matrix *src);

/** matrix_get -- Get element at (i, j). */
double matrix_get(const Matrix *m, int i, int j);

/** matrix_set -- Set element at (i, j). */
void matrix_set(Matrix *m, int i, int j, double value);

/** matrix_add -- C = A + B. */
void matrix_add(const Matrix *A, const Matrix *B, Matrix *C);

/** matrix_sub -- C = A - B. */
void matrix_sub(const Matrix *A, const Matrix *B, Matrix *C);

/** matrix_scale -- B = scalar * A. */
void matrix_scale(const Matrix *A, double scalar, Matrix *B);

/** matrix_multiply -- C = A * B (C must be pre-allocated: rows_A x cols_B). */
void matrix_multiply(const Matrix *A, const Matrix *B, Matrix *C);

/** matrix_vector_multiply -- y = A * x. */
void matrix_vector_multiply(const Matrix *A, const Vector *x, Vector *y);

/** matrix_transpose -- B = A^T. */
void matrix_transpose(const Matrix *A, Matrix *B);

/** matrix_determinant -- Compute det(A) for square matrix (using LU decomposition). */
double matrix_determinant(const Matrix *A);

/** matrix_inverse -- Compute A^{-1} for square non-singular matrix (using Gauss-Jordan). */
bool matrix_inverse(const Matrix *A, Matrix *A_inv);

/** matrix_trace -- Compute trace(A) = sum of diagonal elements. */
double matrix_trace(const Matrix *A);

/** matrix_frobenius_norm -- Compute ||A||_F = sqrt(sum a_ij^2). */
double matrix_frobenius_norm(const Matrix *A);

/** matrix_rank -- Compute numerical rank using SVD-like QR decomposition with tolerance. */
int matrix_rank(const Matrix *A, double tol);

/** matrix_nullspace -- Compute basis for nullspace of A. Returns dimension of nullspace. */
int matrix_nullspace(const Matrix *A, Matrix **basis, double tol);

/** matrix_range -- Compute basis for range/column space of A. Returns rank. */
int matrix_range(const Matrix *A, Matrix **basis, double tol);

/** matrix_equals -- Check element-wise equality within tolerance. */
bool matrix_equals(const Matrix *A, const Matrix *B, double tol);

/** matrix_print -- Print matrix to stdout. */
void matrix_print(const Matrix *m, const char *label);

/** matrix_free -- Free memory. */
void matrix_free(Matrix *m);

/* ============================================================================
 * Vector Field and Scalar Field Operations
 * ============================================================================ */

/** vector_field_create -- Allocate vector field with given eval/jacobian callbacks. */
VectorField *vector_field_create(int dim,
                                  void (*eval)(const Vector*, Vector*),
                                  void (*jacobian)(const Vector*, Matrix*));

/** scalar_field_create -- Allocate scalar field with given eval/gradient callbacks. */
ScalarField *scalar_field_create(int dim,
                                  double (*eval)(const Vector*),
                                  void (*gradient)(const Vector*, Vector*));

/** vector_field_eval -- Evaluate vector field at x. */
void vector_field_eval(const VectorField *vf, const Vector *x, Vector *out);

/** scalar_field_eval -- Evaluate scalar field at x. */
double scalar_field_eval(const ScalarField *sf, const Vector *x);

/** vector_field_free -- Free vector field. */
void vector_field_free(VectorField *vf);

/** scalar_field_free -- Free scalar field. */
void scalar_field_free(ScalarField *sf);

/* ============================================================================
 * Nonlinear System Operations
 * ============================================================================ */

/** nonlinear_system_create -- Allocate control-affine system. */
NonlinearSystem *nonlinear_system_create(int n, int m, int p, const char *name);

/** nonlinear_system_set_drift -- Set drift vector field f. */
void nonlinear_system_set_drift(NonlinearSystem *sys, VectorField *f);

/** nonlinear_system_set_control -- Set i-th control vector field g_i. */
void nonlinear_system_set_control(NonlinearSystem *sys, int i, VectorField *g_i);

/** nonlinear_system_set_output -- Set i-th output map h_i. */
void nonlinear_system_set_output(NonlinearSystem *sys, int i, ScalarField *h_i);

/** nonlinear_system_set_equilibrium -- Set equilibrium point. */
void nonlinear_system_set_equilibrium(NonlinearSystem *sys, const Vector *x_eq);

/** nonlinear_system_eval_dynamics -- Evaluate state derivative x_dot = f(x) + sum g_i(x) u_i. */
void nonlinear_system_eval_dynamics(const NonlinearSystem *sys,
                                     const Vector *x,
                                     const Vector *u,
                                     Vector *x_dot);

/** nonlinear_system_eval_output -- Evaluate output y = h(x). */
void nonlinear_system_eval_output(const NonlinearSystem *sys,
                                   const Vector *x,
                                   Vector *y);

/** nonlinear_system_eval_without_control -- Evaluate f(x) only (no control). */
void nonlinear_system_eval_without_control(const NonlinearSystem *sys,
                                            const Vector *x,
                                            Vector *f_out);

/** nonlinear_system_linearize_at -- Linearize at equilibrium: x' = A dx + B du, y = C dx.
 *  Computes Jacobians A = df/dx|_eq, B = [g_1 ... g_m], C = dh/dx|_eq. */
void nonlinear_system_linearize_at(const NonlinearSystem *sys,
                                    const Vector *x_eq,
                                    Matrix *A,
                                    Matrix *B,
                                    Matrix *C);

/** nonlinear_system_is_equilibrium -- Check if f(x_eq) = 0 within tolerance. */
bool nonlinear_system_is_equilibrium(const NonlinearSystem *sys,
                                      const Vector *x_eq,
                                      double tol);

/** nonlinear_system_free -- Free system and all components. */
void nonlinear_system_free(NonlinearSystem *sys);

/** nonlinear_system_print_info -- Print system summary to stdout. */
void nonlinear_system_print_info(const NonlinearSystem *sys);

#endif /* NONLINEAR_SYSTEM_H */
