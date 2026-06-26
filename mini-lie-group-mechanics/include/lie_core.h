#ifndef LIE_CORE_H
#define LIE_CORE_H

#include <stdbool.h>
#include <stddef.h>
#include <math.h>

/* ==============================================================
 * lie_core.h -- Lie Group and Lie Algebra Core Definitions
 *
 * A Lie group G is a smooth manifold that is also a group, where
 * multiplication and inversion are smooth maps.
 *
 * Lie algebra g = T_e G: the tangent space at identity, equipped
 * with the Lie bracket [.,.] satisfying:
 *   (i)   Bilinearity
 *   (ii)  Antisymmetry: [xi, eta] = -[eta, xi]
 *   (iii) Jacobi identity: [xi,[eta,zeta]] + [eta,[zeta,xi]] + [zeta,[xi,eta]] = 0
 *
 * Key Lie groups in mechanics and control:
 *   SO(3)  -- 3D rotations (rigid body attitude)
 *   SE(3)  -- 3D rigid motions (position + orientation)
 *   SO(2)  -- 2D rotations
 *   SE(2)  -- 2D planar rigid motions
 *   SU(2)  -- Spin group / quaternions (double cover of SO(3))
 *   GL(n)  -- General linear group
 *   SL(n)  -- Special linear group
 *
 * References:
 *   Marsden & Ratiu (1999) Intro to Mechanics and Symmetry, Ch.9
 *   Holm, Schmah, Stoica (2009) Geometric Mechanics and Symmetry
 *   Hall (2015) Lie Groups, Lie Algebras, and Representations
 *   Iserles, Munthe-Kaas, Norsett, Zanna (2000) Lie Group Methods
 *   Murray, Li, Sastry (1994) A Mathematical Intro to Robotic
 *     Manipulation, Appendix A
 * ============================================================== */

/* --- Lie Group Enumeration --- */
typedef enum {
    LIE_GROUP_SO2    = 0,
    LIE_GROUP_SO3    = 1,
    LIE_GROUP_SE2    = 2,
    LIE_GROUP_SE3    = 3,
    LIE_GROUP_SU2    = 4,
    LIE_GROUP_GL3    = 5,
    LIE_GROUP_SL3    = 6,
    LIE_GROUP_CUSTOM = 7
} LieGroupType;

/* --- Lie Algebra Type --- */
typedef enum {
    LIE_ALG_SO2   = 0,
    LIE_ALG_SO3   = 1,
    LIE_ALG_SE2   = 2,
    LIE_ALG_SE3   = 3,
    LIE_ALG_SU2   = 4,
    LIE_ALG_GL3   = 5,
    LIE_ALG_SL3   = 6,
    LIE_ALG_CUSTOM = 7
} LieAlgebraType;

/* --- Matrix (row-major flat storage) --- */
typedef struct {
    double *data;
    int rows;
    int cols;
} LieMatrix;

/* --- Vector --- */
typedef struct {
    double *data;
    int size;
} LieVector;

/* --- Lie Group Element --- */
typedef struct {
    char *name;
    LieGroupType group_type;
    int dim;
    int algebra_dim;
    int mat_size;
    LieMatrix *g;
    LieVector *params;
} LieGroupElement;

/* --- Lie Algebra Element --- */
typedef struct {
    char *name;
    LieAlgebraType algebra_type;
    int dim;
    LieMatrix *mat;
    LieVector *coords;
} LieAlgebraElement;

/* --- Matrix Operations --- */
LieMatrix* lie_matrix_create(int rows, int cols);
void lie_matrix_free(LieMatrix *m);
LieMatrix* lie_matrix_copy(const LieMatrix *src);
void lie_matrix_set(LieMatrix *m, int row, int col, double val);
double lie_matrix_get(const LieMatrix *m, int row, int col);
LieMatrix* lie_matrix_identity(int n);
LieMatrix* lie_matrix_zeros(int rows, int cols);
LieMatrix* lie_matrix_diag(int n, const double *diag);
void lie_matrix_fill(LieMatrix *m, double val);

LieMatrix* lie_matrix_add(const LieMatrix *a, const LieMatrix *b);
LieMatrix* lie_matrix_sub(const LieMatrix *a, const LieMatrix *b);
LieMatrix* lie_matrix_mul(const LieMatrix *a, const LieMatrix *b);
LieMatrix* lie_matrix_scale(const LieMatrix *m, double s);
LieMatrix* lie_matrix_transpose(const LieMatrix *m);

double lie_matrix_trace(const LieMatrix *m);
double lie_matrix_det(const LieMatrix *m);
double lie_matrix_det_2x2(const LieMatrix *m);
double lie_matrix_det_3x3(const LieMatrix *m);
double lie_matrix_det_4x4(const LieMatrix *m);
double lie_matrix_norm_frobenius(const LieMatrix *m);
double lie_matrix_norm_infinity(const LieMatrix *m);
bool lie_matrix_is_identity(const LieMatrix *m, double tol);
bool lie_matrix_is_skew_symmetric(const LieMatrix *m, double tol);
bool lie_matrix_is_orthogonal(const LieMatrix *m, double tol);
bool lie_matrix_is_special_orthogonal(const LieMatrix *m, double tol);
bool lie_matrix_approx_equal(const LieMatrix *a, const LieMatrix *b, double tol);

LieMatrix* lie_matrix_inverse_2x2(const LieMatrix *m);
LieMatrix* lie_matrix_inverse_3x3(const LieMatrix *m);
LieMatrix* lie_matrix_inverse_4x4(const LieMatrix *m);
LieMatrix* lie_matrix_inverse(const LieMatrix *m);
void lie_matrix_solve_linear(LieMatrix *A, double *b, int n, double *x);

LieVector* lie_vector_create(int size);
void lie_vector_free(LieVector *v);
void lie_vector_set(LieVector *v, int idx, double val);
double lie_vector_get(const LieVector *v, int idx);
LieVector* lie_vector_copy(const LieVector *src);
double lie_vector_norm(const LieVector *v);
double lie_vector_norm_sq(const LieVector *v);
double lie_vector_dot(const LieVector *a, const LieVector *b);
LieVector* lie_vector_cross3(const LieVector *a, const LieVector *b);
LieVector* lie_vector_add(const LieVector *a, const LieVector *b);
LieVector* lie_vector_sub(const LieVector *a, const LieVector *b);
LieVector* lie_vector_scale(const LieVector *v, double s);
void lie_vector_normalize(LieVector *v);

LieGroupElement* lie_group_create(LieGroupType gtype, const char *name);
void lie_group_free(LieGroupElement *g);
void lie_group_set_identity(LieGroupElement *g);
LieGroupElement* lie_group_mul(const LieGroupElement *g1, const LieGroupElement *g2);
LieGroupElement* lie_group_inverse(const LieGroupElement *g);
void lie_group_invert_inplace(LieGroupElement *g);
LieVector* lie_group_act_vector(const LieGroupElement *g, const LieVector *x);

LieAlgebraElement* lie_algebra_create(LieAlgebraType atype, const char *name);
void lie_algebra_free(LieAlgebraElement *xi);
LieAlgebraElement* lie_algebra_copy(const LieAlgebraElement *src);
void lie_algebra_set_zero(LieAlgebraElement *xi);

LieAlgebraElement* lie_bracket(const LieAlgebraElement *xi,
                                const LieAlgebraElement *eta);
bool lie_jacobi_identity_check(const LieAlgebraElement *xi,
                                const LieAlgebraElement *eta,
                                const LieAlgebraElement *zeta, double tol);

/* Hat and Vee Maps */
LieMatrix* lie_so3_hat(const LieVector *omega);
LieVector* lie_so3_vee(const LieMatrix *omega_hat);
LieMatrix* lie_se3_hat(const LieVector *xi);
LieVector* lie_se3_vee(const LieMatrix *xi_hat);
LieMatrix* lie_su2_from_vector(const LieVector *v);
LieVector* lie_su2_to_vector(const LieMatrix *u);

/* Exponential Map */
LieMatrix* lie_matrix_exponential(const LieMatrix *A);
LieGroupElement* lie_exp_so3(const LieVector *omega);
LieGroupElement* lie_exp_se3(const LieVector *xi);
LieGroupElement* lie_group_exp(const LieAlgebraElement *xi);

/* Logarithm Map */
LieVector* lie_log_so3(const LieGroupElement *R);
LieVector* lie_log_se3(const LieGroupElement *g);
LieMatrix* lie_matrix_logarithm(const LieMatrix *G);

/* Adjoint Representations */
LieAlgebraElement* lie_Ad(const LieGroupElement *g, const LieAlgebraElement *eta);
LieMatrix* lie_Ad_matrix_so3(const LieMatrix *R);
LieMatrix* lie_Ad_matrix_se3(const LieMatrix *R, const LieVector *p);
LieAlgebraElement* lie_ad(const LieAlgebraElement *xi, const LieAlgebraElement *eta);
LieMatrix* lie_ad_matrix_so3(const LieVector *omega);
LieMatrix* lie_ad_matrix_se3(const LieVector *xi);

/* Coadjoint Representations */
LieVector* lie_Ad_star_so3(const LieMatrix *R, const LieVector *mu);
LieVector* lie_ad_star_so3(const LieVector *omega, const LieVector *mu);

/* Killing Form */
double lie_killing_form(const LieAlgebraElement *xi, const LieAlgebraElement *eta);

/* BCH Formula */
LieAlgebraElement* lie_bch_2nd_order(const LieAlgebraElement *A,
                                      const LieAlgebraElement *B);
LieAlgebraElement* lie_bch_3rd_order(const LieAlgebraElement *A,
                                      const LieAlgebraElement *B);

/* Display */
void lie_matrix_print(const LieMatrix *m, const char *label);
void lie_vector_print(const LieVector *v, const char *label);
void lie_group_print(const LieGroupElement *g);
void lie_algebra_print(const LieAlgebraElement *xi);

#endif /* LIE_CORE_H */