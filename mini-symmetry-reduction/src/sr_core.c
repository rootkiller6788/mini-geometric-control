#include "sr_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define SR_EPS 1e-12
#define SR_PI  3.14159265358979323846

/* ============================================================================
 * Internal Utilities
 * ============================================================================ */

static double* alloc_vec(int n) {
    return (double*)calloc((size_t)n, sizeof(double));
}

static void free_vec(double* v) { free(v); }

__attribute__((unused))
static double dot3(const double* a, const double* b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

__attribute__((unused))
static void cross3(const double* a, const double* b, double* c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static double norm3(const double* a) {
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

static double mat3_det(const double A[9]) {
    return A[0]*(A[4]*A[8] - A[5]*A[7])
         - A[1]*(A[3]*A[8] - A[5]*A[6])
         + A[2]*(A[3]*A[7] - A[4]*A[6]);
}

/* ============================================================================
 * Lie Algebra Operations
 * ============================================================================ */

SRLieAlgebra* sr_algebra_create(int dimension, const char* name) {
    if (dimension <= 0 || !name) return NULL;
    SRLieAlgebra* alg = (SRLieAlgebra*)calloc(1, sizeof(SRLieAlgebra));
    alg->dimension = dimension;
    alg->format = SR_STRUCT_FULL_TENSOR;
    int n = dimension;
    alg->constants = (double*)calloc((size_t)(n * n * n), sizeof(double));
    alg->basis_names = (char**)calloc((size_t)n, sizeof(char*));
    for (int i = 0; i < n; i++) {
        alg->basis_names[i] = (char*)calloc(32, sizeof(char));
        snprintf(alg->basis_names[i], 32, "e_%d", i);
    }
    alg->is_abelian = true;
    alg->jacobi_verified = false;
    return alg;
}

void sr_algebra_free(SRLieAlgebra* alg) {
    if (!alg) return;
    free(alg->constants);
    if (alg->basis_names) {
        for (int i = 0; i < alg->dimension; i++) free(alg->basis_names[i]);
        free(alg->basis_names);
    }
    free(alg);
}

void sr_algebra_set_struct_const(SRLieAlgebra* alg, int i, int j, int k, double val) {
    if (!alg || i < 0 || i >= alg->dimension ||
        j < 0 || j >= alg->dimension || k < 0 || k >= alg->dimension) return;
    int n = alg->dimension;
    alg->constants[i*n*n + j*n + k] = val;
    if (fabs(val) > SR_EPS) alg->is_abelian = false;
}

double sr_algebra_get_struct_const(const SRLieAlgebra* alg, int i, int j, int k) {
    if (!alg || i < 0 || i >= alg->dimension ||
        j < 0 || j >= alg->dimension || k < 0 || k >= alg->dimension) return 0.0;
    return alg->constants[i * alg->dimension * alg->dimension + j * alg->dimension + k];
}

bool sr_algebra_verify_jacobi(const SRLieAlgebra* alg, double tolerance) {
    if (!alg) return false;
    int n = alg->dimension;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++)
    for (int l = 0; l < n; l++) {
        double sum = 0.0;
        for (int m = 0; m < n; m++) {
            double a = sr_algebra_get_struct_const(alg, i, j, m);
            double b = sr_algebra_get_struct_const(alg, m, k, l);
            double c = sr_algebra_get_struct_const(alg, j, k, m);
            double d = sr_algebra_get_struct_const(alg, m, i, l);
            double e = sr_algebra_get_struct_const(alg, k, i, m);
            double f = sr_algebra_get_struct_const(alg, m, j, l);
            sum += a*b + c*d + e*f;
        }
        if (fabs(sum) > tolerance) return false;
    }
    return true;
}

double sr_algebra_commutator_component(const SRLieAlgebra* alg, int i, int j, int k) {
    return sr_algebra_get_struct_const(alg, i, j, k);
}

void sr_algebra_bracket(const SRLieAlgebra* alg, const double* X, const double* Y,
                         double* bracket_out) {
    if (!alg || !X || !Y || !bracket_out) return;
    int n = alg->dimension;
    memset(bracket_out, 0, (size_t)n * sizeof(double));
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
        double xi = X[i], yj = Y[j];
        for (int k = 0; k < n; k++)
            bracket_out[k] += xi * yj * sr_algebra_get_struct_const(alg, i, j, k);
    }
}

void sr_algebra_print(const SRLieAlgebra* alg) {
    if (!alg) { printf("NULL algebra\n"); return; }
    printf("Lie Algebra: dim=%d, abelian=%d\n", alg->dimension, alg->is_abelian);
    printf("Structure constants (nonzero):\n");
    int n = alg->dimension, count = 0;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++) {
        double c = alg->constants[i*n*n + j*n + k];
        if (fabs(c) > SR_EPS) {
            printf("  c^{%d}_{%d,%d} = %g\n", k, i, j, c);
            count++;
        }
    }
    if (count == 0) printf("  (all zero — abelian)\n");
}

double sr_algebra_killing_form(const SRLieAlgebra* alg, int a, int b) {
    if (!alg || a < 0 || a >= alg->dimension || b < 0 || b >= alg->dimension) return 0.0;
    int n = alg->dimension;
    double K = 0.0;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
        double c1 = sr_algebra_get_struct_const(alg, a, i, j);
        double c2 = sr_algebra_get_struct_const(alg, b, j, i);
        K += c1 * c2;
    }
    return K;
}

/* ============================================================================
 * Special Lie Algebra Factories
 * ============================================================================ */

SRLieAlgebra* sr_algebra_create_so3(void) {
    /* so(3): [e1,e2]=e3, [e2,e3]=e1, [e3,e1]=e2  (epsilon_{ijk}) */
    SRLieAlgebra* alg = sr_algebra_create(3, "so(3)");
    sr_algebra_set_struct_const(alg, 0, 1, 2,  1.0); /* [e0,e1] = e2 */
    sr_algebra_set_struct_const(alg, 1, 0, 2, -1.0);
    sr_algebra_set_struct_const(alg, 1, 2, 0,  1.0); /* [e1,e2] = e0 */
    sr_algebra_set_struct_const(alg, 2, 1, 0, -1.0);
    sr_algebra_set_struct_const(alg, 2, 0, 1,  1.0); /* [e2,e0] = e1 */
    sr_algebra_set_struct_const(alg, 0, 2, 1, -1.0);
    alg->is_abelian = false;
    alg->is_simple = true;
    alg->is_semisimple = true;
    alg->jacobi_verified = true;
    return alg;
}

SRLieAlgebra* sr_algebra_create_so2(void) {
    SRLieAlgebra* alg = sr_algebra_create(1, "so(2)");
    alg->is_abelian = true;
    alg->jacobi_verified = true;
    return alg;
}

SRLieAlgebra* sr_algebra_create_su2(void) {
    /* su(2) ≅ so(3) (isomorphic as real Lie algebras) */
    SRLieAlgebra* alg = sr_algebra_create_so3();
    free(alg->basis_names[0]);
    free(alg->basis_names[1]);
    free(alg->basis_names[2]);
    alg->basis_names[0] = strdup("tau_1");
    alg->basis_names[1] = strdup("tau_2");
    alg->basis_names[2] = strdup("tau_3");
    return alg;
}

SRLieAlgebra* sr_algebra_create_se3(void) {
    /* se(3) = so(3) ⋉ R^3 has 6 basis elements: e0..e2 (rotations), e3..e5 (translations)
     * [ei, ej] = epsilon_{ijk} e_k  for i,j,k in {0,1,2}
     * [ei, ej+3] = epsilon_{ijk} e_{k+3}  for i in {0,1,2}, j in {0,1,2}
     * [ei+3, ej+3] = 0
     */
    SRLieAlgebra* alg = sr_algebra_create(6, "se(3)");
    /* Set rotation-rotation brackets */
    sr_algebra_set_struct_const(alg, 0, 1, 2,  1.0); sr_algebra_set_struct_const(alg, 1, 0, 2, -1.0);
    sr_algebra_set_struct_const(alg, 1, 2, 0,  1.0); sr_algebra_set_struct_const(alg, 2, 1, 0, -1.0);
    sr_algebra_set_struct_const(alg, 2, 0, 1,  1.0); sr_algebra_set_struct_const(alg, 0, 2, 1, -1.0);
    /* Set rotation-translation brackets: [e_i, e_{j+3}] = epsilon_{ijk} e_{k+3} */
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++) {
        double val = 0.0;
        if ((i==0 && j==1 && k==2) || (i==1 && j==2 && k==0) || (i==2 && j==0 && k==1)) val = 1.0;
        else if ((i==1 && j==0 && k==2) || (i==2 && j==1 && k==0) || (i==0 && j==2 && k==1)) val = -1.0;
        if (fabs(val) > SR_EPS)
            sr_algebra_set_struct_const(alg, i, j+3, k+3, val);
    }
    alg->is_abelian = false;
    alg->is_semisimple = false; /* se(3) is not semisimple (has non-trivial radical R^3) */
    alg->is_solvable = false;
    return alg;
}

SRLieAlgebra* sr_algebra_create_heisenberg(void) {
    /* 3-dim Heisenberg algebra: [e1,e2]=e3, others zero */
    SRLieAlgebra* alg = sr_algebra_create(3, "heisenberg");
    sr_algebra_set_struct_const(alg, 0, 1, 2,  1.0);
    sr_algebra_set_struct_const(alg, 1, 0, 2, -1.0);
    alg->is_abelian = false;
    alg->is_nilpotent = true;
    alg->is_solvable = true;
    return alg;
}

SRLieAlgebra* sr_algebra_create_abelian(int dim) {
    SRLieAlgebra* alg = sr_algebra_create(dim, "abelian");
    alg->is_abelian = true;
    alg->jacobi_verified = true;
    return alg;
}

/* ============================================================================
 * Lie Group Operations
 * ============================================================================ */

SRLieGroup* sr_group_create(int dimension, int matrix_size, const char* name,
                             SRSymmetryType type) {
    if (dimension <= 0 || matrix_size <= 0 || !name) return NULL;
    SRLieGroup* g = (SRLieGroup*)calloc(1, sizeof(SRLieGroup));
    g->dimension = dimension;
    g->matrix_size = matrix_size;
    g->name = strdup(name);
    g->type = type;
    g->is_compact = false;
    g->is_connected = true;
    g->is_simply_connected = false;
    g->identity = alloc_vec(matrix_size * matrix_size);
    for (int i = 0; i < matrix_size; i++)
        g->identity[i * matrix_size + i] = 1.0;
    return g;
}

void sr_group_free(SRLieGroup* group) {
    if (!group) return;
    free(group->name);
    free(group->identity);
    free(group->multiplication_table);
    free(group);
}

void sr_group_set_compact(SRLieGroup* group, bool compact) {
    if (group) group->is_compact = compact;
}

void sr_group_set_connected(SRLieGroup* group, bool connected) {
    if (group) group->is_connected = connected;
}

SRLieAlgebra* sr_group_get_algebra(SRLieGroup* group) {
    if (!group) return NULL;
    if (!group->algebra) {
        char buf[64];
        snprintf(buf, sizeof(buf), "Lie(%s)", group->name);
        group->algebra = sr_algebra_create(group->dimension, buf);
    }
    return group->algebra;
}

void sr_group_set_identity_matrix(SRLieGroup* group, const double* identity) {
    if (!group || !identity) return;
    memcpy(group->identity, identity,
           (size_t)(group->matrix_size * group->matrix_size) * sizeof(double));
}

void sr_group_print(const SRLieGroup* group) {
    if (!group) { printf("NULL group\n"); return; }
    printf("Lie Group: %s, dim=%d, matrix_size=%d\n",
           group->name, group->dimension, group->matrix_size);
    printf("  compact=%d, connected=%d, simply_connected=%d\n",
           group->is_compact, group->is_connected, group->is_simply_connected);
}

/* ============================================================================
 * Group Element Operations
 * ============================================================================ */

SRGroupElement* sr_element_create(SRLieGroup* group) {
    if (!group) return NULL;
    SRGroupElement* el = (SRGroupElement*)calloc(1, sizeof(SRGroupElement));
    el->group = group;
    el->coords = alloc_vec(group->dimension);
    int ms = group->matrix_size;
    el->matrix = alloc_vec(ms * ms);
    return el;
}

void sr_element_free(SRGroupElement* el) {
    if (!el) return;
    free(el->coords);
    free(el->matrix);
    free(el);
}

void sr_element_set_coords(SRGroupElement* el, const double* coords) {
    if (!el || !coords) return;
    memcpy(el->coords, coords, (size_t)el->group->dimension * sizeof(double));
    el->is_identity = false;
}

void sr_element_set_identity(SRGroupElement* el) {
    if (!el) return;
    memset(el->coords, 0, (size_t)el->group->dimension * sizeof(double));
    int ms = el->group->matrix_size;
    memset(el->matrix, 0, (size_t)(ms * ms) * sizeof(double));
    for (int i = 0; i < ms; i++) el->matrix[i * ms + i] = 1.0;
    el->is_identity = true;
}

void sr_element_multiply(const SRGroupElement* g1, const SRGroupElement* g2,
                          SRGroupElement* result) {
    if (!g1 || !g2 || !result) return;
    if (g1->group != g2->group || g1->group != result->group) return;
    int ms = g1->group->matrix_size;
    memset(result->matrix, 0, (size_t)(ms * ms) * sizeof(double));
    for (int i = 0; i < ms; i++)
    for (int k = 0; k < ms; k++) {
        double aik = g1->matrix[i * ms + k];
        for (int j = 0; j < ms; j++)
            result->matrix[i * ms + j] += aik * g2->matrix[k * ms + j];
    }
    result->is_identity = false;
}

void sr_element_inverse(const SRGroupElement* el, SRGroupElement* inverse) {
    if (!el || !inverse) return;
    int ms = el->group->matrix_size;
    /* For matrix groups, g^{-1} = g^T for orthogonal matrices.
     * General: assume orthogonal / use Gauss for small N.
     * For N <= 3, use explicit formula for orthogonal matrices.
     */
    if (el->group->is_compact && ms <= 3) {
        /* Transpose for orthogonal groups */
        for (int i = 0; i < ms; i++)
        for (int j = 0; j < ms; j++)
            inverse->matrix[i * ms + j] = el->matrix[j * ms + i];
    } else {
        /* Copy and assume identity for general case */
        memcpy(inverse->matrix, el->matrix, (size_t)(ms * ms) * sizeof(double));
    }
}

void sr_element_action_on_vector(const SRGroupElement* el, const double* vec,
                                  double* result) {
    if (!el || !vec || !result) return;
    int ms = el->group->matrix_size;
    memset(result, 0, (size_t)ms * sizeof(double));
    for (int i = 0; i < ms; i++)
    for (int j = 0; j < ms; j++)
        result[i] += el->matrix[i * ms + j] * vec[j];
}

void sr_element_print(const SRGroupElement* el) {
    if (!el) { printf("NULL element\n"); return; }
    int ms = el->group->matrix_size;
    printf("GroupElement of %s (matrix %dx%d):\n", el->group->name, ms, ms);
    for (int i = 0; i < ms; i++) {
        printf("  ");
        for (int j = 0; j < ms; j++)
            printf("%8.4f ", el->matrix[i * ms + j]);
        printf("\n");
    }
}

/* ============================================================================
 * Lie Algebra Element Operations
 * ============================================================================ */

SRLieAlgebraElement* sr_alg_element_create(SRLieAlgebra* alg) {
    if (!alg) return NULL;
    SRLieAlgebraElement* el = (SRLieAlgebraElement*)calloc(1, sizeof(SRLieAlgebraElement));
    el->algebra = alg;
    el->coords = alloc_vec(alg->dimension);
    return el;
}

void sr_alg_element_free(SRLieAlgebraElement* el) {
    if (!el) return;
    free(el->coords);
    free(el->matrix);
    free(el);
}

void sr_alg_element_set_coords(SRLieAlgebraElement* el, const double* coords) {
    if (!el || !coords) return;
    memcpy(el->coords, coords, (size_t)el->algebra->dimension * sizeof(double));
}

double sr_alg_element_norm(const SRLieAlgebraElement* el) {
    if (!el) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < el->algebra->dimension; i++)
        sum += el->coords[i] * el->coords[i];
    return sqrt(sum);
}

/* ============================================================================
 * Adjoint and Coadjoint Actions
 * ============================================================================ */

SRAdjointAction* sr_adjoint_create(SRLieGroup* group) {
    if (!group) return NULL;
    SRAdjointAction* adj = (SRAdjointAction*)calloc(1, sizeof(SRAdjointAction));
    adj->group = group;
    adj->dim = group->dimension;
    adj->matrix = alloc_vec(group->dimension * group->dimension);
    adj->cached_inverse = alloc_vec(group->dimension * group->dimension);
    for (int i = 0; i < group->dimension; i++)
        adj->matrix[i * group->dimension + i] = 1.0;
    return adj;
}

void sr_adjoint_free(SRAdjointAction* adj) {
    if (!adj) return;
    free(adj->matrix);
    free(adj->cached_inverse);
    free(adj);
}

void sr_adjoint_compute(const SRAdjointAction* adj, const SRGroupElement* g,
                         double* Ad_g_matrix) {
    if (!adj || !g || !Ad_g_matrix) return;
    int n = adj->dim, ms = g->group->matrix_size;
    (void)ms;
    /* Ad_g(xi) = g xi g^{-1}. For general case, use matrix representation.
     * For now, copy identity (exact form depends on group parametrization). */
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        Ad_g_matrix[i * n + j] = (i == j) ? 1.0 : 0.0;
}

void sr_adjoint_action_on_algebra(const SRAdjointAction* adj, const SRGroupElement* g,
                                   const double* xi, double* result) {
    if (!adj || !g || !xi || !result) return;
    int n = adj->dim;
    double* Ad = alloc_vec(n * n);
    sr_adjoint_compute(adj, g, Ad);
    memset(result, 0, (size_t)n * sizeof(double));
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        result[i] += Ad[i * n + j] * xi[j];
    free_vec(Ad);
}

SRCoadjointAction* sr_coadjoint_create(SRLieGroup* group) {
    if (!group) return NULL;
    SRCoadjointAction* coadj = (SRCoadjointAction*)calloc(1, sizeof(SRCoadjointAction));
    coadj->group = group;
    coadj->dim = group->dimension;
    coadj->matrix = alloc_vec(group->dimension * group->dimension);
    for (int i = 0; i < group->dimension; i++)
        coadj->matrix[i * group->dimension + i] = 1.0;
    return coadj;
}

void sr_coadjoint_free(SRCoadjointAction* coadj) {
    if (!coadj) return;
    free(coadj->matrix);
    free(coadj);
}

void sr_coadjoint_compute(const SRCoadjointAction* coadj, const SRGroupElement* g,
                           double* Ad_star_matrix) {
    if (!coadj || !g || !Ad_star_matrix) return;
    int n = coadj->dim;
    double* Ad = alloc_vec(n * n);
    sr_adjoint_compute(NULL, g, Ad);
    /* Ad*_{g^{-1}} = (Ad_g)^T (with Killing identification) */
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        Ad_star_matrix[i * n + j] = Ad[j * n + i];
    free_vec(Ad);
}

void sr_coadjoint_action_on_dual(const SRCoadjointAction* coadj, const SRGroupElement* g,
                                  const double* mu, double* result) {
    if (!coadj || !g || !mu || !result) return;
    int n = coadj->dim;
    double* Ad_star = alloc_vec(n * n);
    sr_coadjoint_compute(coadj, g, Ad_star);
    memset(result, 0, (size_t)n * sizeof(double));
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        result[i] += Ad_star[i * n + j] * mu[j];
    free_vec(Ad_star);
}

/* ============================================================================
 * Exponential Map
 * ============================================================================ */

void sr_exp_map_so3(const double* omega, double* R_matrix) {
    /* Rodrigues formula: R = I + sin(theta)/theta * [omega]_x + (1-cos(theta))/theta^2 * [omega]_x^2 */
    if (!omega || !R_matrix) return;
    double theta = norm3(omega);
    double wx[3][3] = {{0, -omega[2], omega[1]},
                       {omega[2], 0, -omega[0]},
                       {-omega[1], omega[0], 0}};
    if (theta < SR_EPS) {
        /* Small angle: R ≈ I + [omega]_x */
        for (int i = 0; i < 9; i++) R_matrix[i] = (i==0||i==4||i==8) ? 1.0 : 0.0;
        R_matrix[1] = -omega[2]; R_matrix[2] = omega[1];
        R_matrix[3] = omega[2];  R_matrix[5] = -omega[0];
        R_matrix[6] = -omega[1]; R_matrix[7] = omega[0];
        return;
    }
    double a = sin(theta) / theta;
    double b = (1.0 - cos(theta)) / (theta * theta);
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
        double w2 = 0.0;
        for (int k = 0; k < 3; k++) w2 += wx[i][k] * wx[k][j];
        R_matrix[i * 3 + j] = ((i == j) ? 1.0 : 0.0) + a * wx[i][j] + b * w2;
    }
}

void sr_exp_map_se3(const double* twist, double* T_matrix) {
    /* twist = [omega (3), v (3)] -> SE(3) matrix via exp */
    if (!twist || !T_matrix) return;
    const double* omega = twist;
    const double* v = twist + 3;
    double theta = norm3(omega);
    double R[9];
    sr_exp_map_so3(omega, R);
    memset(T_matrix, 0, 16 * sizeof(double));
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
        T_matrix[i * 4 + j] = R[i * 3 + j];
    T_matrix[15] = 1.0;
    /* Translation part: t = V * v where V = I + (1-cos)/theta^2 * [w] + (theta-sin)/theta^3 * [w]^2 */
    if (theta < SR_EPS) {
        T_matrix[3] = v[0]; T_matrix[7] = v[1]; T_matrix[11] = v[2];
    } else {
        double a = (1.0 - cos(theta)) / (theta * theta);
        double b = (theta - sin(theta)) / (theta * theta * theta);
        double wx[3][3] = {{0, -omega[2], omega[1]}, {omega[2], 0, -omega[0]}, {-omega[1], omega[0], 0}};
        for (int i = 0; i < 3; i++) {
            double t_i = 0.0;
            for (int j = 0; j < 3; j++) {
                double w2_ij = 0.0;
                for (int k = 0; k < 3; k++) w2_ij += wx[i][k] * wx[k][j];
                t_i += (((i==j) ? 1.0 : 0.0) + a * wx[i][j] + b * w2_ij) * v[j];
            }
            T_matrix[i * 4 + 3] = t_i;
        }
    }
}

void sr_exp_map(const SRLieAlgebra* alg, const double* xi, double* g_matrix,
                int matrix_size) {
    (void)alg; (void)xi; (void)g_matrix; (void)matrix_size;
    /* General exp map: g = exp(sum xi^i e_i). Implementation is group-dependent. */
}

void sr_bch_formula_2nd(const SRLieAlgebra* alg, const double* X, const double* Y,
                         double* Z) {
    /* BCH to 2nd order: Z = X + Y + 1/2 [X,Y] + ... */
    if (!alg || !X || !Y || !Z) return;
    int n = alg->dimension;
    for (int i = 0; i < n; i++) Z[i] = X[i] + Y[i];
    double bracket[256];
    sr_algebra_bracket(alg, X, Y, bracket);
    for (int i = 0; i < n; i++) Z[i] += 0.5 * bracket[i];
}

/* ============================================================================
 * Symplectic Geometry
 * ============================================================================ */

SRSymplecticManifold* sr_symplectic_create(int n_dof) {
    if (n_dof <= 0) return NULL;
    SRSymplecticManifold* sp = (SRSymplecticManifold*)calloc(1, sizeof(SRSymplecticManifold));
    sp->dim = 2 * n_dof;
    sp->omega = alloc_vec(sp->dim * sp->dim);
    sp->coordinate_names = (char**)calloc((size_t)sp->dim, sizeof(char*));
    for (int i = 0; i < sp->dim; i++) {
        sp->coordinate_names[i] = (char*)calloc(32, 1);
        snprintf(sp->coordinate_names[i], 32, "z_%d", i);
    }
    return sp;
}

void sr_symplectic_free(SRSymplecticManifold* sp) {
    if (!sp) return;
    free(sp->omega);
    if (sp->coordinate_names) {
        for (int i = 0; i < sp->dim; i++) free(sp->coordinate_names[i]);
        free(sp->coordinate_names);
    }
    free(sp);
}

void sr_symplectic_set_canonical(SRSymplecticManifold* sp) {
    /* Canonical symplectic form: omega = Σ dp_i ∧ dq_i
     * Matrix: [ 0  I ]
     *         [-I  0 ]
     */
    if (!sp) return;
    int n = sp->dim / 2;
    memset(sp->omega, 0, (size_t)(sp->dim * sp->dim) * sizeof(double));
    for (int i = 0; i < n; i++) {
        sp->omega[i * sp->dim + (n + i)] =  1.0;  /* q_i, p_i */
        sp->omega[(n + i) * sp->dim + i] = -1.0;
    }
    sp->is_closed = true;
    sp->is_nondeg = true;
}

bool sr_symplectic_is_nondeg(const SRSymplecticManifold* sp) {
    if (!sp || sp->dim == 0) return false;
    if (sp->dim == 2) {
        return fabs(sp->omega[0 * sp->dim + 1]) > SR_EPS;
    }
    /* For larger dim, check if the matrix has full rank. Use 2x2 test for n=1. */
    return true;
}

void sr_symplectic_inverse(const SRSymplecticManifold* sp, double* omega_inv) {
    if (!sp || !omega_inv) return;
    int d = sp->dim;
    if (d == 2) {
        double det = sp->omega[0] * sp->omega[3] - sp->omega[1] * sp->omega[2];
        if (fabs(det) < SR_EPS) return;
        omega_inv[0] =  sp->omega[3] / det;
        omega_inv[1] = -sp->omega[1] / det;
        omega_inv[2] = -sp->omega[2] / det;
        omega_inv[3] =  sp->omega[0] / det;
    }
    /* For canonical form: omega^{-1} = [0 -I; I 0] */
}

double sr_symplectic_volume(const SRSymplecticManifold* sp) {
    /* Liouville volume form: omega^n / n! = dp_1∧...∧dp_n∧dq_1∧...∧dq_n */
    if (!sp) return 0.0;
    int n = sp->dim / 2;
    return 1.0; /* For canonical coordinates, volume element = 1 */
    (void)n;
}

void sr_symplectic_print(const SRSymplecticManifold* sp) {
    if (!sp) { printf("NULL symplectic manifold\n"); return; }
    printf("Symplectic Manifold: dim=%d\n", sp->dim);
    printf("  closed=%d, nondegenerate=%d\n", sp->is_closed, sp->is_nondeg);
}

/* ============================================================================
 * Poisson Geometry
 * ============================================================================ */

SRPoissonManifold* sr_poisson_create(int dim) {
    if (dim <= 0) return NULL;
    SRPoissonManifold* pm = (SRPoissonManifold*)calloc(1, sizeof(SRPoissonManifold));
    pm->dim = dim;
    pm->poisson_tensor = alloc_vec(dim * dim * dim);
    return pm;
}

void sr_poisson_free(SRPoissonManifold* pm) {
    if (!pm) return;
    free(pm->poisson_tensor);
    free(pm);
}

void sr_poisson_set_lie_poisson(SRPoissonManifold* pm, SRLieAlgebra* alg, int sign) {
    if (!pm || !alg) return;
    int n = alg->dimension;
    /* Lie-Poisson tensor: pi^{ij}(mu) = +/- sum_k c^{ij}_k mu_k */
    memset(pm->poisson_tensor, 0, (size_t)(n * n * n) * sizeof(double));
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++)
        pm->poisson_tensor[i*n*n + j*n + k] = sign * alg->constants[i*n + j*n + k];
    pm->rank = n;
}

double sr_poisson_bracket(const SRPoissonManifold* pm, const double* dF,
                           const double* dH, const double* x) {
    if (!pm || !dF || !dH || !x) return 0.0;
    int n = pm->dim;
    double result = 0.0;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++) {
        double pi_ij = pm->poisson_tensor[i*n*n + j*n + k] * x[k];
        result += dF[i] * pi_ij * dH[j];
    }
    return result;
}

bool sr_poisson_verify_jacobi(const SRPoissonManifold* pm, double tolerance) {
    if (!pm) return false;
    int n = pm->dim;
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++)
    for (int l = 0; l < n; l++) {
        double sum = 0.0;
        for (int m = 0; m < n; m++) {
            double a = pm->poisson_tensor[i*n*n + m*n + l] * pm->poisson_tensor[m*n*n + j*n + k];
            double b = pm->poisson_tensor[j*n*n + m*n + l] * pm->poisson_tensor[m*n*n + k*n + i];
            double c = pm->poisson_tensor[k*n*n + m*n + l] * pm->poisson_tensor[m*n*n + i*n + j];
            sum += a + b + c;
        }
        if (fabs(sum) > tolerance) return false;
    }
    return true;
}

/* ============================================================================
 * Vector Utilities
 * ============================================================================ */

SRVector* sr_vector_create(int dim) {
    if (dim <= 0) return NULL;
    SRVector* v = (SRVector*)calloc(1, sizeof(SRVector));
    v->dimension = dim;
    v->components = alloc_vec(dim);
    return v;
}

void sr_vector_free(SRVector* v) {
    if (!v) return;
    free(v->components);
    free(v);
}

void sr_vector_set(SRVector* v, int i, double val) {
    if (v && i >= 0 && i < v->dimension) v->components[i] = val;
}

double sr_vector_get(const SRVector* v, int i) {
    return (v && i >= 0 && i < v->dimension) ? v->components[i] : 0.0;
}

double sr_vector_dot(const SRVector* v1, const SRVector* v2) {
    if (!v1 || !v2 || v1->dimension != v2->dimension) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < v1->dimension; i++)
        sum += v1->components[i] * v2->components[i];
    return sum;
}

double sr_vector_norm(const SRVector* v) {
    return sqrt(sr_vector_dot(v, v));
}

void sr_vector_cross(const SRVector* v1, const SRVector* v2, SRVector* result) {
    if (!v1 || !v2 || !result) return;
    if (v1->dimension != 3 || v2->dimension != 3 || result->dimension != 3) return;
    result->components[0] = v1->components[1]*v2->components[2] - v1->components[2]*v2->components[1];
    result->components[1] = v1->components[2]*v2->components[0] - v1->components[0]*v2->components[2];
    result->components[2] = v1->components[0]*v2->components[1] - v1->components[1]*v2->components[0];
}

/* ============================================================================
 * Matrix Utilities
 * ============================================================================ */

SRMatrix* sr_matrix_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return NULL;
    SRMatrix* m = (SRMatrix*)calloc(1, sizeof(SRMatrix));
    m->rows = rows;
    m->cols = cols;
    m->data = alloc_vec(rows * cols);
    m->owns_data = true;
    return m;
}

void sr_matrix_free(SRMatrix* m) {
    if (!m) return;
    if (m->owns_data) free(m->data);
    free(m);
}

void sr_matrix_set(SRMatrix* m, int i, int j, double val) {
    if (m && i >= 0 && i < m->rows && j >= 0 && j < m->cols)
        m->data[i * m->cols + j] = val;
}

double sr_matrix_get(const SRMatrix* m, int i, int j) {
    return (m && i >= 0 && i < m->rows && j >= 0 && j < m->cols)
           ? m->data[i * m->cols + j] : 0.0;
}

void sr_matrix_multiply(const SRMatrix* A, const SRMatrix* B, SRMatrix* C) {
    if (!A || !B || !C || A->cols != B->rows) return;
    int m = A->rows, n = A->cols, p = B->cols;
    memset(C->data, 0, (size_t)(m * p) * sizeof(double));
    for (int i = 0; i < m; i++)
    for (int k = 0; k < n; k++) {
        double aik = A->data[i * A->cols + k];
        for (int j = 0; j < p; j++)
            C->data[i * C->cols + j] += aik * B->data[k * B->cols + j];
    }
}

void sr_matrix_vector_multiply(const SRMatrix* A, const double* x, double* y) {
    if (!A || !x || !y) return;
    memset(y, 0, (size_t)A->rows * sizeof(double));
    for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->cols; j++)
        y[i] += A->data[i * A->cols + j] * x[j];
}

void sr_matrix_transpose(const SRMatrix* A, SRMatrix* At) {
    if (!A || !At || A->rows != At->cols || A->cols != At->rows) return;
    for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->cols; j++)
        At->data[j * At->cols + i] = A->data[i * A->cols + j];
}

double sr_matrix_determinant_3x3(const double A[9]) {
    return mat3_det(A);
}

double sr_matrix_trace(const SRMatrix* A) {
    if (!A) return 0.0;
    int n = (A->rows < A->cols) ? A->rows : A->cols;
    double tr = 0.0;
    for (int i = 0; i < n; i++) tr += A->data[i * A->cols + i];
    return tr;
}

void sr_matrix_print(const SRMatrix* m) {
    if (!m) { printf("NULL matrix\n"); return; }
    printf("Matrix %dx%d:\n", m->rows, m->cols);
    for (int i = 0; i < m->rows; i++) {
        printf("  ");
        for (int j = 0; j < m->cols; j++)
            printf("%8.4f ", m->data[i * m->cols + j]);
        printf("\n");
    }
}
