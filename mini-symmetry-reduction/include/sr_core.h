#ifndef SR_CORE_H
#define SR_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ============================================================================
 * Symmetry Reduction in Geometric Control Theory - Core Types
 *
 * Foundational types for reduction of mechanical and control systems
 * admitting Lie group symmetries.
 *
 * References:
 *   Marsden & Weinstein (1974) "Reduction of symplectic manifolds with symmetry"
 *   Marsden & Ratiu (1999) "Introduction to Mechanics and Symmetry"
 *   Arnold (1989) "Mathematical Methods of Classical Mechanics"
 * ============================================================================ */

/* Group action type */
typedef enum {
    SR_ACTION_FREE        = 0,
    SR_ACTION_PROPER      = 1,
    SR_ACTION_FREE_PROPER = 2,
    SR_ACTION_EFFECTIVE   = 3,
    SR_ACTION_TRANSITIVE  = 4,
    SR_ACTION_HAMILTONIAN = 5
} SRActionType;

/* Reduction method */
typedef enum {
    SR_REDUCTION_MARSDEN_WEINSTEIN = 0,
    SR_REDUCTION_LIE_POISSON       = 1,
    SR_REDUCTION_EULER_POINCARE    = 2,
    SR_REDUCTION_LAGRANGE_POINCARE = 3,
    SR_REDUCTION_COTANGENT_BUNDLE  = 4,
    SR_REDUCTION_SEMIDIRECT        = 5,
    SR_REDUCTION_STAGES            = 6,
    SR_REDUCTION_POISSON           = 7,
    SR_REDUCTION_SINGULAR          = 8
} SRReductionMethod;

/* Symmetry type */
typedef enum {
    SR_SYMMETRY_CONTINUOUS    = 0,
    SR_SYMMETRY_DISCRETE      = 1,
    SR_SYMMETRY_TORAL         = 2,
    SR_SYMMETRY_NONABELIAN    = 3,
    SR_SYMMETRY_SEMIDIRECT    = 4,
    SR_SYMMETRY_GAUGE         = 5,
    SR_SYMMETRY_CONFORMAL     = 6,
    SR_SYMMETRY_SUPERSYMMETRY = 7
} SRSymmetryType;

/* Structure constant storage */
typedef enum {
    SR_STRUCT_FULL_TENSOR = 0,
    SR_STRUCT_SPARSE_CSR  = 1,
    SR_STRUCT_AD_REP      = 2
} SRStructStorage;

/* Momentum map type */
typedef enum {
    SR_MOMENTUM_ANGULAR  = 0,
    SR_MOMENTUM_LINEAR   = 1,
    SR_MOMENTUM_KINETIC  = 2,
    SR_MOMENTUM_QUOTIENT = 3,
    SR_MOMENTUM_IMPULSE  = 4
} SRMomentumType;

/* Basic Numeric Types */
typedef struct { double value; double real; double imag; bool is_complex; } SRComplex;

typedef struct { double* components; int dimension; } SRVector;

typedef struct { double* data; int rows; int cols; bool owns_data; } SRMatrix;

typedef struct { double* data; int dim2; int dim; } SRSymplecticMatrix;

/* Lie algebra g: [e_i, e_j] = sum_k c_{ij}^k e_k */
typedef struct {
    double* constants;
    int dimension;
    char** basis_names;
    SRStructStorage format;
    bool is_nilpotent, is_solvable, is_semisimple, is_simple, is_abelian;
    bool jacobi_verified;
} SRLieAlgebra;

/* Lie group G */
typedef struct {
    int dimension, matrix_size;
    SRLieAlgebra* algebra;
    char* name;
    SRSymmetryType type;
    bool is_compact, is_connected, is_simply_connected;
    double* identity;
    int* multiplication_table;
} SRLieGroup;

/* Group element */
typedef struct {
    SRLieGroup* group;
    double* coords, *matrix;
    bool is_identity;
} SRGroupElement;

/* Lie algebra element */
typedef struct {
    SRLieAlgebra* algebra;
    double* coords, *matrix;
} SRLieAlgebraElement;

/* Adjoint action Ad_g : g -> g */
typedef struct {
    SRLieGroup* group;
    double* matrix, *cached_inverse;
    int dim;
} SRAdjointAction;

/* Coadjoint action Ad*_{g^{-1}} : g* -> g* */
typedef struct {
    SRLieGroup* group;
    double* matrix;
    int dim;
} SRCoadjointAction;

/* Symplectic manifold (P, omega) */
typedef struct {
    int dim;
    double* omega;
    bool is_closed, is_nondeg;
    char** coordinate_names;
} SRSymplecticManifold;

/* Poisson manifold */
typedef struct {
    int dim, rank;
    double* poisson_tensor;
    bool jacobi_verified, is_symplectic;
    SRSymplecticManifold* symplectic_leaf;
} SRPoissonManifold;

/* Lie-Poisson bracket {F,H}(mu) = +/- <mu, [dF, dH]> */
typedef struct {
    SRLieAlgebra* algebra;
    int sign, dim;
    double* structure_constants;
} SRLiePoissonBracket;

/* Momentum map J: P -> g* */
typedef struct {
    SRLieGroup* group;
    SRLieAlgebra* algebra;
    int phase_dim, algebra_dim;
    double* value, *jacobian;
    SRMomentumType type;
    bool is_equivariant, is_conserved;
    double* infinitesimal_generators;
} SRMomentumMap;

/* Marsden-Weinstein reduced space P_mu = J^{-1}(mu) / G_mu */
typedef struct {
    double mu_value;
    SRLieGroup* g_mu;
    int dim_original, dim_level_set, dim_reduced;
    double* reduced_omega, *reduced_hamiltonian, *reduced_coords;
    bool is_regular;
    char* description;
} SRMarsdenWeinsteinSpace;

/* Euler-Poincare reduced system */
typedef struct {
    SRLieAlgebra* algebra;
    SRLieGroup* group;
    int sign, n_config;
    double* configuration, *momentum, *inertia_matrix;
    SRReductionMethod method;
} SREulerPoincareSystem;

/* Lagrange-Poincare reduced system */
typedef struct {
    SRLieGroup* group;
    SRLieAlgebra* algebra;
    double* configuration, *momentum, *group_element;
    double* inertia_matrix, *connection_coeffs;
    int algebra_dim, group_param_dim;
} SRLagrangePoincareSystem;

/* Semidirect product G = S ⋉ V */
typedef struct {
    SRLieGroup* s_subgroup, *v_subgroup;
    int s_dim, v_dim;
    double* rho_matrix, *semidirect_momentum, *lie_algebra_cocycle;
} SRSemidirectProduct;

/* G-invariant Hamiltonian H: P -> R */
typedef struct {
    double (*value)(const double* state, int state_dim, void* params);
    void (*gradient)(const double* state, int state_dim, void* params, double* grad_out);
    void (*hessian)(const double* state, int state_dim, void* params, double* hess_out);
    void* params;
    SRLieGroup* symmetry_group;
    bool is_invariant;
    int state_dim;
} SRInvariantHamiltonian;

/* Reduced control system */
typedef struct {
    int state_dim, control_dim, algebra_dim;
    double* state, *drift, *control_fields;
    SRLieGroup* configuration_group;
    SRLieAlgebraElement* body_velocities;
    double* shape_momentum, *connection_form;
    bool is_driftless, is_reduced;
} SRReducedControlSystem;

/* ============ Core Function Declarations ============ */

/* Lie algebra */
SRLieAlgebra* sr_algebra_create(int dimension, const char* name);
void sr_algebra_free(SRLieAlgebra* alg);
void sr_algebra_set_struct_const(SRLieAlgebra* alg, int i, int j, int k, double val);
double sr_algebra_get_struct_const(const SRLieAlgebra* alg, int i, int j, int k);
bool sr_algebra_verify_jacobi(const SRLieAlgebra* alg, double tolerance);
double sr_algebra_killing_form(const SRLieAlgebra* alg, int i, int j);
double sr_algebra_commutator_component(const SRLieAlgebra* alg, int i, int j, int k);
void sr_algebra_bracket(const SRLieAlgebra* alg, const double* X, const double* Y, double* bracket_out);
void sr_algebra_print(const SRLieAlgebra* alg);

/* Lie group */
SRLieGroup* sr_group_create(int dimension, int matrix_size, const char* name, SRSymmetryType type);
void sr_group_free(SRLieGroup* group);
void sr_group_set_compact(SRLieGroup* group, bool compact);
void sr_group_set_connected(SRLieGroup* group, bool connected);
SRLieAlgebra* sr_group_get_algebra(SRLieGroup* group);
void sr_group_set_identity_matrix(SRLieGroup* group, const double* identity);
void sr_group_print(const SRLieGroup* group);

/* Group element operations */
SRGroupElement* sr_element_create(SRLieGroup* group);
void sr_element_free(SRGroupElement* el);
void sr_element_set_coords(SRGroupElement* el, const double* coords);
void sr_element_set_identity(SRGroupElement* el);
void sr_element_multiply(const SRGroupElement* g1, const SRGroupElement* g2, SRGroupElement* result);
void sr_element_inverse(const SRGroupElement* el, SRGroupElement* inverse);
void sr_element_action_on_vector(const SRGroupElement* el, const double* vec, double* result);
void sr_element_print(const SRGroupElement* el);

/* Lie algebra element operations */
SRLieAlgebraElement* sr_alg_element_create(SRLieAlgebra* alg);
void sr_alg_element_free(SRLieAlgebraElement* el);
void sr_alg_element_set_coords(SRLieAlgebraElement* el, const double* coords);
double sr_alg_element_norm(const SRLieAlgebraElement* el);

/* Adjoint and coadjoint */
SRAdjointAction* sr_adjoint_create(SRLieGroup* group);
void sr_adjoint_free(SRAdjointAction* adj);
void sr_adjoint_compute(const SRAdjointAction* adj, const SRGroupElement* g, double* Ad_g_matrix);
void sr_adjoint_action_on_algebra(const SRAdjointAction* adj, const SRGroupElement* g, const double* xi, double* result);

SRCoadjointAction* sr_coadjoint_create(SRLieGroup* group);
void sr_coadjoint_free(SRCoadjointAction* coadj);
void sr_coadjoint_compute(const SRCoadjointAction* coadj, const SRGroupElement* g, double* Ad_star_matrix);
void sr_coadjoint_action_on_dual(const SRCoadjointAction* coadj, const SRGroupElement* g, const double* mu, double* result);

/* Exponential map */
void sr_exp_map(const SRLieAlgebra* alg, const double* xi, double* g_matrix, int matrix_size);
void sr_exp_map_so3(const double* omega, double* R_matrix);
void sr_exp_map_se3(const double* twist, double* T_matrix);
void sr_bch_formula_2nd(const SRLieAlgebra* alg, const double* X, const double* Y, double* Z);

/* Symplectic geometry */
SRSymplecticManifold* sr_symplectic_create(int n_dof);
void sr_symplectic_free(SRSymplecticManifold* sp);
void sr_symplectic_set_canonical(SRSymplecticManifold* sp);
bool sr_symplectic_is_nondeg(const SRSymplecticManifold* sp);
void sr_symplectic_inverse(const SRSymplecticManifold* sp, double* omega_inv);
double sr_symplectic_volume(const SRSymplecticManifold* sp);
void sr_symplectic_print(const SRSymplecticManifold* sp);

/* Poisson geometry */
SRPoissonManifold* sr_poisson_create(int dim);
void sr_poisson_free(SRPoissonManifold* pm);
void sr_poisson_set_lie_poisson(SRPoissonManifold* pm, SRLieAlgebra* alg, int sign);
double sr_poisson_bracket(const SRPoissonManifold* pm, const double* dF, const double* dH, const double* x);
bool sr_poisson_verify_jacobi(const SRPoissonManifold* pm, double tolerance);

/* Vector utilities */
SRVector* sr_vector_create(int dim);
void sr_vector_free(SRVector* v);
void sr_vector_set(SRVector* v, int i, double val);
double sr_vector_get(const SRVector* v, int i);
double sr_vector_dot(const SRVector* v1, const SRVector* v2);
double sr_vector_norm(const SRVector* v);
void sr_vector_cross(const SRVector* v1, const SRVector* v2, SRVector* result);

/* Matrix utilities */
SRMatrix* sr_matrix_create(int rows, int cols);
void sr_matrix_free(SRMatrix* m);
void sr_matrix_set(SRMatrix* m, int i, int j, double val);
double sr_matrix_get(const SRMatrix* m, int i, int j);
void sr_matrix_multiply(const SRMatrix* A, const SRMatrix* B, SRMatrix* C);
void sr_matrix_vector_multiply(const SRMatrix* A, const double* x, double* y);
void sr_matrix_transpose(const SRMatrix* A, SRMatrix* At);
double sr_matrix_determinant_3x3(const double A[9]);
double sr_matrix_trace(const SRMatrix* A);
void sr_matrix_print(const SRMatrix* m);

/* Special Lie algebra factories */
SRLieAlgebra* sr_algebra_create_so3(void);
SRLieAlgebra* sr_algebra_create_so2(void);
SRLieAlgebra* sr_algebra_create_su2(void);
SRLieAlgebra* sr_algebra_create_se3(void);
SRLieAlgebra* sr_algebra_create_heisenberg(void);
SRLieAlgebra* sr_algebra_create_abelian(int dim);

#endif /* SR_CORE_H */
