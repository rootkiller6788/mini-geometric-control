#ifndef HAMILTONIAN_CONTROL_H
#define HAMILTONIAN_CONTROL_H

#include <stddef.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*=============================================================================
 * mini-hamiltonian-control -- Hamiltonian Control Theory
 *===========================================================================*/

/* L3: Mathematical Structures -- Core types */

typedef struct {
    int      dim2n;
    double **data;
} symplectic_matrix_t;

typedef struct {
    int      n;
    double  *q;
    double  *p;
} phase_point_t;

typedef struct {
    double (*energy)(const double *q, const double *p, int n, void *params);
    void   (*gradient)(const double *q, const double *p, int n,
                       double *grad_q, double *grad_p, void *params);
    void   (*hessian)(const double *q, const double *p, int n,
                      double **hess_qq, double **hess_qp,
                      double **hess_pq, double **hess_pp, void *params);
    void   *params;
    int     n_dof;
} hamiltonian_t;

/* L1: Hamiltonian vector field X_H = J grad(H) */
void hamiltonian_vector_field(const hamiltonian_t *H,
                              const double *q, const double *p,
                              double *dq_dt, double *dp_dt);

void hamiltonian_flow_rk4(const hamiltonian_t *H,
                          double *q, double *p,
                          double dt, int steps);

double hamiltonian_energy_error(const hamiltonian_t *H,
                                const double *q0, const double *p0,
                                const double *q1, const double *p1);

/* L1: Poisson bracket */
double poisson_bracket(const hamiltonian_t *H,
                       double (*F)(const double*, const double*, int, void*),
                       void *params_F,
                       double (*G)(const double*, const double*, int, void*),
                       void *params_G,
                       const double *q, const double *p, int n);

double poisson_jacobi_check(double (*F)(const double*, const double*, int, void*),
                            void *pf,
                            double (*G)(const double*, const double*, int, void*),
                            void *pg,
                            double (*K)(const double*, const double*, int, void*),
                            void *pk,
                            const double *q, const double *p, int n, double h);

/* L1: Canonical transformations */
typedef struct {
    void (*transform)(const double *q, const double *p, int n,
                      double *Q, double *P, void *params);
    void (*inverse)(const double *Q, const double *P, int n,
                    double *q, double *p, void *params);
    void *params;
    int   n_dof;
} canonical_transform_t;

double check_canonical_transform(const canonical_transform_t *T,
                                 const double *q, const double *p, double eps);

/* L3: Generating functions */
typedef enum {
    GEN_F1_Qq, GEN_F2_qP, GEN_F3_pQ, GEN_F4_pP
} gen_func_type_t;

typedef struct {
    gen_func_type_t  type;
    double        (*F)(const double *a, const double *b, int n, void *params);
    void          (*grad_a)(const double *a, const double *b, int n,
                            double *g, void *params);
    void          (*grad_b)(const double *a, const double *b, int n,
                            double *g, void *params);
    void          *params;
    int            n_dof;
} generating_function_t;

int generating_function_transform(const generating_function_t *gf,
                                  const double *old_a, const double *old_b,
                                  double *new_a, double *new_b);

/* L4: Liouville's Theorem */
double phase_volume_element(const phase_point_t *points, int n_points);

double liouville_volume_check(const hamiltonian_t *H,
                              const phase_point_t *cloud, int n_points,
                              double dt, int steps);

/* L5: Action-angle variables */
typedef struct {
    double  action;
    double  angle;
    double  frequency;
    int     converged;
} action_angle_t;

action_angle_t compute_action_angle_separable(const hamiltonian_t *H,
                                              double energy_level,
                                              double *q_turning_points);

/* L5: Normal mode analysis */
typedef struct {
    int      n_modes;
    double  *frequencies;
    double **eigenvectors;
    double   stability;
} normal_modes_t;

normal_modes_t hamiltonian_normal_modes(const hamiltonian_t *H,
                                        const double *q_eq, const double *p_eq);

/* Utility functions */
phase_point_t *phase_point_alloc(int n);
void phase_point_free(phase_point_t *pt);
symplectic_matrix_t *symplectic_matrix_alloc(int n);
void symplectic_matrix_free(symplectic_matrix_t *J);
void symplectic_matrix_standard(symplectic_matrix_t *J);
double **matrix_alloc(int rows, int cols);
void matrix_free(double **m, int rows);
void matvec_mul(double **A, const double *x, double *y, int m, int n);
double vec_norm2(const double *v, int n);
double vec_dot(const double *a, const double *b, int n);
int solve_spd_cholesky(double **A, const double *b, double *x, int n);
int sym_eigen(double **A, int n, double *eigenvalues, double **eigenvectors,
              int max_iter, double tol);

#endif /* HAMILTONIAN_CONTROL_H */
