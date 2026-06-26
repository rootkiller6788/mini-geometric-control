#ifndef POISSON_GEOMETRY_H
#define POISSON_GEOMETRY_H
#include "hamiltonian_control.h"

/*
 * Poisson geometry for Hamiltonian control systems.
 * Generalizes symplectic geometry to degenerate Poisson structures.
 * Reference: Marsden & Ratiu, Introduction to Mechanics and Symmetry
 *            Weinstein, The local structure of Poisson manifolds (1983)
 */

/* L1: Poisson manifold -- bracket defined by Poisson tensor B(z) */
typedef struct {
    int      dim;
    void   (*B_tensor)(const double *z, int dim, double **B_out, void *params);
    void   *params;
} poisson_manifold_t;

/* Evaluate Poisson bracket {F,G}(z) = sum_{ij} B_{ij}(z) dF/dz_i dG/dz_j */
double poisson_bracket_general(const poisson_manifold_t *pm,
                                double (*F)(const double*, int, void*),
                                void *pF,
                                double (*G)(const double*, int, void*),
                                void *pG,
                                const double *z);

/* L1: Lie algebra structure (for Lie-Poisson bracket on g*) */
typedef struct {
    int      dim;
    double  *c_ijk;  /* structure constants, flattened [dim][dim][dim] */
} lie_algebra_t;

/* Lie bracket [xi, eta]_k = sum_{i,j} c_{ij}^k xi_i eta_j */
void lie_bracket(const lie_algebra_t *g,
                 const double *xi, const double *eta,
                 double *bracket_out);

/* Lie-Poisson bracket: {F,G}(mu) = <mu, [dF, dG]> */
double lie_poisson_bracket(const lie_algebra_t *g,
                            double (*F)(const double*, int, void*),
                            void *pF,
                            double (*G)(const double*, int, void*),
                            void *pG,
                            const double *mu);

/* L1: Casimir function -- {C, F} = 0 for all F */
int is_casimir_function(const poisson_manifold_t *pm,
                         double (*C)(const double*, int, void*),
                         void *pC,
                         int n_test_functions,
                         double (**test_fns)(const double*, int, void*),
                         void **test_params,
                         const double *z,
                         double *max_bracket_norm);

/* L2: Symplectic leaves -- dimension = rank(B(z)) */
int symplectic_leaf_dimension(const poisson_manifold_t *pm, const double *z);
int symplectic_leaf_basis(const poisson_manifold_t *pm, const double *z,
                           double **leaf_basis, int *leaf_dim);

/* L3: SO(3) Lie algebra and rigid body dynamics */
typedef struct {
    double  I1, I2, I3;  /* principal moments of inertia */
} rigid_body_params_t;

lie_algebra_t *so3_lie_algebra_create(void);
double rigid_body_hamiltonian(const double *Pi, int dim, void *params);
void rigid_body_gradient(const double *Pi, int dim, double *grad, void *params);
void rigid_body_euler_equations(const rigid_body_params_t *rb,
                                 const double *Pi, double *dPi_dt);

/* L4: Lie-Poisson simulation of reduced dynamics */
int lie_poisson_simulate(const lie_algebra_t *g,
                          double (*h)(const double*, int, void*),
                          void (*grad_h)(const double*, int, double*, void*),
                          void *params,
                          double *mu, double dt, int steps,
                          double **trajectory_out);

/* L5: Energy-Casimir stability method */
typedef struct {
    double  *phi_coeffs;
    int      phi_degree;
    int      n_casimirs;
} energy_casimir_method_t;

int energy_casimir_lyapunov(const lie_algebra_t *g,
                             double (*h)(const double*, int, void*),
                             void (*grad_h)(const double*, int, double*, void*),
                             double (**casimirs)(const double*, int, void*),
                             void (**grad_casimirs)(const double*, int, double*, void*),
                             int n_casimirs,
                             const double *mu_eq,
                             energy_casimir_method_t *ecm,
                             double *V_min_eigenvalue);

/* L6: Rigid body control via momentum wheels */
typedef struct {
    rigid_body_params_t rb;
    double K_p, K_d;
    int    target_axis;
} rigid_body_control_t;

int rigid_body_stabilize(const rigid_body_control_t *ctrl,
                          const double *Pi, const double *Pi_desired,
                          double *u_out);

/* L7: Heavy top (spinning top with gravity) */
typedef struct {
    rigid_body_params_t rb;
    double mass, g, l;
    double Gamma0[3];
} heavy_top_params_t;

void heavy_top_equations(const heavy_top_params_t *ht,
                          const double *Pi, const double *Gamma,
                          double *dPi_dt, double *dGamma_dt);

/* L8: Poisson reduction for controlled systems */
typedef struct {
    lie_algebra_t      *g;
    double (*h_reduced)(const double *mu, const double *u,
                        int dim, int ctrl_dim, void *params);
    void   (*grad_mu_h)(const double *mu, const double *u,
                        int dim, int ctrl_dim, double *grad, void *params);
    void   (*grad_u_h) (const double *mu, const double *u,
                        int dim, int ctrl_dim, double *grad, void *params);
    int     ctrl_dim;
    void   *params;
} reduced_control_system_t;

int poisson_reduction_optimal_control(
    const reduced_control_system_t *rcs,
    const double *mu0, const double *mu_target,
    double tf, int n_steps,
    double **mu_traj, double **u_traj);

void lie_algebra_free(lie_algebra_t *g);

#endif /* POISSON_GEOMETRY_H */
