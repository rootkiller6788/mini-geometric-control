#ifndef FEEDBACK_LINEARIZATION_H
#define FEEDBACK_LINEARIZATION_H

#include <stdbool.h>
#include <stddef.h>
#include "nonlinear_system.h"
#include "lie_derivative.h"
#include "lie_bracket.h"
#include "coordinate_transform.h"

/* ============================================================================
 * Feedback Linearization ? Geometric Control Theory
 *
 * Based on the foundational works of:
 *   Alberto Isidori    (Nonlinear Control Systems, 3rd ed, 1995)
 *   Henk Nijmeijer & Arjan van der Schaft (Nonlinear Dynamical Control Systems, 1990)
 *   Hassan K. Khalil   (Nonlinear Systems, 3rd ed, 2002)
 *   Jean-Jacques Slotine & Weiping Li (Applied Nonlinear Control, 1991)
 *   Roger W. Brockett  (Feedback Invariants for Nonlinear Systems, 1978)
 *   Bronislaw Jakubczyk & Witold Respondek (Linearization of Nonlinear Systems, 1980)
 *   Christopher I. Byrnes & Alberto Isidori (Zero Dynamics, 1984-1991)
 *
 * Core idea: Transform a nonlinear system into an equivalent linear system
 * via nonlinear state feedback and coordinate transformation.
 *
 * Nonlinear system:  x' = f(x) + g(x) u
 *                    y  = h(x)
 *
 * Output y is linearized to rho-th integrator when relative degree rho is well-defined.
 * ============================================================================ */

/* --- Relative Degree Computation Results --- */

typedef struct {
    int rho;                    /* Relative degree: smallest r where L_g L_f^{r-1} h(x) != 0 */
    bool is_well_defined;       /* rho exists and is finite */
    Matrix *decoupling_matrix;  /* Decoupling matrix E(x) for MIMO case */
    double *Lf_powers_h;        /* Array: h, L_f h, ..., L_f^{rho-1} h evaluated at x */
    double Lg_Lf_power;         /* L_g L_f^{rho-1} h(x) -- the high-frequency gain */
} RelativeDegreeResult;

/* --- Input-Output Linearization --- */

typedef struct {
    double alpha;               /* alpha(x) = -L_f^rho h(x) / L_g L_f^{rho-1} h(x) */
    double beta;                /* beta(x) = 1 / L_g L_f^{rho-1} h(x) */
    bool is_feedback_linearizable;
    NonlinearSystem *normalized; /* System in normal form coordinates */
} IOLinearizationResult;

/* --- Input-State Linearization --- */

typedef struct {
    Diffeomorphism *phi;        /* Coordinate transformation z = Phi(x) */
    VectorField *f_transformed; /* f in new coordinates */
    VectorField *g_transformed; /* g in new coordinates */
    bool is_exactly_linearizable;
    int largest_controllable_index; /* dim of largest involutive distribution */
} ISLinearizationResult;

/* --- Zero Dynamics --- */

typedef enum {
    ZERO_DYN_STABLE = 0,        /* Minimum phase: zero dynamics asymptotically stable */
    ZERO_DYN_UNSTABLE = 1,      /* Non-minimum phase */
    ZERO_DYN_MARGINAL = 2,      /* Critical case: marginally stable */
    ZERO_DYN_NONE = 3           /* No zero dynamics (full relative degree = n) */
} ZeroDynamicsStability;

typedef struct {
    int dim_zero_dynamics;      /* n - rho (dimension of zero dynamics) */
    ZeroDynamicsStability stability;
    double *zero_eigenvalues;   /* Eigenvalues of zero-dynamics linearization */
    double lyapunov_exponent;   /* Largest Lyapunov exponent of zero dynamics */
    Vector *zero_manifold;      /* Manifold Z*: {x: h=0, L_f h=0, ..., L_f^{rho-1} h=0} */
} ZeroDynamicsResult;

/* --- Normal Form (Byrnes-Isidori) --- */

typedef struct {
    Vector *xi;                 /* External state xi = (h, L_f h, ..., L_f^{rho-1} h)^T */
    Vector *eta;                /* Internal state eta (zero dynamics coordinates) */
    Diffeomorphism *T;          /* Full transformation T(x) = (xi, eta) */
    VectorField *f_xi;          /* xi' = A*xi + B*[alpha(xi,eta) + beta(xi,eta) u] */
    VectorField *f_eta;         /* eta' = q(xi, eta) -- zero dynamics when xi = 0 */
    VectorField *g_normal;      /* g in normal form */
    double (*h_normal)(const Vector*); /* h(xi) = xi_1 in normal form */
} NormalFormResult;

/* --- MIMO Feedback Linearization --- */

typedef struct {
    int num_inputs;
    int num_outputs;
    int *relative_degrees;      /* Vector relative degree (rho_1, ..., rho_m) */
    bool has_vector_relative_degree;
    bool decoupling_matrix_invertible;
    Matrix *decoupling_matrix;  /* E(x) = [L_{g_j} L_f^{rho_i-1} h_i(x)] */
    Matrix *A;                  /* Decoupled linear system matrix */
    Matrix *B;                  /* Decoupled linear system input matrix */
} MIMOLinearizationResult;

/* ============================================================================
 * Core API -- Feedback Linearization
 * ============================================================================ */

/**
 * compute_relative_degree -- Determine relative degree of SISO nonlinear system.
 *
 * For system x' = f(x) + g(x) u, y = h(x):
 *   rho = min{r : L_g L_f^{r-1} h(x) != 0}
 *
 * Theorem (Isidori, 1995 Sec 4.2): rho is well-defined on an open subset U of R^n
 * iff L_g L_f^k h(x) = 0 for k = 0, ..., rho-2 and L_g L_f^{rho-1} h(x) != 0 on U.
 *
 * @param sys       nonlinear control-affine system
 * @param x         state vector at which to evaluate
 * @param max_order maximum order to check (typically n)
 * @return          relative degree result
 */
RelativeDegreeResult compute_relative_degree(const NonlinearSystem *sys,
                                              const Vector *x,
                                              int max_order);

/**
 * input_output_linearize -- Perform input-output feedback linearization.
 *
 * Compute feedback law u = alpha(x) + beta(x) v that linearizes the I/O map
 * to a chain of rho integrators: y^{(rho)} = v.
 *
 * alpha(x) = -L_f^rho h(x) / L_g L_f^{rho-1} h(x)
 * beta(x)  =  1 / L_g L_f^{rho-1} h(x)
 *
 * Theorem (Isidori, 1995 Sec 4.2): The feedback law exists iff rho is
 * well-defined and beta(x) != 0 on the domain of interest.
 *
 * @param sys  nonlinear system
 * @param x    evaluation point
 * @return     IO linearization result
 */
IOLinearizationResult input_output_linearize(const NonlinearSystem *sys,
                                              const Vector *x);

/**
 * input_state_linearize -- Check full input-state exact linearizability.
 *
 * Theorem (Jakubczyk-Respondek, 1980): A single-input system x'=f(x)+g(x)u
 * is exactly input-state linearizable iff:
 *   1. Span{g, ad_f g, ..., ad_f^{n-1} g} = R^n (controllability)
 *   2. Distribution Delta_k = span{g, ad_f g, ..., ad_f^k g} is involutive
 *      for all k = 0, ..., n-2.
 *
 * This is the geometric condition using Lie brackets and involutivity.
 *
 * @param sys  nonlinear system
 * @return     IS linearization result with coordinate transformation
 */
ISLinearizationResult input_state_linearize(const NonlinearSystem *sys);

/**
 * analyze_zero_dynamics -- Analyze internal/zero dynamics of IO-linearized system.
 *
 * Zero dynamics: eta' = q(0, eta) -- residual dynamics when output forced to zero.
 *
 * Theorem (Byrnes-Isidori, 1991): The closed-loop system is internally stable
 * iff the zero dynamics are asymptotically stable at eta = 0.
 *
 * @param sys  nonlinear system with well-defined relative degree
 * @param x    operating point
 * @return     zero dynamics analysis result
 */
ZeroDynamicsResult analyze_zero_dynamics(const NonlinearSystem *sys,
                                          const Vector *x);

/**
 * compute_normal_form -- Transform system into Byrnes-Isidori normal form.
 *
 * Normal form:
 *   xi'_1 = xi_2
 *   ...
 *   xi'_rho = b(xi,eta) + a(xi,eta) u   (where a != 0)
 *   eta'    = q(xi, eta)                (zero dynamics when xi = 0)
 *   y       = xi_1
 *
 * @param sys  nonlinear system with well-defined relative degree
 * @param x    evaluation point
 * @return     normal form representation
 */
NormalFormResult compute_normal_form(const NonlinearSystem *sys,
                                      const Vector *x);

/**
 * mimo_linearize -- MIMO feedback linearization with vector relative degree.
 *
 * For system with m inputs and m outputs:
 *   y_i^{(rho_i)} = L_f^{rho_i} h_i + sum_j L_{g_j} L_f^{rho_i-1} h_i u_j
 *
 * Theorem (Isidori, 1995 Sec 5.2): MIMO IO linearization is possible iff
 * the decoupling matrix E(x) = [L_{g_j} L_f^{rho_i-1} h_i(x)] is nonsingular.
 *
 * @param sys  MIMO nonlinear system
 * @param x    evaluation point
 * @return     MIMO linearization result
 */
MIMOLinearizationResult mimo_linearize(const NonlinearSystem *sys,
                                        const Vector *x);

/**
 * feedback_law_evaluate -- Compute linearizing feedback u = alpha(x) + beta(x) v.
 *
 * Given a virtual input v (the new input to the linearized system),
 * compute the actual physical input u that achieves linearization.
 *
 * @param io_result  IO linearization result
 * @param v          virtual input
 * @return           physical input u
 */
double feedback_law_evaluate(const IOLinearizationResult *io_result, double v);

/**
 * check_minimum_phase -- Test if system is minimum phase (stable zero dynamics).
 *
 * A system is minimum phase if its zero dynamics are asymptotically stable.
 * For linear systems, all zeros are in the open left-half plane.
 * For nonlinear systems, zero dynamics eta' = q(0,eta) must have eta=0 as
 * asymptotically stable equilibrium.
 *
 * @param zd_result  zero dynamics analysis result
 * @return           true if minimum phase
 */
bool check_minimum_phase(const ZeroDynamicsResult *zd_result);

/**
 * full_state_feedback_linearizable -- Complete check for input-state
 * linearization.
 *
 * Tests both conditions from the Jakubczyk-Respondek theorem:
 *   C1: Controllability rank condition
 *   C2: Involutivity condition
 *
 * @param sys  nonlinear control-affine system
 * @return     true if exactly linearizable by state feedback and diffeomorphism
 */
bool full_state_feedback_linearizable(const NonlinearSystem *sys);

/**
 * compute_decoupling_matrix -- Compute E(x) for MIMO system.
 *
 * E_{ij}(x) = L_{g_j} L_f^{rho_i-1} h_i(x)
 *
 * @param sys              MIMO system
 * @param x                evaluation point
 * @param relative_degrees vector relative degree
 * @return                 decoupling matrix (m x m)
 */
Matrix *compute_decoupling_matrix(const NonlinearSystem *sys,
                                   const Vector *x,
                                   const int *relative_degrees);

/**
 * linearizing_coordinates -- Compute linearizing coordinates for input-state
 * linearization.
 *
 * For a single-input system satisfying the Jakubczyk-Respondek conditions,
 * the linearizing coordinate z = Phi(x) = (T_1(x), ..., T_n(x)) where:
 *   L_g T_1 = 0, ..., L_g T_{n-1} = 0, L_g T_n != 0
 *   T_{k+1} = L_f T_k
 *
 * @param sys  exactly linearizable system
 * @return     diffeomorphism Phi
 */
Diffeomorphism *linearizing_coordinates(const NonlinearSystem *sys);

/**
 * free_output_feedback_linearize -- Linearize with dynamic extension when
 * static feedback fails (no well-defined relative degree).
 *
 * Add integrators before some inputs until relative degree becomes well-defined.
 * This implements the dynamic extension algorithm.
 *
 * @param sys         nonlinear system
 * @param max_extend  maximum number of integrators to add
 * @return            dynamically extended linearization result
 */
IOLinearizationResult free_output_feedback_linearize(const NonlinearSystem *sys,
                                                      int max_extend);

#endif /* FEEDBACK_LINEARIZATION_H */
