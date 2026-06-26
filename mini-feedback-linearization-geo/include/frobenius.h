#ifndef FROBENIUS_H
#define FROBENIUS_H

#include "nonlinear_system.h"
#include "lie_bracket.h"
#include "coordinate_transform.h"

/* ============================================================================
 * Frobenius Theorem ? Involutive Distributions and Complete Integrability
 *
 * A distribution Delta on R^n assigns to each point x a subspace Delta(x)
 * of the tangent space T_x R^n.
 *
 * Delta is INVOLUTIVE if for any two vector fields X, Y in Delta,
 * their Lie bracket [X, Y] is also in Delta.
 *
 * Frobenius Theorem: A nonsingular distribution Delta is completely
 * integrable if and only if it is involutive.
 *
 * Complete integrability means there exists a foliation of R^n by
 * integral submanifolds of dimension = dim(Delta), such that at each
 * point, the tangent space to the leaf equals Delta(x).
 *
 * In control theory, this theorem underpins:
 *   - Input-state feedback linearization (involutivity condition)
 *   - Decomposition into observable/unobservable subspaces
 *   - Exact linearization by dynamic feedback
 *
 * References:
 *   Isidori (1995) Nonlinear Control Systems, Sec 1.5-1.8
 *   Boothby (1986) An Introduction to Differentiable Manifolds, Ch. VI
 *   Warner (1983) Foundations of Differentiable Manifolds and Lie Groups
 *   Frobenius (1877) "Uber das Pfaffsche Problem"
 * ============================================================================ */

/* --- Distribution Type --- */

typedef struct {
    int n;                        /* Ambient space dimension */
    int d;                        /* Distribution dimension (constant rank assumed) */
    VectorField **basis;          /* d basis vector fields spanning Delta at each x */
    int num_basis;
} Distribution;

/* --- Frobenius Result --- */

typedef struct {
    bool is_involutive;           /* Does Delta satisfy the involutivity condition? */
    bool is_integrable;           /* Is Delta completely integrable? (same, by Frobenius) */
    int max_involutive_dim;       /* Largest dim such that sub-distribution is involutive */
    Matrix *integral_manifold;    /* Local representation of integral submanifold */
    Diffeomorphism *foliation;    /* Coordinate chart making Delta = span{d/dz_1, ..., d/dz_d} */
} FrobeniusResult;

/* ============================================================================
 * Distribution API
 * ============================================================================ */

/**
 * distribution_create -- Create a distribution from basis vector fields.
 *
 * @param n           ambient space dimension
 * @param num_basis   number of basis vector fields
 * @param basis       array of vector field pointers (ownership not transferred)
 * @return            allocated distribution
 */
Distribution *distribution_create(int n, int num_basis, VectorField **basis);

/**
 * distribution_add_basis -- Add a vector field to the distribution's spanning set.
 *
 * @param dist  distribution
 * @param vf    vector field to add
 */
void distribution_add_basis(Distribution *dist, VectorField *vf);

/**
 * distribution_eval -- Evaluate all basis vector fields at x.
 *
 * Returns a matrix where column j is basis field j evaluated at x.
 * The distribution Delta(x) is the column space of this matrix.
 *
 * @param dist  distribution
 * @param x     evaluation point
 * @return      n x d matrix
 */
Matrix *distribution_eval(const Distribution *dist, const Vector *x);

/**
 * distribution_dim_at -- Compute pointwise dimension of Delta(x).
 *
 * Uses numerical rank of the evaluation matrix at x.
 *
 * @param dist  distribution
 * @param x     evaluation point
 * @param tol   rank tolerance
 * @return      dim Delta(x)
 */
int distribution_dim_at(const Distribution *dist, const Vector *x, double tol);

/**
 * distribution_free -- Free distribution (does not free basis vector fields).
 *
 * @param dist  distribution to free
 */
void distribution_free(Distribution *dist);

/* ============================================================================
 * Involutivity API
 * ============================================================================ */

/**
 * check_involutivity -- Test if a distribution is involutive.
 *
 * For each pair of basis vector fields X_i, X_j, compute their Lie bracket
 * [X_i, X_j] and check if it lies in the span of Delta(x).
 *
 * Theorem (Frobenius): Involutivity is necessary and sufficient for
 * complete integrability of a non-singular distribution.
 *
 * @param dist  distribution
 * @param x     evaluation point
 * @param tol   numerical tolerance for rank and bracket checks
 * @return      true if involutive at x
 */
bool check_involutivity(const Distribution *dist, const Vector *x, double tol);

/**
 * check_involutivity_global -- Test involutivity on a set of test points.
 *
 * @param dist        distribution
 * @param test_points array of evaluation points
 * @param num_points  number of points
 * @param tol         tolerance
 * @return            true if involutive at all test points
 */
bool check_involutivity_global(const Distribution *dist,
                                Vector **test_points,
                                int num_points,
                                double tol);

/**
 * frobenius_analysis -- Full Frobenius analysis of a distribution.
 *
 *   (1) Check involutivity at point x
 *   (2) If involutive, compute local foliation coordinates
 *   (3) Determine integral manifold structure
 *
 * @param dist  distribution
 * @param x     analysis point
 * @return      Frobenius result
 */
FrobeniusResult frobenius_analysis(const Distribution *dist, const Vector *x);

/**
 * compute_annihilator -- Compute the codistribution Omega = Delta^{perp}.
 *
 * Omega(x) = {w in (R^n)^* : w ? v = 0 for all v in Delta(x)}
 *
 * The annihilator of a distribution is a codistribution (set of covectors)
 * whose dimension is n - dim(Delta).
 *
 * Frobenius (dual formulation): Delta is completely integrable iff
 * d w ? Ideal(Omega) for all w ? Omega.
 *
 * @param dist  distribution
 * @param x     evaluation point
 * @return      n x n-d matrix whose columns span the annihilator at x
 */
Matrix *compute_annihilator(const Distribution *dist, const Vector *x);

/**
 * check_integrability_pfaffian -- Dual Frobenius test using Pfaffian systems.
 *
 * Given a set of 1-forms ?_1, ..., ?_{n-d} that annihilate Delta,
 * check: d?_k ? ?_1 ? ... ? ?_{n-d} = 0 for all k.
 *
 * This is the classical Pfaffian formulation of the Frobenius theorem.
 *
 * @param dist  distribution
 * @param x     evaluation point
 * @param tol   tolerance
 * @return      true if the Pfaffian system is completely integrable
 */
bool check_integrability_pfaffian(const Distribution *dist,
                                   const Vector *x,
                                   double tol);

/**
 * compute_max_involutive_subdist -- Find the largest involutive sub-distribution.
 *
 * For input-state linearization, we need the chain of distributions:
 *   Delta_0 ? Delta_1 ? ... ? Delta_{n-1}
 * where each Delta_k must be involutive and have dimension k+1.
 *
 * If this chain exists, the system is exactly linearizable.
 *
 * @param f     drift vector field
 * @param g     control vector field
 * @param x     evaluation point
 * @return      maximum k such that Delta_k is involutive
 */
int compute_max_involutive_subdist(const VectorField *f,
                                    const VectorField *g,
                                    const Vector *x);

/**
 * build_linearizing_distribution -- Build the chain of distributions for
 * input-state linearization.
 *
 * Delta_0 = span{g}
 * Delta_1 = span{g, ad_f g}
 * ...
 * Delta_{n-1} = span{g, ad_f g, ..., ad_f^{n-1} g}
 *
 * Returns array of n distributions. Caller frees each element.
 *
 * @param f  drift vector field
 * @param g  control vector field
 * @param n  state dimension
 * @return   array of n Distribution pointers
 */
Distribution **build_linearizing_distribution(const VectorField *f,
                                               const VectorField *g,
                                               int n);

#endif /* FROBENIUS_H */
