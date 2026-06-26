#ifndef PARALLEL_TRANSPORT_H
#define PARALLEL_TRANSPORT_H

#include <stdbool.h>
#include "lie_group.h"
#include "principal_bundle.h"

/* ============================================================================
 * Parallel Transport Theory
 *
 * Given a connection on a principal bundle P → M, parallel transport
 * along a curve γ: [0,1] → M is the unique horizontal lift γ̃: [0,1] → P
 * such that π(γ̃(t)) = γ(t) and γ̃'(t) ∈ H_{γ̃(t)}.
 *
 * The parallel transport operator P_γ: P_{γ(0)} → P_{γ(1)} satisfies:
 *   P_γ(p·g) = P_γ(p)·g  (right-equivariance)
 *
 * In terms of the connection 1-form ω, the horizontal lift condition is:
 *   ω(γ̃'(t)) = 0  for all t
 *
 * For a gauge potential A_μ, parallel transport of a section s(x) is
 * governed by the ODE:
 *   d/dt ψ(t) + A_μ(γ(t)) · dγ^μ/dt · ψ(t) = 0
 * where ψ(t) represents the fiber coordinate along the lift.
 *
 * On the lattice, parallel transport along a sequence of links is simply
 * the ordered product of link variables.
 *
 * References:
 *   Kobayashi & Nomizu (1963) Vol I, Ch.II §1-3
 *   Nakahara (2003), §10.1
 *   Marsden & Ratiu (1999), Introduction to Mechanics and Symmetry
 * ============================================================================ */

/* --- Continuous curve in the base manifold --- */
typedef struct {
    int     dim;           /* dimension of the base manifold */
    int     n_points;      /* number of sample points along curve */
    double *points;        /* n_points × dim array of coordinates */
    double  t_start;       /* start parameter value */
    double  t_end;         /* end parameter value */
} BaseCurve;

/* --- Tangent vector at a point --- */
typedef struct {
    int     dim;
    double  v[4];          /* components of tangent vector (max dim 4) */
} TangentVector;

/* --- Parallel transport result along a curve --- */
typedef struct {
    int               n_steps;
    LieGroupElement  *transports; /* transport from start to each step */
    LieGroupElement   total;      /* total transport from start to end */
    bool              success;
} ParallelTransportResult;

/* ============================================================================
 * Base Curve Operations (L2-L3)
 * =========================================================================== */

/** Create a curve sampled at n_points. O(n_points) */
BaseCurve *curve_create(int dim, int n_points);
void curve_free(BaseCurve *curve);

/** Set a point on the curve. O(1) */
void curve_set_point(BaseCurve *curve, int idx, const double *coords);

/** Get a point on the curve. O(1) */
const double *curve_get_point(const BaseCurve *curve, int idx);

/** Get the tangent vector at index i using central differences. O(dim) */
void curve_tangent(const BaseCurve *curve, int idx, TangentVector *out);

/** Create a straight line curve from p0 to p1 with n_points. O(n_points) */
BaseCurve *curve_straight_line(const double *p0, const double *p1,
                                int dim, int n_points);

/** Create a circular curve (great circle on S^2, or circle in R^d). O(n_points) */
BaseCurve *curve_circle(const double *center, double radius,
                         int dim, int n_points);

/* ============================================================================
 * Parallel Transport on the Bundle (L3-L4)
 * =========================================================================== */

/**
 * Parallel transport a fiber element along a continuous curve,
 * using the lattice connection with interpolation.
 *
 * Algorithm: Decompose the curve into lattice links, then
 * multiply link variables in order. Uses nearest-neighbor
 * interpolation to map curve points to lattice sites.
 *
 * Theorem (Kobayashi-Nomizu): Parallel transport exists and is
 * unique for any piecewise-smooth curve with given initial fiber point.
 *
 * Complexity: O(curve_points * dim)
 */
ParallelTransportResult *pt_transport(
    const PrincipalBundle *pb,
    const BaseCurve *curve,
    LieGroupElement initial_fiber);

void pt_result_free(ParallelTransportResult *result);

/**
 * Parallel transport of a tangent vector (associated vector bundle).
 * Uses the connection to transport a vector in the associated
 * vector bundle (e.g., tangent bundle for frame bundles).
 *
 * For SO(3) frame bundle: transports a 3D vector by applying rotation.
 *
 * Complexity: O(curve_points)
 */
void pt_transport_vector(const PrincipalBundle *pb,
                          const BaseCurve *curve,
                          const double *initial_vector,
                          double *final_vector);

/* ============================================================================
 * Geodesic Equation / Horizontal Lifts (L5)
 * =========================================================================== */

/**
 * Compute the horizontal lift of a curve given an initial fiber point.
 * The horizontal lift γ̃ satisfies: (i) π∘γ̃ = γ, (ii) ω(γ̃') = 0.
 *
 * On the lattice, the lift at each step is:
 *   g_{k+1} = g_k · U_k        (right action convention)
 * where U_k is the link variable along the k-th step.
 *
 * Complexity: O(curve_points)
 */
ParallelTransportResult *pt_horizontal_lift(
    const PrincipalBundle *pb,
    const BaseCurve *curve,
    LieGroupElement initial_fiber);

/**
 * Check if a curve in the total space is horizontal.
 * A curve is horizontal if all its tangent vectors are in H_p.
 * On the lattice: check if each step is orthogonal to the vertical
 * direction defined by the connection.
 *
 * Complexity: O(curve_points)
 */
bool pt_is_horizontal(const PrincipalBundle *pb, const BaseCurve *curve);

/* ============================================================================
 * Holonomy Map (L4-L5)
 * =========================================================================== */

/**
 * Compute the holonomy map for a closed curve.
 * For a loop γ based at x, the holonomy is:
 *   hol(γ) = P_γ ∈ G
 * which maps the fiber over x to itself.
 *
 * hol(γ): π^{-1}(x) → π^{-1}(x),  hol(γ)(p) = P_γ(p)
 *
 * Theorem: The set of all holonomies forms a subgroup of G,
 * the holonomy group at x.
 *
 * Complexity: O(curve_points)
 */
LieGroupElement pt_holonomy(const PrincipalBundle *pb, const BaseCurve *loop);

#endif /* PARALLEL_TRANSPORT_H */
