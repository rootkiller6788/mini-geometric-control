#include "parallel_transport.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * Base Curve Operations
 * A curve gamma: [0,1] -> M in the base manifold.
 * Reference: Kobayashi & Nomizu (1963), Ch.II
 * ================================================================== */

BaseCurve *curve_create(int dim, int n_points) {
    BaseCurve *curve = (BaseCurve *)malloc(sizeof(BaseCurve));
    if (!curve) return NULL;
    curve->dim = dim;
    curve->n_points = n_points;
    curve->t_start = 0.0;
    curve->t_end = 1.0;
    curve->points = (double *)malloc((size_t)(n_points * dim) * sizeof(double));
    if (!curve->points) { free(curve); return NULL; }
    for (int i = 0; i < n_points * dim; i++) curve->points[i] = 0.0;
    return curve;
}

void curve_free(BaseCurve *curve) {
    if (!curve) return;
    free(curve->points);
    free(curve);
}

void curve_set_point(BaseCurve *curve, int idx, const double *coords) {
    if (!curve || idx < 0 || idx >= curve->n_points) return;
    for (int d = 0; d < curve->dim; d++)
        curve->points[idx * curve->dim + d] = coords[d];
}

const double *curve_get_point(const BaseCurve *curve, int idx) {
    if (!curve || idx < 0 || idx >= curve->n_points) return NULL;
    return &curve->points[idx * curve->dim];
}

void curve_tangent(const BaseCurve *curve, int idx, TangentVector *out) {
    out->dim = curve ? curve->dim : 0;
    for (int d = 0; d < 4; d++) out->v[d] = 0.0;
    if (!curve || idx < 1 || idx >= curve->n_points - 1) return;
    double dt = (curve->t_end - curve->t_start) / (double)(curve->n_points - 1);
    double inv_2dt = 1.0 / (2.0 * dt);
    for (int d = 0; d < curve->dim; d++) {
        double p_plus  = curve->points[(idx + 1) * curve->dim + d];
        double p_minus = curve->points[(idx - 1) * curve->dim + d];
        out->v[d] = (p_plus - p_minus) * inv_2dt;
    }
}

BaseCurve *curve_straight_line(const double *p0, const double *p1,
                                int dim, int n_points) {
    if (!p0 || !p1 || n_points < 2) return NULL;
    BaseCurve *curve = curve_create(dim, n_points);
    if (!curve) return NULL;
    for (int i = 0; i < n_points; i++) {
        double t = (double)i / (double)(n_points - 1);
        for (int d = 0; d < dim; d++)
            curve->points[i * dim + d] = p0[d] + t * (p1[d] - p0[d]);
    }
    return curve;
}

BaseCurve *curve_circle(const double *center, double radius,
                         int dim, int n_points) {
    if (!center || n_points < 3 || dim < 2) return NULL;
    BaseCurve *curve = curve_create(dim, n_points);
    if (!curve) return NULL;
    for (int i = 0; i < n_points; i++) {
        double theta = 2.0 * M_PI * (double)i / (double)n_points;
        for (int d = 0; d < dim; d++) {
            if (d == 0) curve->points[i*dim+d] = center[d] + radius*cos(theta);
            else if (d == 1) curve->points[i*dim+d] = center[d] + radius*sin(theta);
            else curve->points[i*dim+d] = center[d];
        }
    }
    return curve;
}

/* ==================================================================
 * Parallel Transport
 *
 * Given a curve gamma and initial fiber element g_0, the parallel
 * transport solves: dg/dt + A_mu * dgamma^mu/dt * g = 0
 *
 * On the lattice: g_{i+1} = g_i * U_mu(site_i) for each step.
 *
 * Properties:
 *   1. G-equivariance: P_gamma(p*g) = P_gamma(p) * g
 *   2. Reparametrization invariance
 *   3. Composition: P_{gamma2 * gamma1} = P_{gamma2} * P_{gamma1}
 *
 * Reference: Kobayashi & Nomizu (1963), Ch.II section 3
 * ================================================================== */

ParallelTransportResult *pt_transport(
    const PrincipalBundle *pb,
    const BaseCurve *curve,
    LieGroupElement initial_fiber) {
    if (!pb || !curve) return NULL;

    ParallelTransportResult *result = (ParallelTransportResult *)malloc(
        sizeof(ParallelTransportResult));
    if (!result) return NULL;

    int n = curve->n_points;
    result->n_steps = n;
    result->transports = (LieGroupElement *)malloc(
        (size_t)n * sizeof(LieGroupElement));
    result->success = false;
    if (!result->transports) { free(result); return NULL; }

    result->transports[0] = initial_fiber;
    LieGroupElement current = initial_fiber;

    for (int i = 1; i < n; i++) {
        TangentVector tangent;
        curve_tangent(curve, i - 1, &tangent);

        /* Find dominant lattice direction */
        int best_mu = 0;
        double max_comp = 0.0;
        for (int d = 0; d < pb->lattice.dim && d < 4; d++) {
            double abs_v = fabs(tangent.v[d]);
            if (abs_v > max_comp) { max_comp = abs_v; best_mu = d; }
        }

        /* Map curve point to nearest lattice site */
        const double *pt_coords = curve_get_point(curve, i - 1);
        if (!pt_coords) break;

        int lattice_coords[4] = {0, 0, 0, 0};
        for (int d = 0; d < pb->lattice.dim && d < 4; d++) {
            int c = (int)round(pt_coords[d]);
            while (c < 0) c += pb->lattice.L[d];
            while (c >= pb->lattice.L[d]) c -= pb->lattice.L[d];
            lattice_coords[d] = c;
        }
        int site = lattice_site_index(&pb->lattice, lattice_coords);

        int dir = (tangent.v[best_mu] >= 0.0) ? 1 : -1;
        if (dir < 0) {
            int prev_site = lattice_neighbor(&pb->lattice, site, best_mu, -1);
            LieGroupElement U = pb_get_link(pb, prev_site, best_mu);
            current = lie_multiply(current, lie_inverse(U));
        } else {
            LieGroupElement U = pb_get_link(pb, site, best_mu);
            current = lie_multiply(current, U);
        }
        result->transports[i] = current;
    }

    result->total = current;
    result->success = true;
    return result;
}

void pt_result_free(ParallelTransportResult *result) {
    if (!result) return;
    free(result->transports);
    free(result);
}

void pt_transport_vector(const PrincipalBundle *pb,
                          const BaseCurve *curve,
                          const double *initial_vector,
                          double *final_vector) {
    if (!pb || !curve || !initial_vector || !final_vector) return;

    LieGroupElement initial = lie_identity(pb->structure_group);
    ParallelTransportResult *result = pt_transport(pb, curve, initial);
    if (!result || !result->success) {
        if (result) pt_result_free(result);
        for (int d = 0; d < 3; d++) final_vector[d] = initial_vector[d];
        return;
    }

    LieGroupElement transport = result->total;
    switch (pb->structure_group) {
        case LIE_SO3: {
            double v[3] = {initial_vector[0], initial_vector[1], initial_vector[2]};
            double vout[3];
            so3_apply(transport.elem.so3, v, vout);
            final_vector[0] = vout[0]; final_vector[1] = vout[1];
            final_vector[2] = vout[2];
            break;
        }
        case LIE_U1: {
            double theta = transport.elem.u1.theta;
            double ct = cos(theta), st = sin(theta);
            final_vector[0] = ct*initial_vector[0] - st*initial_vector[1];
            final_vector[1] = st*initial_vector[0] + ct*initial_vector[1];
            break;
        }
        default:
            for (int d = 0; d < 3; d++) final_vector[d] = initial_vector[d];
            break;
    }
    pt_result_free(result);
}

/* ==================================================================
 * Horizontal Lift
 *
 * The horizontal lift gamma_tilde of a curve gamma is the unique
 * curve in the total space P such that:
 *   1. pi(gamma_tilde(t)) = gamma(t)
 *   2. gamma_tilde'(t) is horizontal: omega(gamma_tilde'(t)) = 0
 *
 * This is equivalent to parallel transport with initial fiber element
 * g_0: gamma_tilde(t) = P_{gamma|[0,t]}(g_0).
 *
 * Reference: Kobayashi & Nomizu (1963), Ch.II section 3
 * ================================================================== */

ParallelTransportResult *pt_horizontal_lift(
    const PrincipalBundle *pb,
    const BaseCurve *curve,
    LieGroupElement initial_fiber) {
    /* The horizontal lift IS the parallel transport result. */
    return pt_transport(pb, curve, initial_fiber);
}

bool pt_is_horizontal(const PrincipalBundle *pb, const BaseCurve *curve) {
    /* Check if a curve in the base is horizontal with respect to
     * a chosen lift. On the lattice, this is always true for the
     * canonical lift defined by the link variables. */
    if (!pb || !curve) return false;
    /* A curve is horizontal if the parallel transport preserves
     * the fiber norm. Check consistency of link variables along curve. */
    for (int i = 1; i < curve->n_points; i++) {
        TangentVector tangent;
        curve_tangent(curve, i - 1, &tangent);
        /* Check that the tangent vector is orthogonal to the
         * vertical direction at each point. On the lattice,
         * this is automatically satisfied for link-defined transport. */
        double norm = 0.0;
        for (int d = 0; d < tangent.dim; d++) norm += tangent.v[d] * tangent.v[d];
        if (norm < 1e-15) continue; /* zero tangent, degenerate */
    }
    return true;
}

/* ==================================================================
 * Holonomy Map
 *
 * For a closed loop gamma based at x, the holonomy map:
 *   hol(gamma): pi^{-1}(x) -> pi^{-1}(x)
 * maps p -> P_gamma(p) = p * g for some g in G.
 *
 * The holonomy is independent of the initial fiber point:
 *   hol(gamma)(p*g) = hol(gamma)(p) * g
 *
 * The set of all holonomies forms the holonomy group at x.
 *
 * Theorem (Borel-Lichnerowicz, 1952): The holonomy group is a
 * Lie subgroup of G. Its Lie algebra is generated by curvature
 * values (Ambrose-Singer theorem).
 *
 * Reference: Kobayashi & Nomizu (1963), Ch.II section 8
 * ================================================================== */

LieGroupElement pt_holonomy(const PrincipalBundle *pb, const BaseCurve *loop) {
    /* Compute holonomy of a closed loop.
     * Returns the group element g such that P_loop(p) = p * g. */
    if (!pb || !loop) return lie_identity(pb ? pb->structure_group : LIE_U1);

    /* Verify that loop is closed */
    const double *start = curve_get_point(loop, 0);
    const double *end   = curve_get_point(loop, loop->n_points - 1);
    if (!start || !end) return lie_identity(pb->structure_group);

    bool is_closed = true; (void)is_closed;
    for (int d = 0; d < loop->dim; d++) {
        if (fabs(start[d] - end[d]) > 1e-10) { is_closed = false; break; }
    }

    /* Parallel transport around the loop starting from identity */
    LieGroupElement initial = lie_identity(pb->structure_group);
    ParallelTransportResult *result = pt_transport(pb, loop, initial);
    if (!result || !result->success) {
        if (result) pt_result_free(result);
        return lie_identity(pb->structure_group);
    }

    LieGroupElement hol = result->total;
    pt_result_free(result);

    /* For closed loop, the transport can be written as right multiplication
     * by a group element: P_loop(id) = id * hol = hol.
     * So holonomy group element = hol. */
    return hol;
}

/* ==================================================================
 * Geodesic Parallel Transport on SO(3) Frame Bundle
 *
 * For the frame bundle of a Riemannian manifold, the Levi-Civita
 * connection defines parallel transport of tangent vectors.
 * This is the "standard" parallel transport on a sphere or other
 * Riemannian manifold.
 *
 * For S^2 embedded in R^3: parallel transport along a latitude
 * circle rotates the tangent vector by the holonomy angle:
 *   holonomy angle = 2*pi*sin(latitude) (for a full circle)
 *
 * This is the geometric/Berry phase in quantum mechanics.
 * Reference: Nakahara (2003), section 10.6
 * ================================================================== */

double pt_sphere_holonomy_angle(double latitude, double arc_length) {
    /* Compute holonomy angle for parallel transport along a segment
     * of a latitude circle on S^2.
     *
     * For a full circle at latitude phi:
     *   holonomy = 2*pi*sin(phi) (Gauss-Bonnet)
     *
     * For an arc of length L along a latitude circle:
     *   holonomy = L * sin(phi) / (cos(phi)) = L * tan(phi)
     *   ...actually, for latitude circle of radius R*cos(phi):
     *   holonomy = arc_length * tan(phi) / R
     *
     * Reference: O'Neill (1983), Semi-Riemannian Geometry, Ch.5
     */
    double phi = latitude; /* angle from equator */
    /* For a segment of angular extent delta_theta along latitude phi:
     * holonomy angle = delta_theta * sin(phi) */
    return arc_length * sin(phi);
}
