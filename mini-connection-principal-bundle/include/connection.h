#ifndef CONNECTION_H
#define CONNECTION_H

#include <stdbool.h>
#include "lie_group.h"
#include "principal_bundle.h"

/* ============================================================================
 * Connection Theory on Principal Bundles
 *
 * An Ehresmann connection on a principal G-bundle P is a smooth choice of
 * horizontal subspaces H_p ⊂ T_p P such that:
 *   (1) T_p P = H_p ⊕ V_p (direct sum decomposition)
 *   (2) H_{p·g} = (R_g)_* H_p (right-equivariance)
 *
 * Equivalently, a connection 1-form ω: TP → g satisfying:
 *   (a) ω(A^#) = A for A ∈ g (vertical vectors reproduce algebra element)
 *   (b) R_g^* ω = Ad_{g^{-1}} ◦ ω (equivariance)
 *
 * In lattice gauge theory, the connection is represented by link variables
 * U_μ(x) = P exp(-∫_x^{x+μ̂} A_μ dx^μ), the path-ordered exponential of the
 * gauge potential along a link.
 *
 * The gauge potential A_μ(x) ∈ g in the continuum satisfies:
 *   A_μ → g A_μ g^{-1} + g ∂_μ g^{-1} under gauge transformation.
 *
 * References:
 *   Ehresmann (1951), Colloque de Topologie, Bruxelles
 *   Kobayashi & Nomizu (1963) Vol I, Ch.II §1-2
 *   Wilson (1974), Phys. Rev. D10, 2445
 *   Nakahara (2003), §10.1-10.2
 * ============================================================================ */

/* --- Gauge potential A_μ(x) ∈ g (continuum representation) --- */
typedef struct {
    int               dim;            /* spacetime dimension */
    int               n_points;       /* number of evaluation points */
    double            *coords;        /* n_points × dim array of coordinates */
    LieAlgebraElement *A_values;      /* n_points × dim array of A_μ values */
    LieGroupType       group_type;
} GaugePotential;

/* --- Connection on a lattice (link variable representation) --- */
typedef struct {
    PrincipalBundle   *bundle;        /* underlying principal bundle */
    /* The connection IS the set of link variables in the bundle.
     * Additional fields for connection-specific operations: */
} LatticeConnection;

/* ============================================================================
 * Gauge Potential Operations (L2-L3)
 * =========================================================================== */

/**
 * Allocate gauge potential on a grid of points.
 * A_μ(x) stored as A[idx_pt * dim + mu] ∈ Lie algebra.
 * Complexity: O(n_points * dim)
 */
GaugePotential *gp_create(int dim, int n_points, double *coords,
                          LieGroupType group);

/** Free gauge potential. O(1) */
void gp_free(GaugePotential *gp);

/** Set A_μ at a specific point. O(1) */
void gp_set_A(GaugePotential *gp, int pt, int mu, LieAlgebraElement A);

/** Get A_μ at a specific point. O(1) */
LieAlgebraElement gp_get_A(const GaugePotential *gp, int pt, int mu);

/** Compute the field strength tensor F_{μν} = ∂_μ A_ν - ∂_ν A_μ + [A_μ, A_ν].
 *  Uses central finite differences for derivatives.
 *  Complexity: O(1) per point/direction pair */
LieAlgebraElement gp_field_strength(const GaugePotential *gp, int pt,
                                     int mu, int nu, double dx);

/* ============================================================================
 * Lattice Connection Operations (L2-L3)
 * =========================================================================== */

/**
 * Create a lattice connection (wraps an existing bundle).
 * The connection is represented by the bundle's link variables.
 * Complexity: O(1)
 */
LatticeConnection *lc_create(PrincipalBundle *pb);
void lc_free(LatticeConnection *lc);

/** Get link variable U_μ(x) — delegates to pb_get_link. O(1) */
LieGroupElement lc_get_link(const LatticeConnection *lc, int site, int mu);

/** Set link variable. O(1) */
void lc_set_link(LatticeConnection *lc, int site, int mu, LieGroupElement U);

/**
 * Compute gauge potential A_μ(x) from link variable U_μ(x).
 * A_μ(x) = -log(U_μ(x)) / a  (first-order approximation).
 * For non-abelian groups, log is the principal logarithm.
 * Complexity: O(1) per link
 */
LieAlgebraElement lc_link_to_A(const LatticeConnection *lc, int site, int mu);

/**
 * Compute link variable U_μ(x) from gauge potential A_μ(x).
 * U_μ(x) = exp(-a·A_μ(x)) where a is lattice spacing.
 * Complexity: O(1) per link
 */
LieGroupElement lc_A_to_link(const LatticeConnection *lc,
                              LieAlgebraElement A, double a);

/** Set all link variables from a continuum gauge potential. O(n_links) */
void lc_from_gauge_potential(LatticeConnection *lc, const GaugePotential *gp);

/** Set gauge potential values from link variables. O(n_links) */
void lc_to_gauge_potential(const LatticeConnection *lc, GaugePotential *gp);

/* ============================================================================
 * Connection Properties (L3-L4)
 * =========================================================================== */

/**
 * Verify horizontal distribution properties:
 * (1) Non-degeneracy: H_p ∩ V_p = {0} (encoded in group action)
 * (2) Right-equivariance: verified by gauge covariance of link variables
 * Returns 0 if connection is valid, error code otherwise.
 * Complexity: O(n_links)
 */
int lc_verify_properties(const LatticeConnection *lc);

/** Check if a connection is flat (all plaquettes = identity). O(n_links) */
bool lc_is_flat(const LatticeConnection *lc);

/** Compute the Maurer-Cartan form on the group. O(1) */
LieAlgebraElement lc_maurer_cartan(LieGroupElement g);

/** Compute the connection 1-form evaluated on a tangent vector. O(1) */
LieAlgebraElement lc_omega(const LatticeConnection *lc, int site, int mu);

#endif /* CONNECTION_H */
