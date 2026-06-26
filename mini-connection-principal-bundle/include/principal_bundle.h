#ifndef PRINCIPAL_BUNDLE_H
#define PRINCIPAL_BUNDLE_H

#include <stdbool.h>
#include <stddef.h>
#include "lie_group.h"

/* ============================================================================
 * Principal Bundle Theory
 *
 * A principal G-bundle P over a base manifold M consists of:
 *   - Total space P with free right G-action
 *   - Base manifold M ≅ P/G
 *   - Projection π: P → M
 *   - Local trivializations φ_i: π^{-1}(U_i) → U_i × G
 *
 * In the lattice formulation (Wilson, 1974), the base manifold is discretized
 * as a d-dimensional lattice Λ. The bundle structure is encoded in:
 *   - Link variables U_μ(x) ∈ G for each lattice link (x, x+μ̂)
 *   - Plaquette variables for curvature
 *   - Gauge transformations g(x) ∈ G at each site
 *
 * References:
 *   Kobayashi & Nomizu (1963) Vol I, Ch.II §5
 *   Wilson, Phys. Rev. D10 (1974) 2445
 *   Creutz, Quarks, Gluons and Lattices (1983)
 *   Montvay & Münster, Quantum Fields on a Lattice (1994)
 * ============================================================================ */

/* --- Lattice geometry --- */
typedef struct {
    int    dim;             /* spacetime dimension (typically 2, 3, or 4) */
    int    L[4];            /* lattice size in each dimension */
    int    n_sites;         /* total number of lattice sites = prod(L) */
    int    n_links;         /* total number of directed links */
    double a;               /* lattice spacing */
} LatticeGeometry;

/* --- Transition functions (cocycle condition) --- */
typedef struct {
    int           n_patches;    /* number of open sets in the cover */
    int           n_overlaps;   /* number of non-empty overlaps U_i ∩ U_j */
    LieGroupType  group_type;
    /* g_ij(x): overlap transition function for each overlap pair */
    LieGroupElement *transitions;  /* indexed by overlap index */
    int            *overlap_i;     /* patch i for overlap k */
    int            *overlap_j;     /* patch j for overlap k */
} TransitionFunctions;

/* --- Local section (local trivialization) --- */
typedef struct {
    int            patch_id;
    LieGroupType   group_type;
    /* The local section σ: U_i → P maps base points to bundle points */
    int            n_sample_points;
    double         *base_coords;   /* n_sample_points × base_dim */
    LieGroupElement *fiber_values; /* n_sample_points elements */
} LocalSection;

/* --- Principal bundle (lattice representation) --- */
typedef struct {
    LieGroupType    structure_group;
    LatticeGeometry lattice;
    /* Link variables: U_μ(x) ∈ G for each directed link.
     * Storage: links[l*nsites + x] where l = mu*dir encoding */
    LieGroupElement *link_vars;
    /* Gauge transformation: g(x) ∈ G at each site */
    LieGroupElement *gauge_transform;
    /* Flag indicating whether link variables have been initialized */
    bool             initialized;
} PrincipalBundle;

/* ============================================================================
 * Lattice Geometry (L3)
 * =========================================================================== */

/**
 * Create a d-dimensional lattice with sizes L[0..d-1].
 * Total sites = prod L[i], total links = d * prod L[i].
 * Complexity: O(d)
 */
LatticeGeometry lattice_create(int dim, const int *L, double a);

/** Compute the linear index of a site given coordinates. O(d) */
int lattice_site_index(const LatticeGeometry *geom, const int *coords);

/** Compute coordinates from a linear site index. O(d) */
void lattice_site_coords(const LatticeGeometry *geom, int idx, int *coords);

/** Get neighbor site index in direction mu (+1 forward, -1 backward). O(1) */
int lattice_neighbor(const LatticeGeometry *geom, int site, int mu, int dir);

/** Link index encoding: link = (site_index, mu). O(1) */
int lattice_link_index(const LatticeGeometry *geom, int site, int mu);

/** Check if coordinates are within bounds. O(d) */
bool lattice_in_bounds(const LatticeGeometry *geom, const int *coords);

/* ============================================================================
 * Principal Bundle Operations (L1-L3)
 * =========================================================================== */

/**
 * Allocate and initialize a principal bundle on a lattice.
 * All link variables initialized to identity.
 * Complexity: O(n_links)
 */
PrincipalBundle *pb_create(LieGroupType group, const LatticeGeometry *geom);

/** Free bundle and all associated storage. O(1) */
void pb_free(PrincipalBundle *pb);

/** Get link variable U_μ(x). O(1) */
LieGroupElement pb_get_link(const PrincipalBundle *pb, int site, int mu);

/** Set link variable U_μ(x). O(1) */
void pb_set_link(PrincipalBundle *pb, int site, int mu, LieGroupElement U);

/** Set all link variables to identity (trivial bundle). O(n_links) */
void pb_set_trivial(PrincipalBundle *pb);

/** Copy link variables from src to dest. Both must have identical geometry. O(n_links) */
void pb_copy_links(PrincipalBundle *dest, const PrincipalBundle *src);

/** Check if bundle is trivial (all link vars = identity). O(n_links) */
bool pb_is_trivial(const PrincipalBundle *pb);

/* Gauge transformations */
void pb_gauge_transform(PrincipalBundle *pb);
void pb_set_gauge_at_site(PrincipalBundle *pb, int site, LieGroupElement g);
LieGroupElement pb_get_gauge_at_site(const PrincipalBundle *pb, int site);
void pb_random_gauge_transform(PrincipalBundle *pb, unsigned int seed);

/* Bundle topology */
double pb_winding_number_u1(const PrincipalBundle *pb);

/* ============================================================================
 * Transition Functions (L3)
 * =========================================================================== */

/**
 * Create transition functions for a bundle with n_patches.
 * The cocycle condition: g_ij · g_jk · g_ki = id on U_i ∩ U_j ∩ U_k.
 * Complexity: O(n_overlaps)
 */
TransitionFunctions *tf_create(int n_patches, LieGroupType group);
void tf_free(TransitionFunctions *tf);

/** Set transition function g_ij on overlap U_i ∩ U_j. */
void tf_set(TransitionFunctions *tf, int overlap_idx,
            int i, int j, LieGroupElement g);

/**
 * Verify the cocycle condition for all triple overlaps.
 * Returns true if g_ij(x)·g_jk(x)·g_ki(x) = identity for all triples.
 * Complexity: O(n_overlaps^3)
 */
bool tf_check_cocycle(const TransitionFunctions *tf);

/* ============================================================================
 * Bundle Reconstruction from Transition Functions (L3-L4)
 * =========================================================================== */

/**
 * Reconstruct lattice link variables from transition functions.
 * Given a cover {U_i} and transitions g_ij, construct the bundle by
 * assigning links to patches and using transitions on boundaries.
 * Complexity: O(n_sites * n_patches)
 */
void pb_from_transitions(PrincipalBundle *pb, const TransitionFunctions *tf);

#endif /* PRINCIPAL_BUNDLE_H */
