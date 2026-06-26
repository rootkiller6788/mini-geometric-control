#ifndef HOLONOMY_H
#define HOLONOMY_H

#include <stdbool.h>
#include "lie_group.h"
#include "principal_bundle.h"

/* ============================================================================
 * Holonomy Theory
 *
 * For a connection on a principal bundle P → M, the holonomy group Φ(p)
 * at p ∈ P is the set of all elements g ∈ G obtained by parallel transport
 * around loops based at π(p):
 *
 *   Φ(p) = { g ∈ G | ∃ loop γ based at x=π(p) s.t. P_γ(p) = p·g }
 *
 * The restricted holonomy group Φ^0(p) uses only contractible loops.
 *
 * Ambrose-Singer Theorem (1953): The Lie algebra of the holonomy group
 * Φ(p) is spanned by all values of the curvature Ω(X,Y) and its covariant
 * derivatives at points reachable by parallel transport from p.
 *
 * On the lattice, holonomy along a path γ = (x, μ_1, ..., μ_k) is:
 *   U_γ = U_{μ_k}(x_{k-1}) · ... · U_{μ_2}(x_1) · U_{μ_1}(x_0)
 *
 * References:
 *   Ambrose & Singer (1953), Trans. AMS 75, 428-443
 *   Kobayashi & Nomizu (1963) Vol I, Ch.II §8
 *   Besse, Einstein Manifolds (1987), §10
 * ============================================================================ */

/* --- Path on the lattice --- */
typedef struct {
    int     length;        /* number of links in the path */
    int    *sites;         /* length+1 site indices along path */
    int    *directions;    /* length direction indices (mu_i) for each step */
    bool    is_loop;       /* true if sites[0] == sites[length] */
} LatticePath;

/* --- Holonomy group data --- */
typedef struct {
    LieGroupType     group;
    int              max_generators;
    int              n_generators;
    LieGroupElement *generators;     /* elements spanning the holonomy group */
    int              n_computed;     /* number of distinct holonomy elements found */
    LieGroupElement *elements;       /* computed holonomy elements (loop results) */
} HolonomyGroup;

/* ============================================================================
 * Lattice Path Operations (L3)
 * =========================================================================== */

/**
 * Create a lattice path of given length.
 * Allocates sites array (length+1) and directions array (length).
 * Complexity: O(length)
 */
LatticePath *path_create(int length);
void path_free(LatticePath *path);

/** Set the i-th step of the path. O(1) */
void path_set_step(LatticePath *path, int step, int from_site, int mu);

/** Create a straight path in direction mu of given length. O(length) */
LatticePath *path_straight(const LatticeGeometry *geom, int start_site,
                            int mu, int length);

/** Create a rectangular loop in the μ-ν plane. O(length_mu + length_nu) */
LatticePath *path_rectangle(const LatticeGeometry *geom, int start_site,
                             int mu, int nu, int len_mu, int len_nu);

/** Check if a path is valid (all steps connect). O(length) */
bool path_is_valid(const LatticePath *path, const LatticeGeometry *geom);

/** Copy a path. O(length) */
LatticePath *path_copy(const LatticePath *src);

/** Concatenate two paths (path1 followed by path2). O(l1+l2) */
LatticePath *path_concat(const LatticePath *p1, const LatticePath *p2);

/** Reverse a path. O(length) */
LatticePath *path_reverse(const LatticePath *path);

/* ============================================================================
 * Parallel Transport Along Path (L3-L4)
 * =========================================================================== */

/**
 * Compute holonomy / parallel transport along a path.
 * U_γ = U(path[last]) · ... · U(path[1]) · U(path[0])
 *
 * This is the Wilson line operator. For a closed loop, it gives the
 * holonomy element, which is gauge-covariant:
 *   U_γ → g(start) · U_γ · g(start)^{-1}
 *
 * Complexity: O(path_length)
 */
LieGroupElement pb_parallel_transport_path(const PrincipalBundle *pb,
                                            const LatticePath *path);

/**
 * Check if a loop holonomy is trivial (flat connection along that loop).
 * Complexity: O(path_length)
 */
bool pb_loop_trivial(const PrincipalBundle *pb, const LatticePath *loop);

/* ============================================================================
 * Holonomy Group Computation (L4-L5)
 * =========================================================================== */

/**
 * Allocate and initialize a holonomy group tracker.
 * Complexity: O(1)
 */
HolonomyGroup *holonomy_create(LieGroupType group, int max_generators);

void holonomy_free(HolonomyGroup *hg);

/** Add a computed holonomy element to the group tracker. O(1) */
void holonomy_add_element(HolonomyGroup *hg, LieGroupElement elem);

/** Add a generator spanning element. O(1) */
void holonomy_add_generator(HolonomyGroup *hg, LieGroupElement gen);

/**
 * Compute holonomy group dimension from generators.
 * For U(1): dimension is 0 (flat) or 1 (curved).
 * For SU(2): dimension approximated by number of independent generators.
 * Complexity: O(n_generators)
 */
int holonomy_group_dimension(const HolonomyGroup *hg);

/**
 * Verify Ambrose-Singer: the Lie algebra of holonomy is spanned by
 * curvature values. Check if curvature at each site generates
 * elements in the holonomy algebra.
 * Complexity: O(n_sites * n_generators)
 */
bool holonomy_check_ambrose_singer(const HolonomyGroup *hg,
                                    const PrincipalBundle *pb,
                                    double tolerance);

/* ============================================================================
 * Wilson Loop Computation (L5-L6)
 * =========================================================================== */

/**
 * Compute a rectangular Wilson loop W(R, T) of size R×T.
 * W(R,T) = Tr P exp(i ∮_C A_μ dx^μ) = Tr(U_γ)
 *
 * The Wilson loop is the basic observable in lattice gauge theory.
 * Its expectation value gives the static quark-antiquark potential:
 *   ⟨W(R,T)⟩ ~ exp(-V(R)·T) as T→∞
 *
 * Complexity: O(R+T)
 */
double pb_wilson_loop_rect(const PrincipalBundle *pb, int site,
                            int mu, int nu, int R, int T);

/** Compute Creutz ratio: χ(R,T) = -log(W(R,T)·W(R-1,T-1)/(W(R,T-1)·W(R-1,T))) */
double pb_creutz_ratio(const PrincipalBundle *pb, int site,
                        int mu, int nu, int R, int T);

#endif /* HOLONOMY_H */
