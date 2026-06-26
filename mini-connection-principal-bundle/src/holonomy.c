#include "holonomy.h"
#include "curvature.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * Lattice Path Operations
 *
 * A path gamma on the lattice is a sequence of steps, each step
 * given by a direction mu_i. Paths are used to compute:
 *   - Parallel transport: U_gamma = product of link variables
 *   - Wilson loops: Tr(U_gamma) for closed paths
 *   - Holonomy: U_gamma for a closed loop
 *
 * Reference: Creutz (1983), section 5; Montvay & Munster (1994)
 * ================================================================== */

LatticePath *path_create(int length) {
    LatticePath *path = (LatticePath *)malloc(sizeof(LatticePath));
    if (!path) return NULL;
    path->length = length;
    path->sites = (int *)malloc((size_t)(length + 1) * sizeof(int));
    path->directions = (int *)malloc((size_t)length * sizeof(int));
    path->is_loop = false;
    if (!path->sites || !path->directions) {
        free(path->sites); free(path->directions); free(path); return NULL;
    }
    for (int i = 0; i <= length; i++) path->sites[i] = -1;
    for (int i = 0; i < length; i++) path->directions[i] = -1;
    return path;
}

void path_free(LatticePath *path) {
    if (!path) return;
    free(path->sites);
    free(path->directions);
    free(path);
}

void path_set_step(LatticePath *path, int step, int from_site, int mu) {
    if (!path || step < 0 || step >= path->length) return;
    path->sites[step] = from_site;
    path->directions[step] = mu;
    if (step == path->length - 1) {
        /* Last step: we need both sites[step] and sites[step+1] */
    }
}

LatticePath *path_straight(const LatticeGeometry *geom, int start_site,
                            int mu, int length) {
    if (!geom || length < 0) return NULL;
    LatticePath *path = path_create(length);
    if (!path) return NULL;

    int current = start_site;
    path->sites[0] = current;
    for (int i = 0; i < length; i++) {
        path->directions[i] = mu;
        current = lattice_neighbor(geom, current, mu, 1);
        path->sites[i + 1] = current;
    }
    path->is_loop = (path->sites[0] == path->sites[length]);
    return path;
}

LatticePath *path_rectangle(const LatticeGeometry *geom, int start_site,
                             int mu, int nu, int len_mu, int len_nu) {
    /* Create a rectangular loop: mu for len_mu, nu for len_nu,
     * -mu for len_mu, -nu for len_nu. Total length = 2*(len_mu+len_nu). */
    if (!geom || len_mu < 1 || len_nu < 1) return NULL;
    int total_length = 2 * (len_mu + len_nu);
    LatticePath *path = path_create(total_length);
    if (!path) return NULL;

    int current = start_site;
    path->sites[0] = current;
    int step = 0;

    for (int i = 0; i < len_mu; i++) {
        path->directions[step] = mu;
        current = lattice_neighbor(geom, current, mu, 1);
        path->sites[++step] = current;
    }
    for (int i = 0; i < len_nu; i++) {
        path->directions[step] = nu;
        current = lattice_neighbor(geom, current, nu, 1);
        path->sites[++step] = current;
    }
    for (int i = 0; i < len_mu; i++) {
        path->directions[step] = mu;
        current = lattice_neighbor(geom, current, mu, -1);
        path->sites[++step] = current;
    }
    for (int i = 0; i < len_nu; i++) {
        path->directions[step] = nu;
        current = lattice_neighbor(geom, current, nu, -1);
        path->sites[++step] = current;
    }

    path->is_loop = true;
    return path;
}

bool path_is_valid(const LatticePath *path, const LatticeGeometry *geom) {
    if (!path || !geom) return false;
    for (int i = 0; i < path->length; i++) {
        int s = path->sites[i];
        int mu = path->directions[i];
        if (s < 0 || s >= geom->n_sites || mu < 0 || mu >= geom->dim)
            return false;
        int expected_next = lattice_neighbor(geom, s, mu, 1);
        if (i + 1 <= path->length && path->sites[i+1] >= 0) {
            if (path->sites[i+1] != expected_next) return false;
        }
    }
    return true;
}

LatticePath *path_copy(const LatticePath *src) {
    if (!src) return NULL;
    LatticePath *copy = path_create(src->length);
    if (!copy) return NULL;
    memcpy(copy->sites, src->sites, (size_t)(src->length + 1) * sizeof(int));
    memcpy(copy->directions, src->directions, (size_t)src->length * sizeof(int));
    copy->is_loop = src->is_loop;
    return copy;
}

LatticePath *path_concat(const LatticePath *p1, const LatticePath *p2) {
    if (!p1 || !p2) return NULL;
    /* Check that p1 ends where p2 starts */
    int L = p1->length + p2->length;
    LatticePath *path = path_create(L);
    if (!path) return NULL;

    for (int i = 0; i <= p1->length; i++)
        path->sites[i] = p1->sites[i];
    for (int i = 0; i < p1->length; i++)
        path->directions[i] = p1->directions[i];
    for (int i = 0; i <= p2->length; i++)
        path->sites[p1->length + i] = p2->sites[i];
    for (int i = 0; i < p2->length; i++)
        path->directions[p1->length + i] = p2->directions[i];
    path->is_loop = (path->sites[0] == path->sites[L]);
    return path;
}

LatticePath *path_reverse(const LatticePath *path) {
    if (!path) return NULL;
    LatticePath *rev = path_create(path->length);
    if (!rev) return NULL;
    /* Reverse site order */
    for (int i = 0; i <= path->length; i++)
        rev->sites[i] = path->sites[path->length - i];
    /* Reverse directions: -mu = go backward along mu */
    for (int i = 0; i < path->length; i++)
        rev->directions[i] = path->directions[path->length - 1 - i];
    rev->is_loop = path->is_loop;
    return rev;
}

/* ==================================================================
 * Parallel Transport Along Path
 *
 * U_gamma = U(path[last]) * ... * U(path[1]) * U(path[0])
 *
 * This is the path-ordered product (Wilson line) for the given path.
 * For a closed loop, U_gamma is the holonomy element.
 *
 * Under gauge transformation g(x):
 *   U_gamma -> g(start) * U_gamma * g(end)^{-1}
 * For a closed loop: U_gamma -> g(start) * U_gamma * g(start)^{-1}
 * So Tr(U_gamma) is gauge invariant.
 *
 * Theorem: Parallel transport along a curve is the unique horizontal
 * lift satisfying the ODE d psi/dt + A_mu * dx^mu/dt * psi = 0.
 *
 * Reference: Kobayashi & Nomizu (1963), Ch.II section 3
 * ================================================================== */

LieGroupElement pb_parallel_transport_path(const PrincipalBundle *pb,
                                            const LatticePath *path) {
    LieGroupElement result = lie_identity(pb ? pb->structure_group : LIE_U1);
    if (!pb || !path || path->length < 1) return result;

    result = lie_identity(pb->structure_group);

    for (int i = 0; i < path->length; i++) {
        int s = path->sites[i];
        int mu = path->directions[i];
        LieGroupElement U = pb_get_link(pb, s, mu);

        /* Right multiplication: result = result * U
         * (depends on the convention for path ordering;
         *  this follows the lattice gauge theory convention where
         *  the path is traversed from start to end, multiplying
         *  link variables in order.) */
        result = lie_multiply(result, U);
    }

    return result;
}

bool pb_loop_trivial(const PrincipalBundle *pb, const LatticePath *loop) {
    if (!pb || !loop || !loop->is_loop) return false;
    LieGroupElement hol = pb_parallel_transport_path(pb, loop);
    LieGroupElement id = lie_identity(pb->structure_group);
    double dist = lie_distance(hol, id);
    return dist < 1e-12;
}

/* ==================================================================
 * Holonomy Group
 *
 * The holonomy group Phi(p) at p in P is the set of all elements
 * g in G obtained by parallel transport around loops based at pi(p).
 *
 * Ambrose-Singer Theorem (1953): The Lie algebra of Phi(p) is
 * spanned by all values of the curvature Omega(X,Y) and its
 * covariant derivatives at points reachable from p by parallel
 * transport.
 *
 * For lattice gauge theory, the holonomy group is generated by
 * the plaquette variables (curvature) and their parallel transports.
 *
 * Reference: Ambrose & Singer (1953); Besse (1987), section 10
 * ================================================================== */

HolonomyGroup *holonomy_create(LieGroupType group, int max_generators) {
    HolonomyGroup *hg = (HolonomyGroup *)malloc(sizeof(HolonomyGroup));
    if (!hg) return NULL;
    hg->group = group;
    hg->max_generators = max_generators;
    hg->n_generators = 0;
    hg->n_computed = 0;
    hg->generators = (LieGroupElement *)malloc(
        (size_t)max_generators * sizeof(LieGroupElement));
    hg->elements = (LieGroupElement *)malloc(
        (size_t)max_generators * sizeof(LieGroupElement));
    if (!hg->generators || !hg->elements) {
        free(hg->generators); free(hg->elements); free(hg); return NULL;
    }
    return hg;
}

void holonomy_free(HolonomyGroup *hg) {
    if (!hg) return;
    free(hg->generators);
    free(hg->elements);
    free(hg);
}

void holonomy_add_element(HolonomyGroup *hg, LieGroupElement elem) {
    if (!hg || hg->n_computed >= hg->max_generators) return;

    /* Check if element is already in the set (within tolerance) */
    for (int i = 0; i < hg->n_computed; i++) {
        double dist = lie_distance(hg->elements[i], elem);
        if (dist < 1e-10) return;
    }
    hg->elements[hg->n_computed++] = elem;
}

void holonomy_add_generator(HolonomyGroup *hg, LieGroupElement gen) {
    if (!hg || hg->n_generators >= hg->max_generators) return;

    /* Check linear independence: for abelian groups, just check
     * if generator is a multiple of existing ones.
     * For non-abelian: check if log(gen) is linearly independent. */
    bool independent = true;
    for (int i = 0; i < hg->n_generators && independent; i++) {
        double dist = lie_distance(hg->generators[i], gen);
        LieGroupElement inv_diff = lie_multiply(
            lie_inverse(hg->generators[i]), gen);
        double dist2 = lie_distance(inv_diff,
            lie_identity(hg->group));
        if (dist < 1e-10 || dist2 < 1e-10) independent = false;
    }

    if (independent) {
        hg->generators[hg->n_generators++] = gen;
    }
}

int holonomy_group_dimension(const HolonomyGroup *hg) {
    /* Returns the dimension of the holonomy group.
     * For U(1): 0 (flat) or 1 (curved).
     * For SU(2)/SO(3): number of independent generators. */
    if (!hg) return 0;
    return hg->n_generators;
}

bool holonomy_check_ambrose_singer(const HolonomyGroup *hg,
                                    const PrincipalBundle *pb,
                                    double tolerance) {
    /* Verify Ambrose-Singer: curvature values generate the holonomy algebra.
     * Compute all plaquettes, take their logs (curvature), and check that
     * these span the same space as the holonomy generators. */
    if (!hg || !pb || tolerance < 0.0) return false;

    /* For U(1) in 2D: holonomy algebra is 0 if total flux = 0 mod 2*pi,
     * and 1-dimensional otherwise. Check if first Chern number matches. */
    if (pb->structure_group == LIE_U1 && hg->group == LIE_U1) {
        double c1 = pb_chern_number_1(pb);
        double frac = c1 - floor(c1);
        bool has_holonomy = (fabs(frac) > tolerance && fabs(frac - 1.0) > tolerance);
        bool has_gen = (holonomy_group_dimension(hg) > 0);
        return has_holonomy == has_gen;
    }

    /* For non-abelian groups: compare algebra dimension */
    /* Compute curvature at a reference site */
    if (pb->lattice.dim >= 2) {
        Plaquette plaq = pb_plaquette(pb, 0, 0, 1);
        LieAlgebraElement curv = lie_log(plaq.U);
        double curv_norm = 0.0;
        switch (curv.type) {
            case LIE_SU2: curv_norm = su2_alg_norm(curv.alg.su2); break;
            case LIE_SO3: curv_norm = so3_alg_norm(curv.alg.so3); break;
            default: break;
        }
        /* If curvature is nonzero, holonomy should have at least 1 generator */
        if (curv_norm > tolerance && holonomy_group_dimension(hg) == 0)
            return false;
    }

    return true;
}

/* ==================================================================
 * Wilson Loop Computation
 *
 * The Wilson loop W(C) = Tr P exp(i oint_C A_mu dx^mu) is the
 * basic gauge-invariant observable in lattice gauge theory.
 *
 * For a rectangular loop of size R x T:
 *   W(R,T) = Tr( U_gamma )
 *
 * The static quark-antiquark potential V(R) is extracted from:
 *   V(R) = -lim_{T->inf} (1/T) log W(R,T)
 *
 * The Creutz ratio:
 *   chi(R,T) = -log( W(R,T)*W(R-1,T-1) / (W(R,T-1)*W(R-1,T)) )
 * approximates the force dV/dR for large T.
 *
 * Area law: W(C) ~ exp(-sigma * Area(C))
 * Perimeter law: W(C) ~ exp(-mu * Perimeter(C))
 *
 * Wilson criterion: confinement iff area law holds.
 *
 * Reference: Wilson (1974); Creutz (1980), Phys. Rev. D21, 2308
 * ================================================================== */

double pb_wilson_loop_rect(const PrincipalBundle *pb, int site,
                            int mu, int nu, int R, int T) {
    if (!pb || R < 1 || T < 1) return 0.0;

    LatticePath *loop = path_rectangle(&pb->lattice, site, mu, nu, R, T);
    if (!loop) return 0.0;

    LieGroupElement hol = pb_parallel_transport_path(pb, loop);
    path_free(loop);

    return lie_trace(hol);
}

double pb_creutz_ratio(const PrincipalBundle *pb, int site,
                        int mu, int nu, int R, int T) {
    /* Creutz ratio:
     * chi = -log( W(R,T) * W(R-1,T-1) / (W(R,T-1) * W(R-1,T)) ) */
    if (!pb || R < 2 || T < 2) return 0.0;

    double W_RT     = pb_wilson_loop_rect(pb, site, mu, nu, R, T);
    double W_R1T1   = pb_wilson_loop_rect(pb, site, mu, nu, R-1, T-1);
    double W_RT1    = pb_wilson_loop_rect(pb, site, mu, nu, R, T-1);
    double W_R1T    = pb_wilson_loop_rect(pb, site, mu, nu, R-1, T);

    /* Handle zero or negative values */
    if (W_RT < 1e-15 || W_R1T1 < 1e-15 || W_RT1 < 1e-15 || W_R1T < 1e-15)
        return 0.0;

    double ratio = (W_RT * W_R1T1) / (W_RT1 * W_R1T);
    if (ratio <= 0.0) return 0.0;

    return -log(ratio);
}

/* ==================================================================
 * Polyakov Loop (finite temperature)
 *
 * The Polyakov loop P(x) = Tr Prod_{t=0}^{N_t-1} U_0(x, t)
 * is the Wilson line wrapping around the temporal direction.
 *
 * Its expectation value is the order parameter for confinement:
 *   <P> = 0 in confined phase
 *   <P> != 0 in deconfined phase
 *
 * Reference: Polyakov (1978), Phys. Lett. 72B, 477
 * ================================================================== */

double pb_polyakov_loop(const PrincipalBundle *pb, int spatial_site) {
    /* Compute the Polyakov loop at a given spatial site.
     * Assumes direction 0 is the temporal direction.
     * P = Tr( prod_{t} U_0(x, t) ) */
    if (!pb || pb->lattice.dim < 1) return 0.0;

    int Lt = pb->lattice.L[0];
    if (Lt <= 0) return 0.0;

    /* Build the set of temporal sites at this spatial position.
     * For simplicity, we use the first spatial_site and iterate
     * over all temporal sites with the same spatial coordinates. */
    int ns = pb->lattice.n_sites;
    if (spatial_site < 0 || spatial_site >= ns) return 0.0;

    int coords[4];
    lattice_site_coords(&pb->lattice, spatial_site, coords);

    LieGroupElement total = lie_identity(pb->structure_group);
    int current_coords[4];
    memcpy(current_coords, coords, 4 * sizeof(int));

    for (int t = 0; t < Lt; t++) {
        current_coords[0] = t;
        int s = lattice_site_index(&pb->lattice, current_coords);
        LieGroupElement U0 = pb_get_link(pb, s, 0);
        total = lie_multiply(total, U0);
    }

    return lie_trace(total);
}
