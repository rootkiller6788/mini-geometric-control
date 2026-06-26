#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "principal_bundle.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ==================================================================
 * Lattice Geometry
 *
 * A d-dimensional lattice with sizes L[0..d-1] and spacing a.
 * Sites are labeled by linear index: idx = sum_{mu} coord[mu] * stride[mu].
 * Periodic boundary conditions are used (torus topology).
 *
 * The lattice is the base manifold of the discretized principal bundle.
 * Reference: Montvay & Munster (1994), section 2.1
 * ================================================================== */

LatticeGeometry lattice_create(int dim, const int *L, double a) {
    LatticeGeometry geom;
    geom.dim = dim;
    geom.a = a;
    geom.n_sites = 1;
    for (int mu = 0; mu < dim && mu < 4; mu++) {
        geom.L[mu] = L[mu];
        geom.n_sites *= geom.L[mu];
    }
    for (int mu = dim; mu < 4; mu++) {
        geom.L[mu] = 1;
    }
    geom.n_links = dim * geom.n_sites;
    return geom;
}

int lattice_site_index(const LatticeGeometry *geom, const int *coords) {
    int idx = 0;
    int stride = 1;
    for (int mu = 0; mu < geom->dim; mu++) {
        int c = coords[mu];
        while (c < 0) c += geom->L[mu];
        while (c >= geom->L[mu]) c -= geom->L[mu];
        idx += c * stride;
        stride *= geom->L[mu];
    }
    return idx;
}

void lattice_site_coords(const LatticeGeometry *geom, int idx, int *coords) {
    for (int mu = 0; mu < geom->dim && mu < 4; mu++) {
        int stride = 1;
        for (int nu = 0; nu < mu; nu++) stride *= geom->L[nu];
        coords[mu] = (idx / stride) % geom->L[mu];
    }
    for (int mu = geom->dim; mu < 4; mu++) {
        coords[mu] = 0;
    }
}

int lattice_neighbor(const LatticeGeometry *geom, int site, int mu, int dir) {
    if (mu < 0 || mu >= geom->dim) return site;
    int coords[4];
    lattice_site_coords(geom, site, coords);
    coords[mu] += dir;
    while (coords[mu] < 0) coords[mu] += geom->L[mu];
    while (coords[mu] >= geom->L[mu]) coords[mu] -= geom->L[mu];
    return lattice_site_index(geom, coords);
}

int lattice_link_index(const LatticeGeometry *geom, int site, int mu) {
    return mu * geom->n_sites + site;
}

bool lattice_in_bounds(const LatticeGeometry *geom, const int *coords) {
    for (int mu = 0; mu < geom->dim; mu++) {
        if (coords[mu] < 0 || coords[mu] >= geom->L[mu]) return false;
    }
    return true;
}

/* ==================================================================
 * Principal Bundle Construction
 *
 * A principal G-bundle on a lattice is defined by link variables
 * U_mu(x) in G for each directed link (x, x+mu_hat).
 *
 * The bundle structure is:
 *   - Base: lattice sites
 *   - Fiber: G at each site
 *   - Connection: encoded in link variables
 *
 * Gauge transformations act as:
 *   U_mu(x) -> g(x) * U_mu(x) * g(x+mu)^{-1}
 *
 * The bundle is trivial iff all U_mu(x) can be gauged to identity.
 *
 * Reference: Wilson (1974), Phys. Rev. D10, 2445
 *            Creutz (1983), Quarks, Gluons and Lattices, section 3
 * ================================================================== */

PrincipalBundle *pb_create(LieGroupType group, const LatticeGeometry *geom) {
    PrincipalBundle *pb = (PrincipalBundle *)malloc(sizeof(PrincipalBundle));
    if (!pb) return NULL;

    pb->structure_group = group;
    pb->lattice = *geom;
    pb->initialized = false;

    int n_links = geom->n_links;
    int n_sites = geom->n_sites;
    pb->link_vars = (LieGroupElement *)malloc(
        (size_t)n_links * sizeof(LieGroupElement));
    pb->gauge_transform = (LieGroupElement *)malloc(
        (size_t)n_sites * sizeof(LieGroupElement));

    if (!pb->link_vars || !pb->gauge_transform) {
        free(pb->link_vars);
        free(pb->gauge_transform);
        free(pb);
        return NULL;
    }

    LieGroupElement id = lie_identity(group);
    for (int i = 0; i < n_links; i++) {
        pb->link_vars[i] = id;
    }
    for (int i = 0; i < n_sites; i++) {
        pb->gauge_transform[i] = id;
    }

    pb->initialized = true;
    return pb;
}

void pb_free(PrincipalBundle *pb) {
    if (!pb) return;
    free(pb->link_vars);
    free(pb->gauge_transform);
    free(pb);
}

LieGroupElement pb_get_link(const PrincipalBundle *pb, int site, int mu) {
    if (!pb || site < 0 || site >= pb->lattice.n_sites ||
        mu < 0 || mu >= pb->lattice.dim) {
        return lie_identity(LIE_U1);
    }
    int idx = lattice_link_index(&pb->lattice, site, mu);
    return pb->link_vars[idx];
}

void pb_set_link(PrincipalBundle *pb, int site, int mu, LieGroupElement U) {
    if (!pb || site < 0 || site >= pb->lattice.n_sites ||
        mu < 0 || mu >= pb->lattice.dim) return;
    int idx = lattice_link_index(&pb->lattice, site, mu);
    pb->link_vars[idx] = U;
}

void pb_set_trivial(PrincipalBundle *pb) {
    if (!pb) return;
    LieGroupElement id = lie_identity(pb->structure_group);
    for (int i = 0; i < pb->lattice.n_links; i++) {
        pb->link_vars[i] = id;
    }
}

void pb_copy_links(PrincipalBundle *dest, const PrincipalBundle *src) {
    if (!dest || !src) return;
    if (dest->lattice.n_links != src->lattice.n_links) return;
    size_t bytes = (size_t)src->lattice.n_links * sizeof(LieGroupElement);
    memcpy(dest->link_vars, src->link_vars, bytes);
}

bool pb_is_trivial(const PrincipalBundle *pb) {
    if (!pb) return true;
    LieGroupElement id = lie_identity(pb->structure_group);
    for (int i = 0; i < pb->lattice.n_links; i++) {
        double d = lie_distance(pb->link_vars[i], id);
        if (d > 1e-12) return false;
    }
    return true;
}

/* ==================================================================
 * Transition Functions (Cech cohomology classification)
 *
 * Transition functions g_ij: U_i cap U_j -> G satisfy the cocycle
 * condition: g_ij(x) * g_jk(x) * g_ki(x) = identity.
 *
 * Theorem (Cech): Isomorphism classes of principal G-bundles
 * correspond to H^1(M, G), the first Cech cohomology set.
 * Reference: Kobayashi & Nomizu (1963), Ch.I section 5
 * ================================================================== */

TransitionFunctions *tf_create(int n_patches, LieGroupType group) {
    TransitionFunctions *tf = (TransitionFunctions *)malloc(
        sizeof(TransitionFunctions));
    if (!tf) return NULL;
    tf->n_patches = n_patches;
    tf->n_overlaps = n_patches * n_patches;
    tf->group_type = group;
    tf->transitions = (LieGroupElement *)malloc(
        (size_t)tf->n_overlaps * sizeof(LieGroupElement));
    tf->overlap_i = (int *)malloc((size_t)tf->n_overlaps * sizeof(int));
    tf->overlap_j = (int *)malloc((size_t)tf->n_overlaps * sizeof(int));
    if (!tf->transitions || !tf->overlap_i || !tf->overlap_j) {
        free(tf->transitions); free(tf->overlap_i); free(tf->overlap_j);
        free(tf); return NULL;
    }
    LieGroupElement id = lie_identity(group);
    for (int k = 0; k < tf->n_overlaps; k++) {
        tf->transitions[k] = id;
        tf->overlap_i[k] = -1;
        tf->overlap_j[k] = -1;
    }
    return tf;
}

void tf_free(TransitionFunctions *tf) {
    if (!tf) return;
    free(tf->transitions);
    free(tf->overlap_i);
    free(tf->overlap_j);
    free(tf);
}

void tf_set(TransitionFunctions *tf, int overlap_idx,
            int i, int j, LieGroupElement g) {
    if (!tf || overlap_idx < 0 || overlap_idx >= tf->n_overlaps) return;
    tf->overlap_i[overlap_idx] = i;
    tf->overlap_j[overlap_idx] = j;
    tf->transitions[overlap_idx] = g;
}

bool tf_check_cocycle(const TransitionFunctions *tf) {
    if (!tf || tf->n_patches < 3) return true;
    for (int i = 0; i < tf->n_patches; i++) {
        for (int j = 0; j < tf->n_patches; j++) {
            if (i == j) continue;
            for (int k = 0; k < tf->n_patches; k++) {
                if (i == k || j == k) continue;
                LieGroupElement g_ij = lie_identity(tf->group_type);
                LieGroupElement g_jk = lie_identity(tf->group_type);
                LieGroupElement g_ki = lie_identity(tf->group_type);
                bool f_ij = false, f_jk = false, f_ki = false;
                for (int o = 0; o < tf->n_overlaps; o++) {
                    if (tf->overlap_i[o] == i && tf->overlap_j[o] == j)
                        { g_ij = tf->transitions[o]; f_ij = true; }
                    if (tf->overlap_i[o] == j && tf->overlap_j[o] == k)
                        { g_jk = tf->transitions[o]; f_jk = true; }
                    if (tf->overlap_i[o] == k && tf->overlap_j[o] == i)
                        { g_ki = tf->transitions[o]; f_ki = true; }
                }
                if (f_ij && f_jk && f_ki) {
                    LieGroupElement prod = lie_multiply(
                        lie_multiply(g_ij, g_jk), g_ki);
                    double d = lie_distance(prod,
                        lie_identity(tf->group_type));
                    if (d > 1e-10) return false;
                }
            }
        }
    }
    return true;
}

/* ==================================================================
 * Bundle Reconstruction from Transition Functions
 * Discretizes the clutching construction of principal bundles.
 * Reference: Nash & Sen (1983), Topology and Geometry for Physicists
 * ================================================================== */

void pb_from_transitions(PrincipalBundle *pb, const TransitionFunctions *tf) {
    if (!pb || !tf) return;
    if (pb->structure_group != tf->group_type) return;

    int dim = pb->lattice.dim;
    int ns = pb->lattice.n_sites;
    int np = tf->n_patches;
    if (np < 1) return;

    int *site_patch = (int *)malloc((size_t)ns * sizeof(int));
    if (!site_patch) return;

    int coords[4];
    for (int s = 0; s < ns; s++) {
        lattice_site_coords(&pb->lattice, s, coords);
        int sum = 0;
        for (int mu = 0; mu < dim; mu++) sum += coords[mu];
        site_patch[s] = sum % np;
    }

    LieGroupElement id = lie_identity(pb->structure_group);
    for (int s = 0; s < ns; s++) {
        for (int mu = 0; mu < dim; mu++) {
            int sn = lattice_neighbor(&pb->lattice, s, mu, 1);
            int pi = site_patch[s];
            int pj = site_patch[sn];

            LieGroupElement U = id;
            if (pi != pj) {
                bool found = false;
                for (int o = 0; o < tf->n_overlaps; o++) {
                    if (tf->overlap_i[o] == pi && tf->overlap_j[o] == pj) {
                        U = tf->transitions[o]; found = true; break;
                    }
                }
                if (!found) {
                    for (int o = 0; o < tf->n_overlaps; o++) {
                        if (tf->overlap_i[o] == pj && tf->overlap_j[o] == pi) {
                            U = lie_inverse(tf->transitions[o]); break;
                        }
                    }
                }
            }
            pb_set_link(pb, s, mu, U);
        }
    }
    free(site_patch);
}

/* ==================================================================
 * Gauge Transformation Operations
 *
 * Under g(x) in G at each site:
 *   U_mu(x) -> g(x) * U_mu(x) * g(x+mu)^{-1}
 *
 * Theorem (Elitzur, 1975): Local gauge symmetries cannot be
 * spontaneously broken in lattice gauge theory.
 * Reference: Creutz (1983), section 8
 * ================================================================== */

void pb_gauge_transform(PrincipalBundle *pb) {
    if (!pb) return;
    int dim = pb->lattice.dim;
    int ns = pb->lattice.n_sites;
    for (int s = 0; s < ns; s++) {
        for (int mu = 0; mu < dim; mu++) {
            int sn = lattice_neighbor(&pb->lattice, s, mu, 1);
            LieGroupElement g_s  = pb->gauge_transform[s];
            LieGroupElement g_sn = pb->gauge_transform[sn];
            LieGroupElement U_old = pb_get_link(pb, s, mu);
            LieGroupElement U_new = lie_multiply(
                lie_multiply(g_s, U_old), lie_inverse(g_sn));
            pb_set_link(pb, s, mu, U_new);
        }
    }
}

void pb_set_gauge_at_site(PrincipalBundle *pb, int site, LieGroupElement g) {
    if (!pb || site < 0 || site >= pb->lattice.n_sites) return;
    if (g.type != pb->structure_group) return;
    pb->gauge_transform[site] = g;
}

LieGroupElement pb_get_gauge_at_site(const PrincipalBundle *pb, int site) {
    if (!pb || site < 0 || site >= pb->lattice.n_sites)
        return lie_identity(LIE_U1);
    return pb->gauge_transform[site];
}

void pb_random_gauge_transform(PrincipalBundle *pb, unsigned int seed) {
    if (!pb) return;
    srand(seed);
    int ns = pb->lattice.n_sites;
    LieGroupElement g;
    for (int s = 0; s < ns; s++) {
        switch (pb->structure_group) {
            case LIE_U1: {
                double theta = ((double)rand()/(double)RAND_MAX)*2.0*M_PI;
                g.type = LIE_U1;
                g.elem.u1 = u1_make(theta);
                break;
            }
            case LIE_SU2: {
                double u1 = (double)rand()/(double)RAND_MAX;
                double u2 = (double)rand()/(double)RAND_MAX;
                double u3 = (double)rand()/(double)RAND_MAX;
                double u4 = (double)rand()/(double)RAND_MAX;
                g.type = LIE_SU2;
                g.elem.su2 = su2_make(
                    sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2),
                    sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2),
                    sqrt(-2.0*log(u3))*cos(2.0*M_PI*u4),
                    sqrt(-2.0*log(u3))*sin(2.0*M_PI*u4));
                break;
            }
            case LIE_SO3: {
                double ax=(double)rand()/(double)RAND_MAX-0.5;
                double ay=(double)rand()/(double)RAND_MAX-0.5;
                double az=(double)rand()/(double)RAND_MAX-0.5;
                double angle=((double)rand()/(double)RAND_MAX)*M_PI;
                g.type = LIE_SO3;
                g.elem.so3 = so3_from_axis_angle(ax,ay,az,angle);
                break;
            }
            default:
                g = lie_identity(pb->structure_group);
                break;
        }
        pb->gauge_transform[s] = g;
    }
}

/* ==================================================================
 * Bundle Topology: Check if bundle is nontrivial
 * For U(1) in 2D: compute winding number of link variables
 * ================================================================== */

double pb_winding_number_u1(const PrincipalBundle *pb) {
    /* Compute the U(1) winding number in 2D.
     * This is the magnetic flux through the lattice:
     *   n = (1/2*pi) * sum_{plaquettes} arg(U_{12})
     * For a nontrivial bundle, n != 0.
     * Reference: Creutz (1983), section 10 */
    if (!pb || pb->structure_group != LIE_U1 || pb->lattice.dim < 2)
        return 0.0;

    double total_flux = 0.0;
    int ns = pb->lattice.n_sites;
    int coords[4];

    for (int s = 0; s < ns; s++) {
        lattice_site_coords(&pb->lattice, s, coords);
        /* Only count one plaquette per site in the 0-1 plane */
        LieGroupElement U0 = pb_get_link(pb, s, 0);
        int s1 = lattice_neighbor(&pb->lattice, s, 0, 1);
        LieGroupElement U1 = pb_get_link(pb, s1, 1);
        int s2 = lattice_neighbor(&pb->lattice, s, 1, 1);
        LieGroupElement U0_rev = lie_inverse(pb_get_link(pb, s2, 0));
        LieGroupElement U1_rev = lie_inverse(pb_get_link(pb, s, 1));

        LieGroupElement plaq = lie_multiply(
            lie_multiply(lie_multiply(U0, U1), U0_rev), U1_rev);
        total_flux += u1_log(plaq.elem.u1).a;
    }

    return total_flux / (2.0 * M_PI);
}
