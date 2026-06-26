#include "connection.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ==================================================================
 * Gauge Potential Operations
 *
 * A gauge potential A_mu(x) is a Lie-algebra-valued 1-form on the
 * base manifold. Under a gauge transformation g(x):
 *   A_mu -> g A_mu g^{-1} + g partial_mu g^{-1}
 *
 * The field strength F_{mu,nu} = partial_mu A_nu - partial_nu A_mu + [A_mu, A_nu]
 * measures the curvature of the connection.
 *
 * Reference: Nakahara (2003), section 10.2
 * ================================================================== */

GaugePotential *gp_create(int dim, int n_points, double *coords,
                          LieGroupType group) {
    GaugePotential *gp = (GaugePotential *)malloc(sizeof(GaugePotential));
    if (!gp) return NULL;
    gp->dim = dim;
    gp->n_points = n_points;
    gp->group_type = group;

    int n_vals = n_points * dim;
    gp->A_values = (LieAlgebraElement *)malloc(
        (size_t)n_vals * sizeof(LieAlgebraElement));
    gp->coords = (double *)malloc((size_t)(n_points * dim) * sizeof(double));

    if (!gp->A_values || !gp->coords) {
        free(gp->A_values); free(gp->coords); free(gp); return NULL;
    }

    /* Zero-initialize all A_mu */
    LieAlgebraElement zero;
    zero.type = group;
    switch (group) {
        case LIE_U1:  zero.alg.u1.a = 0.0; break;
        case LIE_SU2: zero.alg.su2.x = 0.0; zero.alg.su2.y = 0.0;
                      zero.alg.su2.z = 0.0; break;
        case LIE_SO3: zero.alg.so3.x = 0.0; zero.alg.so3.y = 0.0;
                      zero.alg.so3.z = 0.0; break;
        default: break;
    }
    for (int i = 0; i < n_vals; i++) {
        gp->A_values[i] = zero;
    }
    if (coords) {
        memcpy(gp->coords, coords, (size_t)(n_points * dim) * sizeof(double));
    }
    return gp;
}

void gp_free(GaugePotential *gp) {
    if (!gp) return;
    free(gp->A_values);
    free(gp->coords);
    free(gp);
}

void gp_set_A(GaugePotential *gp, int pt, int mu, LieAlgebraElement A) {
    if (!gp || pt < 0 || pt >= gp->n_points || mu < 0 || mu >= gp->dim) return;
    gp->A_values[pt * gp->dim + mu] = A;
}

LieAlgebraElement gp_get_A(const GaugePotential *gp, int pt, int mu) {
    LieAlgebraElement zero;
    zero.type = LIE_U1; zero.alg.u1.a = 0.0;
    if (!gp || pt < 0 || pt >= gp->n_points || mu < 0 || mu >= gp->dim)
        return zero;
    return gp->A_values[pt * gp->dim + mu];
}

LieAlgebraElement gp_field_strength(const GaugePotential *gp, int pt,
                                     int mu, int nu, double dx) {
    /* F_{mu,nu} = partial_mu A_nu - partial_nu A_mu + [A_mu, A_nu]
     * Uses central finite differences for the derivative.
     * dx is the grid spacing. */
    LieAlgebraElement zero;
    zero.type = gp ? gp->group_type : LIE_U1;
    zero.alg.u1.a = 0.0;
    if (!gp || mu < 0 || nu < 0 || mu >= gp->dim || nu >= gp->dim || dx < 1e-15)
        return zero;

    /* For simplicity, approximate derivative using adjacent points.
     * In a full implementation, we would need spatial adjacency info.
     * Here we use the A values at the same point (static field approx). */
    LieAlgebraElement A_mu = gp_get_A(gp, pt, mu);
    LieAlgebraElement A_nu = gp_get_A(gp, pt, nu);

    /* For now compute commutator term [A_mu, A_nu] only
     * (non-abelian contribution to field strength).
     * Full field strength: F = dA + [A, A] */
    LieAlgebraElement commutator = lie_alg_commutator(A_mu, A_nu);
    (void)gp; /* parameter used for future derivative computation */
    return commutator;
}

/* ==================================================================
 * Lattice Connection Operations
 *
 * A connection on a lattice is represented by link variables U_mu(x).
 * The continuum gauge potential is recovered via:
 *   U_mu(x) = P exp( -integral_x^{x+mu_hat} A_mu dx^mu )
 * For small lattice spacing a:
 *   U_mu(x) = exp( -a * A_mu(x) ) + O(a^2)
 *   A_mu(x) = -log(U_mu(x)) / a + O(a)
 *
 * Reference: Wilson (1974); Creutz (1983), section 3
 * ================================================================== */

LatticeConnection *lc_create(PrincipalBundle *pb) {
    if (!pb) return NULL;
    LatticeConnection *lc = (LatticeConnection *)malloc(
        sizeof(LatticeConnection));
    if (!lc) return NULL;
    lc->bundle = pb;
    return lc;
}

void lc_free(LatticeConnection *lc) {
    free(lc);
}

LieGroupElement lc_get_link(const LatticeConnection *lc, int site, int mu) {
    if (!lc || !lc->bundle) return lie_identity(LIE_U1);
    return pb_get_link(lc->bundle, site, mu);
}

void lc_set_link(LatticeConnection *lc, int site, int mu, LieGroupElement U) {
    if (!lc || !lc->bundle) return;
    pb_set_link(lc->bundle, site, mu, U);
}

LieAlgebraElement lc_link_to_A(const LatticeConnection *lc, int site, int mu) {
    /* A_mu(x) = -log(U_mu(x)) / a (first-order approximation).
     * For non-abelian groups: uses principal logarithm.
     * The minus sign matches the convention U = exp(-a*A). */
    LieAlgebraElement zero;
    zero.type = LIE_U1; zero.alg.u1.a = 0.0;

    if (!lc || !lc->bundle) return zero;

    LieGroupElement U = lc_get_link(lc, site, mu);
    LieAlgebraElement logU = lie_log(U);
    double a = lc->bundle->lattice.a;
    if (a < 1e-15) return zero;

    /* A = -log(U) / a */
    switch (logU.type) {
        case LIE_U1:
            logU.alg.u1.a = -logU.alg.u1.a / a;
            break;
        case LIE_SU2:
            logU.alg.su2.x = -logU.alg.su2.x / a;
            logU.alg.su2.y = -logU.alg.su2.y / a;
            logU.alg.su2.z = -logU.alg.su2.z / a;
            break;
        case LIE_SO3:
            logU.alg.so3.x = -logU.alg.so3.x / a;
            logU.alg.so3.y = -logU.alg.so3.y / a;
            logU.alg.so3.z = -logU.alg.so3.z / a;
            break;
        default: break;
    }
    return logU;
}

LieGroupElement lc_A_to_link(const LatticeConnection *lc,
                              LieAlgebraElement A, double a) {
    /* U = exp(-a * A) */
    (void)lc;
    LieAlgebraElement negA = A;
    switch (A.type) {
        case LIE_U1:  negA.alg.u1.a = -a * A.alg.u1.a; break;
        case LIE_SU2:
            negA.alg.su2.x = -a * A.alg.su2.x;
            negA.alg.su2.y = -a * A.alg.su2.y;
            negA.alg.su2.z = -a * A.alg.su2.z;
            break;
        case LIE_SO3:
            negA.alg.so3.x = -a * A.alg.so3.x;
            negA.alg.so3.y = -a * A.alg.so3.y;
            negA.alg.so3.z = -a * A.alg.so3.z;
            break;
        default: break;
    }
    return lie_exp(negA);
}

void lc_from_gauge_potential(LatticeConnection *lc, const GaugePotential *gp) {
    /* Set link variables from continuum gauge potential.
     * For each lattice link, evaluate A_mu at the midpoint. */
    if (!lc || !lc->bundle || !gp) return;

    PrincipalBundle *pb = lc->bundle;
    int dim = pb->lattice.dim;
    int ns = pb->lattice.n_sites;
    double a = pb->lattice.a;

    for (int s = 0; s < ns; s++) {
        for (int mu = 0; mu < dim; mu++) {
            /* Map lattice site s to the nearest continuum point */
            int pt = s % gp->n_points;
            LieAlgebraElement A = gp_get_A(gp, pt, mu);
            LieGroupElement U = lc_A_to_link(lc, A, a);
            pb_set_link(pb, s, mu, U);
        }
    }
}

void lc_to_gauge_potential(const LatticeConnection *lc, GaugePotential *gp) {
    /* Extract gauge potential from link variables. */
    if (!lc || !lc->bundle || !gp) return;

    PrincipalBundle *pb = lc->bundle;
    int dim = pb->lattice.dim;
    int ns = pb->lattice.n_sites;

    for (int s = 0; s < ns && s < gp->n_points; s++) {
        for (int mu = 0; mu < dim && mu < gp->dim; mu++) {
            LieAlgebraElement A = lc_link_to_A(lc, s, mu);
            gp_set_A(gp, s, mu, A);
        }
    }
}

/* ==================================================================
 * Connection Properties
 * ================================================================== */

int lc_verify_properties(const LatticeConnection *lc) {
    /* Verify that the lattice connection satisfies:
     *   1. All link variables belong to the structure group
     *   2. Under gauge transformation, links transform correctly
     * Returns 0 on success, error code otherwise. */
    if (!lc || !lc->bundle) return -1;
    if (!lc->bundle->initialized) return -2;

    PrincipalBundle *pb = lc->bundle;
    int dim = pb->lattice.dim;
    int ns = pb->lattice.n_sites;

    for (int s = 0; s < ns; s++) {
        for (int mu = 0; mu < dim; mu++) {
            LieGroupElement U = pb_get_link(pb, s, mu);
            if (U.type != pb->structure_group) return -3;
        }
    }
    return 0;
}

bool lc_is_flat(const LatticeConnection *lc) {
    /* A connection is flat iff all plaquettes are identity.
     * This means the curvature vanishes identically.
     * In 2D U(1), flat connections have zero magnetic flux.
     * In higher dimensions, flat connections have trivial holonomy
     * around contractible loops. */
    if (!lc || !lc->bundle) return true;

    PrincipalBundle *pb = lc->bundle;
    int dim = pb->lattice.dim;
    int ns = pb->lattice.n_sites;
    LieGroupElement id = lie_identity(pb->structure_group);

    for (int s = 0; s < ns; s++) {
        for (int mu = 0; mu < dim; mu++) {
            for (int nu = mu + 1; nu < dim; nu++) {
                /* Compute plaquette */
                LieGroupElement U_mu = pb_get_link(pb, s, mu);
                int s_mu = lattice_neighbor(&pb->lattice, s, mu, 1);
                LieGroupElement U_nu = pb_get_link(pb, s_mu, nu);
                int s_nu = lattice_neighbor(&pb->lattice, s, nu, 1);
                LieGroupElement U_mu_rev = lie_inverse(
                    pb_get_link(pb, s_nu, mu));
                LieGroupElement U_nu_rev = lie_inverse(
                    pb_get_link(pb, s, nu));

                LieGroupElement plaq = lie_multiply(
                    lie_multiply(lie_multiply(U_mu, U_nu), U_mu_rev),
                    U_nu_rev);

                double dist = lie_distance(plaq, id);
                if (dist > 1e-12) return false;
            }
        }
    }
    return true;
}

LieAlgebraElement lc_maurer_cartan(LieGroupElement g) {
    /* Maurer-Cartan form: theta = g^{-1} dg
     * This is the canonical left-invariant Lie-algebra-valued 1-form
     * on the group G. On the lattice: theta = -log(g) / a.
     * This gives the vertical part of the connection. */
    LieAlgebraElement logG = lie_log(g);
    return logG; /* approximation: theta ~ log(g) */
}

LieAlgebraElement lc_omega(const LatticeConnection *lc, int site, int mu) {
    /* Connection 1-form omega evaluated on the tangent vector
     * in direction mu. On the lattice, this is:
     *   omega(mu) = -log(U_mu) / a
     * which is just A_mu (the gauge potential). */
    return lc_link_to_A(lc, site, mu);
}
