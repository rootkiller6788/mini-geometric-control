#include "curvature.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * Curvature / Plaquette Computation
 *
 * The curvature of a connection is measured on the lattice by the
 * plaquette (Wilson loop around an elementary square):
 *
 *   U_{mu,nu}(x) = U_mu(x) * U_nu(x+mu) * U_mu(x+nu)^{-1} * U_nu(x)^{-1}
 *
 * This is the path-ordered exponential of the field strength:
 *   U_{mu,nu}(x) = P exp( i * a^2 * F_{mu,nu}(x) ) + O(a^3)
 *
 * For U(1): U_{mu,nu} = exp( i * a^2 * F_{mu,nu} )
 *   => F_{mu,nu} = -i * log(U_{mu,nu}) / a^2
 *
 * Cartan structure equation: Omega = d omega + (1/2)[omega, omega]
 * In components: F_{mu,nu} = partial_mu A_nu - partial_nu A_mu + [A_mu, A_nu]
 *
 * Reference: Wilson (1974); Nakahara (2003), section 10.3
 * ================================================================== */

Plaquette pb_plaquette(const PrincipalBundle *pb, int site, int mu, int nu) {
    Plaquette plaq;
    plaq.site = site;
    plaq.mu = mu;
    plaq.nu = nu;
    plaq.U = lie_identity(pb ? pb->structure_group : LIE_U1);

    if (!pb || mu < 0 || nu < 0 || mu >= pb->lattice.dim ||
        nu >= pb->lattice.dim || mu == nu) return plaq;

    /* U_{mu,nu}(x) = U_mu(x) * U_nu(x+mu) * U_mu(x+nu)^{-1} * U_nu(x)^{-1} */
    LieGroupElement U_mu = pb_get_link(pb, site, mu);
    int s_mu = lattice_neighbor(&pb->lattice, site, mu, 1);
    LieGroupElement U_nu = pb_get_link(pb, s_mu, nu);
    int s_nu = lattice_neighbor(&pb->lattice, site, nu, 1);
    LieGroupElement U_mu_rev = lie_inverse(pb_get_link(pb, s_nu, mu));
    LieGroupElement U_nu_rev = lie_inverse(pb_get_link(pb, site, nu));

    plaq.U = lie_multiply(
        lie_multiply(lie_multiply(U_mu, U_nu), U_mu_rev), U_nu_rev);
    return plaq;
}

LieAlgebraElement pb_field_strength(const PrincipalBundle *pb, int site,
                                     int mu, int nu) {
    /* F_{mu,nu} = -log(U_{mu,nu}) / a^2 (principal value) */
    LieAlgebraElement zero;
    zero.type = LIE_U1; zero.alg.u1.a = 0.0;

    if (!pb || mu == nu || mu >= pb->lattice.dim || nu >= pb->lattice.dim)
        return zero;

    Plaquette plaq = pb_plaquette(pb, site, mu, nu);
    LieAlgebraElement F = lie_log(plaq.U);
    double a2 = pb->lattice.a * pb->lattice.a;
    if (a2 < 1e-30) return zero;

    /* Convention: U = exp(i * a^2 * F) for U(1)
     *   => F = -i * log(U) / a^2
     * For SU(2)/SO(3): U = exp(i * a^2 * F) similarly */
    switch (F.type) {
        case LIE_U1:
            F.alg.u1.a = -F.alg.u1.a / a2;
            break;
        case LIE_SU2:
            F.alg.su2.x = -F.alg.su2.x / a2;
            F.alg.su2.y = -F.alg.su2.y / a2;
            F.alg.su2.z = -F.alg.su2.z / a2;
            break;
        case LIE_SO3:
            F.alg.so3.x = -F.alg.so3.x / a2;
            F.alg.so3.y = -F.alg.so3.y / a2;
            F.alg.so3.z = -F.alg.so3.z / a2;
            break;
        default: break;
    }
    return F;
}

Plaquette *pb_all_plaquettes(const PrincipalBundle *pb, int *n_plaq) {
    if (!pb || !n_plaq) return NULL;
    int dim = pb->lattice.dim;
    int ns = pb->lattice.n_sites;
    int np = ns * dim * (dim - 1) / 2;
    *n_plaq = np;

    Plaquette *plaqs = (Plaquette *)malloc((size_t)np * sizeof(Plaquette));
    if (!plaqs) return NULL;

    int idx = 0;
    for (int s = 0; s < ns; s++) {
        for (int mu = 0; mu < dim; mu++) {
            for (int nu = mu + 1; nu < dim; nu++) {
                plaqs[idx] = pb_plaquette(pb, s, mu, nu);
                idx++;
            }
        }
    }
    return plaqs;
}

void pb_plaquettes_free(Plaquette *plaqs) {
    free(plaqs);
}

/* ==================================================================
 * Wilson Action
 *
 * S = beta * sum_{x, mu<nu} (1 - (1/N) Re Tr U_{mu,nu}(x))
 *
 * For SU(N): beta = 2N / g^2 (g is the gauge coupling)
 * For U(1):  beta = 1 / g^2
 *
 * Theorem (Wilson, 1974): As a -> 0,
 *   S -> 1/(4g^2) integral d^4x Tr(F_{mu,nu} F^{mu,nu})
 *
 * Reference: Wilson (1974); Creutz (1983), section 3
 * ================================================================== */

double pb_wilson_action_plaq(const PrincipalBundle *pb, int site,
                              int mu, int nu, double beta) {
    if (!pb) return 0.0;
    Plaquette plaq = pb_plaquette(pb, site, mu, nu);
    double tr = lie_trace(plaq.U);
    double N = (double)lie_dimension(pb->structure_group);
    if (N < 1.0) N = 1.0;
    if (pb->structure_group == LIE_U1) N = 1.0; /* U(1) trace is real part */

    /* For U(1): Re Tr U = u1_trace(U) / 2 */
    if (pb->structure_group == LIE_U1) {
        tr = u1_trace(plaq.U.elem.u1) / 2.0; /* Actually just cos(theta) */
    }
    /* For SU(2): tr = su2_trace(U) = 2*a */
    if (pb->structure_group == LIE_SU2) {
        tr = su2_trace(plaq.U.elem.su2);
        N = 2.0;
    }
    /* For SO(3): tr = so3_trace(U) */
    if (pb->structure_group == LIE_SO3) {
        tr = so3_trace(plaq.U.elem.so3);
        N = 3.0;
    }

    return beta * (1.0 - tr / N);
}

double pb_wilson_action_total(const PrincipalBundle *pb, double beta) {
    if (!pb) return 0.0;
    int dim = pb->lattice.dim;
    int ns = pb->lattice.n_sites;
    double total = 0.0;

    for (int s = 0; s < ns; s++) {
        for (int mu = 0; mu < dim; mu++) {
            for (int nu = mu + 1; nu < dim; nu++) {
                total += pb_wilson_action_plaq(pb, s, mu, nu, beta);
            }
        }
    }
    return total;
}

double pb_topological_charge_density(const PrincipalBundle *pb, int site) {
    /* Topological charge density q(x) = (1/32*pi^2) * epsilon_{mu,nu,rho,sigma}
     *   * Tr( F_{mu,nu}(x) * F_{rho,sigma}(x) )
     * For 4D SU(2) gauge theory only.
     * On the lattice: approximated from plaquettes at site x. */
    if (!pb || pb->lattice.dim < 4) return 0.0;
    if (pb->structure_group != LIE_SU2) return 0.0;

    /* Simplified: use plaquette in 0-1 and 2-3 planes */
    LieAlgebraElement F01 = pb_field_strength(pb, site, 0, 1);
    LieAlgebraElement F23 = pb_field_strength(pb, site, 2, 3);

    /* q ~ F_{mu,nu} tilde{F}^{mu,nu} ~ epsilon^{mu,nu,rho,sigma} Tr(F_{mu,nu} F_{rho,sigma}) */
    double q = su2_alg_inner(F01.alg.su2, F23.alg.su2);
    return q / (32.0 * M_PI * M_PI);
}

/* ==================================================================
 * Bianchi Identity
 *
 * D Omega = d Omega + [omega, Omega] = 0
 * In components: D_mu F_{nu,rho} + D_nu F_{rho,mu} + D_rho F_{mu,nu} = 0
 *
 * On the lattice, the Bianchi identity corresponds to the fact that
 * the product of six plaquettes around a 3-cube equals identity.
 * This is a geometric identity, not a dynamical equation.
 *
 * Theorem (Bianchi, 1902): The exterior covariant derivative of
 * the curvature 2-form vanishes identically.
 *
 * Reference: Kobayashi & Nomizu (1963), Ch.II section 5
 * ================================================================== */

bool pb_check_bianchi(const PrincipalBundle *pb, int site,
                      int mu, int nu, int rho, double tolerance) {
    /* Verify D_mu F_{nu,rho} + D_nu F_{rho,mu} + D_rho F_{mu,nu} = 0.
     * On the lattice: construct the 6 plaquettes of a 3-cube and
     * check that their oriented product equals identity.
     *
     * The cube directions are: mu, nu, rho.
     * The six faces are: (mu,nu), (nu,rho), (rho,mu) and their
     * shifted counterparts at x+rho, x+mu, x+nu.
     */
    if (!pb || mu == nu || nu == rho || rho == mu) return true;
    int dim = pb->lattice.dim;
    if (mu >= dim || nu >= dim || rho >= dim) return true;

    /* Three pairs of opposite faces */
    Plaquette p_munu  = pb_plaquette(pb, site, mu, nu);
    int s_rho = lattice_neighbor(&pb->lattice, site, rho, 1);
    Plaquette p_munu_r = pb_plaquette(pb, s_rho, mu, nu);

    Plaquette p_nu_rho  = pb_plaquette(pb, site, nu, rho);
    int s_mu = lattice_neighbor(&pb->lattice, site, mu, 1);
    Plaquette p_nu_rho_m = pb_plaquette(pb, s_mu, nu, rho);

    Plaquette p_rho_mu  = pb_plaquette(pb, site, rho, mu);
    int s_nu = lattice_neighbor(&pb->lattice, site, nu, 1);
    Plaquette p_rho_mu_n = pb_plaquette(pb, s_nu, rho, mu);

    /* Bianchi: product of oriented faces = identity
     * U_{mu,nu}(x) * U_{nu,rho}(x+mu) * U_{rho,mu}(x+nu) *
     * U_{mu,nu}(x+rho)^{-1} * U_{nu,rho}(x)^{-1} * U_{rho,mu}(x+mu)^{-1} = I */
    LieGroupElement prod = lie_multiply(
        lie_multiply(
            lie_multiply(p_munu.U, p_nu_rho_m.U),
            lie_multiply(p_rho_mu_n.U, lie_inverse(p_munu_r.U))),
        lie_multiply(lie_inverse(p_nu_rho.U), lie_inverse(p_rho_mu.U)));

    LieGroupElement id = lie_identity(pb->structure_group);
    double dist = lie_distance(prod, id);
    return dist < tolerance;
}

double pb_bianchi_max_deviation(const PrincipalBundle *pb) {
    /* Compute maximum deviation from Bianchi identity across all 3-cubes.
     * Returns the maximum distance from identity among all 3-cube products. */
    if (!pb) return 0.0;
    int dim = pb->lattice.dim;
    if (dim < 3) return 0.0;
    int ns = pb->lattice.n_sites;

    double max_dev = 0.0;
    for (int s = 0; s < ns; s++) {
        for (int mu = 0; mu < dim; mu++) {
            for (int nu = mu + 1; nu < dim; nu++) {
                for (int rho = nu + 1; rho < dim; rho++) {
                    /* Compute the cube product */
                    Plaquette p_munu = pb_plaquette(pb, s, mu, nu);
                    int s_rho = lattice_neighbor(&pb->lattice, s, rho, 1);
                    Plaquette p_munu_r = pb_plaquette(pb, s_rho, mu, nu);

                    Plaquette p_nu_rho = pb_plaquette(pb, s, nu, rho);
                    int s_mu = lattice_neighbor(&pb->lattice, s, mu, 1);
                    Plaquette p_nu_rho_m = pb_plaquette(pb, s_mu, nu, rho);

                    Plaquette p_rho_mu = pb_plaquette(pb, s, rho, mu);
                    int s_nu = lattice_neighbor(&pb->lattice, s, nu, 1);
                    Plaquette p_rho_mu_n = pb_plaquette(pb, s_nu, rho, mu);

                    LieGroupElement prod = lie_multiply(
                        p_munu.U,
                        lie_multiply(p_nu_rho_m.U,
                            lie_multiply(p_rho_mu_n.U,
                                lie_multiply(lie_inverse(p_munu_r.U),
                                    lie_multiply(lie_inverse(p_nu_rho.U),
                                        lie_inverse(p_rho_mu.U))))));

                    double dev = lie_distance(prod,
                        lie_identity(pb->structure_group));
                    if (dev > max_dev) max_dev = dev;
                }
            }
        }
    }
    return max_dev;
}

/* ==================================================================
 * Chern Numbers / Topological Invariants
 *
 * First Chern number (U(1) in 2D):
 *   c_1 = (1/2*pi) integral d^2x F_{12}(x)
 *   On lattice: c_1 = (1/2*pi) sum_x arg(U_{12}(x))
 *   c_1 is always an integer.
 *
 * Second Chern number / Instanton number (SU(2) in 4D):
 *   c_2 = (1/8*pi^2) integral d^4x Tr(F wedge F)
 *   On lattice: c_2 = (1/16*pi^2) sum epsilon_{mu,nu,rho,sigma} Tr(U_{mu,nu} U_{rho,sigma})
 *
 * Chern-Weil Theorem: Chern numbers are integers independent of
 * the connection. They characterize the bundle topology.
 *
 * Reference: Nakahara (2003), section 11.2; Eguchi, Gilkey & Hanson (1980)
 * ================================================================== */

double pb_chern_number_1(const PrincipalBundle *pb) {
    /* First Chern number for U(1) bundles in 2D.
     * c_1 = (1/2*pi) sum_x F_{12}(x) * a^2
     * On lattice: sum of (1/2*pi) * arg(U_{12}(x)) */
    if (!pb || pb->structure_group != LIE_U1 || pb->lattice.dim < 2)
        return 0.0;

    double total = 0.0;
    int ns = pb->lattice.n_sites;

    for (int s = 0; s < ns; s++) {
        Plaquette plaq = pb_plaquette(pb, s, 0, 1);
        double theta = plaq.U.elem.u1.theta;
        /* Principal value in (-pi, pi] */
        if (theta > M_PI) theta -= 2.0 * M_PI;
        total += theta;
    }

    return total / (2.0 * M_PI);
}

double pb_chern_number_2(const PrincipalBundle *pb) {
    /* Second Chern number for SU(2) bundles in 4D.
     * This is the instanton number (topological charge).
     * Returns an integer for smooth field configurations. */
    if (!pb || pb->structure_group != LIE_SU2 || pb->lattice.dim < 4)
        return 0.0;

    double total = 0.0;
    int ns = pb->lattice.n_sites;

    for (int s = 0; s < ns; s++) {
        /* Compute sum over epsilon tensor of Tr(U_{mu,nu} U_{rho,sigma}) */
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                if (mu == nu) continue;
                for (int rho = 0; rho < 4; rho++) {
                    if (rho == mu || rho == nu) continue;
                    int sigma = 6 - mu - nu - rho; /* 0+1+2+3 = 6 */
                    if (sigma == mu || sigma == nu || sigma == rho) continue;

                    /* epsilon_{mu,nu,rho,sigma}: +1 for even perm, -1 for odd */
                    int perm[4] = {mu, nu, rho, sigma};
                    int inv = 0;
                    for (int i = 0; i < 4; i++)
                        for (int j = i + 1; j < 4; j++)
                            if (perm[i] > perm[j]) inv++;
                    int eps = (inv % 2 == 0) ? 1 : -1;

                    Plaquette p_munu = pb_plaquette(pb, s, mu, nu);
                    Plaquette p_rho_sigma = pb_plaquette(pb, s, rho, sigma);

                    double tr1 = su2_trace(p_munu.U.elem.su2);
                    double tr2 = su2_trace(p_rho_sigma.U.elem.su2);
                    total += eps * tr1 * tr2;
                }
            }
        }
    }

    /* Normalize by 16*pi^2 (each ordered quadruple counted 1/4 factor) */
    return total / (64.0 * M_PI * M_PI);
}
