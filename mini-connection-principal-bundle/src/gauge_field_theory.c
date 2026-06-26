#include "lie_group.h"
#include "principal_bundle.h"
#include "connection.h"
#include "curvature.h"
#include "holonomy.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * Gauge Field Theory Applications
 *
 * This file implements physical applications of principal bundle
 * theory to gauge field theories:
 *   1. Electromagnetism as U(1) gauge theory
 *   2. Yang-Mills theory with SU(2) gauge group
 *   3. Monopole configurations (Dirac monopole)
 *   4. Instanton configurations (Belavin-Polyakov-Schwarz-Tyupkin)
 *
 * Reference: Nakahara (2003), Ch.10-11; Weinberg, QFT Vol.2 (1996)
 * ================================================================== */

/* ==================================================================
 * Electromagnetic Field Configuration (U(1) Gauge Theory)
 *
 * Maxwell's equations in differential form notation:
 *   dF = 0  (Bianchi identity / homogeneous Maxwell)
 *   d*F = *J (inhomogeneous Maxwell)
 *
 * In terms of the gauge potential A_mu:
 *   F_{mu,nu} = partial_mu A_nu - partial_nu A_mu
 *
 * Gauge transformation: A_mu -> A_mu + partial_mu lambda(x)
 * (U(1) is abelian, so no commutator term)
 *
 * Reference: Jackson (1999), Classical Electrodynamics, Ch.11
 * ================================================================== */

typedef struct {
    double E[3];     /* Electric field (F_{0i}) */
    double B[3];     /* Magnetic field (F_{ij}, i<j cyclic) */
    double rho;      /* Charge density */
    double J[3];     /* Current density */
} EMField;

EMField em_field_from_gauge(const GaugePotential *gp, int pt, double dx) {
    /* Extract E and B fields from U(1) gauge potential.
     * E_i = F_{0i}, B_i = (1/2) epsilon_{ijk} F_{jk}
     * Using lattice field strength from plaquettes. */
    EMField f;
    f.rho = 0.0;
    for (int i = 0; i < 3; i++) { f.E[i] = 0.0; f.B[i] = 0.0; f.J[i] = 0.0; }

    if (!gp || gp->dim < 4) return f;

    /* Electric field: E_i = F_{0i} */
    for (int i = 1; i < 4 && i < gp->dim; i++) {
        LieAlgebraElement F0i = gp_field_strength(gp, pt, 0, i, dx);
        f.E[i-1] = F0i.alg.u1.a;
    }

    /* Magnetic field: B_1 = F_{23}, B_2 = F_{31}, B_3 = F_{12} */
    if (gp->dim >= 3) {
        LieAlgebraElement F23 = gp_field_strength(gp, pt, 1, 2, dx);
        f.B[0] = F23.alg.u1.a;
    }
    if (gp->dim >= 4) {
        LieAlgebraElement F31 = gp_field_strength(gp, pt, 2, 0, dx);
        f.B[1] = F31.alg.u1.a;
        LieAlgebraElement F12 = gp_field_strength(gp, pt, 0, 1, dx);
        f.B[2] = F12.alg.u1.a;
    }

    return f;
}

double em_action_total(const PrincipalBundle *pb) {
    /* Compute the Maxwell action S = (1/4) integral d^4x F_{mu,nu} F^{mu,nu}
     * On the lattice: sum over plaquettes of |F|^2 * a^4
     * Using the Wilson action with beta = 1/g^2, the continuum limit is:
     *   S_Wilson -> a^{d-4}/(4g^2) integral F^2 for d dimensions.
     * In 4D: S_Wilson -> (1/4g^2) integral F^2 as a->0. */
    if (!pb || pb->structure_group != LIE_U1) return 0.0;

    double total = 0.0;
    int dim = pb->lattice.dim;
    int ns = pb->lattice.n_sites;
    double a = pb->lattice.a;

    for (int s = 0; s < ns; s++) {
        for (int mu = 0; mu < dim; mu++) {
            for (int nu = mu + 1; nu < dim; nu++) {
                LieAlgebraElement F = pb_field_strength(pb, s, mu, nu);
                double F2 = F.alg.u1.a * F.alg.u1.a;
                /* Volume element: a^4 in 4D */
                double vol = 1.0;
                for (int d = 0; d < dim; d++) vol *= a;
                total += 0.25 * F2 * vol;
            }
        }
    }
    return total;
}

/* ==================================================================
 * Dirac Monopole Configuration (U(1))
 *
 * The Dirac monopole is a U(1) gauge field configuration in 3D
 * with nonzero magnetic charge:
 *   g = (1/4*pi) integral_{S^2} B * dS
 *
 * The gauge potential on two patches (Wu-Yang construction):
 *   North patch: A_phi^N = g * (1 - cos(theta)) / (r * sin(theta))
 *   South patch: A_phi^S = -g * (1 + cos(theta)) / (r * sin(theta))
 *
 * Dirac quantization: e*g = n/2, n in Z
 *
 * Reference: Dirac (1931), Proc. Roy. Soc. A133, 60
 *             Wu & Yang (1975), Phys. Rev. D12, 3845
 * ================================================================== */

void em_set_dirac_monopole(PrincipalBundle *pb, double magnetic_charge) {
    /* Set up a Dirac monopole configuration on a 3D lattice.
     * The monopole sits at the center of the lattice.
     * Link variables are set to represent the Wu-Yang potential. */
    if (!pb || pb->structure_group != LIE_U1 || pb->lattice.dim < 3) return;

    int Lx = pb->lattice.L[0];
    int Ly = pb->lattice.L[1];
    int Lz = pb->lattice.L[2];
    int ns = pb->lattice.n_sites;

    /* Center of the lattice */
    double cx = (double)(Lx - 1) * 0.5;
    double cy = (double)(Ly - 1) * 0.5;
    double cz = (double)(Lz - 1) * 0.5;

    int coords[4];
    for (int s = 0; s < ns; s++) {
        lattice_site_coords(&pb->lattice, s, coords);
        double x = (double)coords[0] - cx;
        double y = (double)coords[1] - cy;
        double z = (double)coords[2] - cz;

        double r = sqrt(x*x + y*y + z*z);
        if (r < 0.5) r = 0.5; /* avoid singularity */

        double cos_theta = z / r;
        double sin_theta = sqrt(fmax(0.0, 1.0 - cos_theta*cos_theta));
        if (sin_theta < 1e-15) sin_theta = 1e-15;
        double phi = atan2(y, x);

        /* Wu-Yang potential: A_phi = g * (1 - cos(theta)) / (r * sin(theta)) */
        double A_phi = magnetic_charge * (1.0 - cos_theta) / (r * sin_theta);

        /* Convert spherical components to Cartesian A_x, A_y, A_z */
        double A_x = -A_phi * sin(phi);
        double A_y =  A_phi * cos(phi);
        double A_z = 0.0;

        /* Set link variables: U_mu = exp(-i * a * A_mu) */
        double a_lat = pb->lattice.a;

        if (pb->lattice.dim >= 1) {
            LieGroupElement Ux;
            Ux.type = LIE_U1;
            Ux.elem.u1 = u1_exp(u1_alg_make(-a_lat * A_x));
            pb_set_link(pb, s, 0, Ux);
        }
        if (pb->lattice.dim >= 2) {
            LieGroupElement Uy;
            Uy.type = LIE_U1;
            Uy.elem.u1 = u1_exp(u1_alg_make(-a_lat * A_y));
            pb_set_link(pb, s, 1, Uy);
        }
        if (pb->lattice.dim >= 3) {
            LieGroupElement Uz;
            Uz.type = LIE_U1;
            Uz.elem.u1 = u1_exp(u1_alg_make(-a_lat * A_z));
            pb_set_link(pb, s, 2, Uz);
        }
    }
}

/* ==================================================================
 * SU(2) Yang-Mills Configuration
 * Yang-Mills theory generalizes Maxwell theory to non-abelian groups.
 * F_{mu,nu}^a = partial_mu A_nu^a - partial_nu A_mu^a + g f^{abc} A_mu^b A_nu^c
 * Reference: Yang & Mills (1954), Phys. Rev. 96, 191
 * ================================================================== */

typedef struct {
    double E[3][3];
    double B[3][3];
} YMField;

YMField ym_field_from_gauge(const PrincipalBundle *pb, int site) {
    YMField f;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            f.E[i][j] = 0.0;
            f.B[i][j] = 0.0;
        }
    if (!pb || pb->structure_group != LIE_SU2 || pb->lattice.dim < 4) return f;

    LieAlgebraElement F01 = pb_field_strength(pb, site, 0, 1);
    LieAlgebraElement F02 = pb_field_strength(pb, site, 0, 2);
    LieAlgebraElement F03 = pb_field_strength(pb, site, 0, 3);

    f.E[0][0] = F01.alg.su2.x; f.E[0][1] = F01.alg.su2.y; f.E[0][2] = F01.alg.su2.z;
    f.E[1][0] = F02.alg.su2.x; f.E[1][1] = F02.alg.su2.y; f.E[1][2] = F02.alg.su2.z;
    f.E[2][0] = F03.alg.su2.x; f.E[2][1] = F03.alg.su2.y; f.E[2][2] = F03.alg.su2.z;

    LieAlgebraElement F23 = pb_field_strength(pb, site, 2, 3);
    LieAlgebraElement F31 = pb_field_strength(pb, site, 3, 1);
    LieAlgebraElement F12 = pb_field_strength(pb, site, 1, 2);

    f.B[0][0] = F23.alg.su2.x; f.B[0][1] = F23.alg.su2.y; f.B[0][2] = F23.alg.su2.z;
    f.B[1][0] = F31.alg.su2.x; f.B[1][1] = F31.alg.su2.y; f.B[1][2] = F31.alg.su2.z;
    f.B[2][0] = F12.alg.su2.x; f.B[2][1] = F12.alg.su2.y; f.B[2][2] = F12.alg.su2.z;

    return f;
}

double ym_energy_density(const YMField *f) {
    if (!f) return 0.0;
    double energy = 0.0;
    for (int i = 0; i < 3; i++)
        for (int a = 0; a < 3; a++)
            energy += 0.5 * (f->E[i][a]*f->E[i][a] + f->B[i][a]*f->B[i][a]);
    return energy;
}

/* ==================================================================
 * SU(2) Instanton (BPST) Configuration
 *
 * The BPST instanton is a classical solution of SU(2) Yang-Mills
 * in 4D Euclidean space with topological charge Q = 1.
 *
 * Regular gauge: A_mu = (2/g) * eta^a_{mu,nu} * x^nu / (x^2 + rho^2)
 * where eta^a_{mu,nu} are the t Hooft symbols and rho is size.
 *
 * Reference: Belavin, Polyakov, Schwarz, Tyupkin (1975), Phys. Lett. 59B, 85
 * ================================================================== */

void ym_set_instanton(PrincipalBundle *pb, double rho) {
    if (!pb || pb->structure_group != LIE_SU2 || pb->lattice.dim < 4) return;

    int L0 = pb->lattice.L[0], L1 = pb->lattice.L[1];
    int L2 = pb->lattice.L[2], L3 = pb->lattice.L[3];
    int ns = pb->lattice.n_sites;
    double a = pb->lattice.a;

    double cx = (double)(L0-1)*0.5, cy = (double)(L1-1)*0.5;
    double cz = (double)(L2-1)*0.5, ct = (double)(L3-1)*0.5;

    int coords[4];
    for (int s = 0; s < ns; s++) {
        lattice_site_coords(&pb->lattice, s, coords);
        double x0 = (double)coords[0] - cx;
        double x1 = (double)coords[1] - cy;
        double x2 = (double)coords[2] - cz;
        double x3 = (double)coords[3] - ct;

        double r2 = x0*x0 + x1*x1 + x2*x2 + x3*x3;
        double denom = r2 + rho*rho;
        if (denom < 1e-15) denom = 1e-15;

        for (int mu = 0; mu < 4 && mu < pb->lattice.dim; mu++) {
            double Amu[3] = {0.0, 0.0, 0.0};
            double x[4] = {x0, x1, x2, x3};

            for (int a = 0; a < 3; a++) {
                double sum = 0.0;
                for (int nu = 0; nu < 4; nu++) {
                    double eta = 0.0;
                    /* t Hooft eta symbols: eta^a_{0,i} = delta^{a,i} */
                    if (mu == 0 && nu == a+1) eta = 1.0;
                    else if (mu == a+1 && nu == 0) eta = -1.0;
                    else {
                        int i = mu - 1, j = nu - 1;
                        if (i >= 0 && i < 3 && j >= 0 && j < 3 && i != j) {
                            int k = 3 - i - j; /* the third index */
                            if (a == k) eta = 1.0;
                            else if (a == j) eta = -1.0;
                        }
                    }
                    sum += eta * x[nu];
                }
                Amu[a] = 2.0 * sum / denom;
            }

            SU2Algebra A_alg;
            A_alg.x = Amu[0]; A_alg.y = Amu[1]; A_alg.z = Amu[2];
            LieGroupElement U;
            U.type = LIE_SU2;
            U.elem.su2 = su2_exp(su2_alg_scale(-a, A_alg));
            pb_set_link(pb, s, mu, U);
        }
    }
}

/* ==================================================================
 * Magnetic Monopole Mass (Bogomolny bound)
 *
 * In the Georgi-Glashow model, SU(2) -> U(1) symmetry breaking
 * produces magnetic monopoles. The Bogomolny bound:
 *   M >= (4*pi/g) * v = g_monopole * v
 *
 * Reference: t Hooft (1974), Nucl. Phys. B79, 276
 *             Polyakov (1974), JETP Lett. 20, 194
 * ================================================================== */

double ym_monopole_mass_estimate(double gauge_coupling, double higgs_vev) {
    return (4.0 * M_PI / gauge_coupling) * higgs_vev;
}

/* ==================================================================
 * Gauge-Invariant Observables
 *
 * Wilson loops, Polyakov loops, topological charge, Creutz ratios.
 * These are the fundamental observables in lattice gauge theory.
 * ================================================================== */

double gf_string_tension_estimate(const PrincipalBundle *pb, int max_R) {
    if (!pb || pb->lattice.dim < 2 || max_R < 3) return 0.0;

    double sum = 0.0;
    int count = 0;
    int step = pb->lattice.n_sites / 10 + 1;

    for (int R = 2; R <= max_R && R < pb->lattice.L[0] / 2; R++) {
        for (int s = 0; s < pb->lattice.n_sites; s += step) {
            double chi = pb_creutz_ratio(pb, s, 0, 1, R, R);
            if (chi > 0.0) {
                sum += chi;
                count++;
            }
        }
    }
    return (count > 0) ? sum / (double)count : 0.0;
}

/* ==================================================================
 * U(1) Confinement Check (compact QED in 3D/4D)
 *
 * Compact U(1) lattice gauge theory exhibits confinement at
 * strong coupling due to monopole condensation (Polyakov 1977).
 *
 * In 3D compact QED, confinement is permanent (no phase transition).
 * In 4D compact QED, there is a confinement-deconfinement phase
 * transition at critical coupling beta_c ~ 1.0.
 *
 * Check: compute Wilson loop and test area vs perimeter law.
 * ================================================================== */

bool gf_is_confining_area_law(const PrincipalBundle *pb, int R_max) {
    /* Test for area law behavior: W(R,T) ~ exp(-sigma * R * T)
     * Check if log(W(R,R)) is linear in R^2 for large R. */
    if (!pb || R_max < 3) return false;

    double prev_logW = 0.0;
    double prev_R2 = 0.0;
    double slope_sum = 0.0;
    int slope_count = 0;

    for (int R = 2; R <= R_max; R++) {
        double avg_logW = 0.0;
        int count = 0;
        int step = pb->lattice.n_sites / 5 + 1;

        for (int s = 0; s < pb->lattice.n_sites; s += step) {
            double W = pb_wilson_loop_rect(pb, s, 0, 1, R, R);
            if (W > 1e-15) {
                avg_logW += log(W);
                count++;
            }
        }

        if (count > 0) {
            avg_logW /= (double)count;
            double R2 = (double)(R * R);
            if (prev_R2 > 0.0) {
                double slope = (avg_logW - prev_logW) / (R2 - prev_R2);
                slope_sum += slope;
                slope_count++;
            }
            prev_logW = avg_logW;
            prev_R2 = R2;
        }
    }

    if (slope_count < 2) return false;
    double avg_slope = slope_sum / (double)slope_count;
    /* Area law: log(W) ~ -sigma * Area, so slope < -epsilon */
    return avg_slope < -0.01;
}
