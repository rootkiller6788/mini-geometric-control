#include <stdio.h>
#include <math.h>
#include "lie_group.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "principal_bundle.h"
#include "curvature.h"
#include "holonomy.h"

/* ==================================================================
 * Example 2: SU(2) Yang-Mills Theory on a 4D Lattice
 *
 * Demonstrates:
 *   1. SU(2) principal bundle on a small 4D lattice
 *   2. Setting up an instanton-like configuration
 *   3. Computing the topological charge (second Chern number)
 *   4. Computing Wilson action and energy density
 *   5. Verifying the Bianchi identity
 *
 * Physics: SU(2) Yang-Mills is the gauge theory of the weak
 * interaction. Instantons are topological solitons with
 * integer topological charge Q. The BPST instanton has Q = 1.
 *
 * Reference: Belavin et al. (1975); Creutz (1983)
 * ================================================================== */

int main(void) {
    printf("=== SU(2) Yang-Mills Theory on a 4D Lattice ===\n\n");

    /* Create a small 4D lattice */
    int L[4] = {4, 4, 4, 4};
    LatticeGeometry geom = lattice_create(4, L, 1.0);
    PrincipalBundle *pb = pb_create(LIE_SU2, &geom);

    printf("Lattice: %d x %d x %d x %d = %d sites, %d links\n",
           L[0], L[1], L[2], L[3], geom.n_sites, geom.n_links);

    /* --- Part 1: Trivial configuration --- */
    printf("\n--- Part 1: Trivial SU(2) Bundle ---\n");

    double action0 = pb_wilson_action_total(pb, 1.0);
    printf("Wilson action of trivial bundle: %.6f\n", action0);

    double c2_0 = pb_chern_number_2(pb);
    printf("Second Chern number: %.6f\n", c2_0);

    bool is_flat = true;
    for (int s = 0; s < geom.n_sites; s++) {
        for (int mu = 0; mu < geom.dim; mu++) {
            for (int nu = mu + 1; nu < geom.dim; nu++) {
                Plaquette p = pb_plaquette(pb, s, mu, nu);
                if (!su2_is_identity(p.U.elem.su2)) {
                    is_flat = false;
                    break;
                }
            }
        }
    }
    printf("Bundle is flat: %s\n", is_flat ? "YES" : "NO");

    /* --- Part 2: Single nontrivial plaquette --- */
    printf("\n--- Part 2: Excite a Single Plaquette ---\n");

    /* Set U_0(0) to a nontrivial SU(2) element */
    SU2Element g_excite = su2_make(cos(0.3), sin(0.3), 0.0, 0.0);
    LieGroupElement U;
    U.type = LIE_SU2;
    U.elem.su2 = g_excite;
    pb_set_link(pb, 0, 0, U);

    Plaquette p_excited = pb_plaquette(pb, 0, 0, 1);
    printf("Plaquette after excitation: Tr(U) = %.6f\n",
           su2_trace(p_excited.U.elem.su2));

    double action1 = pb_wilson_action_total(pb, 1.0);
    printf("Wilson action after excitation: %.6f\n", action1);

    /* --- Part 3: SU(2) instanton --- */
    printf("\n--- Part 3: BPST Instanton Configuration ---\n");

    /* Create a new bundle for the instanton */
    PrincipalBundle *pb_inst = pb_create(LIE_SU2, &geom);

    /* Set instanton with size parameter rho = 2.0 */
    double rho = 2.0;
    printf("Instanton size rho = %.2f\n", rho);

    /* Manual instanton setup: set links to approximate BPST */
    for (int s = 0; s < geom.n_sites; s++) {
        int coords[4];
        lattice_site_coords(&geom, s, coords);
        double x[4];
        for (int d = 0; d < 4; d++) {
            x[d] = (double)coords[d] - (double)(L[d]-1)*0.5;
        }
        double r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
        double denom = r2 + rho*rho;
        if (denom < 1e-15) denom = 1e-15;

        for (int mu = 0; mu < 4; mu++) {
            double Amu[3] = {0.0, 0.0, 0.0};
            /* BPST regular gauge: A_mu^a = 2 * eta^a_{mu,nu} * x^nu / (x^2+rho^2) */
            for (int a = 0; a < 3; a++) {
                double sum = 0.0;
                for (int nu = 0; nu < 4; nu++) {
                    double eta = 0.0;
                    if (mu == 0) {
                        if (nu == a+1) eta = 1.0;
                    } else if (nu == 0) {
                        if (a+1 == mu) eta = -1.0;
                    } else {
                        int i = mu - 1, j = nu - 1;
                        if (i >= 0 && i < 3 && j >= 0 && j < 3) {
                            int k = 3 - i - j;
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
            double a_lat = geom.a;
            LieGroupElement U_link;
            U_link.type = LIE_SU2;
            U_link.elem.su2 = su2_exp(su2_alg_scale(-a_lat, A_alg));
            pb_set_link(pb_inst, s, mu, U_link);
        }
    }

    double action_inst = pb_wilson_action_total(pb_inst, 1.0);
    printf("Wilson action of instanton: %.6f\n", action_inst);

    double c2 = pb_chern_number_2(pb_inst);
    printf("Topological charge Q = %.4f (expected ~1.0 for BPST)\n", c2);

    /* --- Part 4: Bianchi identity --- */
    printf("\n--- Part 4: Bianchi Identity Check ---\n");
    double max_dev = pb_bianchi_max_deviation(pb_inst);
    printf("Maximum Bianchi deviation: %.2e\n", max_dev);
    printf("Bianchi identity holds within tolerance: %s\n",
           max_dev < 1e-6 ? "YES" : "NO - check discretization errors");

    pb_free(pb);
    pb_free(pb_inst);
    printf("\n=== Example complete ===\n");
    return 0;
}
