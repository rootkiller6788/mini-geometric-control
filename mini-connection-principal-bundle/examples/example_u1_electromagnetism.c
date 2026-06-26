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
 * Example 1: U(1) Gauge Theory - Electromagnetism on a Lattice
 *
 * Demonstrates:
 *   1. Creating a U(1) principal bundle on a 2D lattice
 *   2. Setting up a constant magnetic field (homogeneous B)
 *   3. Computing Wilson loops and extracting field strength
 *   4. Computing the first Chern number (magnetic flux quantization)
 *   5. Testing gauge invariance of Wilson loops
 *
 * Physics: A constant magnetic field B in 2D gives:
 *   F_{12} = B (constant)
 *   U_{12}(x) = exp(i * a^2 * B) for each plaquette
 *   Total flux = B * Area = n * (2*pi) (quantized by Dirac)
 *
 * Reference: Wilson (1974); Creutz (1983)
 * ================================================================== */

int main(void) {
    printf("=== U(1) Electromagnetism on a Lattice ===\n\n");

    /* Create a 2D lattice of size 8x8 with spacing a=1.0 */
    int L[2] = {8, 8};
    LatticeGeometry geom = lattice_create(2, L, 1.0);
    PrincipalBundle *pb = pb_create(LIE_U1, &geom);

    printf("Lattice: %d x %d, %d sites, %d links\n",
           geom.L[0], geom.L[1], geom.n_sites, geom.n_links);

    /* --- Part 1: Constant magnetic field --- */
    printf("\n--- Part 1: Constant Magnetic Field ---\n");

    /* Set each plaquette to have flux Phi = B * a^2 */
    double B = 0.3; /* magnetic field strength */
    double flux_per_plaq = B * geom.a * geom.a;

    /* To create constant flux, set link variables U_1(x) to 1,
     * and U_2(x) = exp(i * flux_per_plaq * x_1).
     * This gives U_{12}(x) = U_1(x)*U_2(x+1)*U_1(x+2)^{-1}*U_2(x)^{-1}
     *                     = 1 * exp(i*flux*(x_1+1)) * 1 * exp(-i*flux*x_1)
     *                     = exp(i*flux) */

    int coords[4];
    for (int s = 0; s < geom.n_sites; s++) {
        lattice_site_coords(&geom, s, coords);
        int x1 = coords[0];

        /* U_1(x) = 1 */
        LieGroupElement U1;
        U1.type = LIE_U1;
        U1.elem.u1 = u1_identity();
        pb_set_link(pb, s, 0, U1);

        /* U_2(x) = exp(i * flux * x1) */
        LieGroupElement U2;
        U2.type = LIE_U1;
        U2.elem.u1 = u1_make(flux_per_plaq * (double)x1);
        pb_set_link(pb, s, 1, U2);
    }

    /* Verify constant flux */
    Plaquette p0 = pb_plaquette(pb, 0, 0, 1);
    printf("Plaquette at origin: theta = %.6f (expected %.6f)\n",
           p0.U.elem.u1.theta, flux_per_plaq);

    LieAlgebraElement F = pb_field_strength(pb, 0, 0, 1);
    printf("Field strength F_{12}: %.6f (expected %.6f)\n",
           F.alg.u1.a, B);

    /* --- Part 2: Wilson loops --- */
    printf("\n--- Part 2: Wilson Loops (Area Law) ---\n");

    printf("Rectangular Wilson loops W(R,R):\n");
    for (int R = 1; R <= 4; R++) {
        double W = pb_wilson_loop_rect(pb, 0, 0, 1, R, R);
        double area = (double)(R * R);
        double expected_W = 2.0 * cos(flux_per_plaq * area);
        printf("  W(%d,%d) = %.6f (expected %.6f)\n", R, R, W, expected_W);
    }

    /* Creutz ratios */
    printf("\nCreutz ratios (approximate force):\n");
    for (int R = 2; R <= 4; R++) {
        double chi = pb_creutz_ratio(pb, 0, 0, 1, R, R);
        printf("  chi(%d,%d) = %.6f\n", R, R, chi);
    }

    /* --- Part 3: Topological charge --- */
    printf("\n--- Part 3: Topology (First Chern Number) ---\n");
    double c1 = pb_chern_number_1(pb);
    printf("First Chern number c_1 = %.6f\n", c1);
    printf("Total magnetic flux = %.6f * (2*pi)\n", c1);
    printf("(For constant B=%.3f on %dx%d lattice, expected flux = B*Area/(2pi) = %.3f)\n",
           B, L[0], L[1], B * L[0] * L[1] / (2.0 * M_PI));

    /* --- Part 4: Gauge invariance --- */
    printf("\n--- Part 4: Gauge Invariance ---\n");
    double W_before = pb_wilson_loop_rect(pb, 0, 0, 1, 2, 2);

    /* Apply random gauge transformation */
    pb_random_gauge_transform(pb, 12345);
    pb_gauge_transform(pb);

    double W_after = pb_wilson_loop_rect(pb, 0, 0, 1, 2, 2);
    printf("Wilson loop before gauge: %.6f\n", W_before);
    printf("Wilson loop after gauge:  %.6f\n", W_after);
    printf("Gauge invariance: %s\n",
           fabs(W_before - W_after) < 1e-10 ? "PRESERVED" : "VIOLATED");

    pb_free(pb);
    printf("\n=== Example complete ===\n");
    return 0;
}
