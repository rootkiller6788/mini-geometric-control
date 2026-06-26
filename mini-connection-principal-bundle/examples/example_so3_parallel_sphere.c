#include <stdio.h>
#include <math.h>
#include "lie_group.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "principal_bundle.h"
#include "parallel_transport.h"

/* ==================================================================
 * Example 3: SO(3) Frame Bundle - Parallel Transport on S^2
 *
 * Demonstrates:
 *   1. SO(3) principal bundle as the frame bundle of the sphere
 *   2. Parallel transport of tangent vectors along great circles
 *   3. Holonomy of a latitude circle (geometric phase)
 *   4. SU(2)/SO(3) double cover relationship
 *
 * Mathematics: On S^2, parallel transport along a latitude circle
 * at angle phi (from the equator) rotates tangent vectors by:
 *   holonomy angle = 2*pi*sin(phi)   (for a full circle)
 *
 * This is the Gauss-Bonnet theorem in action: the holonomy equals
 * the integral of curvature over the enclosed area.
 *
 * For a closed geodesic triangle with angles (alpha, beta, gamma):
 *   holonomy = alpha + beta + gamma - pi   (Girard's theorem)
 *
 * Reference: Nakahara (2003), section 10.6; O'Neill (1983), Ch.5
 * ================================================================== */

int main(void) {
    printf("=== SO(3) Frame Bundle: Parallel Transport on S^2 ===\n\n");

    /* --- Part 1: SU(2) / SO(3) double cover --- */
    printf("--- Part 1: SU(2)/SO(3) Double Cover ---\n");

    /* A 90-degree rotation around z-axis */
    SU2Element q = su2_make(cos(M_PI/4.0), 0.0, 0.0, sin(M_PI/4.0));
    SO3Element R = su2_to_so3(q);

    printf("SU(2) quaternion: (%.4f, %.4f, %.4f, %.4f)\n", q.a, q.b, q.c, q.d);
    printf("SO(3) trace = %.4f (expected %.4f for 90deg rot)\n",
           so3_trace(R), 1.0 + 2.0*cos(M_PI/2.0));

    /* Apply rotation to y-axis vector */
    double v[3] = {0.0, 1.0, 0.0};
    double out[3];
    so3_apply(R, v, out);
    printf("R * (0,1,0) = (%.4f, %.4f, %.4f) expected (-1, 0, 0) for 90deg z-rot\n",
           out[0], out[1], out[2]);

    /* Roundtrip: SU(2) -> SO(3) -> SU(2) */
    SU2Element q2 = so3_to_su2(R);
    double dist = su2_distance(q, q2);
    printf("SU2 -> SO3 -> SU2 roundtrip distance: %.2e\n", dist);

    /* --- Part 2: Parallel transport on the sphere --- */
    printf("\n--- Part 2: Parallel Transport on S^2 ---\n");

    /* Create a 3D lattice for the frame bundle */
    int L[3] = {8, 8, 1};
    LatticeGeometry geom = lattice_create(3, L, 1.0);

    /* Create SO(3) principal bundle */
    PrincipalBundle *pb = pb_create(LIE_SO3, &geom);

    /* Set nontrivial connections approximating the Levi-Civita
     * connection on S^2. For demonstration, set links to represent
     * rotation around the z-axis proportional to the y-coordinate
     * (simulating the connection on a sphere). */

    int coords[4];
    for (int s = 0; s < geom.n_sites; s++) {
        lattice_site_coords(&geom, s, coords);
        double x = (double)coords[0] - 3.5; /* center at (3.5, 3.5) */
        double y = (double)coords[1] - 3.5;
        double r_val = sqrt(x*x + y*y); (void)r_val;
        double theta = atan2(y, x);

        /* Connection approximates the spin connection on S^2:
         * omega_phi = cos(phi) where phi is the polar angle.
         * For simplicity, set links with rotation angles proportional
         * to the angular coordinate. */
        double angle = 0.1 * theta; /* rotation proportional to azimuth */

        SO3Element Rx = so3_rodrigues(0, 0, 1, angle);
        LieGroupElement U;
        U.type = LIE_SO3;
        U.elem.so3 = Rx;

        /* Set x-direction links */
        pb_set_link(pb, s, 0, U);

        /* y-direction links: identity for demonstration */
        LieGroupElement Uy;
        Uy.type = LIE_SO3;
        Uy.elem.so3 = so3_identity();
        pb_set_link(pb, s, 1, Uy);
    }

    /* Create a circular curve (approximating a latitude circle) */
    double center[3] = {3.5, 3.5, 0.0};
    BaseCurve *circle = curve_circle(center, 2.5, 3, 20);
    printf("Circular curve: %d points\n", circle->n_points);

    /* Parallel transport around the circle */
    LieGroupElement init = lie_identity(LIE_SO3);
    ParallelTransportResult *pt = pt_transport(pb, circle, init);

    if (pt && pt->success) {
        printf("Parallel transport successful: %d steps\n", pt->n_steps);

        LieGroupElement hol = pt->total;
        double ax, ay, az, angle;
        so3_to_axis_angle(hol.elem.so3, &ax, &ay, &az, &angle);
        printf("Holonomy = rotation by %.4f rad (%.1f deg) around axis (%.3f, %.3f, %.3f)\n",
               angle, angle*180.0/M_PI, ax, ay, az);
    }

    /* --- Part 3: Holonomy as geometric (Berry) phase --- */
    printf("\n--- Part 3: Holonomy / Geometric Phase ---\n");

    /* Transport a vector around the circle */
    double v_init[3] = {1.0, 0.0, 0.0};
    double v_final[3] = {0.0, 0.0, 0.0};
    pt_transport_vector(pb, circle, v_init, v_final);

    printf("Initial vector:  (%.4f, %.4f, %.4f)\n", v_init[0], v_init[1], v_init[2]);
    printf("Transported vector: (%.4f, %.4f, %.4f)\n", v_final[0], v_final[1], v_final[2]);

    double dot = v_init[0]*v_final[0] + v_init[1]*v_final[1] + v_init[2]*v_final[2];
    double rot_angle = acos(fmax(-1.0, fmin(1.0, dot)));
    printf("Rotation angle: %.4f rad = %.2f deg\n", rot_angle, rot_angle*180.0/M_PI);

    /* --- Part 4: Sphere holonomy formula --- */
    printf("\n--- Part 4: Theoretical Holonomy Formula ---\n");
    double phi = 30.0 * M_PI / 180.0; /* latitude 30 degrees */
    printf("For latitude phi = 30 degrees:\n");
    printf("  Full circle holonomy = 2*pi*sin(phi) = %.4f rad = %.2f deg\n",
           2.0*M_PI*sin(phi), 360.0*sin(phi));
    printf("  Half circle holonomy  = pi*sin(phi)    = %.4f rad = %.2f deg\n",
           M_PI*sin(phi), 180.0*sin(phi));

    printf("\n(Gauss-Bonnet: holonomy = integral of Gaussian curvature)\n");
    printf("For S^2: K = 1/R^2, area = 4*pi*R^2, total = 4*pi\n");
    printf("For triangle: holonomy = A/R^2 = spherical excess\n");

    pt_result_free(pt);
    curve_free(circle);
    pb_free(pb);
    printf("\n=== Example complete ===\n");
    return 0;
}
