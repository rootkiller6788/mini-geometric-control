#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nonlinear_system.h"
#include "lie_derivative.h"
#include "lie_bracket.h"
#include "feedback_linearization.h"

int main(void) {
    printf("========================================\n");
    printf("Example: Ball and Beam System\n");
    printf("Full-State Feedback Linearization\n");
    printf("========================================\n\n");

    double ball_mass = 0.11;
    double grav = 9.81;
    double ball_radius = 0.015;
    double beam_inertia = 0.1;
    printf("Parameters:\n");
    printf("  Ball mass = %.3f kg\n", ball_mass);
    printf("  Gravity = %.2f m/s^2\n", grav);
    printf("  Ball radius = %.3f m\n", ball_radius);
    printf("  Beam inertia = %.3f kg*m^2\n", beam_inertia);

    int n = 4, mi = 1, p = 1;
    NonlinearSystem *bb = nonlinear_system_create(n, mi, p, "Ball and Beam");
    VectorField *vf = vector_field_create(n, NULL, NULL);
    VectorField *vg = vector_field_create(n, NULL, NULL);
    ScalarField *vh = scalar_field_create(n, NULL, NULL);
    nonlinear_system_set_drift(bb, vf);
    nonlinear_system_set_control(bb, 0, vg);
    nonlinear_system_set_output(bb, 0, vh);

    Vector *x_eq = vector_create(n);
    vector_set(x_eq, 0.0);
    nonlinear_system_set_equilibrium(bb, x_eq);
    nonlinear_system_print_info(bb);

    printf("\n--- Accessibility Distribution ---\n");
    int rank = controllability_rank(vf, vg, x_eq, 1e-6);
    printf("Controllability rank: %d (need %d)\n", rank, n);

    Matrix *dist = compute_accessibility_distribution(vf, vg, x_eq, n - 1);
    if (dist) {
        printf("Distribution built (%dx%d)\n", dist->rows, dist->cols);
        matrix_free(dist);
    }

    printf("\n--- Relative Degree ---\n");
    RelativeDegreeResult rd = compute_relative_degree(bb, x_eq, n);
    if (rd.is_well_defined) {
        printf("rho = %d\n", rd.rho);
        if (rd.rho == n) printf("Full relative degree!\n");
        IOLinearizationResult io = input_output_linearize(bb, x_eq);
        if (io.is_feedback_linearizable)
            printf("Law: alpha=%f, beta=%f\n", io.alpha, io.beta);
    }

    printf("\n--- Input-State Linearization ---\n");
    printf("Full-state linearizable: %s\n",
           full_state_feedback_linearizable(bb) ? "YES" : "NO");

    ZeroDynamicsResult zd = analyze_zero_dynamics(bb, x_eq);
    printf("\n--- Zero Dynamics: dim=%d, min-phase=%s ---\n",
           zd.dim_zero_dynamics, check_minimum_phase(&zd) ? "YES" : "NO");

    if (rd.Lf_powers_h) free(rd.Lf_powers_h);
    if (zd.zero_eigenvalues) free(zd.zero_eigenvalues);
    if (zd.zero_manifold) vector_free(zd.zero_manifold);
    vector_free(x_eq);
    nonlinear_system_free(bb);

    printf("\n========================================\n");
    printf("Ball and Beam example complete.\n");
    printf("========================================\n");
    return 0;
}
