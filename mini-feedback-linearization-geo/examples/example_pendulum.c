#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nonlinear_system.h"
#include "lie_derivative.h"
#include "lie_bracket.h"
#include "feedback_linearization.h"

int main(void) {
    printf("========================================\n");
    printf("Example: Pendulum on a Cart\n");
    printf("Feedback Linearization Demo\n");
    printf("========================================\n\n");

    int n = 4, m = 1, p = 1;
    NonlinearSystem *pendulum = nonlinear_system_create(n, m, p, "Pendulum on Cart");
    printf("System: %s\n", pendulum->name);
    printf("State dim: %d, Inputs: %d, Outputs: %d\n", n, m, p);

    VectorField *f = vector_field_create(n, NULL, NULL);
    VectorField *g = vector_field_create(n, NULL, NULL);
    ScalarField *h = scalar_field_create(n, NULL, NULL);
    nonlinear_system_set_drift(pendulum, f);
    nonlinear_system_set_control(pendulum, 0, g);
    nonlinear_system_set_output(pendulum, 0, h);

    Vector *x_eq = vector_create(n);
    vector_set(x_eq, 0.0);
    nonlinear_system_set_equilibrium(pendulum, x_eq);
    nonlinear_system_print_info(pendulum);

    printf("\n--- Controllability Analysis ---\n");
    int rank = controllability_rank(f, g, x_eq, 1e-6);
    printf("Accessibility rank at origin: %d (need %d)\n", rank, n);

    printf("\n--- Relative Degree ---\n");
    RelativeDegreeResult rd = compute_relative_degree(pendulum, x_eq, n);
    if (rd.is_well_defined) {
        printf("Relative degree: rho = %d\n", rd.rho);
        printf("L_g L_f^{rho-1} h(x_eq) = %f\n", rd.Lg_Lf_power);
        IOLinearizationResult io = input_output_linearize(pendulum, x_eq);
        if (io.is_feedback_linearizable) {
            printf("\n--- Feedback Law ---\n");
            printf("alpha(x) = %f\n", io.alpha);
            printf("beta(x)  = %f\n", io.beta);
            printf("u = alpha + beta * v\n");
        }
    } else {
        printf("No well-defined relative degree at origin.\n");
    }

    printf("\n--- Input-State Linearization ---\n");
    ISLinearizationResult is = input_state_linearize(pendulum);
    printf("Exactly linearizable: %s\n", is.is_exactly_linearizable ? "YES" : "NO");

    printf("\n--- Zero Dynamics ---\n");
    ZeroDynamicsResult zd = analyze_zero_dynamics(pendulum, x_eq);
    printf("Zero dynamics dim: %d\n", zd.dim_zero_dynamics);
    printf("Stability: %s\n",
        zd.stability == ZERO_DYN_STABLE ? "STABLE" :
        zd.stability == ZERO_DYN_UNSTABLE ? "UNSTABLE" : "OTHER");
    printf("Minimum phase: %s\n", check_minimum_phase(&zd) ? "YES" : "NO");

    if (rd.Lf_powers_h) free(rd.Lf_powers_h);
    if (zd.zero_eigenvalues) free(zd.zero_eigenvalues);
    if (zd.zero_manifold) vector_free(zd.zero_manifold);
    if (is.phi) diffeomorphism_free(is.phi);
    if (is.f_transformed) vector_field_free(is.f_transformed);
    if (is.g_transformed) vector_field_free(is.g_transformed);
    vector_free(x_eq);
    nonlinear_system_free(pendulum);

    printf("\n========================================\n");
    printf("Example complete.\n");
    printf("========================================\n");
    return 0;
}
