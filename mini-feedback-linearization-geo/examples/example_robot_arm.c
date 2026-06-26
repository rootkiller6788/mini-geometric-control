#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nonlinear_system.h"
#include "lie_derivative.h"
#include "lie_bracket.h"
#include "feedback_linearization.h"

int main(void) {
    printf("========================================\n");
    printf("Example: 2-Link Robot Manipulator\n");
    printf("MIMO Feedback Linearization\n");
    printf("========================================\n\n");

    double m1 = 1.0, m2 = 1.0, l1 = 0.5, l2 = 0.5, grav = 9.81;
    printf("Robot Parameters:\n");
    printf("  Link 1: m1=%.1f kg, l1=%.1f m\n", m1, l1);
    printf("  Link 2: m2=%.1f kg, l2=%.1f m\n", m2, l2);
    printf("  Gravity = %.2f m/s^2\n\n", grav);

    int n = 4, mm = 2, p = 2;
    NonlinearSystem *robot = nonlinear_system_create(n, mm, p, "2-Link Robot");
    VectorField *vf = vector_field_create(n, NULL, NULL);
    VectorField *vg1 = vector_field_create(n, NULL, NULL);
    VectorField *vg2 = vector_field_create(n, NULL, NULL);
    ScalarField *vh1 = scalar_field_create(n, NULL, NULL);
    ScalarField *vh2 = scalar_field_create(n, NULL, NULL);
    nonlinear_system_set_drift(robot, vf);
    nonlinear_system_set_control(robot, 0, vg1);
    nonlinear_system_set_control(robot, 1, vg2);
    nonlinear_system_set_output(robot, 0, vh1);
    nonlinear_system_set_output(robot, 1, vh2);

    Vector *x_eq = vector_create(n);
    vector_set(x_eq, 0.0);
    nonlinear_system_set_equilibrium(robot, x_eq);
    nonlinear_system_print_info(robot);

    printf("\n--- MIMO Feedback Linearization ---\n");
    MIMOLinearizationResult mimo = mimo_linearize(robot, x_eq);
    printf("Vector relative degree: %s\n",
           mimo.has_vector_relative_degree ? "YES" : "NO");
    if (mimo.has_vector_relative_degree) {
        for (int i = 0; i < mimo.num_outputs; i++)
            printf("  rho_%d = %d\n", i + 1, mimo.relative_degrees[i]);
        printf("Decoupling matrix invertible: %s\n",
               mimo.decoupling_matrix_invertible ? "YES" : "NO");
    }

    printf("\n--- Computed Torque Control ---\n");
    printf("M(q)*q'' + C*q' + G = tau\n");
    printf("tau = M(q)*v + C*q' + G => q'' = v\n");
    printf("Used in: ABB, KUKA, Fanuc robots\n");

    printf("\n--- SISO Linearization ---\n");
    IOLinearizationResult io = input_output_linearize(robot, x_eq);
    if (io.is_feedback_linearizable)
        printf("alpha=%f, beta=%f\n", io.alpha, io.beta);
    else
        printf("MIMO system - use mimo_linearize()\n");

    if (mimo.relative_degrees) free(mimo.relative_degrees);
    if (mimo.decoupling_matrix) matrix_free(mimo.decoupling_matrix);
    if (mimo.A) matrix_free(mimo.A);
    if (mimo.B) matrix_free(mimo.B);
    vector_free(x_eq);
    nonlinear_system_free(robot);

    printf("\n========================================\n");
    printf("Robot arm example complete.\n");
    printf("========================================\n");
    return 0;
}
