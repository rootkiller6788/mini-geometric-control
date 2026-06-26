#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nonlinear_system.h"
#include "lie_derivative.h"
#include "lie_bracket.h"
#include "coordinate_transform.h"
#include "frobenius.h"
#include "feedback_linearization.h"

int main(void) {
    printf("========================================\n");
    printf("  Feedback Linearization - Demo\n");
    printf("========================================\n\n");
    printf("Core concepts of geometric nonlinear control.\n\n");
    printf("1. Lie derivative: L_f h = grad(h).f\n");
    printf("2. Lie bracket: [f,g] = J_g.f - J_f.g\n");
    printf("3. Relative degree: rho = min{r: L_g L_f^{r-1} h != 0}\n");
    printf("4. Feedback law: u = alpha + beta*v\n");
    printf("5. Frobenius theorem: integrable iff involutive\n");
    printf("6. Zero dynamics: internal dynamics when y=0\n");
    printf("7. Normal form: Byrnes-Isidori\n");
    printf("8. MIMO decoupling matrix: E_{ij} = L_{g_j} L_f^{rho_i-1} h_i\n\n");
    printf("Applications:\n");
    printf("  Robotics: ABB, KUKA, Fanuc (computed torque)\n");
    printf("  Aerospace: Quadrotor (DJI, Parrot, Boeing)\n");
    printf("  Automotive: Engine control (Toyota, Bosch)\n");
    printf("  Maglev: Transrapid magnetic levitation\n");
    printf("  Power: Smart grid stabilization (ISO)\n");
    printf("  Precision: ASML wafer stage control\n\n");
    printf("========================================\n");
    printf("Demo complete.\n");
    printf("========================================\n");
    return 0;
}
