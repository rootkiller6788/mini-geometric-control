/* ============================================================================
 * demo_overview.c — Nonholonomic Planning Concept Demonstrator
 *
 * Interactive overview of key nonholonomic planning concepts:
 *   1. Nonholonomic constraint visualization (unicycle)
 *   2. Lie bracket generation of new motion directions
 *   3. Geometric phase / holonomy
 *   4. Chow-Rashevskii controllability check for all models
 *   5. Growth vector computation
 *
 * This demo serves as the entry point for understanding the module.
 *
 * Reference: Bloch (2003), LaValle (2006).
 * ============================================================================ */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/noplan_core.h"
#include "../include/noplan_geometry.h"
#include "../include/noplan_planning.h"
#include "../include/noplan_models.h"

static void print_separator(const char* title) {
    printf("\n═══════════════════════════════════════════════════\n");
    printf("  %s\n", title);
    printf("═══════════════════════════════════════════════════\n\n");
}

int main(void) {
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║     Nonholonomic Motion Planning — Overview Demo     ║\n");
    printf("║   Geometric Control Theory Meets Path Planning       ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n");

    /* ================================================================
     * Concept 1: What is a nonholonomic constraint?
     * ================================================================ */
    print_separator("1. Nonholonomic Constraints — The Unicycle");

    printf("A nonholonomic constraint is a velocity constraint that\n");
    printf("CANNOT be integrated to a position constraint.\n\n");

    printf("Unicycle example:\n");
    printf("  Constraint: ẋ·sin θ − ẏ·cos θ = 0\n");
    printf("  → No sideways velocity allowed (rolling without slipping)\n");
    printf("  → This constraint is NON-INTEGRABLE → nonholonomic\n\n");

    printf("However, the unicycle CAN reach any (x, y, θ)!\n");
    printf("  → Nonholonomic ≠ Uncontrollable\n");
    printf("  → Lie bracket [g₁, g₂] generates side-slip motion\n\n");

    /* Demonstrate the constraint */
    NonholonomicSystem* uni = noplan_model_unicycle();
    printf("Control vector fields at θ = 0:\n");
    printf("  g₁(0) = (cos 0,  sin 0,  0) = (1, 0, 0)  ← forward\n");
    printf("  g₂(0) = (0,      0,      1) = (0, 0, 1)  ← rotate\n");

    Config q0 = uni->current_q;
    LieBracket lb;
    noplan_lie_bracket(&q0, noplan_unicycle_g1, noplan_unicycle_g2,
                        3, 1e-5, &lb);
    printf("  [g₁,g₂](0) = (%.2f, %.2f, %.2f)          ← side-slip!\n",
           lb.value.v[0], lb.value.v[1], lb.value.v[2]);
    printf("  This 3rd direction is NOT directly actuated\n");
    printf("  but IS reachable via sequential g₁ and g₂.\n");

    /* ================================================================
     * Concept 2: Lie Bracket Motion
     * ================================================================ */
    print_separator("2. Lie Bracket Motion — The Geometric Phase");

    printf("Execute cyclic control: g₁ → g₂ → -g₁ → -g₂\n");
    printf("Each step individually returns to starting position,\n");
    printf("but the composition produces NET MOTION in [g₁,g₂]!\n\n");

    printf("This is the geometric phase / holonomy:\n");
    printf("  • Analogous to the 'falling cat' reorienting in free fall\n");
    printf("  • Same math as Berry phase in quantum mechanics\n");
    printf("  • Enables parallel parking without side-slip actuators\n\n");

    /* Simulate a short Lie bracket cycle */
    noplan_system_reset(uni, &q0);
    double u_cycle[40 * 2];
    for (int i = 0; i < 10; i++) {
        u_cycle[i * 2 + 0] = 0.5; u_cycle[i * 2 + 1] = 0.0;
        u_cycle[(10+i) * 2 + 0] = 0.0; u_cycle[(10+i) * 2 + 1] = 1.0;
        u_cycle[(20+i) * 2 + 0] = -0.5; u_cycle[(20+i) * 2 + 1] = 0.0;
        u_cycle[(30+i) * 2 + 0] = 0.0; u_cycle[(30+i) * 2 + 1] = -1.0;
    }

    Trajectory* cycle_traj = noplan_trajectory_create(41, 2);
    noplan_system_simulate(uni, u_cycle, 40, 0.05, cycle_traj);

    Config q_after = cycle_traj->waypoints[40];
    printf("After one Lie bracket cycle:\n");
    printf("  Δx = %.4f m\n", q_after.q[0] - q0.q[0]);
    printf("  Δy = %.4f m  ← NET side-slip!\n", q_after.q[1] - q0.q[1]);
    printf("  Δθ = %.4f°\n", (q_after.q[2] - q0.q[2]) * 180.0 / M_PI);

    noplan_trajectory_free(cycle_traj);

    /* ================================================================
     * Concept 3: Controllability Analysis
     * ================================================================ */
    print_separator("3. Controllability Analysis — All Canonical Models");

    struct {
        const char* name;
        NonholonomicSystem* sys;
    } models[] = {
        {"Unicycle",           noplan_model_unicycle()},
        {"Car (bicycle)",      noplan_model_car_default()},
        {"Knife Edge",         noplan_model_knife_edge()},
        {"Rolling Ball",       noplan_model_rolling_ball(0.1)},
        {"Snakeboard",         noplan_model_snakeboard_default()},
    };
    int n_models = 5;

    printf("%-20s %4s %4s %7s %8s %6s\n",
           "Model", "n", "m", "LARC", "Brockett", "κ");
    printf("%-20s %4s %4s %7s %8s %6s\n",
           "-----", "-", "-", "----", "--------", "-");

    for (int i = 0; i < n_models; i++) {
        NonholonomicSystem* sys = models[i].sys;
        Config q = sys->current_q;

        bool larc = noplan_larc_holds(sys, &q, 1e-4);
        bool brockett = noplan_brockett_check(sys, &q);
        int deg = noplan_degree_of_nonholonomy(sys, &q, 1e-4);

        printf("%-20s %4d %4d %7s %8s %6d\n",
               models[i].name,
               sys->config_dim,
               sys->n_inputs,
               larc ? "YES" : "NO",
               brockett ? "YES" : "NO",
               deg);

        noplan_system_free(sys);
    }

    printf("\nKey insight column:\n");
    printf("  LARC=YES → System is controllable (can reach any config)\n");
    printf("  Brockett=NO → Cannot use smooth feedback (must plan paths)\n");
    printf("  κ = degree of nonholonomy (higher = harder to plan)\n");

    /* ================================================================
     * Concept 4: Growth Vector
     * ================================================================ */
    print_separator("4. Growth Vectors of Canonical Systems");

    NonholonomicSystem* uni2 = noplan_model_unicycle();
    Config q_uni = uni2->current_q;

    GrowthVector gv;
    noplan_compute_growth_vector(uni2, &q_uni, 3, 1e-5, &gv);

    printf("Unicycle growth vector (κ=%d):\n", gv.kappa);
    for (int k = 1; k <= gv.kappa; k++) {
        printf("  r_%d = dim(Δ_%d) = %d\n", k, k, gv.growth[k]);
    }
    printf("  Δ₁ = span{g₁, g₂}               — base controls\n");
    printf("  Δ₂ = span{g₁, g₂, [g₁,g₂]}       — full tangent space\n");
    printf("  → κ = 2 means we need brackets of depth ≤ 2\n");

    noplan_system_free(uni2);

    /* ================================================================
     * Concept 5: Planning Algorithm Comparison
     * ================================================================ */
    print_separator("5. Planning Algorithm Landscape");

    printf("Algorithm          | Complexity     | Optimal? | Requires\n");
    printf("-------------------|----------------|----------|----------\n");
    printf("RRT                | O(K log K)     | No       | Forward sim\n");
    printf("PRM                | O(N²·steering) | No       | Steering fn\n");
    printf("Sinusoidal         | O(n·m)         | No       | Chained form\n");
    printf("Nilpotent steering | O(2^m·n)       | Local    | Nilp. approx\n");
    printf("RRT-Connect        | O(√K log K)    | No       | Bidirectional\n");
    printf("Lattice A*         | O(V log V)     | Yes*     | Primitives\n");

    printf("\n* Optimal with respect to the lattice discretization.\n");
    printf("  Continuous optimal planning is NP-hard for nonholonomic\n");
    printf("  systems in general (Laumond et al., 1994).\n");

    /* ================================================================
     * Concept 6: Fundamental Theorems
     * ================================================================ */
    print_separator("6. Fundamental Theorems in Nonholonomic Planning");

    printf("1. FROBENIUS THEOREM (1877)\n");
    printf("   Δ involutive ⇔ Δ integrable (holonomic)\n");
    printf("   → Test: compute all [g_i, g_j] and check if in span{g_k}\n\n");

    printf("2. CHOW-RASHEVSKII THEOREM (1938/1939)\n");
    printf("   Lie algebra of {g_i} spans TQ ⇒ system is controllable\n");
    printf("   → Test: compute iterated Lie brackets, check rank\n\n");

    printf("3. BROCKETT'S NECESSARY CONDITION (1983)\n");
    printf("   Smooth feedback stabilization ⇒ (q,u)→f(q,u) surjective\n");
    printf("   → Explains why unicycle cannot be smoothly stabilized\n\n");

    printf("4. SUSSMANN'S ORBIT THEOREM (1987)\n");
    printf("   General necessary condition for STLC (with drift)\n");
    printf("   → Generalizes Chow-Rashevskii to systems with drift\n");

    /* ================================================================
     * Summary
     * ================================================================ */
    print_separator("Summary");

    printf("Nonholonomic motion planning integrates:\n");
    printf("  • Differential geometry (Lie brackets, distributions)\n");
    printf("  • Control theory (controllability, feedback limitations)\n");
    printf("  • Sampling-based planning (RRT, PRM)\n");
    printf("  • Geometric mechanics (SE(2), SE(3), geometric phases)\n\n");

    printf("Core insight:\n");
    printf("  Constraints limit velocities, not positions.\n");
    printf("  Lie brackets create new motion directions.\n");
    printf("  With enough wiggling, we can go anywhere.\n");
    printf("  But smooth feedback won't get us there.\n\n");

    printf("This module implements the complete pipeline from\n");
    printf("theoretical analysis to practical path planning.\n");

    return 0;
}
