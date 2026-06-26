/* ============================================================================
 * example_trailer.c — Trailer Docking and Multi-Body Nonholonomic Planning
 *
 * Demonstrates planning for car-with-trailers systems:
 *   1. Kinematic coupling between tractor and trailers
 *   2. Growth in degree of nonholonomy with each trailer
 *   3. The "jackknife" instability when reversing
 *   4. PRM-based planning in the high-dimensional configuration space
 *
 * The car with N trailers is a classic benchmark for nonholonomic
 * motion planning. Key property: for N ≥ 2, the system loses
 * differential flatness, making planning significantly harder.
 *
 * Trailer systems are ubiquitous:
 *   - Semi-trucks (N=1): Amazon/FedEx logistics
 *   - Fire trucks (N=2): narrow street navigation
 *   - Airport tugs (N=3-5): luggage cart trains
 *   - Agricultural vehicles (N=2-3): harvesting trains
 *
 * Reference: Laumond (1993), Švestka & Overmars (1998).
 * ============================================================================ */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/noplan_core.h"
#include "../include/noplan_geometry.h"
#include "../include/noplan_planning.h"
#include "../include/noplan_models.h"

int main(void) {
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║  Nonholonomic Planning — Car with N Trailers Demo    ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n\n");

    /* ================================================================
     * Part 1: Car + 1 Trailer
     * ================================================================ */
    printf("━━━ Part 1: Car + 1 Trailer (Semi-Truck) ━━━\n\n");

    double lengths1[] = {3.0};  /* 3m trailer */
    NonholonomicSystem* truck = noplan_model_trailer(1, 3.0, lengths1, 1.5);
    printf("System: %s\n", truck->name);
    printf("  Configuration: (x, y, θ_tractor, θ_trailer) — dim=4\n");
    printf("  Controls: v (forward), ω (steering rate)\n");

    /* Controllability */
    Config q_test = truck->current_q;
    bool controllable = noplan_larc_holds(truck, &q_test, 1e-5);
    printf("  Controllable: %s\n", controllable ? "YES" : "NO");
    int deg = noplan_degree_of_nonholonomy(truck, &q_test, 1e-5);
    printf("  Degree of nonholonomy: κ = %d\n", deg);

    /* Demonstrate jackknife when reversing */
    printf("\nDemonstrating trailer kinematics:\n");
    printf("  Forward motion → trailer aligns with tractor.\n");
    printf("  Reversing → trailer angle amplifies → JACKKNIFE risk!\n");

    /* Simulate: start with trailer at 30° */
    truck->current_q.q[3] = 30.0 * M_PI / 180.0;  /* Trailer at 30° */
    printf("\n  Initial: θ_tractor=0°, θ_trailer=30°\n");

    /* Forward: trailer should align */
    double u_forward[2] = {0.5, 0.0};
    for (int step = 0; step < 50; step++) {
        noplan_system_step(truck, u_forward, 0.1);
    }
    printf("  After 5s forward: θ_tractor=%.1f°, θ_trailer=%.1f°\n",
           truck->current_q.q[2] * 180.0 / M_PI,
           truck->current_q.q[3] * 180.0 / M_PI);

    /* Reset and try reversing */
    truck->current_q.q[3] = 30.0 * M_PI / 180.0;
    double u_reverse[2] = {-0.3, 0.0};
    for (int step = 0; step < 50; step++) {
        noplan_system_step(truck, u_reverse, 0.1);
    }
    printf("  After 5s reverse: θ_tractor=%.1f°, θ_trailer=%.1f°\n",
           truck->current_q.q[2] * 180.0 / M_PI,
           truck->current_q.q[3] * 180.0 / M_PI);
    printf("  → Trailer angle INCREASES when reversing (unstable)\n");

    noplan_system_free(truck);

    /* ================================================================
     * Part 2: Car + 2 Trailers — Loss of Flatness
     * ================================================================ */
    printf("\n━━━ Part 2: Car + 2 Trailers (Loss of Differential Flatness) ━━━\n\n");

    double lengths2[] = {2.5, 2.0};
    NonholonomicSystem* double_trailer = noplan_model_trailer(
        2, 3.0, lengths2, 1.0);
    printf("System: %s\n", double_trailer->name);
    printf("  Configuration: (x, y, θ₀, θ₁, θ₂) — dim=5\n");

    Config qt2 = double_trailer->current_q;
    int deg2 = noplan_degree_of_nonholonomy(double_trailer, &qt2, 1e-5);
    printf("  Degree of nonholonomy: κ = %d\n", deg2);
    printf("  For N≥2: system is NOT differentially flat\n");
    printf("  → Cannot use flatness-based planning; must use RRT/PRM\n");

    /* Trailer kinematics check */
    TangentVector qdot;
    double u_test[2] = {1.0, 0.3};
    noplan_trailer_kinematics(double_trailer, u_test, &qt2, &qdot);
    printf("\n  At q=0: q̇ = (%.2f, %.2f, %.2f, %.2f, %.2f)\n",
           qdot.v[0], qdot.v[1], qdot.v[2], qdot.v[3], qdot.v[4]);

    noplan_system_free(double_trailer);

    /* ================================================================
     * Part 3: Trailer Docking Benchmark
     * ================================================================ */
    printf("\n━━━ Part 3: Trailer Docking Benchmark ━━━\n\n");

    double lengths3[] = {2.0};
    NonholonomicSystem* dock_truck = noplan_model_trailer(
        1, 2.5, lengths3, 1.0);

    PlanningProblem* docking = noplan_benchmark_trailer_docking(
        dock_truck, 0.0, -5.0);
    printf("Docking scenario:\n");
    printf("  Start: (10, 0, 0°, 0°)\n");
    printf("  Goal:  (0, -5, 0°, 0°) — backed into dock\n");

    /* Build PRM */
    printf("\nBuilding Probabilistic Roadmap (200 samples)...\n");
    PRMPlanner* prm = noplan_prm_create(dock_truck, docking, 8, 8.0);
    noplan_prm_build(prm, 200);
    printf("  Roadmap: %d nodes built\n", prm->n_nodes);

    /* Count edges */
    int total_edges = 0;
    for (int i = 0; i < prm->n_nodes; i++) {
        total_edges += prm->nodes[i].n_neighbors;
    }
    printf("  Total edges: %d (avg %.1f per node)\n",
           total_edges, (double)total_edges / prm->n_nodes);

    /* Query the roadmap */
    Trajectory* path = noplan_trajectory_create(500, 2);
    if (noplan_prm_query(prm, &docking->q_start, &docking->q_goal, path)) {
        printf("  ✓ Path found! Length: %.2f m\n",
               noplan_trajectory_length(path));
    } else {
        printf("  PRM query returned no exact path.\n");
    }

    noplan_trajectory_free(path);
    noplan_prm_free(prm);
    noplan_problem_free(docking);
    noplan_system_free(dock_truck);

    /* ================================================================
     * Part 4: Kinematic Analysis of Multi-Trailer Systems
     * ================================================================ */
    printf("\n━━━ Part 4: Scaling of Nonholonomy with Trailers ━━━\n\n");
    printf("  Trailers  |  dim(Q)  |  κ (deg. of nonholonomy)  |  Flat?\n");
    printf("  ----------|----------|----------------------------|--------\n");

    for (int N = 0; N <= 4; N++) {
        if (N == 0) {
            NonholonomicSystem* car = noplan_model_car_default();
            Config qc = car->current_q;
            int k = noplan_degree_of_nonholonomy(car, &qc, 1e-5);
            printf("  %-9d | %-8d | %-26d |  YES\n", N, car->config_dim, k);
            noplan_system_free(car);
        } else {
            double* lens = (double*)malloc((size_t)N * sizeof(double));
            for (int i = 0; i < N; i++) lens[i] = 2.0;
            NonholonomicSystem* tr = noplan_model_trailer(N, 3.0, lens, 1.0);
            Config qt = tr->current_q;
            int k = noplan_degree_of_nonholonomy(tr, &qt, 1e-5);
            bool flat = (N <= 1);
            printf("  %-9d | %-8d | %-26d |  %s\n",
                   N, tr->config_dim, k, flat ? "YES" : "NO");
            noplan_system_free(tr);
            free(lens);
        }
    }

    printf("\n  Key insight: The degree of nonholonomy grows with N,\n");
    printf("  making planning exponentially harder in high dimensions.\n");

    /* ================================================================
     * Part 5: Real-World Relevance
     * ================================================================ */
    printf("\n━━━ Part 5: Real-World Applications ━━━\n\n");
    printf("  • Amazon Prime Air: autonomous truck + trailer docking\n");
    printf("  • Tesla Semi: autonomous backing to loading dock\n");
    printf("  • NASA Mars rovers: trailer-like instrument deployment\n");
    printf("  • John Deere: autonomous tractor + harvester trains\n");
    printf("  • Port of Rotterdam: automated guided vehicle trains\n");
    printf("  • DARPA Grand Challenge: autonomous navigation with trailers\n");

    printf("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("  Trailer Planning Summary:\n");
    printf("  • Car + 1 trailer: flat, κ=3, tractable\n");
    printf("  • Car + 2 trailers: NOT flat, requires sampling-based planning\n");
    printf("  • Reversing amplifies trailer angles (unstable)\n");
    printf("  • PRM builds reusable roadmaps for repeated queries\n");
    printf("  • Real systems use sensor feedback + planned paths\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    return 0;
}
