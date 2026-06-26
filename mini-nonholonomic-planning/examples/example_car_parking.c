/* ============================================================================
 * example_car_parking.c — Car-Like Robot Parallel Parking
 *
 * Demonstrates nonholonomic motion planning for a car-like vehicle:
 *   1. Parallel parking using Reeds-Shepp maneuvers
 *   2. Garage parking (perpendicular)
 *   3. RRT planning with car kinematics
 *   4. Shortcut path optimization
 *
 * The car model (bicycle kinematics) with bounded steering angle
 * is the canonical example of a nonholonomic system that requires
 * back-and-forth maneuvers for tight parking — a task that autonomous
 * vehicles (Tesla, Waymo, etc.) must solve in the real world.
 *
 * Key constraints:
 *   - Minimum turning radius R_min = L / tan(φ_max)
 *   - Cannot move sideways → parallel parking requires multiple maneuvers
 *   - Curvature bounded by steering limits
 *
 * Reference: Reeds & Shepp (1990), Laumond et al. (1994).
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
    printf("║   Nonholonomic Motion Planning — Car Parking Demo    ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n\n");

    /* ================================================================
     * Part 1: Car Model Properties
     * ================================================================ */
    printf("━━━ Part 1: Car Model ━━━\n\n");

    NonholonomicSystem* car = noplan_model_car_default();
    noplan_model_print_kinematics(car);

    /* Compute the minimum turning radius */
    double L = 2.5, phi_max = 0.610865;
    double R_min = L / tan(phi_max);
    printf("  Wheelbase L = %.1f m\n", L);
    printf("  Max steering angle φ_max = %.1f°\n", phi_max * 180.0 / M_PI);
    printf("  Min turning radius R_min = %.2f m\n", R_min);

    /* Controllability check */
    Config q_test = car->current_q;
    bool controllable = noplan_larc_holds(car, &q_test, 1e-5);
    printf("  Controllable: %s\n", controllable ? "YES ✓" : "NO");
    if (controllable) {
        int deg = noplan_degree_of_nonholonomy(car, &q_test, 1e-5);
        printf("  Degree of nonholonomy: κ = %d\n", deg);
    }

    /* ================================================================
     * Part 2: Dubins vs Reeds-Shepp Path Lengths
     * ================================================================ */
    printf("\n━━━ Part 2: Dubins vs Reeds-Shepp ━━━\n\n");

    /* Test configurations */
    struct { double sx, sy, st, gx, gy, gt; const char* desc; } cases[] = {
        {0, 0, 0, 5, 0, 0,              "Straight forward"},
        {0, 0, 0, 0, 3, 0,              "Lateral displacement (parallel park)"},
        {0, 0, 0, -3, 0, M_PI,          "Reverse direction"},
        {0, 0, 0, 5, 5, M_PI/2,         "Arbitrary pose"},
    };

    printf("Path length comparison (R_min = %.2f m):\n\n", R_min);
    printf("  %-30s %12s %12s\n", "Scenario", "Dubins", "Reeds-Shepp");
    printf("  %-30s %12s %12s\n", "--------", "------", "-----------");

    for (int i = 0; i < 4; i++) {
        Config qs = noplan_config_create(3, 'S');
        Config qg = noplan_config_create(3, 'S');
        qs.q[0] = cases[i].sx; qs.q[1] = cases[i].sy; qs.q[2] = cases[i].st;
        qg.q[0] = cases[i].gx; qg.q[1] = cases[i].gy; qg.q[2] = cases[i].gt;

        double dubins_len = noplan_dubins_length(&qs, &qg, R_min);
        double rs_len = noplan_reeds_shepp_length(&qs, &qg, R_min);

        printf("  %-30s %10.2f m %10.2f m\n", cases[i].desc,
               dubins_len, rs_len);
    }

    /* ================================================================
     * Part 3: Parallel Parking via RRT
     * ================================================================ */
    printf("\n━━━ Part 3: Parallel Parking Planning ━━━\n\n");

    PlanningProblem* park = noplan_benchmark_parallel_parking(car, 2.5, 5.0, 2.0);
    printf("Parallel parking scenario:\n");
    printf("  Start: (0, 0, 0°) — aligned with curb\n");
    printf("  Goal:  (0, %.1f, 0°) — in parking spot\n", 2.5);
    printf("  Obstacles: 2 bounding cars\n");

    /* Attempt RRT planning */
    RRTPlanner* rrt = noplan_rrt_create(car, park, 0.3);
    rrt->goal_tolerance = 0.3;
    noplan_rrt_set_goal_bias(rrt, 0.15);

    printf("\nPlanning with RRT (max 2000 iterations)...\n");
    Trajectory* plan = noplan_trajectory_create(1000, 2);

    bool success = noplan_rrt_plan(rrt, 2000, plan);

    if (success) {
        printf("  ✓ Parking maneuver found!\n");
        printf("  Tree size: %d nodes\n", rrt->n_nodes);
        printf("  Path waypoints: %d\n", plan->n_waypoints);
        printf("  Estimated maneuver time: %.1f s\n", plan->total_time);

        /* Check for collisions */
        bool collision = noplan_collision_check_trajectory(park, plan);
        printf("  Collision-free: %s\n", collision ? "NO ✗" : "YES ✓");

        /* Show maneuver phases */
        printf("\n  Maneuver phases (every 25%% of path):\n");
        int n = plan->n_waypoints;
        for (int pct = 0; pct <= 100; pct += 25) {
            int idx = n * pct / 100;
            if (idx >= n) idx = n - 1;
            Config q = plan->waypoints[idx];
            printf("    %3d%%: (%.2f, %.2f, %.1f°)\n",
                   pct, q.q[0], q.q[1], q.q[2] * 180.0 / M_PI);
        }
    } else {
        printf("  RRT did not find exact path to goal.\n");
        printf("  Tree explored %d nodes.\n", rrt->n_nodes);

        /* Show best attempt */
        double best_d = INFINITY;
        int best_i = 0;
        for (int i = 0; i < rrt->n_nodes; i++) {
            double d = noplan_config_distance_se2(
                &rrt->nodes[i].q, &park->q_goal);
            if (d < best_d) { best_d = d; best_i = i; }
        }
        printf("  Closest approach: %.3f m at node %d\n", best_d, best_i);
        Config q_best = rrt->nodes[best_i].q;
        printf("    (%.2f, %.2f, %.1f°)\n",
               q_best.q[0], q_best.q[1], q_best.q[2] * 180.0 / M_PI);
    }

    noplan_trajectory_free(plan);
    noplan_rrt_free(rrt);
    noplan_problem_free(park);

    /* ================================================================
     * Part 4: Garage Parking
     * ================================================================ */
    printf("\n━━━ Part 4: Garage Parking ━━━\n\n");

    PlanningProblem* garage = noplan_benchmark_garage_parking(car, 3.0);
    printf("Garage parking scenario:\n");
    printf("  Start: (-5, 0, 90°) — approaching garage\n");
    printf("  Goal:  (0, 3, 90°) — inside garage\n");

    RRTPlanner* rrt2 = noplan_rrt_create(car, garage, 0.2);
    rrt2->goal_tolerance = 0.4;
    noplan_rrt_set_goal_bias(rrt2, 0.2);

    printf("Planning with RRT (max 1500 iterations)...\n");
    Trajectory* plan2 = noplan_trajectory_create(800, 2);

    bool success2 = noplan_rrt_plan(rrt2, 1500, plan2);
    if (success2) {
        printf("  ✓ Garage parking maneuver found!\n");
        printf("  Path waypoints: %d, Length: %.2f m\n",
               plan2->n_waypoints, noplan_trajectory_length(plan2));
    } else {
        printf("  RRT did not reach goal (explored %d nodes).\n", rrt2->n_nodes);
    }

    noplan_trajectory_free(plan2);
    noplan_rrt_free(rrt2);
    noplan_problem_free(garage);

    /* ================================================================
     * Part 5: Path Optimization (Shortcut)
     * ================================================================ */
    printf("\n━━━ Part 5: Path Optimization ━━━\n\n");

    /* Quick RRT to get a rough path, then optimize */
    Config qs2 = noplan_config_create(3, 'S');
    Config qg2 = noplan_config_create(3, 'S');
    qg2.q[0] = 4.0; qg2.q[1] = 3.0; qg2.q[2] = M_PI / 6;
    noplan_system_reset(car, &qs2);

    PlanningProblem* prob3 = noplan_problem_create(car, &qs2, &qg2);
    RRTPlanner* rrt3 = noplan_rrt_create(car, prob3, 0.25);
    rrt3->goal_tolerance = 0.3;

    Trajectory* raw_path = noplan_trajectory_create(600, 2);
    if (noplan_rrt_plan(rrt3, 800, raw_path)) {
        printf("Raw RRT path:  %d waypoints, length=%.2f m\n",
               raw_path->n_waypoints, noplan_trajectory_length(raw_path));

        /* Apply shortcut smoothing */
        noplan_path_shortcut(raw_path, car, prob3, 50);

        printf("Optimized path: %d waypoints, length=%.2f m\n",
               raw_path->n_waypoints, noplan_trajectory_length(raw_path));
    } else {
        printf("Could not find initial path for optimization demo.\n");
    }

    noplan_trajectory_free(raw_path);
    noplan_rrt_free(rrt3);
    noplan_problem_free(prob3);

    /* ================================================================
     * Summary
     * ================================================================ */
    printf("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("  Auto-Parking Summary:\n");
    printf("  • Car model: L=%.1fm, φ_max=%d°, R_min=%.1fm\n",
           L, (int)(phi_max * 180.0 / M_PI), R_min);
    printf("  • Reeds-Shepp allows reverse → shorter parking paths\n");
    printf("  • RRT finds nonholonomic parking maneuvers\n");
    printf("  • Tesla Autopark / Waymo use similar principles\n");
    printf("  • Path optimization reduces waypoint count\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    noplan_system_free(car);
    return 0;
}
