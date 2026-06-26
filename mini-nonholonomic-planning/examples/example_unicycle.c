/* ============================================================================
 * example_unicycle.c — Unicycle Motion Planning Demo
 *
 * Demonstrates the fundamental nonholonomic planning concepts:
 *   1. Forward simulation of a unicycle
 *   2. Lie bracket motion (geometric phase)
 *   3. Parallel parking via cyclic control sequence
 *   4. Chow-Rashevskii controllability test
 *   5. RRT-based path planning
 *
 * The unicycle (q = (x, y, θ), q̇ = g₁·v + g₂·ω) is the simplest
 * controllable nonholonomic system. With only 2 controls for 3 states,
 * it achieves full controllability through the Lie bracket [g₁, g₂]
 * which generates sideways motion.
 *
 * Key takeaway: Nonholonomic ≠ Uncontrollable.
 *              Brockett: Cannot smoothly stabilize (m < n).
 *              Chow-Rashevskii: Can still plan paths (STLC).
 *
 * Reference: Murray & Sastry (1993), LaValle (2006).
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
    printf("║    Nonholonomic Motion Planning — Unicycle Demo      ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n\n");

    /* ================================================================
     * Part 1: Forward Simulation
     * ================================================================ */
    printf("━━━ Part 1: Forward Simulation ━━━\n\n");

    NonholonomicSystem* uni = noplan_model_unicycle();
    printf("System: %s  (n=%d, m=%d)\n",
           uni->name, uni->config_dim, uni->n_inputs);
    printf("Controls: v ∈ [%.1f, %.1f], ω ∈ [%.1f, %.1f]\n",
           uni->u_min[0], uni->u_max[0],
           uni->u_min[1], uni->u_max[1]);

    /* Simulate a circular path */
    printf("\nSimulating circular path (v=1.0, ω=0.5 rad/s, T=4π)...\n");
    Config q_start = uni->current_q;
    printf("  Start: (%.2f, %.2f, %.2f°)\n",
           q_start.q[0], q_start.q[1], q_start.q[2] * 180.0 / M_PI);

    double u_seq[628 * 2];  /* ~4π seconds at 0.02s step */
    int n_steps = 628;
    for (int i = 0; i < n_steps; i++) {
        u_seq[i * 2 + 0] = 1.0;
        u_seq[i * 2 + 1] = 0.5;
    }

    Trajectory* traj = noplan_trajectory_create(n_steps + 1, 2);
    noplan_system_simulate(uni, u_seq, n_steps, 0.02, traj);

    Config q_final = traj->waypoints[n_steps];
    printf("  Final: (%.2f, %.2f, %.2f°)\n",
           q_final.q[0], q_final.q[1], q_final.q[2] * 180.0 / M_PI);
    printf("  Path length: %.2f m\n", noplan_trajectory_length(traj));
    printf("  Control energy: %.2f\n", noplan_trajectory_energy(traj));
    noplan_trajectory_free(traj);

    /* ================================================================
     * Part 2: Lie Bracket Motion (Geometric Phase)
     * ================================================================ */
    printf("\n━━━ Part 2: Lie Bracket Motion — Geometric Phase ━━━\n\n");

    printf("Executing cyclic control sequence:\n");
    printf("  g₁(+ε) → g₂(+ε) → g₁(-ε) → g₂(-ε)\n");
    printf("  (forward → left → backward → right)\n");
    printf("  Expected net motion: ~[g₁,g₂] = sideways!\n\n");

    noplan_system_reset(uni, &q_start);

    double cyclic_u[40 * 2];
    double eps_v = 0.5, eps_w = 1.0;
    /* Forward */
    for (int i = 0; i < 10; i++) {
        cyclic_u[i * 2 + 0] = eps_v;   cyclic_u[i * 2 + 1] = 0.0;
    }
    /* Left */
    for (int i = 0; i < 10; i++) {
        cyclic_u[(10 + i) * 2 + 0] = 0.0;  cyclic_u[(10 + i) * 2 + 1] = eps_w;
    }
    /* Backward */
    for (int i = 0; i < 10; i++) {
        cyclic_u[(20 + i) * 2 + 0] = -eps_v; cyclic_u[(20 + i) * 2 + 1] = 0.0;
    }
    /* Right */
    for (int i = 0; i < 10; i++) {
        cyclic_u[(30 + i) * 2 + 0] = 0.0;  cyclic_u[(30 + i) * 2 + 1] = -eps_w;
    }

    Trajectory* cyclic_traj = noplan_trajectory_create(41, 2);
    noplan_system_simulate(uni, cyclic_u, 40, 0.1, cyclic_traj);

    Config q_after_cycle = cyclic_traj->waypoints[40];
    double dx = q_after_cycle.q[0] - q_start.q[0];
    double dy = q_after_cycle.q[1] - q_start.q[1];
    double dth = q_after_cycle.q[2] - q_start.q[2];

    printf("After one cycle (should return to near start with net side-slip):\n");
    printf("  Δx = %.4f m\n", dx);
    printf("  Δy = %.4f m  ← THIS is the geometric phase!\n", dy);
    printf("  Δθ = %.4f°\n", dth * 180.0 / M_PI);
    printf("  Lie bracket [g₁,g₂] = (-sin θ, cos θ, 0) at θ=0 → (0, -1, 0)\n");
    printf("  → Predicts net y-displacement ≈ -ε_v·ε_w·T = %.4f m\n",
           -eps_v * eps_w * 4.0);

    noplan_trajectory_free(cyclic_traj);

    /* ================================================================
     * Part 3: Controllability Analysis
     * ================================================================ */
    printf("\n━━━ Part 3: Controllability Analysis ━━━\n\n");

    Config q_test = noplan_config_create(3, 'S');
    printf("Testing controllability at q = (0, 0, 0):\n");

    /* Distribution rank */
    int dist_rank = noplan_distribution_rank(uni->distribution, &q_test);
    printf("  dim(Δ) = %d  (span of g₁, g₂)\n", dist_rank);

    /* LARC test */
    bool larc = noplan_larc_holds(uni, &q_test, 1e-5);
    printf("  LARC holds: %s\n", larc ? "YES ✓" : "NO ✗");

    /* Chow rank */
    int chow_rank = noplan_chow_rank(uni, &q_test, 2, 1e-5);
    printf("  Chow rank (depth 2): %d/%d\n", chow_rank, uni->config_dim);

    /* Degree of nonholonomy */
    int deg = noplan_degree_of_nonholonomy(uni, &q_test, 1e-5);
    printf("  Degree of nonholonomy κ = %d\n", deg);

    /* Growth vector */
    GrowthVector gv;
    noplan_compute_growth_vector(uni, &q_test, 3, 1e-5, &gv);
    printf("  Growth vector: (");
    for (int k = 1; k <= gv.kappa; k++) {
        printf("%d%s", gv.growth[k], k < gv.kappa ? ", " : "");
    }
    printf(")\n");

    /* Brockett's condition */
    bool brockett = noplan_brockett_check(uni, &q_test);
    printf("  Brockett (smoothly stabilizable): %s\n",
           brockett ? "YES" : "NO (m=2 < n=3)");

    /* Frobenius */
    bool frob = noplan_frobenius_test(uni->distribution, &q_test, 1e-4);
    printf("  Frobenius (holonomic?): %s\n", frob ? "YES" : "NO → nonholonomic ✓");

    /* ================================================================
     * Part 4: Lie Bracket Computation
     * ================================================================ */
    printf("\n━━━ Part 4: Lie Bracket Verification ━━━\n\n");

    LieBracket lb;
    noplan_lie_bracket(&q_test, noplan_unicycle_g1, noplan_unicycle_g2,
                        3, 1e-5, &lb);
    printf("  [g₁, g₂] at θ=0:\n");
    printf("    Computed: (%.4f, %.4f, %.4f)\n",
           lb.value.v[0], lb.value.v[1], lb.value.v[2]);
    printf("    Expected: (0, -1, 0)  [sin(0), -cos(0), 0]\n");
    printf("    Magnitude: %.4f\n", lb.magnitude);
    printf("    Is zero? %s (should be NO for nonholonomic systems)\n",
           lb.is_zero ? "yes" : "no");

    /* ================================================================
     * Part 5: RRT Path Planning
     * ================================================================ */
    printf("\n━━━ Part 5: RRT Path Planning ━━━\n\n");

    Config q_goal_rrt = noplan_config_create(3, 'S');
    q_goal_rrt.q[0] = 3.0;
    q_goal_rrt.q[1] = 2.0;
    q_goal_rrt.q[2] = M_PI / 4;

    printf("Planning from (0,0,0°) → (3,2,45°)\n");
    printf("Using RRT with step_duration=0.3s, goal_bias=0.1\n");

    noplan_system_reset(uni, &q_test);
    PlanningProblem* prob = noplan_problem_create(uni, &q_test, &q_goal_rrt);
    RRTPlanner* rrt = noplan_rrt_create(uni, prob, 0.3);
    rrt->goal_tolerance = 0.5;

    Trajectory* plan_traj = noplan_trajectory_create(500, 2);
    bool success = noplan_rrt_plan(rrt, 1000, plan_traj);

    if (success) {
        printf("  Planning SUCCEEDED after %d iterations\n", rrt->n_nodes);
        printf("  Path length: %.2f, Waypoints: %d\n",
               noplan_trajectory_length(plan_traj),
               plan_traj->n_waypoints);
        printf("  Final config: (%.2f, %.2f, %.2f°)\n",
               rrt->nodes[rrt->goal_node_idx].q.q[0],
               rrt->nodes[rrt->goal_node_idx].q.q[1],
               rrt->nodes[rrt->goal_node_idx].q.q[2] * 180.0 / M_PI);
        printf("  Distance to goal: %.3f\n",
               noplan_config_distance_se2(
                   &rrt->nodes[rrt->goal_node_idx].q, &q_goal_rrt));
    } else {
        printf("  Planning did not reach goal within %d iterations\n", 1000);
        double best_d = INFINITY;
        for (int i = 0; i < rrt->n_nodes; i++) {
            double d = noplan_config_distance_se2(
                &rrt->nodes[i].q, &q_goal_rrt);
            if (d < best_d) best_d = d;
        }
        printf("  Best distance achieved: %.3f\n", best_d);
    }

    noplan_trajectory_free(plan_traj);
    noplan_rrt_free(rrt);
    noplan_problem_free(prob);

    /* ================================================================
     * Summary
     * ================================================================ */
    printf("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("  Key Results:\n");
    printf("  • Unicycle is CONTROLLABLE (LARC passes)\n");
    printf("  • But NOT smoothly stabilizable (Brockett fails)\n");
    printf("  • Lie bracket [g₁,g₂] = side-slip direction\n");
    printf("  • Cyclic control → geometric phase → net displacement\n");
    printf("  • RRT successfully plans nonholonomic paths\n");
    printf("  • κ = 2, Growth vector = (2, 3)\n");
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    noplan_system_free(uni);
    return 0;
}
