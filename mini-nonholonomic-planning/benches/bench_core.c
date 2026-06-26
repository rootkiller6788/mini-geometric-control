/* ============================================================================
 * bench_core.c — Performance Benchmarks
 *
 * Measures execution time for key operations:
 *   - Lie bracket computation throughput
 *   - RRT extension rate (nodes/second)
 *   - PRM construction rate (edges/second)
 *   - Chow-Rashevskii rank computation scaling
 *   - Distribution involutivity check
 *
 * Provides timing data for complexity analysis.
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "../include/noplan_core.h"
#include "../include/noplan_geometry.h"
#include "../include/noplan_planning.h"
#include "../include/noplan_models.h"

static double get_time_ms(void) {
    return (double)clock() / CLOCKS_PER_SEC * 1000.0;
}

#define BENCH(name, iterations, block) do { \
    double t0 = get_time_ms(); \
    for (int _i = 0; _i < (iterations); _i++) { block } \
    double t1 = get_time_ms(); \
    double elapsed = t1 - t0; \
    printf("  %-30s %8.2f ms  (%8.1f ops/s)\n", \
           name, elapsed, (iterations) * 1000.0 / elapsed); \
} while(0)

int main(void) {
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║   Nonholonomic Planning — Performance Benchmarks     ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n\n");

    /* ================================================================
     * Lie Bracket Computation
     * ================================================================ */
    printf("━━━ Lie Bracket Computation ━━━\n\n");

    Config q = noplan_config_create(3, 'S');
    int n_bracket = 5000;

    BENCH("[g₁,g₂] (finite diff)", n_bracket, {
        LieBracket lb;
        noplan_lie_bracket(&q, noplan_unicycle_g1, noplan_unicycle_g2,
                           3, 1e-5, &lb);
    });

    /* ================================================================
     * Distribution Rank
     * ================================================================ */
    printf("\n━━━ Distribution Rank ━━━\n\n");

    NonholonomicSystem* uni = noplan_model_unicycle();
    int n_rank = 10000;

    BENCH("Distribution rank (MGS)", n_rank, {
        noplan_distribution_rank(uni->distribution, &q);
    });

    /* ================================================================
     * QR Rank Decomposition
     * ================================================================ */
    printf("\n━━━ QR Rank Decomposition ━━━\n\n");

    double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};  /* Rank-2 approx */
    int n_qr = 20000;

    BENCH("QR rank (3x3)", n_qr, {
        noplan_qr_rank(A, 3, 3, 1e-8);
    });

    /* ================================================================
     * Configuration Distance
     * ================================================================ */
    printf("\n━━━ Configuration Distance ━━━\n\n");

    Config a = noplan_config_create(3, 'E');
    Config b = noplan_config_create(3, 'E');
    b.q[0] = 5.0; b.q[1] = 5.0;
    int n_dist = 500000;

    BENCH("Euclidean distance (3D)", n_dist, {
        noplan_config_distance(&a, &b);
    });

    BENCH("SE(2) distance", n_dist, {
        noplan_config_distance_se2(&a, &b);
    });

    /* ================================================================
     * SE(2) Operations
     * ================================================================ */
    printf("\n━━━ SE(2) Lie Group ━━━\n\n");

    SE2 g = se2_create(1.0, 2.0, M_PI / 3);
    SE2 h = se2_create(0.5, -0.5, M_PI / 6);
    int n_se2 = 100000;

    BENCH("SE(2) compose", n_se2, {
        se2_compose(&g, &h);
    });

    BENCH("SE(2) inverse", n_se2, {
        se2_inverse(&g);
    });

    BENCH("SE(2) exp", n_se2, {
        SE2 out;
        se2_exp(1.0, 0.5, 0.3, &out);
    });

    /* ================================================================
     * RRT Extension
     * ================================================================ */
    printf("\n━━━ RRT Planning ━━━\n\n");

    Config q_goal = noplan_config_create(3, 'S');
    q_goal.q[0] = 3.0; q_goal.q[1] = 2.0;

    noplan_system_reset(uni, &q);
    PlanningProblem* prob = noplan_problem_create(uni, &q, &q_goal);
    RRTPlanner* rrt = noplan_rrt_create(uni, prob, 0.3);

    int n_extend = 100;
    double t0_ext = get_time_ms();
    for (int i = 0; i < n_extend; i++) {
        noplan_rrt_extend(rrt);
    }
    double t1_ext = get_time_ms();
    printf("  RRT: %d extensions in %.2f ms (%.1f nodes/s, %d total nodes)\n",
           n_extend, t1_ext - t0_ext,
           n_extend * 1000.0 / (t1_ext - t0_ext), rrt->n_nodes);

    noplan_rrt_free(rrt);
    noplan_problem_free(prob);

    /* ================================================================
     * Chow-Rashevskii Test
     * ================================================================ */
    printf("\n━━━ Chow-Rashevskii Test ━━━\n\n");

    int n_chow = 500;
    BENCH("LARC test (unicycle, depth 2)", n_chow, {
        noplan_larc_holds(uni, &q, 1e-5);
    });

    /* ================================================================
     * System Simulation
     * ================================================================ */
    printf("\n━━━ System Simulation ━━━\n\n");

    double u[2] = {1.0, 0.5};
    int n_step = 100000;

    noplan_system_reset(uni, &q);
    BENCH("Euler step (unicycle)", n_step, {
        noplan_system_step(uni, u, 0.01);
    });

    /* ================================================================
     * Summary
     * ================================================================ */
    printf("\n━━━ Summary ━━━\n\n");
    printf("  All benchmarks complete.\n");
    printf("  Key bottlenecks:\n");
    printf("    - Lie bracket (finite diff): O(n²) per evaluation\n");
    printf("    - RRT extension: O(n_nodes·n + n_obstacles)\n");
    printf("    - Chow-Rashevskii: exponential in max_depth\n");
    printf("    - PRM construction: O(N²·steering_cost)\n");

    noplan_system_free(uni);
    return 0;
}
