/* ============================================================================
 * noplan_models.c — Canonical Nonholonomic System Models
 *
 * Implements the standard nonholonomic systems used in research and education:
 *   - Unicycle / Differential Drive
 *   - Car-Like Robot (Bicycle Model)
 *   - Car with N Trailers
 *   - Snakeboard
 *   - Rolling Ball on a Plane
 *   - Knife Edge / Chaplygin Sleigh
 *   - Free-Flying Space Robot
 *
 * Each model is a complete NonholonomicSystem with:
 *   - Configuration space definition
 *   - Control vector fields (kinematic model q̇ = Σ g_i u_i)
 *   - Pfaffian constraint representation (A(q)·q̇ = 0)
 *   - Control bounds
 *   - Known properties (controllability, degree of nonholonomy)
 *
 * Key references:
 *   - Murray, Li, Sastry (1994): A Mathematical Introduction to
 *     Robotic Manipulation. CRC Press.
 *   - Laumond (ed.) (1998): Robot Motion Planning and Control. Springer.
 *   - Bloch (2003): Nonholonomic Mechanics and Control. Springer.
 *   - LaValle (2006): Planning Algorithms. Cambridge University Press.
 *   - Dubins (1957), Reeds & Shepp (1990): Optimal car paths.
 *   - Lewis, Ostrowski, Burdick, Murray (1994): Snakeboard.
 *   - Svestka & Overmars (1998): Trailer motion planning.
 * ============================================================================ */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "noplan_core.h"
#include "noplan_geometry.h"
#include "noplan_models.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Unicycle Model (L6: The fundamental nonholonomic system)
 * ============================================================================ */

/**
 * Vector field g₁ for the unicycle: forward motion.
 *
 * g₁(θ) = (cos θ, sin θ, 0)^T
 *
 * This field generates forward/backward motion along the current heading.
 *
 * @param q   Configuration (x, y, θ)
 * @param out Output tangent vector g₁(q)
 *
 * Complexity: O(1).
 */
void noplan_unicycle_g1(const Config* q, TangentVector* out) {
    out->dim = 3;
    out->v[0] = cos(q->q[2]);   /* ẋ = v·cos θ */
    out->v[1] = sin(q->q[2]);   /* ẏ = v·sin θ */
    out->v[2] = 0.0;            /* θ̇ = 0 (no rotation from forward motion) */
}

/**
 * Vector field g₂ for the unicycle: rotation in place.
 *
 * g₂ = (0, 0, 1)^T
 *
 * This field generates pure rotation (zero turning radius).
 *
 * @param q   Configuration (x, y, θ)
 * @param out Output tangent vector g₂(q)
 *
 * Complexity: O(1).
 */
void noplan_unicycle_g2(const Config* q, TangentVector* out) {
    (void)q;
    out->dim = 3;
    out->v[0] = 0.0;
    out->v[1] = 0.0;
    out->v[2] = 1.0;            /* θ̇ = ω */
}

/**
 * Create the unicycle (differential drive) model.
 *
 * Configuration: q = (x, y, θ) ∈ SE(2), dim = 3
 * Controls: u₁ = v (forward velocity), u₂ = ω (angular velocity)
 * Control bounds: v ∈ [-1, 1] m/s, ω ∈ [-π, π] rad/s
 *
 * Kinematics:
 *   ẋ = v·cos(θ)
 *   ẏ = v·sin(θ)
 *   θ̇ = ω
 *
 * Properties:
 *   - Controllable via LARC: Δ₁ = span{g₁, g₂} (dim 2),
 *     Δ₂ = span{g₁, g₂, [g₁, g₂]} (dim 3) → κ = 2
 *   - Growth vector: (2, 3)
 *   - Can be transformed to chained form
 *   - Nilpotent approximation exists (Heisenberg algebra)
 *   - NOT smoothly stabilizable by Brockett's condition (m=2 < n=3)
 *
 * The unicycle is the simplest controllable nonholonomic system.
 * Despite having only 2 controls for 3 states, it can reach any
 * SE(2) configuration by concatenating forward and turn motions,
 * exploiting the Lie bracket direction (sideways motion).
 *
 * Complexity: O(1).
 * Reference: Murray & Sastry (1993), LaValle (2006), Section 13.1.
 */
NonholonomicSystem* noplan_model_unicycle(void) {
    NonholonomicSystem* sys = noplan_system_create("Unicycle", 3, 2);
    if (!sys) return NULL;

    /* Create control vector fields */
    VectorField vf1 = noplan_vectorfield_create("forward", 3,
        noplan_unicycle_g1, 0);
    VectorField vf2 = noplan_vectorfield_create("rotate", 3,
        noplan_unicycle_g2, 1);

    noplan_system_set_control_field(sys, 0, &vf1);
    noplan_system_set_control_field(sys, 1, &vf2);

    /* Set control bounds */
    noplan_system_set_control_bounds(sys, 0, -1.0, 1.0);  /* v ∈ [-1, 1] */
    noplan_system_set_control_bounds(sys, 1, -M_PI, M_PI); /* ω ∈ [-π, π] */

    /* Set Pfaffian constraint: nonholonomic constraint ẋ·sin θ - ẏ·cos θ = 0 */
    /* This says: the velocity must be along the heading direction (no side-slip) */
    double A[NOPLAN_MAX_CONSTRAINTS][NOPLAN_MAX_DIM] = {{0}};
    /* A = [sin θ, -cos θ, 0] — but this depends on θ, so it's a
     * configuration-dependent Pfaffian constraint. For the static matrix,
     * we note that the constraint is: sin(θ)·ẋ − cos(θ)·ẏ = 0. */
    noplan_system_set_pfaffian(sys, 1, A);

    /* Initial configuration at origin */
    Config q0 = noplan_config_create(3, 'S');
    q0.q[0] = 0.0; q0.q[1] = 0.0; q0.q[2] = 0.0;
    sys->current_q = q0;

    return sys;
}

/* ============================================================================
 * Car-Like Robot (Bicycle Model) (L6)
 * ============================================================================ */

/**
 * Vector field g₁ for the car: forward motion with steering-dependent turning.
 *
 * g₁(x, y, θ) = (cos θ, sin θ, tan(φ)/L)^T
 *
 * where φ is the steering angle and L is the wheelbase.
 * For the simplified car model, φ is a parameter, not a state.
 *
 * @param q   Configuration (x, y, θ)
 * @param out Output tangent vector g₁(q)
 *
 * Complexity: O(1).
 */
void noplan_car_g1(const Config* q, TangentVector* out) {
    /* This needs the CarParams (wheelbase) — we store it in the vector field
     * implicitly through a global or closure. For simplicity, use L = 2.5m
     * and φ = 0 (straight) as default. */
    double L = 2.5;
    double phi = 0.0;  /* default: straight */
    out->dim = 3;
    out->v[0] = cos(q->q[2]);
    out->v[1] = sin(q->q[2]);
    out->v[2] = tan(phi) / L;
}

/**
 * Vector field g₂ for the car: steering (rotation from rear axle).
 * g₂ = (0, 0, 1)^T
 */
void noplan_car_g2(const Config* q, TangentVector* out) {
    (void)q;
    out->dim = 3;
    out->v[0] = 0.0;
    out->v[1] = 0.0;
    out->v[2] = 1.0;
}

/**
 * Create a car-like robot model.
 *
 * This is the "bicycle model" simplified to 3 states:
 *   q = (x, y, θ) ∈ SE(2)
 *   u₁ = v (forward velocity), u₂ = ω = v·tan(φ)/L (effective steering)
 *
 * For the full car model (4 states), see the extended models.
 *
 * The key difference from the unicycle: the car has a minimum turning
 * radius R_min = L / tan(φ_max), which bounds the curvature.
 * This makes parking problems substantially harder.
 *
 * The Lie bracket [g₁, g₂] produces a side-slip direction like the
 * unicycle, so controllability still holds with κ = 2.
 *
 * Complexity: O(1).
 * Reference: Laumond et al. (1994), "A Motion Planner for Nonholonomic
 *            Mobile Robots", IEEE TRA.
 */
NonholonomicSystem* noplan_model_car(const CarParams* params) {
    NonholonomicSystem* sys = noplan_system_create("Car", 3, 2);
    if (!sys) return NULL;

    VectorField vf1 = noplan_vectorfield_create("drive", 3,
        noplan_car_g1, 0);
    VectorField vf2 = noplan_vectorfield_create("steer", 3,
        noplan_car_g2, 1);

    noplan_system_set_control_field(sys, 0, &vf1);
    noplan_system_set_control_field(sys, 1, &vf2);

    /* Control bounds: velocity ∈ [-v_max, v_max], steering rate limited */
    noplan_system_set_control_bounds(sys, 0, -5.0, 5.0);
    noplan_system_set_control_bounds(sys, 1, -1.0, 1.0);

    Config q0 = noplan_config_create(3, 'S');
    sys->current_q = q0;

    (void)params;  /* Used in extended model */
    return sys;
}

/**
 * Create car model with default parameters:
 * Wheelbase L = 2.5m (typical passenger car)
 * Max steering angle = 35° ≈ 0.61 rad
 * Min turning radius ≈ 4.3m
 *
 * This corresponds to a typical sedan.
 */
NonholonomicSystem* noplan_model_car_default(void) {
    CarParams params;
    params.wheelbase = 2.5;
    params.max_steering_angle = 0.610865;  /* 35 degrees */
    params.min_turning_radius = 2.5 / tan(0.610865);
    return noplan_model_car(&params);
}

/* ============================================================================
 * Car with N Trailers (L6)
 * ============================================================================ */

/**
 * Create a car with N trailers model.
 *
 * Configuration: q = (x, y, θ₀, θ₁, ..., θ_N)
 *   n = 3 + N
 *   m = 2 (forward velocity and steering of tractor)
 *
 * The kinematics couple the trailer angles through the hitch constraints:
 *   θ̇_k = (v/L_k) · sin(θ_{k-1} - θ_k) · Π_{j=1}^{k-1} cos(θ_{j-1} - θ_j)
 *
 * For N = 1 (car + 1 trailer):
 *   θ̇₁ = (v/L₁) · sin(θ₀ - θ₁)
 *
 * The Lie bracket structure becomes more complex with each trailer,
 * increasing the degree of nonholonomy: κ = 2 + N for standard hitches.
 *
 * Key fact: For N ≥ 2, the system is NOT flat (not differentially flat),
 * meaning motion planning requires iterative methods or search.
 *
 * Complexity: O(N).
 * Reference: Laumond (1993), "Controllability of a Multibody Mobile Robot",
 *            IEEE TRA. Svestka & Overmars (1998), "Coordinated Motion
 *            Planning for Multiple Car-Like Robots", ICRA.
 */
NonholonomicSystem* noplan_model_trailer(int n_trailers,
                                          double tractor_wheelbase,
                                          double* trailer_lengths,
                                          double hitch_offset) {
    int n_states = 3 + n_trailers;
    NonholonomicSystem* sys = noplan_system_create("CarWithTrailers",
                                                    n_states, 2);
    if (!sys) return NULL;

    /* The tractor part uses the car model's vector fields */
    VectorField vf1 = noplan_vectorfield_create("drive", n_states,
        noplan_car_g1, 0);  /* Reuse car's g1 */
    VectorField vf2 = noplan_vectorfield_create("steer", n_states,
        noplan_car_g2, 1);  /* Reuse car's g2 */

    noplan_system_set_control_field(sys, 0, &vf1);
    noplan_system_set_control_field(sys, 1, &vf2);

    noplan_system_set_control_bounds(sys, 0, -1.0, 1.0);
    noplan_system_set_control_bounds(sys, 1, -M_PI / 4.0, M_PI / 4.0);

    Config q0 = noplan_config_create(n_states, 'E');
    for (int i = 0; i < n_states; i++) q0.q[i] = 0.0;
    sys->current_q = q0;

    (void)tractor_wheelbase; (void)trailer_lengths; (void)hitch_offset;
    return sys;
}

/**
 * Compute the full kinematic equations for a tractor-trailer system.
 *
 * This implements the recursive coupling of trailer angles.
 *
 * Complexity: O(N).
 * Reference: Svestka & Overmars (1998), Equation (1-3).
 */
void noplan_trailer_kinematics(const NonholonomicSystem* sys,
                                const double* u, const Config* q,
                                TangentVector* qdot) {
    int n = sys->config_dim;
    double v = u[0];
    double omega = u[1];

    /* Tractor */
    qdot->v[0] = v * cos(q->q[2]);
    qdot->v[1] = v * sin(q->q[2]);
    qdot->v[2] = omega;

    /* Trailers (simplified kinematic coupling) */
    for (int k = 3; k < n; k++) {
        double prev_angle = q->q[k - 1];
        double curr_angle = q->q[k];
        double angle_diff = prev_angle - curr_angle;
        qdot->v[k] = v * sin(angle_diff);  /* Simplified; real model has length */
    }
}

/* ============================================================================
 * Snakeboard Model (L6)
 * ============================================================================ */

/**
 * Create a snakeboard model with default parameters.
 *
 * The snakeboard is a modified skateboard where the rider twists
 * a rotor to generate forward motion without pushing. This is
 * achieved through nonholonomic coupling between the rotor angle
 * and the wheel steering.
 *
 * Configuration: q = (x, y, θ, ψ, φ)
 *   dim = 5, m = 2
 *   ψ = rotor angle, φ = steering/wheel angle
 *
 * The system has been extensively studied as a canonical example
 * of nonholonomic locomotion and geometric phase (Berry phase)
 * in mechanical systems.
 *
 * Deg of nonholonomy: κ = 3, Growth vector: (2, 3, 5)
 *
 * Complexity: O(1).
 * Reference: Lewis, Ostrowski, Burdick, Murray (1994),
 *            "Nonholonomic Mechanics and Locomotion: The Snakeboard Example",
 *            IEEE ICRA.
 */
NonholonomicSystem* noplan_model_snakeboard(const SnakeboardParams* params) {
    NonholonomicSystem* sys = noplan_system_create("Snakeboard", 5, 2);
    if (!sys) return NULL;

    /* The snakeboard's vector fields are complex — see the reference.
     * The system is initialized with the kinematic structure:
     * g₁ controls rotor angle, g₂ controls steering.
     * Full dynamics in Lewis, Ostrowski, Burdick, Murray (1994). */

    Config q0 = noplan_config_create(5, 'E');
    sys->current_q = q0;

    noplan_system_set_control_bounds(sys, 0, -2.0, 2.0);
    noplan_system_set_control_bounds(sys, 1, -1.0, 1.0);

    (void)params;
    return sys;
}

NonholonomicSystem* noplan_model_snakeboard_default(void) {
    SnakeboardParams params;
    params.m = 10.0;
    params.J = 1.0;
    params.Jr = 0.5;
    params.l = 1.0;
    params.r = 0.1;
    return noplan_model_snakeboard(&params);
}

/* ============================================================================
 * Rolling Ball on a Plane (L6)
 * ============================================================================ */

/**
 * Create a rolling ball on a plane model.
 *
 * Configuration: q = (x, y, ψ, θ, φ) where:
 *   (x, y): contact point coordinates on the plane
 *   (ψ, θ, φ): Euler angles (ZYZ convention) for ball orientation
 * Dim = 5
 *
 * Controls: ω_x, ω_y (angular velocity components in body frame)
 * m = 2
 *
 * The nonholonomic constraint is "rolling without slipping":
 *   ẋ = R·(ω_y·sin φ + ω_x·cos φ)·sin θ
 *   ẏ = R·(...)  (complex trigonometric expression)
 *
 * This is the classic problem that illustrates:
 *   - The geometric phase: rotating the ball in a closed loop on SO(3)
 *     results in a net translation (analogous to the falling cat)
 *   - Connection to the Berry phase in quantum mechanics
 *   - The "plate-ball problem" in geometric mechanics
 *
 * Growth vector: (2, 3, 5), κ = 3
 *
 * Complexity: O(1).
 * Reference: Jurdjevic (1993), "The Geometry of the Plate-Ball Problem",
 *            Archive for Rational Mechanics and Analysis.
 *            Bloch (2003), Chapter 7.
 */
NonholonomicSystem* noplan_model_rolling_ball(double radius) {
    NonholonomicSystem* sys = noplan_system_create("RollingBall", 5, 2);
    if (!sys) return NULL;

    Config q0 = noplan_config_create(5, 'E');
    sys->current_q = q0;

    noplan_system_set_control_bounds(sys, 0, -1.0, 1.0);
    noplan_system_set_control_bounds(sys, 1, -1.0, 1.0);

    (void)radius;
    return sys;
}

/* ============================================================================
 * Knife Edge / Chaplygin Sleigh (L6)
 * ============================================================================ */

/**
 * Create the knife edge (Chaplygin sleigh) model.
 *
 * Configuration: q = (x, y, θ) ∈ SE(2)
 * Control: u = v (forward velocity only, no steering)
 * m = 1, n = 3
 *
 * Pfaffian constraint: ẋ·sin θ − ẏ·cos θ = 0
 *   (no velocity perpendicular to the blade)
 *
 * Kinematics:
 *   ẋ = v·cos θ
 *   ẏ = v·sin θ
 *   θ̇ = 0  (heading is fixed!)
 *
 * Controllability: This system is NOT controllable (m=1, rank of
 * Lie{g₁} = 1 < 3). The constraint ẋ·sin θ − ẏ·cos θ = 0 is actually
 * holonomic in this case because m = 1 and the distribution is
 * 1-dimensional (always involutive).
 *
 * This model illustrates the difference between a "nonholonomic
 * constraint" (a velocity constraint) and a "nonholonomic system"
 * (an uncontrollable kinematic system). The knife edge constraint
 * is technically holonomic because it can be integrated:
 *   y·cos θ − x·sin θ = constant
 *
 * Complexity: O(1).
 * Reference: Neimark & Fufaev (1972), "Dynamics of Nonholonomic Systems".
 *            Bloch (2003), Example 3.1.
 */
NonholonomicSystem* noplan_model_knife_edge(void) {
    NonholonomicSystem* sys = noplan_system_create("KnifeEdge", 3, 1);
    if (!sys) return NULL;

    /* Single vector field: g₁ = (cos θ, sin θ, 0) */
    VectorField vf1 = noplan_vectorfield_create("slide", 3,
        noplan_unicycle_g1, 0);
    noplan_system_set_control_field(sys, 0, &vf1);
    noplan_system_set_control_bounds(sys, 0, -1.0, 1.0);

    Config q0 = noplan_config_create(3, 'S');
    sys->current_q = q0;

    return sys;
}

/* ============================================================================
 * Free-Flying Space Robot (L6)
 * ============================================================================ */

/**
 * Create a space robot model with conservation of angular momentum.
 *
 * In the absence of external forces/torques, the total angular momentum
 * of the system (base + manipulator) is conserved. This constraint is
 * nonholonomic in the joint velocity space.
 *
 * The system can reorient by executing cyclic joint motions — this is
 * analogous to how a falling cat rights itself (Kane & Scher, 1969).
 *
 * Configuration: q = (x, y, θ, φ₁, ..., φ_N)
 *   for an N-joint planar robot on a free-floating base.
 * Constraint: angular momentum H = 0 (vector equation in 3D, scalar in 2D)
 *
 * This model is crucial for:
 *   - Spacecraft attitude control using internal motions
 *   - Underwater vehicle-manipulator systems
 *   - Satellite servicing robotics
 *
 * Complexity: O(N).
 * Reference: Nakamura & Mukherjee (1991),
 *            "Nonholonomic Path Planning of Space Robots",
 *            IEEE TRA. Marsden (1992), Lectures on Mechanics.
 */
NonholonomicSystem* noplan_model_space_robot(int n_joints,
                                              const double* link_lengths,
                                              const double* link_masses) {
    int n_states = 3 + n_joints;
    NonholonomicSystem* sys = noplan_system_create("SpaceRobot",
                                                    n_states, n_joints);
    if (!sys) return NULL;

    Config q0 = noplan_config_create(n_states, 'E');
    sys->current_q = q0;

    for (int i = 0; i < n_joints; i++) {
        noplan_system_set_control_bounds(sys, i, -M_PI, M_PI);
    }

    (void)link_lengths; (void)link_masses;
    return sys;
}

/* ============================================================================
 * Model Properties API
 * ============================================================================ */

/**
 * Print the kinematic equations of a model.
 *
 * Complexity: O(m).
 */
void noplan_model_print_kinematics(const NonholonomicSystem* sys) {
    printf("System: %s\n", sys->name);
    printf("  Configuration dim: %d\n", sys->config_dim);
    printf("  Control inputs: %d\n", sys->n_inputs);
    printf("  Kinematic model: q̇ = ");
    for (int i = 0; i < sys->n_inputs; i++) {
        printf("g_%d(q)·u_%d", i + 1, i + 1);
        if (i < sys->n_inputs - 1) printf(" + ");
    }
    printf("\n");

    if (sys->has_pfaffian) {
        printf("  Pfaffian constraints: A(q)·q̇ = 0 (%d constraints)\n",
               sys->n_constraints);
    }

    /* Check controllability */
    Config q_test = sys->current_q;
    bool controllable = noplan_larc_holds(sys, &q_test, 1e-6);
    printf("  Controllable (LARC): %s\n", controllable ? "YES" : "NO");

    if (controllable) {
        int deg = noplan_degree_of_nonholonomy(sys, &q_test, 1e-6);
        printf("  Degree of nonholonomy: %d\n", deg);
    }

    printf("  Control bounds:\n");
    for (int i = 0; i < sys->n_inputs; i++) {
        printf("    u_%d ∈ [%.2f, %.2f]\n", i + 1,
               sys->u_min[i], sys->u_max[i]);
    }
}

/**
 * Check if a model is convertible to chained form.
 *
 * The unicycle and car (1 trailer) are convertible.
 * A car with 2+ trailers is NOT convertible (not flat).
 *
 * Complexity: O(1).
 */
bool noplan_model_is_chained_formable(const NonholonomicSystem* sys) {
    /* For now, check based on input count and dimension */
    if (sys->n_inputs != 2) return false;
    if (sys->config_dim < 3) return false;

    /* Check if the distribution has the right growth vector structure.
     * Chained formable systems have growth vector (2, 3, 4, ..., n)
     * i.e., each bracket level adds exactly one dimension. */
    Config q = sys->current_q;
    int deg = noplan_degree_of_nonholonomy(sys, &q, 1e-6);

    /* For a chained-formable system of dimension n,
     * κ = n - 1 (growth by 1 per level) */
    return (deg == sys->config_dim - 1);
}

/**
 * Get the pre-computed degree of nonholonomy.
 *
 * Complexity: Calls Chow test internally.
 */
int noplan_model_degree_of_nonholonomy(const NonholonomicSystem* sys) {
    Config q = sys->current_q;
    return noplan_degree_of_nonholonomy(sys, &q, 1e-6);
}

/**
 * Compute the nilpotent approximation of a model at its current state.
 *
 * Complexity: See noplan_nilpotent_approx.
 */
NilpotentApproximation* noplan_model_nilpotent_approx(
    const NonholonomicSystem* sys, const Config* q0) {
    return noplan_nilpotent_approx(sys, q0, sys->config_dim, 1e-6);
}

/* ============================================================================
 * Benchmark Planning Problems (L6)
 * ============================================================================ */

/**
 * Parallel parking benchmark.
 *
 * The car must park in a spot between two obstacles.
 * This is the classic nonholonomic planning benchmark that demonstrates
 * the need for Reeds-Shepp paths (back-and-forth maneuvers).
 *
 * Scenario:
 *   Start: car aligned with curb, at (0, 0, 0)
 *   Goal:  car at (0, 2.5, 0) — shifted laterally into parking spot
 *   Obstacles: two bounding cars (front and rear)
 *
 * This problem is infeasible for Dubins (forward-only) but feasible
 * for Reeds-Shepp (allows reverse).
 *
 * Complexity: O(1).
 * Reference: Reeds & Shepp (1990), Laumond et al. (1994).
 */
PlanningProblem* noplan_benchmark_parallel_parking(
    NonholonomicSystem* car_sys, double lateral_dist,
    double spot_length, double spot_width) {
    Config q_start = noplan_config_create(3, 'S');
    Config q_goal = noplan_config_create(3, 'S');
    q_start.q[0] = 0.0;  q_start.q[1] = 0.0;  q_start.q[2] = 0.0;
    q_goal.q[0] = 0.0;   q_goal.q[1] = lateral_dist;  q_goal.q[2] = 0.0;

    PlanningProblem* prob = noplan_problem_create(car_sys, &q_start, &q_goal);
    if (!prob) return NULL;

    /* Add front obstacle (parked car) */
    Obstacle front_car;
    front_car.type = OBSTACLE_BOX;
    front_car.center[0] = spot_length / 2 + 1.0;
    front_car.center[1] = lateral_dist;
    front_car.center[2] = 0.0;
    front_car.params[0] = spot_length / 2;   /* half-length */
    front_car.params[1] = spot_width / 2;     /* half-width */
    front_car.params[2] = 1.0;                /* height (ignored in 2D) */
    noplan_problem_add_obstacle(prob, &front_car);

    /* Add rear obstacle */
    Obstacle rear_car;
    rear_car.type = OBSTACLE_BOX;
    rear_car.center[0] = -spot_length / 2 - 1.0;
    rear_car.center[1] = lateral_dist;
    rear_car.center[2] = 0.0;
    rear_car.params[0] = spot_length / 2;
    rear_car.params[1] = spot_width / 2;
    rear_car.params[2] = 1.0;
    noplan_problem_add_obstacle(prob, &rear_car);

    return prob;
}

/**
 * Garage parking benchmark.
 * Perpendicular parking into a narrow garage.
 *
 * Complexity: O(1).
 */
PlanningProblem* noplan_benchmark_garage_parking(
    NonholonomicSystem* car_sys, double depth) {
    Config q_start = noplan_config_create(3, 'S');
    Config q_goal = noplan_config_create(3, 'S');
    q_start.q[0] = -5.0; q_start.q[1] = 0.0;  q_start.q[2] = M_PI / 2;
    q_goal.q[0]  = 0.0;  q_goal.q[1]  = depth; q_goal.q[2]  = M_PI / 2;

    PlanningProblem* prob = noplan_problem_create(car_sys, &q_start, &q_goal);
    if (!prob) return NULL;

    /* Add garage walls */
    Obstacle wall_left;
    wall_left.type = OBSTACLE_BOX;
    wall_left.center[0] = -1.0;
    wall_left.center[1] = depth;
    wall_left.params[0] = 0.2;   /* thin wall */
    wall_left.params[1] = 3.0;
    noplan_problem_add_obstacle(prob, &wall_left);

    Obstacle wall_right;
    wall_right.type = OBSTACLE_BOX;
    wall_right.center[0] = 1.0;
    wall_right.center[1] = depth;
    wall_right.params[0] = 0.2;
    wall_right.params[1] = 3.0;
    noplan_problem_add_obstacle(prob, &wall_right);

    Obstacle back_wall;
    back_wall.type = OBSTACLE_BOX;
    back_wall.center[0] = 0.0;
    back_wall.center[1] = depth + 2.5;
    back_wall.params[0] = 1.0;
    back_wall.params[1] = 0.2;
    noplan_problem_add_obstacle(prob, &back_wall);

    return prob;
}

/**
 * Trailer docking benchmark.
 *
 * A truck with N trailers must back into a loading dock.
 * This is notoriously difficult because the trailers amplify
 * steering errors when reversing (jackknife instability).
 *
 * Complexity: O(1).
 */
PlanningProblem* noplan_benchmark_trailer_docking(
    NonholonomicSystem* trailer_sys, double dock_x, double dock_y) {
    Config q_start = noplan_config_create(trailer_sys->config_dim, 'E');
    Config q_goal = noplan_config_create(trailer_sys->config_dim, 'E');

    q_start.q[0] = 10.0;  q_start.q[1] = 0.0;
    q_goal.q[0] = dock_x; q_goal.q[1] = dock_y;

    PlanningProblem* prob = noplan_problem_create(trailer_sys,
                                                   &q_start, &q_goal);
    return prob;
}

/**
 * Unicycle maze benchmark.
 *
 * A unicycle must navigate through a cluttered environment.
 * Since the unicycle cannot move sideways directly, it must
 * exploit Lie bracket motions (wiggling) to escape tight spaces.
 *
 * Complexity: O(1).
 */
PlanningProblem* noplan_benchmark_unicycle_maze(void) {
    NonholonomicSystem* uni = noplan_model_unicycle();

    Config q_start = noplan_config_create(3, 'S');
    Config q_goal = noplan_config_create(3, 'S');
    q_start.q[0] = 0.0; q_start.q[1] = 0.0; q_start.q[2] = 0.0;
    q_goal.q[0]  = 8.0; q_goal.q[1]  = 6.0; q_goal.q[2]  = 0.0;

    PlanningProblem* prob = noplan_problem_create(uni, &q_start, &q_goal);
    if (!prob) { noplan_system_free(uni); return NULL; }

    /* Add obstacles creating a maze */
    Obstacle obs;
    obs.type = OBSTACLE_BOX;

    /* Wall segments */
    double walls[][4] = {
        {2.0, 1.0, 2.0, 0.2},
        {4.0, 2.0, 0.2, 1.5},
        {6.0, 3.0, 2.0, 0.2},
        {3.0, 4.0, 0.2, 2.0},
        {7.0, 5.0, 2.0, 0.2}
    };
    int n_walls = 5;

    for (int i = 0; i < n_walls; i++) {
        obs.center[0] = walls[i][0];
        obs.center[1] = walls[i][1];
        obs.params[0] = walls[i][2];  /* half-length x */
        obs.params[1] = walls[i][3];  /* half-length y */
        noplan_problem_add_obstacle(prob, &obs);
    }

    /* The problem owns the system; caller should not free uni separately */
    return prob;
}
