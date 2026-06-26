#ifndef NOPLAN_MODELS_H
#define NOPLAN_MODELS_H

#include "noplan_core.h"

/* ============================================================================
 * Nonholonomic Motion Planning — Canonical Models
 *
 * L6: Canonical Nonholonomic Problems
 *
 *   Each model is a complete specification of a nonholonomic system
 *   including configuration space, constraint distribution, vector fields,
 *   and known properties (controllability, chained-form convertibility, etc.)
 *
 *   Models follow the kinematic form: q̇ = Σ g_i(q) u_i
 * ============================================================================ */

/* ---------------------------------------------------------------------------
 * L6: Unicycle (Rolling Disk / Differential Drive)
 *
 * Configuration: q = (x, y, θ) ∈ SE(2)
 * Controls: u₁ = forward velocity v, u₂ = angular velocity ω
 *
 * Kinematics:
 *   ẋ = v · cos(θ)
 *   ẏ = v · sin(θ)
 *   θ̇ = ω
 *
 * Vector fields: g₁ = (cos θ, sin θ, 0), g₂ = (0, 0, 1)
 * Lie bracket:   g₃ = [g₁, g₂] = (-sin θ, cos θ, 0)
 *
 * Controllability: Rank(g₁, g₂, g₃) = 3 = n  →  LARC holds
 * Degree of nonholonomy: κ = 2 (need brackets of depth ≤ 2)
 * Growth vector: (2, 3)
 *
 * This is the simplest controllable nonholonomic system.
 * It models differential-drive robots, the Segway, and spherical rolling.
 * --------------------------------------------------------------------------- */
NonholonomicSystem* noplan_model_unicycle(void);

/* Gives the vector fields for the unicycle at configuration q */
void noplan_unicycle_g1(const Config* q, TangentVector* out);
void noplan_unicycle_g2(const Config* q, TangentVector* out);

/* ---------------------------------------------------------------------------
 * L6: Car-Like Robot (Bicycle Model / Simple Car / Dubins Car)
 *
 * Configuration: q = (x, y, θ) ∈ SE(2)  [simplified, no steering angle state]
 * Controls: u₁ = forward velocity v, u₂ = steering rate ω
 *
 * This is the "simple car" with bounded curvature:
 *   ẋ = v · cos(θ)
 *   ẏ = v · sin(θ)
 *   θ̇ = v · tan(φ) / L  where φ is steering angle, L is wheelbase
 *
 * Vector fields: g₁ = (cos θ, sin θ, tan(φ)/L), g₂ = (0, 0, 1)
 *
 * Lie bracket: g₃ = [g₁, g₂] = (sin θ/L·sec²(φ), -cos θ/L·sec²(φ), 0)
 *
 * For the Reeds-Shepp car, v ∈ {-1, +1} (forward/backward)
 * For the Dubins car, v ∈ {+1} only (forward only)
 *
 * This is the canonical model for autonomous vehicle parking.
 * --------------------------------------------------------------------------- */
typedef struct {
    double wheelbase;                 /* L: distance between front and rear axles */
    double max_steering_angle;        /* φ_max: maximum steering angle */
    double min_turning_radius;        /* L / tan(φ_max) */
} CarParams;

NonholonomicSystem* noplan_model_car(const CarParams* params);
NonholonomicSystem* noplan_model_car_default(void);  /* L=2.5m, φ_max=35° */

void noplan_car_g1(const Config* q, TangentVector* out);
void noplan_car_g2(const Config* q, TangentVector* out);

/* ---------------------------------------------------------------------------
 * L6: Car with N Trailers (Fire Truck / Truck-Trailer)
 *
 * Configuration: q = (x, y, θ₀, θ₁, ..., θ_N)
 *   where θ₀ is tractor heading, θ_i are trailer headings
 * Dim = 3 + N (one angle per trailer)
 * Controls: u₁ = forward velocity of tractor, u₂ = steering angle of tractor
 *
 * Constraint: each trailer is hitched to the previous vehicle.
 * The hitch point velocity must be compatible with the trailer orientation.
 *
 * This is the classic "nonholonomic motion planning" benchmark.
 * For N ≥ 2, the system is not convertible to chained form.
 * Controllability still holds via the kinematic model.
 *
 * Reference: Laumond (1993), Švestka & Overmars (1998)
 * --------------------------------------------------------------------------- */
NonholonomicSystem* noplan_model_trailer(int n_trailers,
                                          double tractor_wheelbase,
                                          double* trailer_lengths,
                                          double hitch_offset);

void noplan_trailer_kinematics(const NonholonomicSystem* sys,
                                const double* u, const Config* q,
                                TangentVector* qdot);

/* ---------------------------------------------------------------------------
 * L6: Snakeboard
 *
 * Configuration: q = (x, y, θ, ψ, φ)
 *   (x, y): position of center of mass
 *   θ: orientation of the board
 *   ψ: rotor angle (controls coupling between wheels)
 *   φ: wheel angle (steering)
 *
 * Controls: u₁ = ψ̇ (rotor velocity), u₂ = φ̇ (steering velocity)
 * Dim = 5, m = 2
 *
 * The snakeboard moves by exploiting the kinematic coupling between
 * rotor oscillation and wheel steering — a nonholonomic gait.
 * Degree of nonholonomy: κ = 3
 * Growth vector: (2, 3, 5)
 *
 * Reference: Lewis, Ostrowski, Burdick, Murray (1994)
 *            "Nonholonomic Mechanics and Locomotion: The Snakeboard Example"
 * --------------------------------------------------------------------------- */
typedef struct {
    double m;                         /* Total mass */
    double J;                         /* Moment of inertia */
    double Jr;                        /* Rotor inertia */
    double l;                         /* Half-length between wheel trucks */
    double r;                         /* Wheel radius */
} SnakeboardParams;

NonholonomicSystem* noplan_model_snakeboard(const SnakeboardParams* params);
NonholonomicSystem* noplan_model_snakeboard_default(void);

/* ---------------------------------------------------------------------------
 * L6: Rolling Ball on a Plane (Plate-Ball Problem)
 *
 * Configuration: q = (x, y, ψ, θ, φ) where
 *   (x, y): contact point on the plane
 *   (ψ, θ, φ): Euler angles of the ball
 * Dim = 5 (or 8 including angular velocities)
 *
 * The rolling without slipping constraint:
 *   ẋ = R·(ω_y·sin(φ) + ω_x·cos(φ))... (nonholonomic)
 *
 * This is a classic problem in geometric mechanics.
 * The ball-plate system is kinematic (no slipping) and exhibits
 * nonholonomic behavior related to the Berry phase.
 *
 * Controllability holds via the Lie algebra so(3) of rotations.
 * Growth vector: (2, 3, 5)
 *
 * Reference: Jurdjevic (1993), Bloch (2003)
 * --------------------------------------------------------------------------- */
NonholonomicSystem* noplan_model_rolling_ball(double radius);

/* ---------------------------------------------------------------------------
 * L6: Knife Edge / Sleigh (Chaplygin Sleigh)
 *
 * Configuration: q = (x, y, θ) ∈ SE(2)
 * Constraint: ẋ·sin(θ) − ẏ·cos(θ) = 0  (no sideways motion)
 * Controls: u = forward velocity
 *
 * Pfaffian constraint: A(q)·q̇ = 0 where A = [sin θ, -cos θ, 0]
 * This is a 1-control 3-state system. m=1 < n=3, but controllable?
 * No — this system is NOT controllable because m=1, and the Lie
 * algebra generated by a single vector field is 1-dimensional.
 *
 * This illustrates that m=2 is typically needed for controllability
 * of a 3D system (unless the system has special structure).
 *
 * Reference: Neimark & Fufaev (1972), Bloch (2003)
 * --------------------------------------------------------------------------- */
NonholonomicSystem* noplan_model_knife_edge(void);

/* ---------------------------------------------------------------------------
 * L6: Free-Flying Space Robot (Angular Momentum Conservation)
 *
 * Configuration: q = (x, y, θ, φ₁, ..., φ_N) for N joints
 *   (x, y, θ): base position and orientation
 *   φ_i: joint angles
 *
 * Constraint: total angular momentum = 0 (nonholonomic if base uncontrolled)
 * This is a nonholonomic constraint because it depends on velocities.
 *
 * The space robot can reorient itself by cyclic joint motions —
 * this is the "falling cat" problem (Kane & Scher, 1969).
 *
 * Reference: Nakamura & Mukherjee (1991), Marsden (1992)
 * --------------------------------------------------------------------------- */
NonholonomicSystem* noplan_model_space_robot(int n_joints,
                                              const double* link_lengths,
                                              const double* link_masses);

/* ---------------------------------------------------------------------------
 * Model Properties API
 * --------------------------------------------------------------------------- */

/**
 * Print the kinematic equations of a model in human-readable form.
 */
void    noplan_model_print_kinematics(const NonholonomicSystem* sys);

/**
 * Check if a model is in chained form (or can be converted).
 */
bool    noplan_model_is_chained_formable(const NonholonomicSystem* sys);

/**
 * Get the degree of nonholonomy for a pre-computed model.
 */
int     noplan_model_degree_of_nonholonomy(const NonholonomicSystem* sys);

/**
 * Compute the nilpotent approximation of a standard model at a point.
 */
NilpotentApproximation* noplan_model_nilpotent_approx(
    const NonholonomicSystem* sys, const Config* q0);

/* ---------------------------------------------------------------------------
 * L6: Canonical Parking / Benchmark Problems
 * --------------------------------------------------------------------------- */

/**
 * Parallel parking problem: car initially parallel to curb,
 * must park between two obstacles. Standard benchmark.
 *
 * Initial: (0, 0, 0)
 * Goal:    (0, d, 0) where d is lateral displacement
 *
 * For Reeds-Shepp: d = 2.0m, parking spot length = 1.5× car length
 */
PlanningProblem* noplan_benchmark_parallel_parking(
    NonholonomicSystem* car_sys, double lateral_dist,
    double spot_length, double spot_width);

/**
 * Garage parking (perpendicular): car must enter a narrow garage.
 * Initial: (0, 0, π/2) facing garage
 * Goal:    (0, d, π/2) inside the garage
 */
PlanningProblem* noplan_benchmark_garage_parking(
    NonholonomicSystem* car_sys, double depth);

/**
 * Trailer docking: truck with N trailers must back into a loading dock.
 * Classical hard problem for N ≥ 2.
 */
PlanningProblem* noplan_benchmark_trailer_docking(
    NonholonomicSystem* trailer_sys, double dock_x, double dock_y);

/**
 * Unicycle maze: navigate through a cluttered environment.
 * Must exploit side-slip (Lie bracket motion) to escape tight spaces.
 */
PlanningProblem* noplan_benchmark_unicycle_maze(void);

#endif /* NOPLAN_MODELS_H */
