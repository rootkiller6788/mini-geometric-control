#ifndef LIE_MECHANICS_H
#define LIE_MECHANICS_H

#include "lie_core.h"
#include "lie_actions.h"
#include "lie_reduction.h"
#include "lie_integration.h"

/* ==============================================================
 * lie_mechanics.h -- Mechanical Systems on Lie Groups
 *
 * Pre-built mechanical system models on Lie groups:
 * rigid bodies, tops, pendulums, multi-body systems,
 * and controlled vehicles (spacecraft, quadrotors, AUVs).
 *
 * References:
 *   Bloch (2003) Nonholonomic Mechanics and Control
 *   Bullo & Lewis (2004) Geometric Control of Mechanical Systems
 *   Murray, Li, Sastry (1994) A Mathematical Introduction to
 *     Robotic Manipulation
 *   Marsden & Ratiu (1999) Introduction to Mechanics and Symmetry
 * ============================================================== */

/* --- Double Pendulum on SO(3) x SO(3) ---
 * Two rigid bodies coupled by spherical joints.
 * Lagrangian: L = 1/2 Omega1^T I1 Omega1 + 1/2 Omega2^T I2 Omega2
 *                 - m1 g h1(z) - m2 g h2(z)
 *
 * Ref: Marsden & Scheurle (1993) ZAMP 44:17-43
 */
typedef struct {
    LieGroupElement *R1;
    LieGroupElement *R2;
    LieVector *Omega1;
    LieVector *Omega2;
    LieInertiaTensor *I1;
    LieInertiaTensor *I2;
    double m1, m2;
    double l1, l2;
    double g;
    LieVector *Gamma1;
    LieVector *Gamma2;
} LieDoublePendulum;

LieDoublePendulum* lie_double_pendulum_create(double m1, double m2,
                                                double l1, double l2);
void lie_double_pendulum_free(LieDoublePendulum *dp);
void lie_double_pendulum_rhs(const LieDoublePendulum *dp,
                              LieVector *dOmega1, LieVector *dOmega2);

/* --- Free Rigid Body (Force-Free Top) ---
 * Euler equations with zero torque:
 *   I1 dOmega1/dt = (I2 - I3) Omega2 Omega3
 *   I2 dOmega2/dt = (I3 - I1) Omega3 Omega1
 *   I3 dOmega3/dt = (I1 - I2) Omega1 Omega2
 *
 * Integrable. Constants: energy T, angular momentum ||Pi||^2.
 *
 * Ref: Landau & Lifshitz (1976) Mechanics, SS37
 *      Arnold (1989) Mathematical Methods of Classical Mechanics
 */
typedef struct {
    LieMatrix *R;
    LieVector *Omega;
    LieInertiaTensor *I;
    double energy;
    double momentum_sq;
} LieFreeTop;

LieFreeTop* lie_free_top_create(double I1, double I2, double I3);
void lie_free_top_free(LieFreeTop *top);
void lie_free_top_init(LieFreeTop *top, double Omega1, double Omega2, double Omega3);
void lie_free_top_rhs(const LieFreeTop *top, LieVector *dOmega_dt);
void lie_free_top_analytical(const LieFreeTop *top, double t,
                              LieVector *Omega_out);

/* --- Satellite Attitude Dynamics ---
 * Model of a rigid spacecraft with control and disturbance torques.
 *
 * Ref: Wertz (1978) Spacecraft Attitude Determination and Control
 *      Wie (2008) Space Vehicle Dynamics and Control
 */
typedef struct {
    LieMatrix *R;
    LieVector *Omega;
    LieInertiaTensor *I;
    LieVector *tau_control;
    LieVector *tau_disturbance;
    double mu_earth;
    double orbit_radius;
    double n_orbit;
} LieSatelliteState;

LieSatelliteState* lie_satellite_create(double Ixx, double Iyy, double Izz);
void lie_satellite_free(LieSatelliteState *sat);
void lie_gravity_gradient_torque(const LieSatelliteState *sat,
                                  const LieVector *o3_body, LieVector *tau_gg);
void lie_satellite_rhs(const LieSatelliteState *sat, double t,
                        LieVector *dOmega_dt);

/* --- Quadrotor on SE(3) ---
 * State: (p, R, v, Omega) in SE(3) x R^6
 *
 * Dynamics: p_dot=v, v_dot=g*e3+(f/m)*R*e3, R_dot=R*Omega^hat,
 *           Omega_dot=I^{-1}(tau - Omega x I Omega)
 *
 * Ref: Lee, Leok, McClamroch (2010) CDC 2010
 *      Mellinger & Kumar (2011) ICRA 2011
 */
typedef struct {
    LieGroupElement *g;
    LieVector *p;
    LieMatrix *R;
    LieVector *v;
    LieVector *Omega;
    double mass;
    LieInertiaTensor *I;
    double thrust;
    LieVector *tau;
} LieQuadrotorState;

LieQuadrotorState* lie_quadrotor_create(double mass, double Ixx, double Iyy, double Izz);
void lie_quadrotor_free(LieQuadrotorState *qr);
void lie_quadrotor_rhs(const LieQuadrotorState *qr,
                        LieVector *p_dot, LieVector *v_dot,
                        LieVector *Omega_dot);
void lie_quadrotor_geometric_controller(
    const LieQuadrotorState *qr,
    const LieVector *p_d, const LieVector *v_d, const LieVector *a_d,
    const LieMatrix *R_d, const LieVector *Omega_d,
    double *thrust_out, LieVector *tau_out,
    double kp, double kd, double kR, double kOmega);

/* --- Satellite Trajectory --- */
typedef struct {
    int n_steps;
    double *t;
    LieGroupElement **attitude;
    LieVector **omega;
} LieSatelliteTrajectory;

LieSatelliteTrajectory* lie_satellite_simulate(
    const LieSatelliteState *sat, LieRHSFunction f,
    double t0, double tf, int n_steps);
void lie_satellite_trajectory_free(LieSatelliteTrajectory *traj);

/* --- Underwater Vehicle on SE(3) ---
 * Rigid body in fluid with added mass, damping, restoring forces.
 *
 * Ref: Fossen (2011) Handbook of Marine Craft Hydrodynamics
 */
typedef struct {
    LieMatrix *R;
    LieVector *p;
    LieVector *v;
    LieVector *Omega;
    LieMatrix *M_rb;
    LieMatrix *M_added;
    LieMatrix *D;
    LieVector *g_restore;
    double rho_water;
    double volume;
} LieUnderwaterVehicle;

LieUnderwaterVehicle* lie_underwater_vehicle_create(void);
void lie_underwater_vehicle_free(LieUnderwaterVehicle *uv);
void lie_underwater_vehicle_rhs(const LieUnderwaterVehicle *uv,
                                 const LieVector *thrust,
                                 LieVector *v_dot, LieVector *Omega_dot);

#endif /* LIE_MECHANICS_H */