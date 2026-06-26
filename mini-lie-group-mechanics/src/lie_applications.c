/* ==============================================================
 * lie_applications.c -- Real-World Applications of Lie Group Mechanics
 *
 * L7 Applications demonstrating the use of Lie group methods for:
 *   - Apollo Lunar Module attitude dynamics
 *   - SpaceX rocket stage separation on SE(3)
 *   - Quadrotor trajectory generation (NASA Mars helicopter)
 *   - GPS satellite constellation orbit modeling
 *   - Tesla vehicle inertial navigation
 *
 * Each application demonstrates a distinct use case of geometric
 * mechanics on Lie groups with real-world parameters.
 *
 * References:
 *   Wertz (1978) Spacecraft Attitude Determination and Control
 *   Lee, Leok, McClamroch (2010) Geometric Tracking Control on SE(3)
 *   Fossen (2011) Handbook of Marine Craft Hydrodynamics
 *   Murray, Li, Sastry (1994) A Mathematical Introduction to
 *     Robotic Manipulation
 * ============================================================== */

#include "lie_mechanics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==============================================================
 * Apollo Lunar Module Attitude Simulation
 *
 * The Apollo LM (1969) used a cold-gas reaction control system.
 * Inertia properties (approximate, kg*m^2):
 *   Ixx = 15000, Iyy = 35000, Izz = 25000
 *
 * The LM operated in 1/6 Earth gravity during descent.
 *
 * Ref: NASA Technical Note D-6846 "Apollo Experience Report -
 *      Mission Planning for Lunar Module Descent and Ascent"
 * ============================================================== */
void lie_apollo_lm_attitude_demo(void) {
    printf("\n=== Apollo Lunar Module Attitude Dynamics ===\n");
    printf("Simulating LM with Ixx=15000, Iyy=35000, Izz=25000 kg*m^2\n");

    LieSatelliteState *lm = lie_satellite_create(15000.0, 35000.0, 25000.0);
    /* Lunar orbit: radius ~1838 km (Moon radius), GM_moon = 4902.8 km^3/s^2 */
    lm->mu_earth = 4902.8;
    lm->orbit_radius = 1850.0;
    lm->n_orbit = sqrt(lm->mu_earth / (lm->orbit_radius*lm->orbit_radius*lm->orbit_radius));

    /* Initial angular velocity: 5 deg/s about z-axis */
    double omega0 = 5.0 * M_PI / 180.0;
    lm->Omega->data[2] = omega0;
    lm->Omega->data[0] = 0.01; /* small perturbation */

    /* Gravity gradient torque at 20 deg pitch */
    LieVector o3_body = {.data=(double[3]){sin(20.0*M_PI/180.0), 0, cos(20.0*M_PI/180.0)},.size=3};
    LieVector tau_gg; tau_gg.size = 3; tau_gg.data = (double[3]){0};
    lie_gravity_gradient_torque(lm, &o3_body, &tau_gg);

    printf("Gravity gradient torque (N*m): ");
    lie_vector_print(&tau_gg, "");

    /* Precession period */
    double I_ratio = (lm->I->inertia->data[4] - lm->I->inertia->data[0])
                    / lm->I->inertia->data[0];
    printf("Inertia ratio (Iyy-Ixx)/Ixx = %.4f\n", I_ratio);
    printf("Expected precession ~ %.1f seconds per revolution\n",
           2.0*M_PI/(omega0 * I_ratio));

    lie_satellite_free(lm);
}

/* ==============================================================
 * SpaceX Falcon 9 Stage Separation
 *
 * During stage separation, the second stage separates from the
 * first stage with a relative velocity. The dynamics on SE(3)
 * capture both attitude and translation.
 *
 * Falcon 9 parameters (approximate):
 *   Stage 2 mass: ~111,500 kg (fueled)
 *   Inertia: Ixx=50000, Iyy=400000, Izz=400000 kg*m^2
 *
 * Ref: SpaceX Falcon 9 User's Guide (2020)
 * ============================================================== */
void lie_spacex_separation_demo(void) {
    printf("\n=== SpaceX Falcon 9 Stage Separation (SE(3) model) ===\n");

    /* Create quadrotor state to represent the upper stage */
    LieQuadrotorState *stage2 = lie_quadrotor_create(111500.0, 50000.0, 400000.0, 400000.0);

    /* Initial separation: 1 m/s relative velocity along body z */
    stage2->v->data[2] = 1.0;
    /* Small attitude perturbation */
    LieVector omega0 = {.data=(double[3]){0.0, 0.01, 0.0},.size=3};
    memcpy(stage2->Omega->data, omega0.data, 3*sizeof(double));

    /* Simulate one time step */
    LieVector p_dot, v_dot, Omega_dot;
    p_dot.size=3; p_dot.data=(double[3]){0};
    v_dot.size=3; v_dot.data=(double[3]){0};
    Omega_dot.size=3; Omega_dot.data=(double[3]){0};

    stage2->thrust = 0.0;  /* coast phase */
    lie_quadrotor_rhs(stage2, &p_dot, &v_dot, &Omega_dot);

    printf("Position derivative: "); lie_vector_print(&p_dot, "");
    printf("Velocity derivative: "); lie_vector_print(&v_dot, "");
    printf("Omega derivative: "); lie_vector_print(&Omega_dot, "");

    lie_quadrotor_free(stage2);
}

/* ==============================================================
 * NASA Ingenuity Mars Helicopter
 *
 * Ingenuity (2020) is a 1.8 kg coaxial helicopter operating on Mars.
 * It uses a swashplate-less design with two counter-rotating rotors.
 *
 * Mars gravity: 3.721 m/s^2
 * Rotor speed: ~2400 RPM
 *
 * The attitude dynamics on SO(3) combined with the reduced-gravity
 * environment make geometric control essential.
 *
 * Ref: Balaram et al. (2018) "Mars Helicopter Technology
 *      Demonstrator", AIAA Atmospheric Flight Mechanics Conference
 * ============================================================== */
void lie_mars_helicopter_demo(void) {
    printf("\n=== NASA Ingenuity Mars Helicopter Dynamics ===\n");
    printf("Mass: 1.8 kg, Mars gravity: 3.721 m/s^2\n");

    LieQuadrotorState *heli = lie_quadrotor_create(1.8, 0.02, 0.02, 0.04);

    /* Hover condition on Mars */
    heli->thrust = 1.8 * 3.721; /* 6.698 N */
    /* Initial attitude: small tilt */
    LieVector omega0 = {.data=(double[3]){0.05, 0.0, 2400.0*2.0*M_PI/60.0},.size=3};
    memcpy(heli->Omega->data, omega0.data, 3*sizeof(double));

    LieVector p_dot, v_dot, Omega_dot;
    p_dot.size=3; p_dot.data=(double[3]){0};
    v_dot.size=3; v_dot.data=(double[3]){0};
    Omega_dot.size=3; Omega_dot.data=(double[3]){0};

    lie_quadrotor_rhs(heli, &p_dot, &v_dot, &Omega_dot);
    printf("Thrust: %.3f N\n", heli->thrust);
    printf("Angular acceleration: "); lie_vector_print(&Omega_dot, "");

    lie_quadrotor_free(heli);
}

/* ==============================================================
 * GPS Satellite Constellation
 *
 * GPS satellites orbit at ~20,200 km altitude (MEO).
 * Each satellite must maintain precise attitude for antenna pointing.
 *
 * Orbital period: ~12 hours
 * Gravity gradient torque is significant at MEO due to large solar panels.
 *
 * Typical GPS Block IIF satellite:
 *   Mass: ~1,630 kg
 *   Inertia: Ixx=1500, Iyy=4000, Izz=3500 kg*m^2 (approx)
 *
 * Ref: GPS Directorate (2020) "GPS Constellation Status"
 * ============================================================== */
void lie_gps_satellite_demo(void) {
    printf("\n=== GPS Satellite Attitude at MEO ===\n");
    printf("GPS Block IIF: 20,200 km altitude\n");

    LieSatelliteState *gps = lie_satellite_create(1500.0, 4000.0, 3500.0);
    gps->orbit_radius = 26561.0; /* Earth radius + 20200 km */
    gps->mu_earth = 398600.4418;
    gps->n_orbit = sqrt(gps->mu_earth / pow(gps->orbit_radius, 3.0));

    /* Orbit period */
    double T_orbit = 2.0 * M_PI / gps->n_orbit;
    printf("Orbital period: %.1f minutes\n", T_orbit / 60.0);

    /* Gravity gradient torque with nadir pointing (o3_body = [0,0,1]) */
    LieVector o3_nadir = {.data=(double[3]){0,0,1},.size=3};
    LieVector tau_gg; tau_gg.size=3; tau_gg.data=(double[3]){0};
    lie_gravity_gradient_torque(gps, &o3_nadir, &tau_gg);
    printf("Gravity gradient torque (N*m/km^2): ");
    lie_vector_print(&tau_gg, "");

    /* Free precession period */
    double Ix = gps->I->inertia->data[0];
    double Iy = gps->I->inertia->data[4];
    double Iz = gps->I->inertia->data[8];
    gps->Omega->data[2] = 0.1; /* rad/s */
    double Pi_z = Iz * gps->Omega->data[2];
    double T_precess = 2.0*M_PI * Iy / (Pi_z * (Iy - Ix));
    printf("Free precession period: %.2f s\n", fabs(T_precess));

    lie_satellite_free(gps);
}

/* ==============================================================
 * Tesla Vehicle Inertial Navigation
 *
 * A Tesla Model 3 uses IMU (Inertial Measurement Unit) for
 * dead-reckoning navigation. The attitude is updated on SO(3)
 * using gyroscope measurements.
 *
 * Typical automotive IMU: 6-axis (3 accel + 3 gyro)
 * Gyro bias stability: ~5 deg/hr
 *
 * Ref: Titterton & Weston (2004) Strapdown Inertial Navigation
 *      Technology, 2nd ed.
 * ============================================================== */
void lie_tesla_imu_demo(void) {
    printf("\n=== Tesla Vehicle IMU Attitude Propagation on SO(3) ===\n");

    /* Simulate a 90-degree turn over 3 seconds */
    LieGroupElement *R0 = lie_group_create(LIE_GROUP_SO3, "initial");
    LieVector omega; omega.size=3; omega.data=(double[3]){0};
    omega.data[2] = (M_PI/2.0) / 3.0; /* rad/s about z-axis */

    printf("Input angular velocity (rad/s): ");
    lie_vector_print(&omega, "");

    /* Compute attitude after 3 seconds: R = exp(omega * 3) */
    LieVector *omega_t = lie_vector_scale(&omega, 3.0);
    LieGroupElement *R_final = lie_exp_so3(omega_t);

    printf("Final rotation matrix (90 deg about z):\n");
    lie_matrix_print(R_final->g, "  R");

    /* Verify: R should be [[0, -1, 0], [1, 0, 0], [0, 0, 1]] */
    printf("Expected: R_z(90 deg) = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]\n");

    lie_group_free(R0); lie_vector_free(omega_t); lie_group_free(R_final);
}

/* ==============================================================
 * Main Application Demo
 * ============================================================== */

int main(void) {
    printf("============================================\n");
    printf("  Lie Group Mechanics -- Applications Demo\n");
    printf("============================================\n");

    lie_apollo_lm_attitude_demo();
    lie_spacex_separation_demo();
    lie_mars_helicopter_demo();
    lie_gps_satellite_demo();
    lie_tesla_imu_demo();

    printf("\n=== All application demos completed ===\n");
    return 0;
}

