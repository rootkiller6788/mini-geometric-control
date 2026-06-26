/* ==============================================================
 * example_quadrotor.c -- Quadrotor Control on SE(3)
 *
 * Demonstrates geometric tracking control of a quadrotor UAV
 * on the SE(3) configuration manifold.
 *
 * The quadrotor tracks a circular trajectory in the x-y plane
 * while maintaining level attitude.
 *
 * Ref: Lee, Leok, McClamroch (2010) "Geometric Tracking Control
 *      of a Quadrotor UAV on SE(3)", CDC 2010
 * ============================================================== */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lie_core.h"
#include "lie_mechanics.h"

int main(void) {
    printf("=== Quadrotor Geometric Tracking on SE(3) ===\n\n");

    /* Create a 0.5 kg quadrotor */
    LieQuadrotorState *qr = lie_quadrotor_create(0.5, 0.0023, 0.0023, 0.0046);
    printf("Mass: 0.5 kg\n");
    printf("Inertia: Ixx=0.0023, Iyy=0.0023, Izz=0.0046 kg*m^2\n");

    /* Initial state: at origin, hovering */
    qr->thrust = 0.5 * 9.81;  /* hover thrust */

    /* Desired trajectory: circle in x-y plane at 1m height */
    double radius = 2.0;
    double omega_des = 1.0;  /* rad/s */

    double h = 0.01;
    int steps = 300;
    printf("Simulating %d steps with dt=%.3f\n", steps, h);
    printf("Desired: circle radius=%.1fm, omega=%.1f rad/s\n\n", radius, omega_des);

    /* Controller gains */
    double kp = 4.0, kd = 3.0, kR = 8.0, kOmega = 4.0;

    for (int k = 0; k < steps; k++) {
        double t = k * h;

        /* Desired state at time t */
        LieVector p_d, v_d, a_d;
        p_d.size = 3; v_d.size = 3; a_d.size = 3;
        double pd[3] = {radius*cos(omega_des*t), radius*sin(omega_des*t), 1.0};
        double vd[3] = {-radius*omega_des*sin(omega_des*t),
                         radius*omega_des*cos(omega_des*t), 0.0};
        double ad[3] = {-radius*omega_des*omega_des*cos(omega_des*t),
                        -radius*omega_des*omega_des*sin(omega_des*t), 0.0};
        p_d.data = pd; v_d.data = vd; a_d.data = ad;

        /* Desired attitude: body z aligned with thrust direction */
        LieMatrix *R_d = lie_matrix_identity(3);

        /* Compute geometric controller */
        LieVector tau_out;
        tau_out.size = 3; double td[3] = {0,0,0}; tau_out.data = td;
        double thrust_cmd;

        lie_quadrotor_geometric_controller(qr, &p_d, &v_d, &a_d, R_d,
            &(LieVector){.data=(double[3]){0,0,0},.size=3},
            &thrust_cmd, &tau_out, kp, kd, kR, kOmega);

        qr->thrust = thrust_cmd;
        memcpy(qr->tau->data, tau_out.data, 3*sizeof(double));

        /* Simulate one step (Euler on vector part + exp map on rotation) */
        LieVector p_dot, v_dot, Omega_dot;
        p_dot.size=3; p_dot.data=(double[3]){0};
        v_dot.size=3; v_dot.data=(double[3]){0};
        Omega_dot.size=3; Omega_dot.data=(double[3]){0};

        lie_quadrotor_rhs(qr, &p_dot, &v_dot, &Omega_dot);

        /* Update translation */
        for (int i = 0; i < 3; i++) qr->p->data[i] += h * v_dot.data[i];
        for (int i = 0; i < 3; i++) qr->v->data[i] += h * v_dot.data[i];

        /* Update rotation: R <- R * exp(h*Omega) */
        for (int i = 0; i < 3; i++) qr->Omega->data[i] += h * Omega_dot.data[i];
        LieVector *hOmega = lie_vector_scale(qr->Omega, h);
        LieGroupElement *exp_hO = lie_exp_so3(hOmega);
        LieMatrix *Rnew = lie_matrix_mul(qr->R, exp_hO->g);
        memcpy(qr->R->data, Rnew->data, 9*sizeof(double));
        lie_vector_free(hOmega); lie_group_free(exp_hO); lie_matrix_free(Rnew);

        lie_matrix_free(R_d);

        if (k % 50 == 0) {
            printf("  t=%.2f: pos=(%6.2f,%6.2f,%6.2f) "
                   "thrust=%.3fN\n",
                   t, qr->p->data[0], qr->p->data[1], qr->p->data[2],
                   qr->thrust);
        }
    }

    printf("\nFinal position: (%.3f, %.3f, %.3f)\n",
           qr->p->data[0], qr->p->data[1], qr->p->data[2]);
    printf("Final velocity: (%.3f, %.3f, %.3f)\n",
           qr->v->data[0], qr->v->data[1], qr->v->data[2]);

    lie_quadrotor_free(qr);
    printf("\n=== Example completed ===\n");
    return 0;
}

