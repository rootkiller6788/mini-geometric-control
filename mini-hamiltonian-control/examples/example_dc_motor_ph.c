/* example_dc_motor_ph.c -- DC motor as port-Hamiltonian system */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hamiltonian_control.h"
#include "port_hamiltonian.h"

int main(void) {
    printf("=== DC Motor: Port-Hamiltonian Model ===\n");

    dc_motor_params_t params = {0.01, 1.0, 0.001, 0.0001, 0.01};
    port_hamiltonian_system_t *motor = dc_motor_create(&params);

    double x[2] = {0.0, 0.0}; /* phi=0, p=0 */
    double u[2] = {12.0, 0.0}; /* V=12V, tau_load=0 */
    double y[2];
    double dt = 0.001;

    printf("Parameters: L=%.3fH, R=%.1fOhm, J=%.4f, b=%.4f, K=%.3f\n",
           params.L, params.R, params.J_m, params.b, params.K);
    printf("Input: V=12V, tau_load=0\n\n");

    for (int k = 0; k <= 100; k++) {
        if (k % 20 == 0) {
            double i = y[0]; /* current */
            double w = y[1]; /* angular velocity */
            printf("t=%.3f: i=%.4f A, omega=%.2f rad/s, phi=%.4f Wb, p=%.4f kg.m^2/s\n",
                   k*dt, i, w, x[0], x[1]);
        }
        ph_symplectic_step(motor, u, x, dt, y);
    }

    port_hamiltonian_system_free(motor);
    printf("\nSteady-state: i ~ V/R = %.1f A, omega ~ V/K = %.0f rad/s\n",
           12.0/params.R, 12.0/params.K);
    return 0;
}