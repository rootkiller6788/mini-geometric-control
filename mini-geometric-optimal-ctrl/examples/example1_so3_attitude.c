#include "../include/geo_optctrl_core.h"
#include "../include/lie_group_methods.h"
#include "../include/symplectic_geometry.h"
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== Example 1: SO(3) Attitude Maneuver ===

");
    double omega[3] = {0.0, 0.0, M_PI/2.0};
    double R[3][3];
    geo_so3_exp(omega, R);
    printf("Rotation by 90 deg about Z axis:
");
    for (int i=0;i<3;i++) {
        printf("  [  0.0000   0.0000   0.0000]
", R[i][0], R[i][1], R[i][2]);
    }
    double R_recover[3][3];
    double omega2[3] = {0, 0, -M_PI/2.0};
    geo_so3_exp(omega2, R_recover);
    printf("
Verify: R(-90) * R(90) should be identity:
");
    double result[3][3] = {{0}};
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++)
        result[i][j] += R_recover[i][k] * R[k][j];
    for (int i=0;i<3;i++)
        printf("  [  0.0000   0.0000   0.0000]
", result[i][0], result[i][1], result[i][2]);
    printf("
Euler-Poincare rigid body dynamics:
");
    double I[3] = {2.0, 3.0, 4.0};
    double Omega[3] = {1.0, 0.5, 0.2};
    double tau[3] = {0.0, 0.0, 0.0};
    double dOmega[3];
    geo_euler_poincare_rigid_body(I, Omega, tau, dOmega);
    printf("  I=(0.0,0.0,0.0) Omega=(0.0,0.0,0.0)
", I[0],I[1],I[2],Omega[0],Omega[1],Omega[2]);
    printf("  dOmega/dt = (0.0000, 0.0000, 0.0000)
", dOmega[0], dOmega[1], dOmega[2]);
    printf("
Example complete.
");
    return 0;
}
