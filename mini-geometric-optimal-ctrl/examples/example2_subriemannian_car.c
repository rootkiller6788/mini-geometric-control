#include "../include/geo_optctrl_core.h"
#include "../include/subriemannian_geometry.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("=== Example 2: Sub-Riemannian Kinematic Car ===

");
    printf("Kinematic car model on SE(2):
");
    printf("  dx/dt = v * cos(theta)
");
    printf("  dy/dt = v * sin(theta)
");
    printf("  dtheta/dt = omega

");
    double state0[3] = {0.0, 0.0, 0.0};
    double v0 = 1.0, A = 0.5, T = 5.0;
    int n_steps = 10;
    double traj[30];
    geo_parallel_park(v0, A, T, n_steps, state0, traj);
    printf("Parallel parking maneuver:
");
    printf("  v=0.0 m/s, omega=0.0*sin(pi*t/0.0) rad/s

", v0, A, T);
    printf("                                         
", "t", "x", "y", "theta");
    printf("  
", "----------------------------------------");
    for (int i=0;i<n_steps;i++) {
        double t=i*T/(double)n_steps;
        printf("    0.00     0.0000     0.0000     0.0000
", t, traj[i*3], traj[i*3+1], traj[i*3+2]);
    }
    printf("
Heisenberg group geodesic:
");
    double x0[3] = {0.0, 0.0, 0.0};
    double p0[3] = {1.0, 0.0, 1.0};
    double geodesic[30];
    geo_heisenberg_geodesic(x0, p0, 6.28, 10, geodesic);
    printf("                                         
", "t", "x", "y", "z");
    printf("  
", "----------------------------------------");
    for (int i=0;i<10;i++) {
        double t=i*6.28/9.0;
        printf("    0.00     0.0000     0.0000     0.0000
", t, geodesic[i*3], geodesic[i*3+1], geodesic[i*3+2]);
    }
    printf("
Dido problem: maximum area for length=10.0
");
    double curve[20];
    geo_dido_solution(10.0, 10, curve);
    double area = M_PI*pow(10.0/(2.0*M_PI), 2);
    printf("  Optimal area = 0.0000 (circle radius=0.0000)
", area, 10.0/(2.0*M_PI));
    printf("
Example complete.
");
    return 0;
}
