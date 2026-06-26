#include "../include/geo_optctrl_core.h"
#include "../include/lie_group_methods.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("===============================
");
    printf(" Geometric Optimal Control Demo
");
    printf("===============================

");
    printf("1. SO(3) Exponential Map (Rodrigues):
");
    double omegas[4][3] = {{1,0,0},{0,1,0},{0,0,1},{0.5,0.3,0.1}};
    for (int k=0;k<4;k++) {
        double R[3][3];
        geo_so3_exp(omegas[k], R);
        printf("   omega=(0.0,0.0,0.0):
", omegas[k][0], omegas[k][1], omegas[k][2]);
        for (int i=0;i<3;i++)
            printf("     [ 0.0000  0.0000  0.0000]
", R[i][0],R[i][1],R[i][2]);
    }
    printf("
2. Euler-Poincare Rigid Body:
");
    printf("   I=diag(2,3,4), free rotation with Omega=(1,0.5,0.2)
");
    double I[3]={2,3,4}, O[3]={1,0.5,0.2}, tau[3]={0,0,0}, dO[3];
    geo_euler_poincare_rigid_body(I,O,tau,dO);
    printf("   dOmega/dt = (0.0000, 0.0000, 0.0000)
", dO[0],dO[1],dO[2]);
    printf("
3. Arnold Cat Map (chaos on torus):
");
    double xi=0.2,yi=0.3;
    for (int n=0;n<8;n++) {
        geo_arnold_cat_map(&xi,&yi,1);
        printf("   n=0: (0.000000, 0.000000)
", n+1, xi, yi);
    }
    printf("
Demo complete.
");
    return 0;
}
