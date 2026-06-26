#include "../include/lie_group_methods.h"
#include "../include/symplectic_geometry.h"
#include <stdio.h>
#include <time.h>

int main(void) {
    printf("Benchmark: SO(3) exponential map
");
    clock_t start = clock();
    double omega[3] = {0.1, 0.2, 0.3};
    double R[3][3];
    for (int i=0;i<100000;i++) {
        omega[0]+=0.00001;
        geo_so3_exp(omega, R);
    }
    clock_t end = clock();
    double elapsed = (double)(end-start)/CLOCKS_PER_SEC;
    printf("  100000 exp calls: 0.0000 sec (0.0 us/call)
", elapsed, elapsed*10.0);
    printf("Benchmark: Euler-Poincare dynamics
");
    start = clock();
    double I[3]={2,3,4},O[3]={1,0.5,0.2},t[3]={0,0,0},dO[3];
    for (int i=0;i<1000000;i++) geo_euler_poincare_rigid_body(I,O,t,dO);
    end = clock();
    elapsed = (double)(end-start)/CLOCKS_PER_SEC;
    printf("  1000000 dynamics calls: 0.0000 sec (0.0 ns/call)
", elapsed, elapsed*1000.0);
    return 0;
}
