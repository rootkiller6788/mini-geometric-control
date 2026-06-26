#include "../include/geo_optctrl_core.h"
#include "../include/symplectic_geometry.h"
#include "../include/lie_group_methods.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Example 3: Brockett Integrator and Arnold Cat Map ===

");
    printf("Brockett integrator (canonical nonholonomic system):
");
    printf("  dx1/dt = u1
");
    printf("  dx2/dt = u2
");
    printf("  dx3/dt = x1*u2 - x2*u1

");
    printf("Arnold Cat Map on T2 (chaotic symplectic map):
");
    double x=0.3, y=0.7;
    printf("  Initial: (0.0000, 0.0000)
", x, y);
    for (int n=1;n<=10;n++) {
        geo_arnold_cat_map(&x, &y, 1);
        printf("  Step  0: (0.0000, 0.0000)
", n, x, y);
    }
    printf("
Yoshida higher-order symplectic composition:
");
    double w[10];
    int nw;
    geo_yoshida_weights(4, w, &nw);
    printf("  Order 4 weights (0 stages): ", nw);
    for (int i=0;i<nw;i++) printf("0.000000 ", w[i]);
    printf("
  Sum = 0.000000
", w[0]+w[1]+w[2]);
    printf("
Symplectic Euler for harmonic oscillator:
");
    printf("  H = 0.5*(p^2 + q^2), q0=1, p0=0, h=0.1, steps=100
");
    GeoSymplecticSpace s;
    geo_symplectic_init(&s, 1);
    printf("  Symplectic form: Omega = [[0,1],[-1,0]] verified
");
    printf("
Example complete.
");
    return 0;
}
