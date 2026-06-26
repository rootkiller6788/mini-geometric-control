import os

base = r"F:\nano-everything\mini-complex-control-theory\21. mini-geometric-control\mini-feedback-linearization-geo"
os.chdir(base)

src_fl = r"""#include "feedback_linearization.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

#define DEFAULT_TOL 1e-8

RelativeDegreeResult compute_relative_degree(const NonlinearSystem *sys,
                                              const Vector *x,
                                              int max_order) {
    RelativeDegreeResult result;
    memset(&result, 0, sizeof(result));
    result.rho = -1;
    result.is_well_defined = false;
    if (!sys || !x || max_order <= 0) return result;
    if (sys->m != 1 || sys->p != 1) return result;
    VectorField *g = sys->g[0];
    VectorField *f = sys->f;
    ScalarField *h = sys->h[0];
    if (!f || !g || !h) return result;
    for (int k = 0; k < max_order; k++) {
        double LgLfk_h = lie_derivative_mixed(g, f, h, x, k);
        if (fabs(LgLfk_h) > DEFAULT_TOL) {
            result.rho = k + 1;
            result.is_well_defined = true;
            result.Lg_Lf_power = LgLfk_h;
            break;
        }
    }
    if (result.is_well_defined && result.rho > 0) {
        result.Lf_powers_h = compute_all_lie_derivatives(f, h, x, result.rho);
    }
    return result;
}
"""
with open(os.path.join(base, "src", "feedback_linearization.c"), "w") as f:
    f.write(src_fl)
print("feedback_linearization.c part 1 written")
