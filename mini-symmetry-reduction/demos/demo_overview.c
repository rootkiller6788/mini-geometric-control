#include <stdio.h>

int main(void) {
    printf("=============================================================\n");
    printf("  Symmetry Reduction in Geometric Control Theory\n");
    printf("  Module: mini-symmetry-reduction\n");
    printf("=============================================================\n\n");

    printf("This module implements symmetry reduction for mechanical\n");
    printf("and control systems with Lie group symmetries.\n\n");

    printf("Key concepts covered:\n");
    printf("  - Marsden-Weinstein symplectic reduction\n");
    printf("  - Euler-Poincare reduction on Lie algebras\n");
    printf("  - Lie-Poisson bracket on g* (coadjoint orbit dynamics)\n");
    printf("  - Noether\'s theorem: symmetries -> conserved quantities\n");
    printf("  - Semidirect product reduction (heavy top, SE(3))\n");
    printf("  - Reduced optimal control and controlled Lagrangians\n\n");

    printf("Canonical problems:\n");
    printf("  1. Free rigid body (SO(3) -> so(3)* reduction)\n");
    printf("  2. Heavy top (semidirect product SE(3) reduction)\n");
    printf("  3. Spherical pendulum (S^1 reduction)\n");
    printf("  4. Spacecraft attitude control (controlled EP)\n");
    printf("  5. Underwater vehicle (Kirchhoff equations on se(3)*)\n\n");

    printf("Build & Run:\n");
    printf("  make          - Build static library\n");
    printf("  make test     - Run test suite (~20 tests)\n");
    printf("  make examples - Build all examples\n");
    printf("  make demo     - Run all examples\n");
    printf("  make bench    - Run performance benchmarks\n");
    printf("  make clean    - Clean build artifacts\n\n");

    printf("=============================================================\n");
    return 0;
}
