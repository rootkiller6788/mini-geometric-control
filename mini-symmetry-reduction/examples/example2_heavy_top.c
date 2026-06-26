#include "sr_core.h"
#include "sr_reduction.h"
#include "sr_poisson.h"
#include "sr_dynamics.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    printf("============================================================\n");
    printf("  Example 2: Heavy Top — Semidirect Product Reduction\n");
    printf("  S^1 symmetry + SO(3) ⋉ R^3 semidirect product\n");
    printf("============================================================\n\n");

    /* Heavy top coordinates on se(3)*: (Pi, Gamma) */
    double Pi[3] = {0.5, 0.3, 1.0};
    double Gamma[3] = {0.2, 0.1, 0.96};

    /* Normalize Gamma to unit (gravity direction) */
    double nrm = sqrt(Gamma[0]*Gamma[0] + Gamma[1]*Gamma[1] + Gamma[2]*Gamma[2]);
    Gamma[0] /= nrm; Gamma[1] /= nrm; Gamma[2] /= nrm;

    printf("Angular momentum: Pi = (%.3f, %.3f, %.3f)\n", Pi[0], Pi[1], Pi[2]);
    printf("Gravity direction in body frame: Gamma = (%.3f, %.3f, %.3f)\n",
           Gamma[0], Gamma[1], Gamma[2]);

    /* Casimirs */
    double c1 = sr_heavy_top_casimir(Pi, Gamma);
    double c2 = sr_heavy_top_casimir2(Pi, Gamma);

    printf("\nCasimir functions:\n");
    printf("  C1 = Pi · Gamma = %.4f (vertical angular momentum)\n", c1);
    printf("  C2 = ||Gamma||^2 = %.4f (gravity vector magnitude)\n", c2);

    /* Heavy top Lie-Poisson bracket test */
    printf("\nLie-Poisson bracket structure on se(3)*:\n");
    printf("  {Pi_i, Pi_j} = -epsilon_{ijk} Pi_k\n");
    printf("  {Pi_i, Gamma_j} = -epsilon_{ijk} Gamma_k\n");
    printf("  {Gamma_i, Gamma_j} = 0\n");

    /* Test bracket computation */
    double nabla_Pi_F[3] = {1,0,0};
    double nabla_Gamma_F[3] = {0,0,0};
    double nabla_Pi_H[3] = {0,1,0};
    double nabla_Gamma_H[3] = {0,0,0};
    double bracket_val;
    sr_heavy_top_bracket(Pi, Gamma, nabla_Pi_F, nabla_Gamma_F,
                          nabla_Pi_H, nabla_Gamma_H, &bracket_val);
    printf("\n  {Pi_1, Pi_2} = %.4f (expected -Pi_3 = %.4f)\n", bracket_val, -Pi[2]);

    /* Energy-momentum analysis */
    printf("\n--- Semidirect Product Reduction ---\n");
    printf("Heavy top has S^1 symmetry (rotation about gravity direction).\n");
    printf("Reduced phase space: se(3)* coadjoint orbit.\n");
    printf("Semidirect structure: SE(3) = SO(3) ⋉ R^3\n");
    printf("Reduction produces 4-dimensional symplectic leaves (2 Casimirs).\n");

    printf("\n--- Physical Interpretation ---\n");
    printf("Pi · Gamma = const: angular momentum component along gravity axis\n");
    printf("||Gamma|| = 1: unit vector along gravity in body frame\n");
    printf("Hamiltonian: H = 1/2 Pi·I^{-1}·Pi + mgl Gamma·chi\n");
    printf("where chi is the center of mass direction in body frame.\n");

    printf("\n============================================================\n");
    printf("  Example 2 Complete\n");
    printf("============================================================\n");
    return 0;
}
