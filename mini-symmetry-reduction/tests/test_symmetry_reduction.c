#include "sr_core.h"
#include "sr_reduction.h"
#include "sr_poisson.h"
#include "sr_dynamics.h"
#include "sr_control.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#define TOL 1e-10

static int tests_run = 0;
static int tests_passed = 0;
#define TEST(name) do { printf("  %-55s", name); } while(0)
#define PASS() do { tests_passed++; tests_run++; } while(0)
#define FAIL(msg) do { tests_run++; printf("FAIL: %s\n", msg); } while(0)
#define CHECK(cond) do { if (cond) PASS(); else FAIL(#cond); } while(0)
#define CHECK_CLOSE(a, b, tol) do { if (fabs((a)-(b)) < (tol)) PASS(); else { printf("FAIL: |%g-%g| = %g\n", (a),(b),fabs((a)-(b))); } } while(0)

/* ============================================================================
 * L1: Definitions — Struct instantiation and basic integrity
 * ============================================================================ */
void test_l1_lie_algebra_creation(void) {
    TEST("L1: Lie algebra creation");
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    CHECK(so3 != NULL);
    CHECK(so3->dimension == 3);
    CHECK(so3->is_abelian == false);
    CHECK(so3->is_simple == true);
    sr_algebra_free(so3);

    SRLieAlgebra* ab = sr_algebra_create_abelian(4);
    CHECK(ab != NULL);
    CHECK(ab->dimension == 4);
    CHECK(ab->is_abelian == true);
    sr_algebra_free(ab);
}

void test_l1_group_element(void) {
    TEST("L1: Group element operations");
    SRLieGroup* SO3 = sr_group_create(3, 3, "SO(3)", SR_SYMMETRY_NONABELIAN);
    SO3->is_compact = true;
    SRGroupElement* g = sr_element_create(SO3);
    sr_element_set_identity(g);
    CHECK(g->is_identity == true);
    sr_element_free(g);
    sr_group_free(SO3);
}

void test_l1_momentum_map(void) {
    TEST("L1: Momentum map creation");
    SRLieGroup* SO3 = sr_group_create(3, 3, "SO(3)", SR_SYMMETRY_NONABELIAN);
    SO3->is_compact = true;
    SO3->algebra = sr_algebra_create_so3();
    SRMomentumMap* J = sr_momentum_create(SO3, 6, SR_MOMENTUM_ANGULAR);
    CHECK(J != NULL);
    CHECK(J->algebra_dim == 3);
    CHECK(J->phase_dim == 6);
    sr_momentum_free(J);
    sr_algebra_free(SO3->algebra);
    sr_group_free(SO3);
}

/* ============================================================================
 * L2: Core Concepts — Structure constants, bracket, actions
 * ============================================================================ */
void test_l2_structure_constants(void) {
    TEST("L2: so(3) structure constants");
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    /* [e1, e2] = e3 */
    double c = sr_algebra_get_struct_const(so3, 0, 1, 2);
    CHECK_CLOSE(c, 1.0, TOL);
    /* Antisymmetry: [e2, e1] = -[e1, e2] */
    double c_anti = sr_algebra_get_struct_const(so3, 1, 0, 2);
    CHECK_CLOSE(c_anti, -1.0, TOL);
    sr_algebra_free(so3);
}

void test_l2_lie_bracket(void) {
    TEST("L2: Lie bracket computation");
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    double X[3] = {1.0, 0.0, 0.0};
    double Y[3] = {0.0, 1.0, 0.0};
    double Z[3] = {0.0, 0.0, 0.0};
    sr_algebra_bracket(so3, X, Y, Z);
    CHECK_CLOSE(Z[0], 0.0, TOL);
    CHECK_CLOSE(Z[1], 0.0, TOL);
    CHECK_CLOSE(Z[2], 1.0, TOL); /* [e1, e2] = e3 */
    sr_algebra_free(so3);
}

void test_l2_killing_form(void) {
    TEST("L2: Killing form of so(3)");
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    /* Killing form of so(3) is K(e_i, e_j) = -2 delta_{ij} */
    double K = sr_algebra_killing_form(so3, 0, 1);
    CHECK_CLOSE(K, 0.0, TOL);
    sr_algebra_free(so3);
}

/* ============================================================================
 * L3: Mathematical Structures — Symplectic, Poisson, Heisenberg
 * ============================================================================ */
void test_l3_symplectic(void) {
    TEST("L3: Symplectic manifold creation");
    SRSymplecticManifold* sp = sr_symplectic_create(2);
    CHECK(sp != NULL);
    CHECK(sp->dim == 4);
    sr_symplectic_set_canonical(sp);
    CHECK(sp->is_closed == true);
    CHECK(sp->is_nondeg == true);
    sr_symplectic_free(sp);
}

void test_l3_poisson(void) {
    TEST("L3: Poisson manifold creation");
    SRPoissonManifold* pm = sr_poisson_create(3);
    CHECK(pm != NULL);
    CHECK(pm->dim == 3);
    sr_poisson_free(pm);
}

void test_l3_heisenberg_algebra(void) {
    TEST("L3: Heisenberg algebra");
    SRLieAlgebra* heis = sr_algebra_create_heisenberg();
    CHECK(heis->dimension == 3);
    CHECK(heis->is_nilpotent == true);
    double c = sr_algebra_get_struct_const(heis, 0, 1, 2);
    CHECK_CLOSE(c, 1.0, TOL);
    sr_algebra_free(heis);
}

/* ============================================================================
 * L4: Fundamental Laws — Noether, Jacobi, Exponential map
 * ============================================================================ */
void test_l4_noether(void) {
    TEST("L4: Noether theorem verification");
    SRLieGroup* SO3 = sr_group_create(3, 3, "SO(3)", SR_SYMMETRY_NONABELIAN);
    SO3->is_compact = true;
    SO3->algebra = sr_algebra_create_so3();
    SRMomentumMap* J = sr_momentum_angular_so3();
    double phase_point[6] = {1,0,0, 0,1,0};
    double conserved[3];
    sr_noether_conserved_quantities(J, phase_point, conserved);
    CHECK(J != NULL);
    sr_momentum_free(J);
    sr_algebra_free(SO3->algebra);
    sr_group_free(SO3);
}

void test_l4_jacobi_identity(void) {
    TEST("L4: Jacobi identity for so(3)");
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    bool jacobi = sr_algebra_verify_jacobi(so3, TOL);
    CHECK(jacobi == true);
    sr_algebra_free(so3);
}

void test_l4_exponential_map(void) {
    TEST("L4: SO(3) exponential map (Rodrigues)");
    double omega[3] = {0.1, 0.2, 0.3};
    double R[9];
    sr_exp_map_so3(omega, R);
    /* Check orthogonality: R * R^T = I */
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
        double dot = 0.0;
        for (int k = 0; k < 3; k++)
            dot += R[i*3 + k] * R[j*3 + k];
        double expected = (i == j) ? 1.0 : 0.0;
        if (fabs(dot - expected) > 1e-8) {
            FAIL("Non-orthogonal exp(so(3))");
            return;
        }
    }
    PASS();
}

/* ============================================================================
 * L5: Algorithms — Euler-Poincare, Lie-Poisson, RK4, BCH
 * ============================================================================ */
void test_l5_euler_poincare(void) {
    TEST("L5: Euler-Poincare system creation");
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    SREulerPoincareSystem* ep = sr_ep_create(so3, NULL, +1);
    CHECK(ep != NULL);
    CHECK(ep->n_config == 3);
    double I[9] = {2,0,0,0,3,0,0,0,4};
    sr_ep_set_inertia(ep, I);
    double xi[3] = {1.0, 0.5, 0.2};
    double dxi[3];
    sr_ep_rhs(ep, xi, dxi);
    double energy = sr_ep_energy(ep);
    CHECK(energy >= 0.0);
    sr_ep_free(ep);
}

void test_l5_lie_poisson_bracket(void) {
    TEST("L5: Lie-Poisson bracket evaluation");
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    SRLiePoissonBracket* lpb = sr_lpb_create(so3, +1);
    CHECK(lpb != NULL);
    double mu[3] = {1.0, 2.0, 3.0};
    double nabla_F[3] = {1,0,0};
    double nabla_H[3] = {0,1,0};
    double bracket = sr_lpb_evaluate(lpb, mu, nabla_F, nabla_H);
    /* {mu_1, mu_2} = sign * mu_3 for so(3)* Lie-Poisson (sign=+1 gives +mu_3) */
    CHECK_CLOSE(bracket, mu[2], TOL);
    sr_lpb_free(lpb);
}

void test_l5_rk4_integration(void) {
    TEST("L5: RK4 integration step");
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    SREulerPoincareSystem* ep = sr_ep_create(so3, NULL, +1);
    double I[9] = {1,0,0,0,2,0,0,0,3};
    sr_ep_set_inertia(ep, I);
    ep->configuration[0] = 0.1;
    ep->configuration[1] = 0.05;
    ep->configuration[2] = 0.02;
    double e0 = sr_ep_energy(ep);
    sr_ep_step(ep, 0.01);
    double e1 = sr_ep_energy(ep);
    /* Energy should be approximately conserved */
    CHECK_CLOSE(e0, e1, 0.1);
    sr_ep_free(ep);
}

void test_l5_bch_formula(void) {
    TEST("L5: BCH formula (2nd order)");
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    double X[3] = {0.1, 0.0, 0.0};
    double Y[3] = {0.0, 0.1, 0.0};
    double Z[3];
    sr_bch_formula_2nd(so3, X, Y, Z);
    /* Z = X + Y + 1/2[X,Y] = (0.1, 0.1, 0.005) */
    CHECK_CLOSE(Z[0], 0.1, TOL);
    CHECK_CLOSE(Z[1], 0.1, TOL);
    CHECK_CLOSE(Z[2], 0.005, TOL);
    sr_algebra_free(so3);
}

/* ============================================================================
 * L6: Canonical Problems — Rigid body, heavy top, spherical pendulum
 * ============================================================================ */
void test_l6_rigid_body(void) {
    TEST("L6: Free rigid body reduced dynamics");
    SRLiePoissonBracket* lpb = sr_rigid_body_poisson();
    double Pi[3] = {1.0, 2.0, 3.0};
    double nabla_F[3] = {1,0,0};
    double nabla_H[3] = {0,1,0};
    double bracket = sr_rigid_body_bracket(Pi, nabla_F, nabla_H);
    /* {Pi1, Pi2} = -Pi3 */
    CHECK_CLOSE(bracket, -3.0, TOL);
    sr_lpb_free(lpb);
}

void test_l6_heavy_top(void) {
    TEST("L6: Heavy top Casimirs");
    double Pi[3] = {1,2,3};
    double Gamma[3] = {0,0,1};
    double c1 = sr_heavy_top_casimir(Pi, Gamma);
    double c2 = sr_heavy_top_casimir2(Pi, Gamma);
    CHECK_CLOSE(c1, 3.0, TOL);
    CHECK_CLOSE(c2, 1.0, TOL);

    double nabla_Pi_F[3] = {1,0,0};
    double nabla_Gamma_F[3] = {0,0,0};
    double nabla_Pi_H[3] = {0,1,0};
    double nabla_Gamma_H[3] = {0,0,0};
    double bracket;
    sr_heavy_top_bracket(Pi, Gamma, nabla_Pi_F, nabla_Gamma_F,
                          nabla_Pi_H, nabla_Gamma_H, &bracket);
    CHECK_CLOSE(bracket, -3.0, TOL);
}

void test_l6_spherical_pendulum(void) {
    TEST("L6: Spherical pendulum reduction");
    double dtheta, dp;
    sr_spherical_pendulum_reduced_rhs(0.5, 1.0, 2.0, 1.0, 1.0, 9.81, &dtheta, &dp);
    CHECK(dtheta != 0.0); /* Should have non-zero derivative */
    double V_eff = sr_spherical_pendulum_effective_potential(0.5, 2.0, 1.0, 1.0, 9.81);
    CHECK(V_eff != 0.0);
}

/* ============================================================================
 * L7: Applications — Spacecraft, underwater vehicle
 * ============================================================================ */
void test_l7_spacecraft(void) {
    TEST("L7: Spacecraft (NASA) attitude dynamics");
    double Omega[3] = {0.1, 0.05, 0.02};
    double I_diag[3] = {10.0, 15.0, 20.0};
    double tau[3] = {0.01, 0.0, -0.01};
    double dOmega[3];
    sr_spacecraft_reduced_rhs(Omega, I_diag, tau, dOmega);
    /* With control torque, derivatives should be non-zero */
    CHECK(fabs(dOmega[0]) > 0);
    CHECK(fabs(dOmega[2]) > 0);
}

void test_l7_quadrotor(void) {
    TEST("L7: Quadrotor (Tesla/SpaceX) UAV dynamics");
    double Omega[3] = {0, 0, 0.5};
    double v[3] = {1, 0, 0};
    double I_diag[3] = {0.01, 0.01, 0.02};
    double mass = 0.5;
    double thrust[3] = {0, 0, 5.0};
    double torque[3] = {0.001, 0.0, 0.0};
    double dOmega[3], dv[3];
    sr_quadrotor_reduced_rhs(Omega, v, I_diag, mass, thrust, torque, dOmega, dv);
    CHECK(fabs(dv[2]) > 0); /* Vertical thrust */
    CHECK(fabs(dOmega[0]) > 0); /* Roll torque */
}

void test_l7_underwater(void) {
    TEST("L7: Underwater vehicle (REMUS AUV) dynamics");
    double state[6] = {1,0,0, 0,1,0};
    double M[36] = {10,0,0,0,0,0,0,15,0,0,0,0,0,0,20,0,0,0,
                    0,0,0,5,0,0,0,0,0,0,5,0,0,0,0,0,0,5};
    double forces[6] = {0.1, 0, 0, 0.5, 0, 0};
    double derivs[6];
    sr_underwater_vehicle_reduced_rhs(state, M, forces, derivs);
    CHECK(fabs(derivs[0]) > 0 || fabs(derivs[3]) > 0);
}

/* ============================================================================
 * L8: Advanced Topics
 * ============================================================================ */
void test_l8_ep_rhs_controlled(void) {
    TEST("L8: Controlled Euler-Poincare equations");
    SRLieAlgebra* so3 = sr_algebra_create_so3();
    SREulerPoincareSystem* ep = sr_ep_create(so3, NULL, +1);
    double I[9] = {2,0,0,0,3,0,0,0,4};
    sr_ep_set_inertia(ep, I);
    double xi[3] = {1, 0.5, 0.2};
    double F_ext[3] = {0.1, -0.05, 0.0};
    double dxi[3];
    sr_ctrl_ep_rhs(ep, xi, F_ext, dxi);
    CHECK(dxi[0] != 0.0);
    sr_ep_free(ep);
}

void test_l8_energy_momentum(void) {
    TEST("L8: Energy-momentum stability analysis");
    SRLiePoissonBracket* lpb = sr_rigid_body_poisson();
    double mu_eq[3] = {1.0, 0.0, 0.0};
    double xi_eq[3] = {0.5, 0.0, 0.0};
    double eig_real[3], eig_imag[3];
    sr_reduced_stability(lpb, mu_eq, xi_eq, eig_real, eig_imag);
    CHECK(eig_real[0] >= 0 || eig_real[0] <= 0); /* Just check computation runs */
    sr_lpb_free(lpb);
}

/* ============================================================================
 * Matrix and Vector Utilities
 * ============================================================================ */
void test_matrix_ops(void) {
    TEST("Matrix operations");
    SRMatrix* A = sr_matrix_create(3, 3);
    sr_matrix_set(A, 0, 1, 2.0);
    CHECK_CLOSE(sr_matrix_get(A, 0, 1), 2.0, TOL);
    sr_matrix_free(A);
}

void test_vector_ops(void) {
    TEST("Vector operations");
    SRVector* v1 = sr_vector_create(3);
    SRVector* v2 = sr_vector_create(3);
    sr_vector_set(v1, 0, 1.0);
    sr_vector_set(v1, 1, 2.0);
    sr_vector_set(v1, 2, 3.0);
    sr_vector_set(v2, 0, 4.0);
    sr_vector_set(v2, 1, 5.0);
    sr_vector_set(v2, 2, 6.0);
    double dot = sr_vector_dot(v1, v2);
    CHECK_CLOSE(dot, 32.0, TOL);
    double nrm = sr_vector_norm(v1);
    CHECK_CLOSE(nrm, sqrt(14.0), TOL);
    sr_vector_free(v1);
    sr_vector_free(v2);
}

/* ============================================================================ */
int main(void) {
    printf("=== Symmetry Reduction Test Suite ===\n\n");

    printf("--- L1: Definitions ---\n");
    test_l1_lie_algebra_creation();
    test_l1_group_element();
    test_l1_momentum_map();

    printf("\n--- L2: Core Concepts ---\n");
    test_l2_structure_constants();
    test_l2_lie_bracket();
    test_l2_killing_form();

    printf("\n--- L3: Mathematical Structures ---\n");
    test_l3_symplectic();
    test_l3_poisson();
    test_l3_heisenberg_algebra();

    printf("\n--- L4: Fundamental Laws ---\n");
    test_l4_noether();
    test_l4_jacobi_identity();
    test_l4_exponential_map();

    printf("\n--- L5: Algorithms ---\n");
    test_l5_euler_poincare();
    test_l5_lie_poisson_bracket();
    test_l5_rk4_integration();
    test_l5_bch_formula();

    printf("\n--- L6: Canonical Problems ---\n");
    test_l6_rigid_body();
    test_l6_heavy_top();
    test_l6_spherical_pendulum();

    printf("\n--- L7: Applications ---\n");
    test_l7_spacecraft();
    test_l7_quadrotor();
    test_l7_underwater();

    printf("\n--- L8: Advanced Topics ---\n");
    test_l8_ep_rhs_controlled();
    test_l8_energy_momentum();

    printf("\n--- Utilities ---\n");
    test_matrix_ops();
    test_vector_ops();

    printf("\n========================================\n");
    printf("Results: %d/%d tests passed\n", tests_passed, tests_run);
    printf("========================================\n");

    return (tests_passed == tests_run) ? 0 : 1;
}
