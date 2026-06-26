/* ==============================================================
 * test_lie.c -- Test Suite for mini-lie-group-mechanics
 *
 * Tests all core APIs: matrix ops, vector ops, Lie group/algebra
 * operations, exponential/log maps, bracket, BCH, rigid body,
 * heavy top, satellite, quadrotor, and geometric integration.
 *
 * All tests use standard assert().
 * ============================================================== */

#include "lie_core.h"
#include "lie_actions.h"
#include "lie_reduction.h"
#include "lie_integration.h"
#include "lie_mechanics.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define TEST_TOL 1e-10

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_run = 0;
static int tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  TEST %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASSED\n"); } while(0)

/* ---- Matrix Tests ---- */
static void test_matrix_create_free(void) {
    TEST("matrix_create_free");
    LieMatrix *m = lie_matrix_create(3, 3);
    assert(m != NULL);
    assert(m->rows == 3 && m->cols == 3);
    lie_matrix_free(m);
    PASS();
}

static void test_matrix_set_get(void) {
    TEST("matrix_set_get");
    LieMatrix *m = lie_matrix_create(2, 2);
    lie_matrix_set(m, 0, 1, 3.14);
    double v = lie_matrix_get(m, 0, 1);
    assert(fabs(v - 3.14) < TEST_TOL);
    lie_matrix_free(m);
    PASS();
}

static void test_matrix_identity(void) {
    TEST("matrix_identity");
    LieMatrix *I = lie_matrix_identity(3);
    assert(lie_matrix_is_identity(I, TEST_TOL));
    lie_matrix_free(I);
    PASS();
}

static void test_matrix_mul(void) {
    TEST("matrix_mul");
    LieMatrix *A = lie_matrix_create(2, 2);
    LieMatrix *B = lie_matrix_create(2, 2);
    lie_matrix_set(A, 0, 0, 1); lie_matrix_set(A, 0, 1, 2);
    lie_matrix_set(A, 1, 0, 3); lie_matrix_set(A, 1, 1, 4);
    lie_matrix_set(B, 0, 0, 5); lie_matrix_set(B, 0, 1, 6);
    lie_matrix_set(B, 1, 0, 7); lie_matrix_set(B, 1, 1, 8);
    LieMatrix *C = lie_matrix_mul(A, B);
    assert(fabs(lie_matrix_get(C, 0, 0) - 19) < TEST_TOL);
    assert(fabs(lie_matrix_get(C, 0, 1) - 22) < TEST_TOL);
    assert(fabs(lie_matrix_get(C, 1, 0) - 43) < TEST_TOL);
    assert(fabs(lie_matrix_get(C, 1, 1) - 50) < TEST_TOL);
    lie_matrix_free(A); lie_matrix_free(B); lie_matrix_free(C);
    PASS();
}

static void test_matrix_transpose(void) {
    TEST("matrix_transpose");
    LieMatrix *A = lie_matrix_create(2, 3);
    lie_matrix_set(A, 0, 0, 1); lie_matrix_set(A, 0, 1, 2); lie_matrix_set(A, 0, 2, 3);
    lie_matrix_set(A, 1, 0, 4); lie_matrix_set(A, 1, 1, 5); lie_matrix_set(A, 1, 2, 6);
    LieMatrix *At = lie_matrix_transpose(A);
    assert(At->rows == 3 && At->cols == 2);
    assert(fabs(lie_matrix_get(At, 0, 0) - 1) < TEST_TOL);
    assert(fabs(lie_matrix_get(At, 2, 1) - 6) < TEST_TOL);
    lie_matrix_free(A); lie_matrix_free(At);
    PASS();
}

static void test_matrix_det(void) {
    TEST("matrix_det");
    LieMatrix *A = lie_matrix_create(3, 3);
    lie_matrix_set(A, 0, 0, 1); lie_matrix_set(A, 0, 1, 0); lie_matrix_set(A, 0, 2, 0);
    lie_matrix_set(A, 1, 0, 0); lie_matrix_set(A, 1, 1, 2); lie_matrix_set(A, 1, 2, 0);
    lie_matrix_set(A, 2, 0, 0); lie_matrix_set(A, 2, 1, 0); lie_matrix_set(A, 2, 2, 3);
    double d = lie_matrix_det(A);
    assert(fabs(d - 6.0) < TEST_TOL);
    lie_matrix_free(A);
    PASS();
}

static void test_matrix_inverse(void) {
    TEST("matrix_inverse");
    LieMatrix *A = lie_matrix_create(2, 2);
    lie_matrix_set(A, 0, 0, 4); lie_matrix_set(A, 0, 1, 7);
    lie_matrix_set(A, 1, 0, 2); lie_matrix_set(A, 1, 1, 6);
    LieMatrix *Ainv = lie_matrix_inverse(A);
    LieMatrix *I = lie_matrix_mul(A, Ainv);
    assert(lie_matrix_is_identity(I, TEST_TOL));
    lie_matrix_free(A); lie_matrix_free(Ainv); lie_matrix_free(I);
    PASS();
}

static void test_matrix_orthogonal_check(void) {
    TEST("matrix_orthogonal_check");
    LieMatrix *I = lie_matrix_identity(3);
    assert(lie_matrix_is_special_orthogonal(I, TEST_TOL));
    lie_matrix_free(I);
    PASS();
}

/* ---- Vector Tests ---- */
static void test_vector_cross(void) {
    TEST("vector_cross3");
    LieVector *a = lie_vector_create(3);
    LieVector *b = lie_vector_create(3);
    a->data[0] = 1; a->data[1] = 0; a->data[2] = 0;
    b->data[0] = 0; b->data[1] = 1; b->data[2] = 0;
    LieVector *c = lie_vector_cross3(a, b);
    assert(fabs(c->data[2] - 1.0) < TEST_TOL);
    lie_vector_free(a); lie_vector_free(b); lie_vector_free(c);
    PASS();
}

static void test_vector_norm(void) {
    TEST("vector_norm");
    LieVector *v = lie_vector_create(3);
    v->data[0] = 3; v->data[1] = 4; v->data[2] = 0;
    double n = lie_vector_norm(v);
    assert(fabs(n - 5.0) < TEST_TOL);
    lie_vector_free(v);
    PASS();
}

/* ---- Lie Group Tests ---- */
static void test_group_create_identity(void) {
    TEST("group_create_identity");
    LieGroupElement *g = lie_group_create(LIE_GROUP_SO3, "test");
    assert(g != NULL);
    assert(lie_matrix_is_identity(g->g, TEST_TOL));
    lie_group_free(g);
    PASS();
}

static void test_group_mul(void) {
    TEST("group_mul");
    LieGroupElement *g1 = lie_group_create(LIE_GROUP_SO3, "g1");
    LieGroupElement *g2 = lie_group_create(LIE_GROUP_SO3, "g2");
    LieGroupElement *prod = lie_group_mul(g1, g2);
    assert(prod != NULL);
    assert(lie_matrix_is_identity(prod->g, TEST_TOL));
    lie_group_free(g1); lie_group_free(g2); lie_group_free(prod);
    PASS();
}

/* ---- SO(3) exp/log Tests ---- */
static void test_exp_so3_identity(void) {
    TEST("exp_so3_zero");
    LieVector *omega = lie_vector_create(3);
    LieGroupElement *R = lie_exp_so3(omega);
    assert(lie_matrix_is_identity(R->g, TEST_TOL));
    lie_vector_free(omega); lie_group_free(R);
    PASS();
}

static void test_exp_so3_90deg(void) {
    TEST("exp_so3_90deg_z");
    LieVector *omega = lie_vector_create(3);
    omega->data[2] = M_PI / 2.0;
    LieGroupElement *R = lie_exp_so3(omega);
    /* R should be rotation by 90 deg about z */
    assert(fabs(lie_matrix_get(R->g, 0, 0) - 0.0) < TEST_TOL);
    assert(fabs(lie_matrix_get(R->g, 0, 1) + 1.0) < TEST_TOL);
    assert(fabs(lie_matrix_get(R->g, 1, 0) - 1.0) < TEST_TOL);
    assert(fabs(lie_matrix_get(R->g, 2, 2) - 1.0) < TEST_TOL);
    assert(lie_matrix_is_special_orthogonal(R->g, TEST_TOL));
    lie_vector_free(omega); lie_group_free(R);
    PASS();
}

static void test_exp_log_roundtrip(void) {
    TEST("exp_log_so3_roundtrip");
    LieVector *omega = lie_vector_create(3);
    omega->data[0] = 0.5; omega->data[1] = 0.3; omega->data[2] = -0.2;
    LieGroupElement *R = lie_exp_so3(omega);
    LieVector *omega2 = lie_log_so3(R);
    assert(fabs(omega2->data[0] - omega->data[0]) < TEST_TOL);
    assert(fabs(omega2->data[1] - omega->data[1]) < TEST_TOL);
    assert(fabs(omega2->data[2] - omega->data[2]) < TEST_TOL);
    lie_vector_free(omega); lie_vector_free(omega2); lie_group_free(R);
    PASS();
}

/* ---- Lie Algebra Tests ---- */
static void test_bracket_so3(void) {
    TEST("bracket_so3");
    LieAlgebraElement *xi = lie_algebra_create(LIE_ALG_SO3, "xi");
    LieAlgebraElement *eta = lie_algebra_create(LIE_ALG_SO3, "eta");
    xi->coords->data[0] = 1; eta->coords->data[1] = 1;
    LieAlgebraElement *b = lie_bracket(xi, eta);
    assert(fabs(b->coords->data[2] - 1.0) < TEST_TOL);
    lie_algebra_free(xi); lie_algebra_free(eta); lie_algebra_free(b);
    PASS();
}

static void test_jacobi_identity(void) {
    TEST("jacobi_identity_so3");
    LieAlgebraElement *a = lie_algebra_create(LIE_ALG_SO3, "a");
    LieAlgebraElement *b = lie_algebra_create(LIE_ALG_SO3, "b");
    LieAlgebraElement *c = lie_algebra_create(LIE_ALG_SO3, "c");
    a->coords->data[0]=1;a->coords->data[1]=2;a->coords->data[2]=3;
    b->coords->data[0]=4;b->coords->data[1]=5;b->coords->data[2]=6;
    c->coords->data[0]=7;c->coords->data[1]=8;c->coords->data[2]=9;
    assert(lie_jacobi_identity_check(a, b, c, TEST_TOL));
    lie_algebra_free(a); lie_algebra_free(b); lie_algebra_free(c);
    PASS();
}

static void test_killing_form(void) {
    TEST("killing_form_so3");
    LieAlgebraElement *a = lie_algebra_create(LIE_ALG_SO3, "a");
    a->coords->data[0] = 1; a->coords->data[1] = 2; a->coords->data[2] = 3;
    LieAlgebraElement *b = lie_algebra_create(LIE_ALG_SO3, "b");
    b->coords->data[0] = 4; b->coords->data[1] = 5; b->coords->data[2] = 6;
    double k = lie_killing_form(a, b);
    double expected = -2.0 * (1.0*4.0 + 2.0*5.0 + 3.0*6.0);
    assert(fabs(k - expected) < TEST_TOL);
    lie_algebra_free(a); lie_algebra_free(b);
    PASS();
}

/* ---- Hat/Vee Tests ---- */
static void test_hat_vee_roundtrip(void) {
    TEST("hat_vee_so3_roundtrip");
    LieVector *omega = lie_vector_create(3);
    omega->data[0]=1; omega->data[1]=2; omega->data[2]=3;
    LieMatrix *hat = lie_so3_hat(omega);
    assert(lie_matrix_is_skew_symmetric(hat, TEST_TOL));
    LieVector *omega2 = lie_so3_vee(hat);
    assert(fabs(omega2->data[0]-1)<TEST_TOL);
    assert(fabs(omega2->data[1]-2)<TEST_TOL);
    assert(fabs(omega2->data[2]-3)<TEST_TOL);
    lie_vector_free(omega); lie_matrix_free(hat); lie_vector_free(omega2);
    PASS();
}

/* ---- BCH Formula Tests ---- */
static void test_bch_2nd_order(void) {
    TEST("bch_2nd_order");
    LieAlgebraElement *A = lie_algebra_create(LIE_ALG_SO3, "A");
    LieAlgebraElement *B = lie_algebra_create(LIE_ALG_SO3, "B");
    A->coords->data[0] = 0.1; B->coords->data[1] = 0.2;
    LieAlgebraElement *bch = lie_bch_2nd_order(A, B);
    assert(bch != NULL);
    assert(fabs(bch->coords->data[0] - 0.1) < TEST_TOL);
    assert(fabs(bch->coords->data[1] - 0.2) < TEST_TOL);
    assert(fabs(bch->coords->data[2] - 0.01) < 1e-3); /* 0.5 * [A,B]_z = 0.5*0.1*0.2 = 0.01 */
    lie_algebra_free(A); lie_algebra_free(B); lie_algebra_free(bch);
    PASS();
}

/* ---- Rigid Body Tests ---- */
static void test_rigid_body_create(void) {
    TEST("rigid_body_create");
    LieRigidBodyState *rb = lie_rigid_body_create(1, 2, 3, 0, 0, 0);
    assert(rb != NULL);
    assert(fabs(rb->I_body->inertia->data[0] - 1.0) < TEST_TOL);
    lie_rigid_body_free(rb);
    PASS();
}

static void test_rigid_body_rhs_no_torque(void) {
    TEST("rigid_body_rhs_zero_torque");
    LieRigidBodyState *rb = lie_rigid_body_create(1, 2, 3, 0, 0, 0);
    rb->Omega->data[0] = 1.0;
    LieVector *Pi_dot = lie_vector_create(3);
    LieMatrix *R_dot = lie_matrix_create(3, 3);
    lie_rigid_body_rhs(rb, Pi_dot, R_dot);
    /* For Omega = (1,0,0), Pi = (1,0,0), PixO = (0,0,0) */
    assert(fabs(Pi_dot->data[0]) < TEST_TOL);
    lie_vector_free(Pi_dot); lie_matrix_free(R_dot);
    lie_rigid_body_free(rb);
    PASS();
}

static void test_rigid_body_energy(void) {
    TEST("rigid_body_energy");
    LieRigidBodyState *rb = lie_rigid_body_create(1, 2, 3, 0, 0, 0);
    rb->Omega->data[0] = 2.0;
    double T = lie_rigid_body_kinetic_energy(rb);
    assert(fabs(T - 2.0) < TEST_TOL); /* 0.5 * 1 * 2^2 = 2 */
    lie_rigid_body_free(rb);
    PASS();
}

/* ---- Heavy Top Tests ---- */
static void test_heavy_top_energy(void) {
    TEST("heavy_top_energy");
    LieHeavyTopState *ht = lie_heavy_top_create(1, 1, 2, 1, 1);
    ht->Omega->data[2] = 2.0;
    double E = lie_heavy_top_energy(ht);
    assert(E > 0);
    lie_heavy_top_free(ht);
    PASS();
}

/* ---- Inertia Tensor Tests ---- */
static void test_inertia_apply(void) {
    TEST("inertia_apply");
    LieInertiaTensor *I = lie_inertia_create(2, 3, 4, 0, 0, 0);
    LieVector *omega = lie_vector_create(3);
    omega->data[0] = 1;
    LieVector *L = lie_inertia_apply(I, omega);
    assert(fabs(L->data[0] - 2.0) < TEST_TOL);
    lie_vector_free(omega); lie_vector_free(L); lie_inertia_free(I);
    PASS();
}

/* ---- SO(3) Adjoint Tests ---- */
static void test_Ad_so3_identity(void) {
    TEST("Ad_identity_so3");
    LieMatrix *R = lie_matrix_identity(3);
    LieMatrix *Ad = lie_Ad_matrix_so3(R);
    assert(lie_matrix_is_identity(Ad, TEST_TOL));
    lie_matrix_free(R); lie_matrix_free(Ad);
    PASS();
}

static void test_ad_so3(void) {
    TEST("ad_so3");
    LieVector *omega = lie_vector_create(3);
    omega->data[0] = 1;
    LieMatrix *ad = lie_ad_matrix_so3(omega);
    assert(lie_matrix_is_skew_symmetric(ad, TEST_TOL));
    lie_vector_free(omega); lie_matrix_free(ad);
    PASS();
}

/* ---- SE(3) Tests ---- */
static void test_exp_log_se3_roundtrip(void) {
    TEST("exp_log_se3_roundtrip");
    LieVector *xi = lie_vector_create(6);
    xi->data[2] = 1.0;  /* translation along z */
    xi->data[3] = 0.2;   /* rotation about x */
    LieGroupElement *g = lie_exp_se3(xi);
    assert(g != NULL);
    LieVector *xi2 = lie_log_se3(g);
    assert(fabs(xi2->data[2] - 1.0) < 1e-6);
    lie_vector_free(xi); lie_group_free(g); lie_vector_free(xi2);
    PASS();
}

/* ---- Momentum Map Test ---- */
static void test_angular_momentum(void) {
    TEST("angular_momentum");
    LieVector *r = lie_vector_create(3);
    LieVector *p = lie_vector_create(3);
    r->data[0] = 1; p->data[1] = 1;
    LieVector *L = lie_angular_momentum_compute(r, p);
    assert(fabs(L->data[2] - 1.0) < TEST_TOL);
    lie_vector_free(r); lie_vector_free(p); lie_vector_free(L);
    PASS();
}

/* ---- Free Top Tests ---- */
static void test_free_top_init(void) {
    TEST("free_top_init");
    LieFreeTop *top = lie_free_top_create(1, 2, 3);
    lie_free_top_init(top, 1, 0, 0);
    assert(fabs(top->energy - 0.5) < TEST_TOL);
    assert(fabs(top->momentum_sq - 1.0) < TEST_TOL);
    lie_free_top_free(top);
    PASS();
}

/* ---- Satellite Tests ---- */
static void test_satellite_create(void) {
    TEST("satellite_create");
    LieSatelliteState *sat = lie_satellite_create(100, 200, 300);
    assert(sat != NULL);
    assert(sat->n_orbit > 0);
    lie_satellite_free(sat);
    PASS();
}

/* ---- Quadrotor Tests ---- */
static void test_quadrotor_create(void) {
    TEST("quadrotor_create");
    LieQuadrotorState *qr = lie_quadrotor_create(1.0, 0.01, 0.01, 0.02);
    assert(qr != NULL);
    assert(qr->mass == 1.0);
    lie_quadrotor_free(qr);
    PASS();
}

/* ---- Geodesic Tests ---- */
static void test_geodesic_distance(void) {
    TEST("geodesic_distance_so3");
    LieMatrix *R1 = lie_matrix_identity(3);
    /* R2 = 90-degree rotation about z */
    LieMatrix *R2 = lie_matrix_create(3, 3);
    lie_matrix_set(R2, 0, 0, 0); lie_matrix_set(R2, 0, 1, -1);
    lie_matrix_set(R2, 1, 0, 1); lie_matrix_set(R2, 1, 1, 0);
    lie_matrix_set(R2, 2, 2, 1);
    double dist = lie_geodesic_distance_so3(R1, R2);
    assert(fabs(dist - M_PI/2.0) < TEST_TOL);
    lie_matrix_free(R1); lie_matrix_free(R2);
    PASS();
}

/* ---- Manifold State Tests ---- */
static void test_manifold_state(void) {
    TEST("manifold_state");
    LieManifoldState *s = lie_manifold_state_create(3, "test");
    assert(s != NULL);
    assert(s->manifold_dim == 3);
    lie_manifold_state_free(s);
    PASS();
}

/* ---- Transform Tests ---- */
static void test_transform_compose(void) {
    TEST("transform_compose");
    LieRigidTransform *tf1 = lie_rigid_transform_create(3);
    LieRigidTransform *tf2 = lie_rigid_transform_create(3);
    tf1->translation->data[0] = 1.0;
    tf2->translation->data[1] = 2.0;
    LieRigidTransform *tf3 = lie_transform_compose(tf1, tf2);
    assert(fabs(tf3->translation->data[0] - 1.0) < TEST_TOL);
    assert(fabs(tf3->translation->data[1] - 2.0) < TEST_TOL);
    lie_rigid_transform_free(tf1); lie_rigid_transform_free(tf2); lie_rigid_transform_free(tf3);
    PASS();
}

/* ============================================================== */

int main(void) {
    printf("============================================\n");
    printf("  mini-lie-group-mechanics Test Suite\n");
    printf("============================================\n");

    test_matrix_create_free();
    test_matrix_set_get();
    test_matrix_identity();
    test_matrix_mul();
    test_matrix_transpose();
    test_matrix_det();
    test_matrix_inverse();
    test_matrix_orthogonal_check();

    test_vector_cross();
    test_vector_norm();

    test_group_create_identity();
    test_group_mul();

    test_exp_so3_identity();
    test_exp_so3_90deg();
    test_exp_log_roundtrip();

    test_bracket_so3();
    test_jacobi_identity();
    test_killing_form();

    test_hat_vee_roundtrip();
    test_bch_2nd_order();

    test_rigid_body_create();
    test_rigid_body_rhs_no_torque();
    test_rigid_body_energy();

    test_heavy_top_energy();
    test_inertia_apply();
    test_Ad_so3_identity();
    test_ad_so3();
    test_exp_log_se3_roundtrip();
    test_angular_momentum();
    test_free_top_init();
    test_satellite_create();
    test_quadrotor_create();
    test_geodesic_distance();
    test_manifold_state();
    test_transform_compose();

    printf("\n============================================\n");
    printf("  Results: %d/%d tests passed\n", tests_passed, tests_run);
    printf("============================================\n");

    return (tests_passed == tests_run) ? 0 : 1;
}

