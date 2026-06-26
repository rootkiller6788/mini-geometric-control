#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "lie_group.h"
#include "principal_bundle.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "connection.h"
#include "curvature.h"
#include "holonomy.h"
#include "parallel_transport.h"

static int tests_run = 0;
static int tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)
#define CHECK(cond, msg) do { if (cond) PASS(); else FAIL(msg); } while(0)

/* L1: U(1) group operations */
static void test_u1_basics(void) {
    TEST("U1 identity");
    U1Element e = u1_identity();
    CHECK(u1_is_identity(e), "identity not recognized");

    TEST("U1 multiply");
    U1Element a = u1_make(1.0);
    U1Element b = u1_make(2.0);
    U1Element c = u1_multiply(a, b);
    CHECK(fabs(c.theta - 3.0) < 1e-10, "multiplication wrong");

    TEST("U1 inverse");
    U1Element inv = u1_inverse(a);
    U1Element prod = u1_multiply(a, inv);
    CHECK(u1_is_identity(prod), "a * a^{-1} != identity");

    TEST("U1 distance");
    double d = u1_distance(u1_make(0.0), u1_make(0.1));
    CHECK(fabs(d - 0.1) < 1e-10, "distance wrong");

    TEST("U1 trace");
    double tr = u1_trace(u1_make(0.0));
    CHECK(fabs(tr - 2.0) < 1e-10, "trace of identity should be 2");

    TEST("U1 exp/log");
    U1Algebra X = u1_alg_make(1.5);
    U1Element g = u1_exp(X);
    U1Algebra Y = u1_log(g);
    CHECK(fabs(Y.a - 1.5) < 1e-10, "exp/log roundtrip failed");

    TEST("U1 power");
    U1Element p = u1_power(u1_make(0.5), 3);
    CHECK(fabs(p.theta - 1.5) < 1e-10, "power wrong");
}

/* L1: SU(2) group operations */
static void test_su2_basics(void) {
    TEST("SU2 identity");
    SU2Element e = su2_identity();
    CHECK(su2_is_identity(e), "identity not recognized");

    TEST("SU2 multiply associativity");
    SU2Element a = su2_make(0.5, 0.5, 0.5, 0.5);
    SU2Element b = su2_make(0.6, 0.0, 0.8, 0.0);
    SU2Element c = su2_make(0.0, 1.0, 0.0, 0.0);
    SU2Element ab_c = su2_multiply(su2_multiply(a, b), c);
    SU2Element a_bc = su2_multiply(a, su2_multiply(b, c));
    double d = su2_distance(ab_c, a_bc);
    CHECK(d < 1e-6, "associativity failed");

    TEST("SU2 inverse");
    SU2Element inv = su2_inverse(a);
    SU2Element prod = su2_multiply(a, inv);
    CHECK(su2_is_identity(prod), "a * a^{-1} != identity");

    TEST("SU2 exp/log roundtrip");
    SU2Algebra X = su2_alg_make(1.0, 0.5, -0.3);
    SU2Element g = su2_exp(X);
    SU2Algebra Y = su2_log(g);
    double dn = su2_alg_norm(su2_alg_add(Y, su2_alg_scale(-1.0, X)));
    CHECK(dn < 1e-8, "exp/log roundtrip failed");

    TEST("SU2 commutator Jacobi");
    SU2Algebra X1 = su2_alg_make(1, 0, 0);
    SU2Algebra X2 = su2_alg_make(0, 1, 0);
    SU2Algebra X3 = su2_alg_make(0, 0, 1);
    SU2Algebra t1 = su2_alg_commutator(X1, su2_alg_commutator(X2, X3));
    SU2Algebra t2 = su2_alg_commutator(X2, su2_alg_commutator(X3, X1));
    SU2Algebra t3 = su2_alg_commutator(X3, su2_alg_commutator(X1, X2));
    SU2Algebra jacobi = su2_alg_add(su2_alg_add(t1, t2), t3);
    CHECK(su2_alg_norm(jacobi) < 1e-10, "Jacobi identity failed");
}

/* L1: SO(3) group operations */
static void test_so3_basics(void) {
    TEST("SO3 identity");
    SO3Element I = so3_identity();
    CHECK(fabs(so3_trace(I) - 3.0) < 1e-10, "identity trace != 3");

    TEST("SO3 Rodrigues 90deg");
    SO3Element Rx = so3_rodrigues(1, 0, 0, M_PI/2.0);
    double v[3] = {0, 1, 0};
    double out[3];
    so3_apply(Rx, v, out);
    CHECK(fabs(out[2] - 1.0) < 1e-10 && fabs(out[0]) < 1e-10,
          "90deg rotation wrong");

    TEST("SO3 multiply");
    SO3Element R1 = so3_from_axis_angle(0, 0, 1, 1.0);
    SO3Element R2 = so3_from_axis_angle(0, 0, 1, -1.0);
    SO3Element R = so3_multiply(R1, R2);
    CHECK(so3_distance(R, I) < 1e-10, "R * R^{-1} != I");

    TEST("SO3 exp/log");
    SO3Algebra X = so3_alg_make(0.3, -0.4, 0.5);
    SO3Element R3 = so3_exp(X);
    SO3Algebra Y = so3_log(R3);
    CHECK(so3_alg_norm(so3_alg_commutator(X, so3_alg_scale(-1.0, Y))) < 1e-10
          || fabs(so3_alg_norm(X) - so3_alg_norm(Y)) < 1e-8,
          "SO3 exp/log roundtrip");
}

/* L3: Lattice geometry */
static void test_lattice_geometry(void) {
    TEST("Lattice creation");
    int L[3] = {4, 4, 4};
    LatticeGeometry geom = lattice_create(3, L, 1.0);
    CHECK(geom.n_sites == 64 && geom.n_links == 192, "lattice size wrong");

    TEST("Site index/coords roundtrip");
    int coords[4] = {2, 3, 1, 0};
    int idx = lattice_site_index(&geom, coords);
    int coords2[4];
    lattice_site_coords(&geom, idx, coords2);
    CHECK(coords[0]==coords2[0] && coords[1]==coords2[1] && coords[2]==coords2[2],
          "coords roundtrip failed");

    TEST("Neighbor");
    int nb = lattice_neighbor(&geom, idx, 0, 1);
    int coords_nb[4];
    lattice_site_coords(&geom, nb, coords_nb);
    CHECK(coords_nb[0]==3 && coords_nb[1]==3 && coords_nb[2]==1,
          "neighbor wrong");
}

/* L3: Principal bundle */
static void test_principal_bundle(void) {
    int L[2] = {4, 4};
    LatticeGeometry geom = lattice_create(2, L, 1.0);

    TEST("Bundle creation");
    PrincipalBundle *pb = pb_create(LIE_U1, &geom);
    CHECK(pb != NULL && pb->initialized, "bundle creation failed");

    TEST("Bundle trivial check");
    CHECK(pb_is_trivial(pb), "new bundle should be trivial");

    TEST("Set/get link");
    LieGroupElement U;
    U.type = LIE_U1;
    U.elem.u1 = u1_make(0.5);
    pb_set_link(pb, 0, 0, U);
    LieGroupElement V = pb_get_link(pb, 0, 0);
    CHECK(fabs(V.elem.u1.theta - 0.5) < 1e-10, "set/get link failed");

    pb_free(pb);
}

/* L4: Plaquette and curvature */
static void test_curvature(void) {
    int L[2] = {4, 4};
    LatticeGeometry geom = lattice_create(2, L, 1.0);
    PrincipalBundle *pb = pb_create(LIE_U1, &geom);

    TEST("Plaquette of trivial bundle");
    Plaquette p = pb_plaquette(pb, 0, 0, 1);
    CHECK(u1_is_identity(p.U.elem.u1), "trivial plaquette != identity");

    TEST("Field strength of trivial bundle");
    LieAlgebraElement F = pb_field_strength(pb, 0, 0, 1);
    CHECK(fabs(F.alg.u1.a) < 1e-10, "F != 0 for trivial bundle");

    /* Set consistent link variables for a constant flux configuration.
     * U_0(x) = 1, U_1(x) = exp(i * flux * x_0).
     * This gives U_01(x) = exp(i * flux) for all x. */
    double flux = 0.5;
    for (int s = 0; s < geom.n_sites; s++) {
        int coords_test[4];
        lattice_site_coords(&geom, s, coords_test);
        LieGroupElement U0; U0.type = LIE_U1;
        U0.elem.u1 = u1_identity();
        pb_set_link(pb, s, 0, U0);
        LieGroupElement U1; U1.type = LIE_U1;
        U1.elem.u1 = u1_make(flux * (double)coords_test[0]);
        pb_set_link(pb, s, 1, U1);
    }

    TEST("Plaquette constant flux");
    Plaquette p2 = pb_plaquette(pb, 0, 0, 1);
    double theta = p2.U.elem.u1.theta;
    /* For constant flux, U_01 should equal exp(i*flux) */
    double expected = flux;
    if (expected > M_PI) expected -= 2.0*M_PI;
    CHECK(fabs(theta - flux) < 1e-10 || fabs(theta - flux - 2.0*M_PI) < 1e-10 ||
          fabs(theta - flux + 2.0*M_PI) < 1e-10,
          "plaquette should have constant flux");

    TEST("Wilson action non-negative");
    double action = pb_wilson_action_total(pb, 1.0);
    CHECK(action >= 0.0, "Wilson action negative");

    TEST("Bianchi identity trivial 2D");
    bool bianchi = pb_check_bianchi(pb, 0, 0, 1, 0, 1e-10);
    CHECK(bianchi, "Bianchi should hold");

    TEST("Chern number trivial bundle");
    double c1 = pb_chern_number_1(pb);
    CHECK(fabs(c1) < 1e-10, "Chern number of trivial bundle != 0");

    pb_free(pb);
}

/* L5: Parallel transport and holonomy */
static void test_holonomy(void) {
    int L[2] = {4, 4};
    LatticeGeometry geom = lattice_create(2, L, 1.0);
    PrincipalBundle *pb = pb_create(LIE_U1, &geom);

    TEST("Path creation");
    LatticePath *path = path_straight(&geom, 0, 0, 3);
    CHECK(path != NULL && path->length == 3, "path creation failed");

    TEST("Path validity");
    CHECK(path_is_valid(path, &geom), "path should be valid");

    TEST("Parallel transport trivial");
    LieGroupElement transport = pb_parallel_transport_path(pb, path);
    CHECK(u1_is_identity(transport.elem.u1), "trivial transport != identity");

    path_free(path);

    TEST("Rectangle path");
    LatticePath *rect = path_rectangle(&geom, 0, 0, 1, 2, 2);
    CHECK(rect != NULL && rect->is_loop, "rectangle path failed");

    TEST("Wilson loop trivial");
    double W = pb_wilson_loop_rect(pb, 0, 0, 1, 2, 2);
    CHECK(fabs(W - 2.0) < 1e-10, "Wilson loop of trivial bundle != N");

    path_free(rect);

    TEST("Creutz ratio trivial");
    double chi = pb_creutz_ratio(pb, 0, 0, 1, 2, 2);
    CHECK(fabs(chi) < 1e-10, "Creutz ratio of trivial != 0");

    pb_free(pb);
}

/* L5: Gauge transformation */
static void test_gauge_transform(void) {
    int L[2] = {4, 4};
    LatticeGeometry geom = lattice_create(2, L, 1.0);
    PrincipalBundle *pb = pb_create(LIE_U1, &geom);

    /* Set nontrivial links first */
    LieGroupElement U;
    U.type = LIE_U1;
    for (int s = 0; s < geom.n_sites; s++) {
        for (int mu = 0; mu < geom.dim; mu++) {
            U.elem.u1 = u1_make(0.1 * (s + mu));
            pb_set_link(pb, s, mu, U);
        }
    }

    double action_before = pb_wilson_action_total(pb, 1.0);

    TEST("Gauge invariance of Wilson action");
    pb_random_gauge_transform(pb, 42);
    pb_gauge_transform(pb);
    double action_after = pb_wilson_action_total(pb, 1.0);
    CHECK(fabs(action_before - action_after) < 1e-6,
          "Wilson action not gauge invariant!");

    pb_free(pb);
}

/* L2: Connection properties */
static void test_connection(void) {
    int L[3] = {4, 4, 4};
    LatticeGeometry geom = lattice_create(3, L, 1.0);
    PrincipalBundle *pb = pb_create(LIE_U1, &geom);
    LatticeConnection *lc = lc_create(pb);

    TEST("Connection creation");
    CHECK(lc != NULL, "connection creation failed");

    TEST("Connection flat check");
    CHECK(lc_is_flat(lc), "trivial connection should be flat");

    TEST("Connection verify");
    CHECK(lc_verify_properties(lc) == 0, "connection properties invalid");

    lc_free(lc);
    pb_free(pb);
}

/* L8: SU(2) <-> SO(3) double cover */
static void test_double_cover(void) {
    TEST("SU2 -> SO3 -> SU2 roundtrip");
    SU2Element q = su2_make(0.5, 0.5, 0.5, 0.5);
    SO3Element R = su2_to_so3(q);
    SU2Element q2 = so3_to_su2(R);
    double d = su2_distance(q, q2);
    /* Since SU(2)->SO(3) is 2:1, q and -q map to same R */
    SU2Element q_neg = su2_make(-q.a, -q.b, -q.c, -q.d);
    double d_neg = su2_distance(q, q_neg);
    CHECK(d < 1e-10 || d_neg < 1e-10, "double cover roundtrip failed");

    TEST("SU2 -> SO3 preserves multiplication");
    SU2Element qa = su2_make(0.8, 0.0, 0.6, 0.0);
    SU2Element qb = su2_make(0.6, 0.0, 0.0, 0.8);
    SO3Element Ra = su2_to_so3(qa);
    SO3Element Rb = su2_to_so3(qb);
    SO3Element Rab = so3_multiply(Ra, Rb);
    SO3Element R_direct = su2_to_so3(su2_multiply(qa, qb));
    CHECK(so3_distance(Rab, R_direct) < 1e-10, "homomorphism property failed");
}

int main(void) {
    printf("=== Connection on Principal Bundle Test Suite ===\n\n");

    printf("L1: Lie Group Operations\n");
    test_u1_basics();
    test_su2_basics();
    test_so3_basics();

    printf("\nL3: Lattice & Bundle\n");
    test_lattice_geometry();
    test_principal_bundle();

    printf("\nL4: Curvature & Plaquettes\n");
    test_curvature();

    printf("\nL5: Holonomy & Wilson Loops\n");
    test_holonomy();

    printf("\nL5: Gauge Transformations\n");
    test_gauge_transform();

    printf("\nL2: Connection Properties\n");
    test_connection();

    printf("\nL8: SU(2)/SO(3) Double Cover\n");
    test_double_cover();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
