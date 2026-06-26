#include "../include/geo_optctrl_core.h"
#include "../include/symplectic_geometry.h"
#include "../include/lie_group_methods.h"
#include "../include/pontryagin_maximum.h"
#include "../include/subriemannian_geometry.h"
#include "../include/geometric_integrators.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) printf("  TEST: %s ... ", name)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define CHECK(cond, msg) do { if (cond) PASS(); else FAIL(msg); } while(0)

void test_lie_bracket_so3(void) {
    TEST("so(3) Lie bracket = cross product");
    double X[3][3] = {{0,-3,2},{3,0,-1},{-2,1,0}};
    double Y[3][3] = {{0,-1,0},{1,0,0},{0,0,0}};
    double Z[3][3];
    geo_so3_bracket(X, Y, Z);
    double omega[3];
    geo_so3_vee(Z, omega);
    CHECK(fabs(omega[0]-2.0)<1e-10 && fabs(omega[1]+1.0)<1e-10 && fabs(omega[2]-0.0)<1e-10, "bracket mismatch");
}

void test_so3_exp_log(void) {
    TEST("SO(3) exp and log are inverses");
    double omega[3] = {0.3, -0.5, 0.2};
    double R[3][3];
    geo_so3_exp(omega, R);
    double omega2[3];
    int ret = geo_so3_log(R, omega2);
    CHECK(ret==0, "log returned error");
    CHECK(fabs(omega[0]-omega2[0])<1e-10, "exp/log roundtrip x");
}

void test_symplectic_form(void) {
    TEST("Symplectic form canonical structure");
    GeoSymplecticSpace s;
    geo_symplectic_init(&s, 3);
    CHECK(s.n==3 && s.dim_2n==6, "dimension");
    CHECK(s.matrix[0][3]==1.0, "upper right");
    CHECK(s.matrix[3][0]==-1.0, "lower left");
}

void test_liouville_form(void) {
    TEST("Liouville 1-form on T*Q");
    double q[2]={1,2}, p[2]={3,4}, theta[4];
    geo_liouville_one_form(q, p, theta, 2);
    CHECK(theta[0]==0.0 && theta[2]==3.0 && theta[3]==4.0, "Liouville form");
}

void test_euler_poincare(void) {
    TEST("Euler-Poincare rigid body");
    double I[3]={2,3,4}, O[3]={1,0.5,0.2}, t[3]={0,0,0}, dO[3];
    geo_euler_poincare_rigid_body(I, O, t, dO);
    CHECK(fabs(dO[0]+0.05)<1e-10, "Euler-Poincare dO0");
}

void test_so3_kinematics(void) {
    TEST("SO(3) kinematics Rdot = R*Omega_hat");
    double R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    double O[3] = {0,0,1}, dR[3][3];
    geo_so3_kinematics(R, O, dR);
    CHECK(fabs(dR[0][1]+1.0)<1e-10 && fabs(dR[1][0]-1.0)<1e-10, "kinematics");
}

void test_canonical_symplectic(void) {
    TEST("Canonical symplectic form matrix");
    double m[GEO_MAX_DIM][GEO_MAX_DIM];
    geo_canonical_symplectic_form(3, m);
    CHECK(m[0][3]==1.0 && m[3][0]==-1.0, "canonical form");
}

void test_heisenberg(void) {
    TEST("Heisenberg group [X,Y]=Z");
    double x[3]={0,0,0}, Xv[3], Yv[3];
    geo_heisenberg_X(x, Xv, NULL);
    geo_heisenberg_Y(x, Yv, NULL);
    CHECK(Xv[0]==1.0 && Yv[1]==1.0, "Heisenberg fields");
    GeoVectorField vfX={3,0,geo_heisenberg_X,NULL};
    GeoVectorField vfY={3,0,geo_heisenberg_Y,NULL};
    double br[3];
    geo_lie_bracket(&vfX, &vfY, x, br, 3, 1e-6);
    CHECK(fabs(br[2]-1.0)<1e-2, "[X,Y] != d/dz");
}

void test_arnold_cat_map(void) {
    TEST("Arnold cat map bounds");
    double x=0.3, y=0.7;
    geo_arnold_cat_map(&x, &y, 1);
    CHECK(x>=0.0 && x<1.0 && y>=0.0 && y<1.0, "cat map");
}

void test_yoshida_weights(void) {
    TEST("Yoshida weights order 4");
    double w[10]; int nw;
    geo_yoshida_weights(4, w, &nw);
    CHECK(nw==3, "Yoshida nw");
    CHECK(fabs(w[0]+w[1]+w[2]-1.0)<1e-10, "Yoshida sum");
}

void test_legendre_transform(void) {
    TEST("Legendre transform mechanical");
    double q[2]={0,0}, v[2]={3,4}, m=2.0, p[2], H;
    geo_legendre_transform(q,v,&m,NULL,NULL,p,&H,2);
    CHECK(fabs(p[0]-6.0)<1e-10, "momentum");
    CHECK(fabs(H-25.0)<1e-10, "energy");
}

void test_gauss_legendre(void) {
    TEST("Gauss-Legendre Butcher tableau");
    GeoButcherTableau t;
    geo_gauss_legendre_tableau(&t, 1);
    CHECK(t.stages==1 && t.c[0]==0.5 && t.b[0]==1.0, "1-stage GL");
    geo_gauss_legendre_tableau(&t, 2);
    CHECK(t.stages==2 && t.b[0]==0.5, "2-stage GL");
}

void test_product_manifold(void) {
    TEST("Product manifold embed/project");
    double x1[3]={1,2,3}, x2[2]={4,5}, r[5], r1[3], r2[2];
    geo_product_embed(x1,3,x2,2,r);
    CHECK(r[0]==1&&r[3]==4, "embed");
    geo_product_project(r,5,r1,3,r2,2);
    CHECK(r1[0]==1&&r2[0]==4, "project");
}

void test_so3_jacobi(void) {
    TEST("Jacobi identity on so(3)");
    double X[3][3]={{0,-1,0},{1,0,0},{0,0,0}};
    double Y[3][3]={{0,0,-1},{0,0,0},{1,0,0}};
    double Z[3][3]={{0,0,0},{0,0,-1},{0,1,0}};
    int ok = geo_so3_jacobi_identity_check(X,Y,Z,1e-10);
    CHECK(ok==1, "Jacobi identity");
}

void test_dido(void) {
    TEST("Dido problem circle");
    double L=10.0, pts[40];
    geo_dido_solution(L, 20, pts);
    double R=L/(2.0*3.1415926535);
    CHECK(fabs(sqrt(pts[0]*pts[0]+pts[1]*pts[1])-R)<1e-10, "Dido");
}

int main(void) {
    printf("=== Geometric Optimal Control Test Suite ===\n\n");
    test_lie_bracket_so3();
    test_so3_exp_log();
    test_symplectic_form();
    test_liouville_form();
    test_euler_poincare();
    test_so3_kinematics();
    test_canonical_symplectic();
    test_heisenberg();
    test_arnold_cat_map();
    test_yoshida_weights();
    test_legendre_transform();
    test_gauss_legendre();
    test_product_manifold();
    test_so3_jacobi();
    test_dido();
    printf("\n=== Results: %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}