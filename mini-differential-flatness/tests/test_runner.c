#include "flatness_core.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "flatness_applications.h"

int main(void) {
  int ok = 0, fail = 0;
  printf("=== Tests ===\n\n");

  { FlatPolynomial p; flat_poly_init(&p, 3);
    if (p.num_terms == 0 && p.num_vars == 3) { printf("PASS: poly_init\n"); ok++; }
    else { printf("FAIL: poly_init\n"); fail++; } }

  { FlatPolynomial p; flat_poly_init(&p, 1);
    int pw2[20] = {2}, pw0[20] = {0};
    flat_poly_add_term(&p, 1.0, pw2); flat_poly_add_term(&p, -4.0, pw0);
    double x[] = {3.0}; double v = flat_poly_evaluate(&p, x);
    if (fabs(v - 5.0) < 1e-10) { printf("PASS: poly_eval\n"); ok++; }
    else { printf("FAIL: poly_eval got %f\n", v); fail++; } }

  { FlatVectorField vf; flat_vf_init(&vf, 3, "test");
    if (vf.dimension == 3) { printf("PASS: vf_init\n"); ok++; }
    else { printf("FAIL: vf_init\n"); fail++; } }

  { FlatControlAffineSystem s; int r = flat_make_quadrotor_system(&s);
    if (r == 0 && s.state_dim == 12) { printf("PASS: quadrotor\n"); ok++; }
    else { printf("FAIL: quadrotor\n"); fail++; } }

  { FlatControlAffineSystem s; int r = flat_make_crane_system(&s);
    if (r == 0) { printf("PASS: crane\n"); ok++; }
    else { printf("FAIL: crane\n"); fail++; } }

  { double A[4] = {0,1,0,0}, B[2] = {0,1};
    if (flat_pbh_test(A, B, 2, 1) == 1) { printf("PASS: pbh\n"); ok++; }
    else { printf("FAIL: pbh\n"); fail++; } }

  { double M[6] = {1,0,0,0,1,0}; int r = flat_matrix_rank(M, 2, 3, 1e-10);
    if (r == 2) { printf("PASS: rank\n"); ok++; }
    else { printf("FAIL: rank got %d\n", r); fail++; } }

  { double A[4] = {4,2,2,3}, L[4]; int r = flat_cholesky_decomposition(A, 2, L);
    if (r == 0) { printf("PASS: cholesky\n"); ok++; }
    else { printf("FAIL: cholesky\n"); fail++; } }

  { double A[4] = {0,1,-2,-3}, lr[2], li[2]; int r = flat_eigenvalues_2x2(A, lr, li);
    if (r == 0) { printf("PASS: eigenvalues\n"); ok++; }
    else { printf("FAIL: eigenvalues\n"); fail++; } }

  { QuadrotorParams pr = {1,0.25,0.01,0.01,0.02,1,0.1}; double u[]={10,0,0,0},w[4];
    int r = flat_quadrotor_rotor_speeds(&pr,u,w);
    if (r == 0) { printf("PASS: rotor\n"); ok++; }
    else { printf("FAIL: rotor\n"); fail++; } }

  printf("\n=== %d passed, %d failed ===\n", ok, fail);
  return fail > 0 ? 1 : 0;
}
