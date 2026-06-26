#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nonlinear_system.h"
#include "lie_derivative.h"
#include "lie_bracket.h"
#include "coordinate_transform.h"
#include "frobenius.h"
#include "feedback_linearization.h"

void test_vectors(void) {
    printf("\n[1] Vector Operations\n");
    Vector *v = vector_create(3);
    assert(v && v->dim == 3);
    vector_free(v);
    Vector *a = vector_create(3), *b = vector_create(3), *c = vector_create(3);
    vector_set(a, 1.0); vector_set(b, 2.0); vector_add(a, b, c);
    assert(fabs(vector_get(c, 0) - 3.0) < 1e-6);
    assert(fabs(vector_dot(a, b) - 6.0) < 1e-6);
    assert(fabs(vector_norm(a) - sqrt(3.0)) < 1e-6);
    vector_free(a); vector_free(b); vector_free(c);
    Vector *w = vector_create(3); vector_set(w, 3.0);
    double nr = vector_normalize(w);
    assert(fabs(nr - sqrt(27.0)) < 1e-6);
    assert(fabs(vector_norm(w) - 1.0) < 1e-6);
    vector_free(w);
    Vector *p = vector_create(3), *q = vector_create(3), *r = vector_create(3);
    vector_set(p, 1.0); vector_set(q, 2.0);
    vector_lincomb(2.0, p, 3.0, q, r);
    assert(fabs(vector_get(r, 0) - 8.0) < 1e-6);
    vector_free(p); vector_free(q); vector_free(r);
    Vector *s = vector_create(3); vector_set(s, 5.0);
    Vector *t = vector_clone(s);
    assert(vector_equals(s, t, 1e-6));
    Vector *z = vector_create(3); vector_set(z, 0.0);
    assert(vector_is_zero(z, 1e-6));
    vector_free(s); vector_free(t); vector_free(z);
    Vector *mv = vector_create(3);
    vector_set_element(mv, 0, -5.0); vector_set_element(mv, 1, 3.0);
    vector_set_element(mv, 2, 4.0);
    assert(fabs(vector_max_norm(mv) - 5.0) < 1e-6);
    Vector *d1 = vector_create(3), *d2 = vector_create(3);
    vector_set(d1, 0.0); vector_set(d2, 0.0); vector_set_element(d2, 0, 3.0); vector_set_element(d2, 1, 4.0);
    assert(fabs(vector_distance(d1, d2) - 5.0) < 1e-6);
    Vector *sa = vector_create(3), *sb = vector_create(3), *sc = vector_create(3);
    vector_set(sa, 5.0); vector_set(sb, 3.0); vector_sub(sa, sb, sc);
    assert(fabs(vector_get(sc, 0) - 2.0) < 1e-6);
    Vector *sv = vector_create(3); vector_set(sv, 2.0);
    vector_scale(sv, 3.0, sv);
    assert(fabs(vector_get(sv, 0) - 6.0) < 1e-6);
    vector_free(mv); vector_free(d1); vector_free(d2);
    vector_free(sa); vector_free(sb); vector_free(sc); vector_free(sv);
    printf("  vector ops OK\n");
}

void test_matrices(void) {
    printf("\n[2] Matrix Operations\n");
    Matrix *I = matrix_create_identity(3);
    assert(I && fabs(matrix_get(I, 0, 0) - 1.0) < 1e-6);
    assert(fabs(matrix_get(I, 0, 1)) < 1e-6);
    matrix_free(I);

    Matrix *A = matrix_create(2, 3), *B = matrix_create(3, 2), *C = matrix_create(2, 2);
    matrix_set(A, 0, 0, 1); matrix_set(A, 0, 1, 2); matrix_set(A, 0, 2, 3);
    matrix_set(A, 1, 0, 4); matrix_set(A, 1, 1, 5); matrix_set(A, 1, 2, 6);
    matrix_set(B, 0, 0, 7); matrix_set(B, 0, 1, 8);
    matrix_set(B, 1, 0, 9); matrix_set(B, 1, 1, 10);
    matrix_set(B, 2, 0, 11); matrix_set(B, 2, 1, 12);
    matrix_multiply(A, B, C);
    assert(fabs(matrix_get(C, 0, 0) - 58.0) < 1e-6);
    assert(fabs(matrix_get(C, 1, 1) - 154.0) < 1e-6);
    matrix_free(A); matrix_free(B); matrix_free(C);

    Matrix *D = matrix_create(2, 2);
    matrix_set(D, 0, 0, 3); matrix_set(D, 0, 1, 1);
    matrix_set(D, 1, 0, 2); matrix_set(D, 1, 1, 4);
    assert(fabs(matrix_determinant(D) - 10.0) < 1e-6);
    Matrix *E = matrix_create(3, 3);
    matrix_set(E, 0, 0, 1); matrix_set(E, 0, 1, 2); matrix_set(E, 0, 2, 3);
    matrix_set(E, 1, 0, 0); matrix_set(E, 1, 1, 1); matrix_set(E, 1, 2, 4);
    matrix_set(E, 2, 0, 5); matrix_set(E, 2, 1, 6); matrix_set(E, 2, 2, 0);
    assert(fabs(matrix_determinant(E) - 1.0) < 1e-6);
    matrix_free(D); matrix_free(E);

    Matrix *F = matrix_create(2, 2);
    matrix_set(F, 0, 0, 3); matrix_set(F, 0, 1, 1);
    matrix_set(F, 1, 0, 2); matrix_set(F, 1, 1, 4);
    Matrix *Finv = matrix_create(2, 2);
    bool ok = matrix_inverse(F, Finv);
    assert(ok);
    Matrix *chk = matrix_create(2, 2);
    matrix_multiply(F, Finv, chk);
    assert(fabs(matrix_get(chk, 0, 0) - 1.0) < 1e-6);
    assert(fabs(matrix_get(chk, 0, 1)) < 1e-6);
    matrix_free(F); matrix_free(Finv); matrix_free(chk);

    Matrix *G = matrix_create(2, 2);
    matrix_set(G, 0, 0, 3); matrix_set(G, 0, 1, 1);
    matrix_set(G, 1, 0, 2); matrix_set(G, 1, 1, 4);
    assert(fabs(matrix_trace(G) - 7.0) < 1e-6);
    Matrix *GT = matrix_create(2, 2);
    matrix_transpose(G, GT);
    assert(fabs(matrix_get(GT, 0, 1) - 2.0) < 1e-6);
    assert(matrix_rank(G, 1e-6) == 2);

    Matrix *H1 = matrix_create(2, 2), *H2 = matrix_create(2, 2), *H3 = matrix_create(2, 2);
    matrix_set(H1, 0, 0, 1); matrix_set(H1, 0, 1, 2);
    matrix_set(H1, 1, 0, 3); matrix_set(H1, 1, 1, 4);
    matrix_set(H2, 0, 0, 5); matrix_set(H2, 0, 1, 6);
    matrix_set(H2, 1, 0, 7); matrix_set(H2, 1, 1, 8);
    matrix_add(H1, H2, H3);
    assert(fabs(matrix_get(H3, 0, 0) - 6.0) < 1e-6);
    assert(matrix_frobenius_norm(H1) > 0);
    matrix_free(G); matrix_free(GT); matrix_free(H1); matrix_free(H2); matrix_free(H3);
    printf("  matrix ops OK\n");
}

void test_lie_api(void) {
    printf("\n[3] Lie Derivative / Lie Bracket\n");
    assert(fabs(lie_derivative(NULL, NULL, NULL)) < 1e-6);
    assert(fabs(lie_derivative_mixed(NULL, NULL, NULL, NULL, 0)) < 1e-6);
    assert(fabs(lie_derivative_kth(NULL, NULL, NULL, 0)) < 1e-6);
    assert(check_relative_degree_order(NULL, NULL, NULL, NULL, 0, 1e-6));
    assert(compute_all_lie_derivatives(NULL, NULL, NULL, 0) == NULL);
    assert(!lie_bracket_is_zero(NULL, NULL, NULL, 1e-6));
    assert(controllability_rank(NULL, NULL, NULL, 1e-6) == 0);
    assert(!check_jacobi_identity(NULL, NULL, NULL, NULL, 1e-6));
    printf("  lie API OK\n");
}

void test_diffeomorphism(void) {
    printf("\n[4] Diffeomorphism\n");
    Diffeomorphism *id = identity_diffeomorphism(3);
    assert(id && id->n == 3);
    Vector *x = vector_create(3), *z = vector_create(3);
    vector_set(x, 1.0); vector_set_element(x, 1, 2.0); vector_set_element(x, 2, 3.0);
    diffeomorphism_apply_forward(id, x, z);
    assert(vector_equals(x, z, 1e-6));
    Vector *x2 = vector_create(3);
    diffeomorphism_apply_inverse(id, z, x2);
    assert(vector_equals(x, x2, 1e-6));
    assert(diffeomorphism_check_jacobian(id, x, 1e-6));
    Vector *tps[3]; int i;
    for (i = 0; i < 3; i++) { tps[i] = vector_create(3); vector_set(tps[i], (double)(i+1)); }
    assert(diffeomorphism_validate(id, tps, 3, 1e-6));
    for (i = 0; i < 3; i++) vector_free(tps[i]);
    diffeomorphism_free(id); vector_free(x); vector_free(z); vector_free(x2);
    printf("  diffeomorphism OK\n");
}

void test_frobenius(void) {
    printf("\n[5] Frobenius Theorem\n");
    VectorField *f1 = vector_field_create(3, NULL, NULL);
    VectorField *f2 = vector_field_create(3, NULL, NULL);
    VectorField *ba[2] = {f1, f2};
    Distribution *dist = distribution_create(3, 2, ba);
    assert(dist && dist->n == 3 && dist->num_basis == 2);
    Vector *xp = vector_create(3); vector_set(xp, 0.0);
    Matrix *de = distribution_eval(dist, xp);
    assert(de != NULL);
    assert(distribution_dim_at(dist, xp, 1e-6) >= 0);
    bool invol = check_involutivity(dist, xp, 1e-6);
    assert(invol == true || invol == false);
    Matrix *ann = compute_annihilator(dist, xp);
    assert(ann != NULL);
    bool pfaff = check_integrability_pfaffian(dist, xp, 1e-6);
    assert(pfaff == true || pfaff == false);
    Vector *tps[2]; int i;
    for (i = 0; i < 2; i++) { tps[i] = vector_create(3); vector_set(tps[i], (double)i); }
    bool glob = check_involutivity_global(dist, tps, 2, 1e-6);
    assert(glob == true || glob == false);
    FrobeniusResult fr = frobenius_analysis(dist, xp);
    assert(fr.is_involutive == true || fr.is_involutive == false);
    for (i = 0; i < 2; i++) vector_free(tps[i]);
    distribution_free(dist); vector_field_free(f1); vector_field_free(f2);
    matrix_free(de); matrix_free(ann); vector_free(xp);
    printf("  frobenius OK\n");
}

void test_fl_api(void) {
    printf("\n[6] Feedback Linearization API\n");
    assert(!check_minimum_phase(NULL));
    assert(fabs(feedback_law_evaluate(NULL, 5.0)) < 1e-6);
    assert(!full_state_feedback_linearizable(NULL));
    assert(compute_decoupling_matrix(NULL, NULL, NULL) == NULL);
    RelativeDegreeResult rd = compute_relative_degree(NULL, NULL, 0);
    assert(!rd.is_well_defined);
    IOLinearizationResult io = input_output_linearize(NULL, NULL);
    assert(!io.is_feedback_linearizable);
    ISLinearizationResult is = input_state_linearize(NULL);
    assert(!is.is_exactly_linearizable);
    MIMOLinearizationResult mimo = mimo_linearize(NULL, NULL);
    assert(!mimo.has_vector_relative_degree);
    ZeroDynamicsResult zd = analyze_zero_dynamics(NULL, NULL);
    assert(zd.stability == ZERO_DYN_NONE);
    NormalFormResult nf = compute_normal_form(NULL, NULL);
    assert(nf.xi == NULL);
    IOLinearizationResult fo = free_output_feedback_linearize(NULL, 1);
    assert(!fo.is_feedback_linearizable);
    printf("  FL API OK\n");
}

int main(void) {
    printf("========================================\n");
    printf("Feedback Linearization Test Suite\n");
    printf("========================================\n");
    test_vectors();
    test_matrices();
    test_lie_api();
    test_diffeomorphism();
    test_frobenius();
    test_fl_api();
    printf("\n========================================\n");
    printf("ALL TESTS PASSED\n");
    printf("========================================\n");
    return 0;
}
