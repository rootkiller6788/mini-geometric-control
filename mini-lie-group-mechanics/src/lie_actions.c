/* ==============================================================
 * lie_actions.c -- Group Actions on Manifolds
 *
 * Implements: manifold states, rigid transforms, momentum maps,
 * infinitesimal generators, orbits, left/right translations,
 * Maurer-Cartan form, coadjoint orbits, locked inertia tensor.
 *
 * References:
 *   Marsden & Ratiu (1999) Ch.4, 9, 11, 13
 *   Abraham & Marsden (1978) Foundations of Mechanics
 * ============================================================== */

#include "lie_actions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/* --- Manifold State --- */
LieManifoldState* lie_manifold_state_create(int dim, const char *name) {
    assert(dim > 0);
    LieManifoldState *s = (LieManifoldState*)malloc(sizeof(LieManifoldState));
    assert(s);
    s->name = name ? strdup(name) : strdup("unnamed");
    s->manifold_dim = dim;
    s->position = lie_vector_create(dim);
    s->velocity = lie_vector_create(dim);
    return s;
}
void lie_manifold_state_free(LieManifoldState *s) {
    if (s) { free(s->name); lie_vector_free(s->position); lie_vector_free(s->velocity); free(s); }
}

/* --- Rigid Transform (R, t): x -> R x + t --- */
LieRigidTransform* lie_rigid_transform_create(int dim) {
    assert(dim >= 2 && dim <= 3);
    LieRigidTransform *tf = (LieRigidTransform*)malloc(sizeof(LieRigidTransform));
    assert(tf);
    tf->rotation = lie_matrix_identity(dim);
    tf->translation = lie_vector_create(dim);
    return tf;
}
void lie_rigid_transform_free(LieRigidTransform *tf) {
    if (tf) { lie_matrix_free(tf->rotation); lie_vector_free(tf->translation); free(tf); }
}

LieVector* lie_transform_point(const LieRigidTransform *tf, const LieVector *p) {
    assert(tf && p && p->size == tf->rotation->cols);
    int n = p->size;
    LieVector *result = lie_vector_create(n);
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++) s += tf->rotation->data[i*n+j] * p->data[j];
        result->data[i] = s + tf->translation->data[i];
    }
    return result;
}

LieRigidTransform* lie_transform_compose(const LieRigidTransform *tf1, const LieRigidTransform *tf2) {
    assert(tf1 && tf2 && tf1->rotation->rows == tf2->rotation->rows);
    int n = tf1->rotation->rows;
    LieRigidTransform *tf = lie_rigid_transform_create(n);
    lie_matrix_free(tf->rotation);
    tf->rotation = lie_matrix_mul(tf1->rotation, tf2->rotation);
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++) s += tf1->rotation->data[i*n+j] * tf2->translation->data[j];
        tf->translation->data[i] = s + tf1->translation->data[i];
    }
    return tf;
}

LieRigidTransform* lie_transform_inverse(const LieRigidTransform *tf) {
    assert(tf);
    int n = tf->rotation->rows;
    LieRigidTransform *inv = lie_rigid_transform_create(n);
    lie_matrix_free(inv->rotation);
    inv->rotation = lie_matrix_transpose(tf->rotation);
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++) s -= inv->rotation->data[i*n+j] * tf->translation->data[j];
        inv->translation->data[i] = s;
    }
    return inv;
}

LieRigidTransform* lie_se3_to_transform(const LieGroupElement *g) {
    assert(g && g->group_type == LIE_GROUP_SE3);
    LieRigidTransform *tf = lie_rigid_transform_create(3);
    lie_matrix_free(tf->rotation);
    tf->rotation = lie_matrix_create(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) tf->rotation->data[i*3+j] = g->g->data[i*4+j];
        tf->translation->data[i] = g->g->data[i*4+3];
    }
    return tf;
}

LieGroupElement* lie_transform_to_se3(const LieRigidTransform *tf) {
    assert(tf && tf->rotation->rows == 3);
    LieGroupElement *g = lie_group_create(LIE_GROUP_SE3, "from_tf");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) g->g->data[i*4+j] = tf->rotation->data[i*3+j];
        g->g->data[i*4+3] = tf->translation->data[i];
    }
    for (int j = 0; j < 3; j++) g->g->data[3*4+j] = 0.0;
    g->g->data[15] = 1.0;
    return g;
}

/* --- Momentum Map --- */
LieMomentum* lie_momentum_create(void) {
    LieMomentum *mom = (LieMomentum*)malloc(sizeof(LieMomentum));
    assert(mom);
    mom->linear_momentum = lie_vector_create(3);
    mom->angular_momentum = lie_vector_create(3);
    return mom;
}
void lie_momentum_free(LieMomentum *mom) {
    if (mom) { lie_vector_free(mom->linear_momentum); lie_vector_free(mom->angular_momentum); free(mom); }
}

LieVector* lie_angular_momentum_compute(const LieVector *r, const LieVector *p) {
    assert(r && p && r->size == 3 && p->size == 3);
    return lie_vector_cross3(r, p);  /* L = r x p */
}

LieVector* lie_momentum_map_so3_position(const LieVector *q, const LieVector *p) {
    /* J(q,p) = q x p (angular momentum for SO(3) action on T*R^3) */
    return lie_angular_momentum_compute(q, p);
}

/* --- Infinitesimal Generator --- */
LieVector* lie_infinitesimal_generator_so3(const LieVector *omega, const LieVector *x) {
    /* xi_{R^3}(x) = omega x x for SO(3) acting on R^3 by rotation */
    return lie_vector_cross3(omega, x);
}

/* --- Orbit Check --- */
bool lie_same_so3_orbit(const LieVector *x1, const LieVector *x2, double tol) {
    assert(x1 && x2);
    /* Two points are on the same SO(3) orbit iff they have the same norm */
    double n1 = lie_vector_norm(x1);
    double n2 = lie_vector_norm(x2);
    return fabs(n1 - n2) < tol;
}

/* --- Left/Right Translations --- */
LieGroupElement* lie_left_translate(const LieGroupElement *g, const LieGroupElement *h) {
    return lie_group_mul(g, h);  /* L_g(h) = g h */
}
LieGroupElement* lie_right_translate(const LieGroupElement *g, const LieGroupElement *h) {
    return lie_group_mul(h, g);  /* R_g(h) = h g */
}

LieMatrix* lie_tangent_left_translate(const LieGroupElement *g, const LieMatrix *xi_matrix) {
    /* T_e L_g(xi) = g * xi */
    return lie_matrix_mul(g->g, xi_matrix);
}

/* --- Maurer-Cartan Form --- */
LieMatrix* lie_maurer_cartan_form(const LieGroupElement *g, const LieMatrix *tangent_vector) {
    /* omega_g(v) = g^{-1} * v */
    LieMatrix *ginv = lie_matrix_inverse(g->g);
    LieMatrix *result = lie_matrix_mul(ginv, tangent_vector);
    lie_matrix_free(ginv);
    return result;
}

/* --- Coadjoint Orbit Radius --- */
double lie_coadjoint_orbit_radius_so3(const LieVector *mu) {
    /* Coadjoint orbits of SO(3) are spheres: ||mu|| = const */
    return lie_vector_norm(mu);
}

/* --- Locked Inertia Tensor --- */
LieInertiaTensor* lie_inertia_create(double Ixx, double Iyy, double Izz,
                                      double Ixy, double Ixz, double Iyz) {
    LieInertiaTensor *I = (LieInertiaTensor*)malloc(sizeof(LieInertiaTensor));
    assert(I);
    I->inertia = lie_matrix_create(3, 3);
    lie_matrix_set(I->inertia, 0, 0, Ixx);
    lie_matrix_set(I->inertia, 0, 1, Ixy);
    lie_matrix_set(I->inertia, 0, 2, Ixz);
    lie_matrix_set(I->inertia, 1, 0, Ixy);
    lie_matrix_set(I->inertia, 1, 1, Iyy);
    lie_matrix_set(I->inertia, 1, 2, Iyz);
    lie_matrix_set(I->inertia, 2, 0, Ixz);
    lie_matrix_set(I->inertia, 2, 1, Iyz);
    lie_matrix_set(I->inertia, 2, 2, Izz);
    return I;
}
void lie_inertia_free(LieInertiaTensor *I) {
    if (I) { lie_matrix_free(I->inertia); free(I); }
}

LieVector* lie_inertia_apply(const LieInertiaTensor *I, const LieVector *omega) {
    assert(I && omega && omega->size == 3);
    LieVector *result = lie_vector_create(3);
    result->data[0] = I->inertia->data[0]*omega->data[0] + I->inertia->data[1]*omega->data[1] + I->inertia->data[2]*omega->data[2];
    result->data[1] = I->inertia->data[3]*omega->data[0] + I->inertia->data[4]*omega->data[1] + I->inertia->data[5]*omega->data[2];
    result->data[2] = I->inertia->data[6]*omega->data[0] + I->inertia->data[7]*omega->data[1] + I->inertia->data[8]*omega->data[2];
    return result;
}

LieInertiaTensor* lie_inertia_inverse(const LieInertiaTensor *I) {
    assert(I);
    LieMatrix *Iinv = lie_matrix_inverse_3x3(I->inertia);
    LieInertiaTensor *inv = (LieInertiaTensor*)malloc(sizeof(LieInertiaTensor));
    assert(inv);
    inv->inertia = Iinv;
    return inv;
}

LieMatrix* lie_inertia_to_matrix(const LieInertiaTensor *I) {
    assert(I);
    return lie_matrix_copy(I->inertia);
}

