/* ==============================================================
 * lie_reduction.c -- Lie-Poisson and Euler-Poincare Reduction
 *
 * Implements: Reduced Lagrangian, Euler-Poincare equations,
 * Lie-Poisson bracket, Rigid body on SO(3), Heavy top dynamics,
 * Noether conserved quantities, SE(3) semidirect product reduction.
 *
 * Key theorems:
 *   Euler-Poincare: d/dt(dl/dxi) = ad*_xi(dl/dxi)  (Marsden & Ratiu SS13.5)
 *   Lie-Poisson: dF/dt = {F, H} on g*  (Marsden & Ratiu SS13.7)
 *   Euler rigid body: I*Omega_dot = I*Omega x Omega + tau
 *   Heavy top: Pi_dot = Pi x Omega + mg*l*Gamma x chi
 *
 * References:
 *   Marsden & Ratiu (1999) Introduction to Mechanics and Symmetry
 *   Arnold (1989) Mathematical Methods of Classical Mechanics
 *   Holm, Marsden, Ratiu (1998) Phys. Rev. Lett.
 * ============================================================== */

#include "lie_reduction.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/* --- Reduced Lagrangian --- */
LieReducedLagrangian* lie_reduced_lagrangian_create(
    const LieInertiaTensor *I, const char *name) {
    assert(I);
    LieReducedLagrangian *l = (LieReducedLagrangian*)malloc(sizeof(LieReducedLagrangian));
    assert(l);
    l->name = name ? strdup(name) : strdup("unnamed");
    l->locked_inertia = (LieInertiaTensor*)I;
    l->potential_offset = 0.0;
    return l;
}
void lie_reduced_lagrangian_free(LieReducedLagrangian *l) {
    if (l) { free(l->name); free(l); }
}

double lie_reduced_lagrangian_eval(const LieReducedLagrangian *l, const LieVector *xi) {
    assert(l && xi);
    LieVector *Ixi = lie_inertia_apply(l->locked_inertia, xi);
    double val = 0.5 * lie_vector_dot(xi, Ixi);
    lie_vector_free(Ixi);
    return val;
}

LieVector* lie_reduced_lagrangian_deriv(const LieReducedLagrangian *l, const LieVector *xi) {
    /* delta l / delta xi = I * xi (the body angular momentum) */
    return lie_inertia_apply(l->locked_inertia, xi);
}

/* --- Euler-Poincare RHS --- */
void lie_euler_poincare_rhs(const LieReducedLagrangian *l,
                             const LieVector *xi,
                             const LieVector *external_force,
                             LieVector *dxi_dt) {
    assert(l && xi && dxi_dt);
    /* d/dt (delta l / delta xi) = ad*_xi (delta l / delta xi) + F_ext */
    LieVector *mu = lie_reduced_lagrangian_deriv(l, xi);  /* delta l / delta xi */
    LieVector *adstar_mu = lie_ad_star_so3(xi, mu);

    /* dxi/dt = I^{-1} (ad*_xi(mu) + F_ext) */
    LieInertiaTensor *Iinv = lie_inertia_inverse(l->locked_inertia);
    LieVector *rhs = lie_vector_add(adstar_mu, external_force);
    LieVector *result = lie_inertia_apply(Iinv, rhs);

    for (int i = 0; i < xi->size; i++) dxi_dt->data[i] = result->data[i];

    lie_vector_free(mu); lie_vector_free(adstar_mu);
    lie_vector_free(rhs); lie_vector_free(result);
    lie_inertia_free(Iinv);
}

/* --- Lie-Poisson Bracket for SO(3) --- */
double lie_lie_poisson_bracket_so3(const LieVector *mu,
    void (*grad_F)(const LieVector*, LieVector*),
    void (*grad_H)(const LieVector*, LieVector*),
    const LieVector *mu_current) {
    (void)mu; /* mu parameter passed for interface: {F, H} evaluated at mu */
    /* {F, H}(mu) = - mu . (nabla F x nabla H) */
    LieVector dF = {.size=3}; dF.data = (double[3]){0};
    LieVector dH = {.size=3}; dH.data = (double[3]){0};
    grad_F(mu_current, &dF);
    grad_H(mu_current, &dH);
    LieVector *cross = lie_vector_cross3(&dF, &dH);
    double bracket = -lie_vector_dot(mu_current, cross);
    lie_vector_free(cross);
    return bracket;
}

/* --- Rigid Body on SO(3) --- */
LieRigidBodyState* lie_rigid_body_create(double Ixx, double Iyy, double Izz,
                                          double Ixy, double Ixz, double Iyz) {
    LieRigidBodyState *rb = (LieRigidBodyState*)malloc(sizeof(LieRigidBodyState));
    assert(rb);
    rb->I_body = lie_inertia_create(Ixx, Iyy, Izz, Ixy, Ixz, Iyz);
    rb->R = lie_matrix_identity(3);
    rb->Pi = lie_vector_create(3);
    rb->Omega = lie_vector_create(3);
    rb->tau = lie_vector_create(3);
    return rb;
}
void lie_rigid_body_free(LieRigidBodyState *rb) {
    if (rb) {
        lie_inertia_free(rb->I_body); lie_matrix_free(rb->R);
        lie_vector_free(rb->Pi); lie_vector_free(rb->Omega);
        lie_vector_free(rb->tau); free(rb);
    }
}
void lie_rigid_body_reset(LieRigidBodyState *rb) {
    assert(rb);
    for (int i = 0; i < 9; i++) rb->R->data[i] = (i % 4 == 0) ? 1.0 : 0.0;
    for (int i = 0; i < 3; i++) {
        rb->Pi->data[i] = 0.0; rb->Omega->data[i] = 0.0; rb->tau->data[i] = 0.0;
    }
}

void lie_rigid_body_rhs(const LieRigidBodyState *rb,
                         LieVector *Pi_dot, LieMatrix *R_dot) {
    assert(rb && Pi_dot && R_dot);
    /* Pi_dot = Pi x Omega + tau */
    LieVector *Pi_x_Omega = lie_vector_cross3(rb->Pi, rb->Omega);
    for (int i = 0; i < 3; i++) Pi_dot->data[i] = Pi_x_Omega->data[i] + rb->tau->data[i];
    lie_vector_free(Pi_x_Omega);

    /* R_dot = R * Omega_hat */
    LieMatrix *Omega_hat = lie_so3_hat(rb->Omega);
    LieMatrix *ROmega = lie_matrix_mul(rb->R, Omega_hat);
    memcpy(R_dot->data, ROmega->data, 9 * sizeof(double));
    lie_matrix_free(Omega_hat); lie_matrix_free(ROmega);
}

double lie_rigid_body_kinetic_energy(const LieRigidBodyState *rb) {
    assert(rb);
    LieVector *IOmega = lie_inertia_apply(rb->I_body, rb->Omega);
    double T = 0.5 * lie_vector_dot(rb->Omega, IOmega);
    lie_vector_free(IOmega);
    return T;
}

/* --- Heavy Top --- */
LieHeavyTopState* lie_heavy_top_create(double Ixx, double Iyy, double Izz,
                                        double mass, double com_dist) {
    LieHeavyTopState *ht = (LieHeavyTopState*)malloc(sizeof(LieHeavyTopState));
    assert(ht);
    ht->I_body = lie_inertia_create(Ixx, Iyy, Izz, 0.0, 0.0, 0.0);
    ht->R = lie_matrix_identity(3);
    ht->Pi = lie_vector_create(3);
    ht->Omega = lie_vector_create(3);
    ht->Gamma = lie_vector_create(3);
    ht->chi = lie_vector_create(3);
    ht->mass = mass;
    ht->g_accel = 9.81;
    ht->com_distance = com_dist;
    /* Default: gravity direction is z, COM along z in body frame */
    ht->Gamma->data[2] = 1.0;
    ht->chi->data[2] = 1.0;
    return ht;
}
void lie_heavy_top_free(LieHeavyTopState *ht) {
    if (ht) {
        lie_inertia_free(ht->I_body); lie_matrix_free(ht->R);
        lie_vector_free(ht->Pi); lie_vector_free(ht->Omega);
        lie_vector_free(ht->Gamma); lie_vector_free(ht->chi);
        free(ht);
    }
}
void lie_heavy_top_reset(LieHeavyTopState *ht) {
    assert(ht);
    for (int i = 0; i < 9; i++) ht->R->data[i] = (i % 4 == 0) ? 1.0 : 0.0;
    for (int i = 0; i < 3; i++) {
        ht->Pi->data[i] = 0.0; ht->Omega->data[i] = 0.0;
    }
    ht->Gamma->data[0] = 0.0; ht->Gamma->data[1] = 0.0; ht->Gamma->data[2] = 1.0;
}

void lie_heavy_top_rhs(const LieHeavyTopState *ht,
                        LieVector *Pi_dot, LieVector *Gamma_dot, LieMatrix *R_dot) {
    assert(ht && Pi_dot && Gamma_dot && R_dot);
    /* Pi_dot = Pi x Omega + m*g*l * Gamma x chi */
    LieVector *Pi_x_Omega = lie_vector_cross3(ht->Pi, ht->Omega);
    LieVector *Gamma_x_chi = lie_vector_cross3(ht->Gamma, ht->chi);
    double mg = ht->mass * ht->g_accel * ht->com_distance;
    for (int i = 0; i < 3; i++)
        Pi_dot->data[i] = Pi_x_Omega->data[i] + mg * Gamma_x_chi->data[i];
    lie_vector_free(Pi_x_Omega); lie_vector_free(Gamma_x_chi);

    /* Gamma_dot = Gamma x Omega */
    LieVector *Gamma_x_Omega = lie_vector_cross3(ht->Gamma, ht->Omega);
    for (int i = 0; i < 3; i++) Gamma_dot->data[i] = Gamma_x_Omega->data[i];
    lie_vector_free(Gamma_x_Omega);

    /* R_dot = R * Omega_hat */
    LieMatrix *Oh = lie_so3_hat(ht->Omega);
    LieMatrix *ROh = lie_matrix_mul(ht->R, Oh);
    memcpy(R_dot->data, ROh->data, 9 * sizeof(double));
    lie_matrix_free(Oh); lie_matrix_free(ROh);
}

double lie_heavy_top_energy(const LieHeavyTopState *ht) {
    assert(ht);
    LieVector *IOmega = lie_inertia_apply(ht->I_body, ht->Omega);
    double T = 0.5 * lie_vector_dot(ht->Omega, IOmega);
    lie_vector_free(IOmega);
    double V = ht->mass * ht->g_accel * ht->com_distance * ht->Gamma->data[2];
    return T + V;
}

/* --- Noether Conserved Quantity --- */
double lie_noether_conserved_quantity(const LieReducedLagrangian *l, const LieVector *xi) {
    /* For a left-invariant Lagrangian: J(xi) = <delta l / delta xi, xi>
     * This is the Noether conserved quantity associated with left symmetry. */
    LieVector *mu = lie_reduced_lagrangian_deriv(l, xi);
    double J = lie_vector_dot(mu, xi);
    lie_vector_free(mu);
    return J;
}

/* --- SE(3) Momentum --- */
LieSE3Momentum* lie_se3_momentum_create(void) {
    LieSE3Momentum *mom = (LieSE3Momentum*)malloc(sizeof(LieSE3Momentum));
    assert(mom);
    mom->Pi = lie_vector_create(3);
    mom->P = lie_vector_create(3);
    return mom;
}
void lie_se3_momentum_free(LieSE3Momentum *mom) {
    if (mom) { lie_vector_free(mom->Pi); lie_vector_free(mom->P); free(mom); }
}

void lie_se3_lie_poisson_rhs(const LieSE3Momentum *mu,
                              const LieVector *body_velocity,
                              const LieVector *body_angular_velocity,
                              LieVector *dmu_dt) {
    assert(mu && body_velocity && body_angular_velocity && dmu_dt);
    /* SE(3) Lie-Poisson bracket for semidirect product:
     * dPi/dt = Pi x Omega + P x v
     * dP/dt = P x Omega */
    LieVector *PixOmega = lie_vector_cross3(mu->Pi, body_angular_velocity);
    LieVector *Pxv = lie_vector_cross3(mu->P, body_velocity);
    LieVector *PxOmega = lie_vector_cross3(mu->P, body_angular_velocity);

    for (int i = 0; i < 3; i++) {
        dmu_dt->data[i] = PxOmega->data[i];
        dmu_dt->data[i+3] = PixOmega->data[i] + Pxv->data[i];
    }

    lie_vector_free(PixOmega); lie_vector_free(Pxv); lie_vector_free(PxOmega);
}

