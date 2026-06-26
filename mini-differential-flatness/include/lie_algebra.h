/**
 * @file lie_algebra.h
 * @brief Lie algebraic tools for differential flatness analysis
 *
 * L1: Lie derivative L_f h = grad(h)*f, Lie bracket [f,g] = dg*f - df*g
 * L3: Distribution Delta = span{f1,...,fk}, involutive distribution
 * L4: Frobenius theorem: integrable iff involutive
 *
 * Reference: Isidori (1995) Nonlinear Control Systems, Springer.
 * Reference: Nijmeijer & van der Schaft (1990) Nonlinear Dynamical Control Systems
 */

#ifndef LIE_ALGEBRA_H
#define LIE_ALGEBRA_H

#include "flatness_core.h"

/* Lie derivative: L_f h = grad(h)*f */
int lie_derivative(const FlatPolynomial *h, const FlatVectorField *f,
                   int x_dim, FlatPolynomial *out);

/* Lie derivative of vector field: L_f g */
int lie_derivative_vector_field(const FlatVectorField *g,
                                 const FlatVectorField *f, FlatVectorField *out);

/* k-th iterated Lie derivative: L_f^k h */
int lie_derivative_iterated(const FlatPolynomial *h, const FlatVectorField *f,
                             int x_dim, int k, FlatPolynomial *out);

/* Lie bracket: [f, g] = dg*f - df*g */
int lie_bracket(const FlatVectorField *f, const FlatVectorField *g,
                FlatVectorField *out);

/* Iterated Lie bracket: ad_f^k g */
int lie_bracket_iterated(const FlatVectorField *f, const FlatVectorField *g,
                          int k, FlatVectorField *out);

/* Lie bracket of indexed fields */
int lie_bracket_indexed(const FlatVectorField *fields, int n, int i, int j,
                         FlatVectorField *out);

/* Distribution structure */
typedef struct { FlatVectorField *basis; int dim; int max_dim; } FlatDistribution;

int  flat_distribution_init(FlatDistribution *dist, int max_dim);
void flat_distribution_free(FlatDistribution *dist);
int  flat_distribution_add(FlatDistribution *dist, const FlatVectorField *vf);
int  flat_distribution_contains(const FlatDistribution *dist, const FlatVectorField *vf);
int  flat_distribution_is_involutive(const FlatDistribution *dist);
int  flat_distribution_involutive_closure(FlatDistribution *dist);
int  flat_distribution_sequence(const FlatControlAffineSystem *sys,
                                 FlatDistribution *seq, int max_steps);
int  flat_check_accessibility_rank(const FlatControlAffineSystem *sys);

/* Differential 1-form */
typedef struct { FlatPolynomial coefficients[20]; int dim; } FlatOneForm;

int flat_exterior_derivative(const FlatPolynomial *h, int n, FlatOneForm *out);
int flat_wedge_product(const FlatOneForm *omega, const FlatOneForm *eta,
                       int n, double *out_coeffs);

#endif /* LIE_ALGEBRA_H */
