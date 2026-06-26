#ifndef COORDINATE_TRANSFORM_H
#define COORDINATE_TRANSFORM_H

#include "nonlinear_system.h"

/* ============================================================================
 * Coordinate Transformations and Diffeomorphisms
 *
 * A diffeomorphism is a smooth, invertible map with smooth inverse.
 *
 *   Phi: R^n -> R^n,  z = Phi(x)
 *
 * such that det(dPhi/dx) != 0 (Jacobian is non-singular).
 *
 * In geometric control, coordinate transformations are used to:
 *   1. Transform nonlinear systems to normal forms
 *   2. Find linearizing coordinates
 *   3. Decompose systems into external/internal dynamics
 *
 * The transformation of a vector field f under diffeomorphism Phi:
 *   f_tilde(z) = [dPhi/dx ? f(x)] evaluated at x = Phi^{-1}(z)
 *
 * References:
 *   Isidori (1995) Nonlinear Control Systems, Appendix A
 *   Boothby (1986) An Introduction to Differentiable Manifolds
 *   Lee (2012) Introduction to Smooth Manifolds, 2nd ed.
 * ============================================================================ */

/* --- Diffeomorphism Type --- */

typedef struct {
    int n;                            /* Dimension */
    /* Forward map: z = Phi(x) */
    void (*forward)(const Vector *x, Vector *z);
    /* Inverse map: x = Phi^{-1}(z) */
    void (*inverse)(const Vector *z, Vector *x);
    /* Jacobian of forward map: J_Phi = dPhi/dx at x */
    void (*jacobian)(const Vector *x, Matrix *J);
    /* Jacobian of inverse: d(Phi^{-1})/dz at z */
    void (*inv_jacobian)(const Vector *z, Matrix *J_inv);
    /* Check if Phi is globally defined (vs local only) */
    bool is_global;
} Diffeomorphism;

/* --- Coordinate Transformation Results --- */

typedef struct {
    Diffeomorphism *transform;        /* The coordinate transformation */
    NonlinearSystem *original;        /* Original system */
    NonlinearSystem *transformed;     /* System in new coordinates */
    bool is_valid;                    /* Transform is a valid diffeomorphism */
} CoordinateTransformResult;

/* ============================================================================
 * Diffeomorphism API
 * ============================================================================ */

/**
 * diffeomorphism_create -- Allocate a diffeomorphism with given callbacks.
 *
 * @param n             dimension
 * @param forward       z = Phi(x)
 * @param inverse       x = Phi^{-1}(z)
 * @param jacobian      J = dPhi/dx
 * @param inv_jacobian  J^{-1} = dPhi^{-1}/dz
 * @return              allocated diffeomorphism
 */
Diffeomorphism *diffeomorphism_create(int n,
                                       void (*forward)(const Vector*, Vector*),
                                       void (*inverse)(const Vector*, Vector*),
                                       void (*jacobian)(const Vector*, Matrix*),
                                       void (*inv_jacobian)(const Vector*, Matrix*));

/**
 * diffeomorphism_apply_forward -- Apply forward map z = Phi(x).
 *
 * @param phi  diffeomorphism
 * @param x    original coordinates
 * @param z    transformed coordinates (output)
 */
void diffeomorphism_apply_forward(const Diffeomorphism *phi,
                                   const Vector *x,
                                   Vector *z);

/**
 * diffeomorphism_apply_inverse -- Apply inverse map x = Phi^{-1}(z).
 *
 * @param phi  diffeomorphism
 * @param z    transformed coordinates
 * @param x    original coordinates (output)
 */
void diffeomorphism_apply_inverse(const Diffeomorphism *phi,
                                   const Vector *z,
                                   Vector *x);

/**
 * diffeomorphism_check_jacobian -- Verify det(dPhi/dx) != 0 at point x.
 *
 * A diffeomorphism requires non-singular Jacobian everywhere.
 *
 * @param phi  diffeomorphism
 * @param x    check point
 * @param tol  numerical tolerance
 * @return     true if Jacobian is non-singular
 */
bool diffeomorphism_check_jacobian(const Diffeomorphism *phi,
                                    const Vector *x,
                                    double tol);

/**
 * diffeomorphism_validate -- Verify Phi ? Phi^{-1} = id within tolerance.
 *
 * Checks the round-trip property on a grid of test points.
 *
 * @param phi        diffeomorphism
 * @param test_points array of test vectors
 * @param num_points number of test points
 * @param tol        tolerance
 * @return           true if round-trip holds for all test points
 */
bool diffeomorphism_validate(const Diffeomorphism *phi,
                              Vector **test_points,
                              int num_points,
                              double tol);

/**
 * diffeomorphism_free -- Free diffeomorphism.
 *
 * Does not free the callback functions themselves (caller's responsibility).
 *
 * @param phi  diffeomorphism to free
 */
void diffeomorphism_free(Diffeomorphism *phi);

/* ============================================================================
 * System Transformation API
 * ============================================================================ */

/**
 * transform_system -- Apply diffeomorphism to a nonlinear system.
 *
 * Given system x' = f(x) + sum g_i(x) u_i, y = h(x)
 * and diffeomorphism z = Phi(x):
 *
 * The transformed system is:
 *   z' = f_tilde(z) + sum g_i_tilde(z) u_i
 *   y  = h_tilde(z)
 *
 * where:
 *   f_tilde(z) = (dPhi/dx) ? f(x) |_{x=Phi^{-1}(z)}
 *   g_i_tilde(z) = (dPhi/dx) ? g_i(x) |_{x=Phi^{-1}(z)}
 *   h_tilde(z) = h(x) |_{x=Phi^{-1}(z)}
 *
 * @param sys   original system
 * @param phi   diffeomorphism
 * @return      transformed system (caller must free)
 */
NonlinearSystem *transform_system(const NonlinearSystem *sys,
                                   const Diffeomorphism *phi);

/**
 * transform_vector_field -- Transform a single vector field under diffeomorphism.
 *
 * f_tilde(z) = (dPhi/dx) ? f(x) |_{x=Phi^{-1}(z)}
 *
 * Creates a new vector field with callbacks that evaluate the transformed
 * field using Phi and its Jacobian.
 *
 * @param f    original vector field
 * @param phi  diffeomorphism
 * @return     transformed vector field
 */
VectorField *transform_vector_field(const VectorField *f,
                                     const Diffeomorphism *phi);

/**
 * transform_scalar_field -- Transform a scalar field under diffeomorphism.
 *
 * h_tilde(z) = h(x) |_{x=Phi^{-1}(z)}
 *
 * @param h    original scalar field
 * @param phi  diffeomorphism
 * @return     transformed scalar field
 */
ScalarField *transform_scalar_field(const ScalarField *h,
                                     const Diffeomorphism *phi);

/**
 * identity_diffeomorphism -- Create identity map Phi(x) = x.
 *
 * @param n  dimension
 * @return   identity diffeomorphism
 */
Diffeomorphism *identity_diffeomorphism(int n);

/**
 * linear_diffeomorphism -- Create linear map Phi(x) = A x.
 *
 * Requires A to be non-singular. Inverse is x = A^{-1} z.
 *
 * @param A  non-singular n x n matrix
 * @return   linear diffeomorphism
 */
Diffeomorphism *linear_diffeomorphism(const Matrix *A);

/**
 * compose_diffeomorphisms -- Compose Phi_2 ? Phi_1: (Phi_2?Phi_1)(x) = Phi_2(Phi_1(x)).
 *
 * The composition of two diffeomorphisms is a diffeomorphism.
 *
 * @param phi1  inner diffeomorphism
 * @param phi2  outer diffeomorphism
 * @return      composed diffeomorphism (must have same dimension)
 */
Diffeomorphism *compose_diffeomorphisms(const Diffeomorphism *phi1,
                                         const Diffeomorphism *phi2);

/**
 * inverse_diffeomorphism -- Create the inverse diffeomorphism.
 *
 * Given Phi, returns Phi^{-1} by swapping forward and inverse callbacks.
 *
 * @param phi  diffeomorphism
 * @return     inverse diffeomorphism
 */
Diffeomorphism *inverse_diffeomorphism(const Diffeomorphism *phi);

#endif /* COORDINATE_TRANSFORM_H */
