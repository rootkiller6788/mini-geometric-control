/**
 * lie_group_methods.h ? Lie Group Methods for Geometric Optimal Control
 *
 * References:
 *   Hall (2015) "Lie Groups, Lie Algebras, and Representations"
 *   Iserles, Munthe-Kaas, Norsett, Zanna (2000) "Lie-Group Methods"
 *   Marsden & Ratiu (1999) "Introduction to Mechanics and Symmetry"
 *
 * Covers:
 *   - Lie groups SO(3), SE(3) as matrix groups
 *   - Lie algebras so(3), se(3) with commutator bracket
 *   - Exponential and logarithm maps
 *   - Adjoint representations (Ad, ad)
 *   - Group actions and reduced equations
 *   - Euler-Poincare equations for rigid bodies
 *
 * Knowledge coverage: L1 Definitions, L3 Math Structures, L5 Algorithms
 */

#ifndef LIE_GROUP_METHODS_H
#define LIE_GROUP_METHODS_H

#include "geo_optctrl_core.h"

/* -----------------------------------------------------------------------------
 * L1: Lie algebra element (skew-symmetric matrices)
 * ----------------------------------------------------------------------------- */

/**
 * so(3) ? Lie algebra of SO(3), the rotation group.
 *
 * Elements: 3?3 skew-symmetric matrices parametrized by ? ? ?? via
 * the hat map ^: ?? ? so(3):
 *   ?? = [[0, -??, ??],
 *         [??, 0, -??],
 *         [-??, ??, 0]]
 *
 * The Lie bracket is the matrix commutator:
 *   [???, ???] = ?????? - ?????? = (?? ? ??)^
 *
 * This encodes the cross product in ?? as a Lie bracket.
 */

/** Map from ?? to so(3): hat map ? ? ?? */
void geo_so3_hat(const double omega[3], double hat[3][3]);

/** Map from so(3) to ??: vee map ?? ? ? */
void geo_so3_vee(const double hat[3][3], double omega[3]);

/** so(3) Lie bracket [X, Y] = XY - YX */
void geo_so3_bracket(const double X[3][3], const double Y[3][3],
                     double result[3][3]);

/**
 * se(3) ? Lie algebra of SE(3), rigid body motions.
 *
 * Elements: pairs (?, v) ? ?? ? ?? represented as 4?4 matrices:
 *   ?? = [[??, v],
 *         [0,  0]]
 *
 * Lie bracket: [(??,v?), (??,v?)] = (?????, ???v? - ???v?)
 */

/** Construct se(3) element from (omega, v) in ?? ? ?? */
void geo_se3_hat(const double omega[3], const double v[3],
                 double hat[4][4]);

/** Extract (omega, v) from se(3) matrix */
void geo_se3_vee(const double hat[4][4], double omega[3], double v[3]);

/** se(3) Lie bracket */
void geo_se3_bracket(const double X[4][4], const double Y[4][4],
                     double result[4][4]);

/* -----------------------------------------------------------------------------
 * L3: Exponential map exp: g ? G
 * ----------------------------------------------------------------------------- */

/**
 * SO(3) exponential map: Rodrigues formula.
 *
 * exp(??) = I? + (sin ?)/? ? ?? + (1-cos ?)/?? ? ???
 * where ? = ||?||.
 *
 * This maps angular velocity ? to the rotation matrix R = exp(??).
 * The exponential map is surjective but not injective:
 * exp(??) = exp(?'?) iff ? and ?' differ by 2?k in direction.
 *
 * Complexity: O(1) ? closed form via Rodrigues
 *
 * @param omega  angular velocity vector in ??
 * @param R      output rotation matrix (3?3, column-major for storage: R[i][j])
 */
void geo_so3_exp(const double omega[3], double R[3][3]);

/**
 * Matrix exponential for general n?n matrices via scaling-and-squaring.
 *
 * Uses the algorithm of Higham (2009):
 *   exp(A) = (exp(A/2^k))^{2^k}
 *
 * @param A    n?n matrix (row-major: A[i][j])
 * @param M    output exp(A) (n?n)
 * @param n    matrix dimension
 * @param tol  tolerance (typically 1e-14)
 */
void geo_matrix_exp(const double *A, double *M, int n, double tol);

/* -----------------------------------------------------------------------------
 * L3: Logarithm map log: G ? g (inverse of exp near identity)
 * ----------------------------------------------------------------------------- */

/**
 * SO(3) logarithm map: inverse of Rodrigues formula.
 *
 * For R = exp(??) with ||?|| < ?:
 *   ? = arccos((tr(R) - 1)/2)
 *   ?? = ?/(2 sin ?) ? (R - R?)
 *
 * @param R      rotation matrix (3?3)
 * @param omega  output angular velocity (valid when near identity)
 * @return       0 on success, -1 if near singularity (? close to ?)
 */
int geo_so3_log(const double R[3][3], double omega[3]);

/**
 * Matrix logarithm for general matrices (principal logarithm).
 *
 * Uses inverse scaling-and-squaring with Schur decomposition.
 * Valid when A has no eigenvalues on the negative real axis.
 *
 * @param A  input matrix (n?n)
 * @param L  output log(A) (n?n)
 * @param n  dimension
 * @return   0 on success, -1 if eigenvalues on branch cut
 */
int geo_matrix_log(const double *A, double *L, int n);

/* -----------------------------------------------------------------------------
 * L2: Adjoint representations
 * ----------------------------------------------------------------------------- */

/**
 * Adjoint map Ad_g: g ? g, the derivative of conjugation:
 *
 *   Ad_g(?) = d/dt|_{t=0} (g exp(t?) g^{-1}) = g ? g^{-1}
 *
 * For SO(3): Ad_R(??) = R ?? R? = (R?)^
 *
 * For SE(3): Ad_{(R,t)}(?,v) = (R?, -R??t + Rv)
 */
void geo_so3_Ad(const double R[3][3], const double eta[3][3],
                double result[3][3]);

/**
 * Adjoint map on the Lie algebra ad_X(Y) = [X, Y].
 *
 * ad: g ? gl(g), ad_X(Y) = [X, Y]
 *
 * This is the derivative of Ad at the identity:
 *   ad_X = d/dt|_{t=0} Ad_{exp(tX)}
 *
 * For so(3): ad_?(?) = ???? - ???? = (? ? ?)^
 *
 * For se(3): ad_{(?,v)}(?,u) = (???, ??u - ??v)
 */
void geo_so3_ad(const double X[3][3], const double Y[3][3],
                double result[3][3]);

/**
 * Jacobi identity for Lie bracket on so(3) verification:
 *   [X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0
 *
 * @param X,Y,Z  so(3) matrices
 * @param tol    tolerance for numerical check
 * @return       1 if holds, 0 otherwise
 */
int geo_so3_jacobi_identity_check(const double X[3][3],
                                  const double Y[3][3],
                                  const double Z[3][3], double tol);

/* -----------------------------------------------------------------------------
 * L5: Lie group integration methods
 * ----------------------------------------------------------------------------- */

/**
 * RK-MK (Runge-Kutta-Munthe-Kaas) integrator for ODEs on Lie groups.
 *
 * For ? = A??(A) where A ? G and ?(A) ? g, the RK-MK method
 * integrates in the Lie algebra and maps to the group via exp.
 *
 * This integrator automatically stays on the group manifold,
 * preserving group structure to machine precision.
 *
 * Reference: Munthe-Kaas (1999) "High order Runge-Kutta methods on manifolds"
 *
 * @param group_dim    Lie group dimension (3 for SO(3), 6 for SE(3))
 * @param omega        function computing Lie algebra element ?(A)
 * @param omega_ctx    context for omega
 * @param A0           initial group element (flattened)
 * @param t            time to integrate
 * @param n_steps      number of steps
 * @param A_out        output group element
 * @param use_exp      function pointer to exponential map
 */
typedef void (*GeoGroupExp)(const double *algebra, double *group);

void geo_rk_munthe_kaas(int group_dim,
                        void (*omega)(const double *A, void *ctx, double *omega_out),
                        void *omega_ctx,
                        const double *A0, double t, int n_steps,
                        double *A_out, GeoGroupExp use_exp);

/**
 * Commutator-free Lie group integrator (Crouch-Grossman type).
 *
 * Uses products of exponentials without commutator terms.
 * Order 2: A_{n+1} = A_n ? exp(h??(A_n + h/2?A_n??(A_n)))
 *
 * @param group_dim    Lie group dimension
 * @param rhs          right hand side ?: G ? g
 * @param rhs_ctx      context
 * @param A0           initial group element
 * @param h            step size
 * @param n_steps      number of steps
 * @param A_out        output
 * @param use_exp      exponential map
 */
void geo_crouch_grossman_2(int group_dim,
                           void (*rhs)(const double *A, void *ctx, double *omega),
                           void *rhs_ctx,
                           const double *A0, double h, int n_steps,
                           double *A_out, GeoGroupExp use_exp);

/* -----------------------------------------------------------------------------
 * L3: Euler-Poincare equations for rigid body
 * ----------------------------------------------------------------------------- */

/**
 * Euler-Poincare equations on so(3)* for a free rigid body:
 *
 *   I??? + ? ? (I??) = ?
 *
 * where:
 *   I = diag(I?, I?, I?)  moment of inertia tensor
 *   ? ? ??                body angular velocity
 *   ? ? ??                external torque
 *
 * This is the reduced form of Euler's rigid body equations after
 * symmetry reduction by SO(3) (i.e., working on the Lie algebra instead
 * of the full configuration space).
 *
 * Reference: Marsden & Ratiu (1999), Ch. 13
 *
 * @param I        principal moments of inertia [I?, I?, I?]
 * @param Omega    body angular velocity [??, ??, ??]
 * @param tau      external torque [??, ??, ??]
 * @param dOmega   output d?/dt [???, ???, ???]
 */
void geo_euler_poincare_rigid_body(const double I[3],
                                   const double Omega[3],
                                   const double tau[3],
                                   double dOmega[3]);

/**
 * Full rigid body kinematics on SO(3):
 *
 *   ? = R???
 *
 * where R ? SO(3) is the orientation matrix and ?? ? so(3).
 *
 * @param R      rotation matrix (3?3)
 * @param Omega  body angular velocity (3-vector)
 * @param dR     output ? (3?3)
 */
void geo_so3_kinematics(const double R[3][3], const double Omega[3],
                        double dR[3][3]);

/**
 * Cat map (Arnold's cat map) on the torus T? = ??/??.
 *
 * A classic example of a hyperbolic (chaotic) dynamical system on a
 * compact manifold. The cat map is given by:
 *   (x_{n+1}, y_{n+1}) = ((2x_n + y_n) mod 1, (x_n + y_n) mod 1)
 *
 * This is area-preserving (symplectic with respect to dx?dy) and
 * demonstrates Anosov diffeomorphism properties.
 *
 * Reference: Arnold & Avez (1968) "Ergodic Problems of Classical Mechanics"
 */
void geo_arnold_cat_map(double *x, double *y, int n_iter);

#endif /* LIE_GROUP_METHODS_H */
