#ifndef LIE_GROUP_H
#define LIE_GROUP_H

#include <stdbool.h>
#include <stddef.h>

/* ============================================================================
 * Lie Group Theory for Principal Bundles
 *
 * Core Lie groups used in gauge theory and geometric control:
 *   U(1)  — circle group, abelian (electromagnetism / QED)
 *   SU(2) — special unitary group, non-abelian (weak interaction / Yang-Mills)
 *   SO(3) — rotation group (frame bundles / geometric control)
 *
 * References:
 *   Kobayashi & Nomizu, Foundations of Differential Geometry, Vol. I (1963)
 *   Nakahara, Geometry, Topology and Physics (2003)
 *   Hall, Lie Groups, Lie Algebras, and Representations (2015)
 * ============================================================================ */

/* --- Group type enumeration --- */
typedef enum {
    LIE_U1  = 1,
    LIE_SU2 = 2,
    LIE_SO3 = 3,
    LIE_R   = 4,
    LIE_RN  = 5
} LieGroupType;

/* --- U(1) group: angle representation in [0, 2π) --- */
typedef struct {
    double theta;          /* angle parameter */
} U1Element;

typedef struct {
    double a;              /* Lie algebra element: i*a, a ∈ R */
} U1Algebra;

/* --- SU(2) group: quaternion representation (a + bi + cj + dk) --- */
typedef struct {
    double a, b, c, d;     /* a^2 + b^2 + c^2 + d^2 = 1 */
} SU2Element;

typedef struct {
    double x, y, z;        /* su(2) ≅ R^3 via Pauli matrices */
} SU2Algebra;

/* --- SO(3) group: 3×3 rotation matrix --- */
typedef struct {
    double m[3][3];        /* orthogonal matrix, det = +1 */
} SO3Element;

typedef struct {
    double x, y, z;        /* so(3) ≅ R^3, cross-product matrix */
} SO3Algebra;

/* --- Generic Lie group element (tagged union) --- */
typedef struct {
    LieGroupType type;
    union {
        U1Element   u1;
        SU2Element  su2;
        SO3Element  so3;
    } elem;
} LieGroupElement;

/* --- Generic Lie algebra element (tagged union) --- */
typedef struct {
    LieGroupType type;
    union {
        U1Algebra  u1;
        SU2Algebra su2;
        SO3Algebra so3;
    } alg;
} LieAlgebraElement;

/* ============================================================================
 * U(1) Group Operations (abelian, L1-L2)
 * =========================================================================== */

/**
 * Construct U(1) element from angle.
 * Theorem source: Standard parametrization of S^1 ≅ U(1).
 * Complexity: O(1)
 */
U1Element u1_make(double theta);
U1Element u1_identity(void);
U1Element u1_multiply(U1Element g, U1Element h);
U1Element u1_inverse(U1Element g);
U1Element u1_power(U1Element g, int n);
bool      u1_is_identity(U1Element g);
double    u1_distance(U1Element g, U1Element h);
double    u1_trace(U1Element g);

/* U(1) Lie algebra */
U1Algebra u1_alg_make(double a);
U1Element u1_exp(U1Algebra X);
U1Algebra u1_log(U1Element g);
U1Algebra u1_alg_add(U1Algebra X, U1Algebra Y);
U1Algebra u1_alg_scale(double s, U1Algebra X);
double    u1_alg_norm(U1Algebra X);

/* ============================================================================
 * SU(2) Group Operations (non-abelian, L1-L2)
 * =========================================================================== */

/**
 * Construct SU(2) element from quaternion components.
 * Automatically normalizes to unit quaternion.
 * Theorem source: SU(2) ≅ S^3 via quaternion representation.
 * Complexity: O(1)
 */
SU2Element su2_make(double a, double b, double c, double d);
SU2Element su2_identity(void);
SU2Element su2_multiply(SU2Element g, SU2Element h);
SU2Element su2_inverse(SU2Element g);
SU2Algebra su2_adjoint_action(SU2Element g, SU2Algebra X);
bool      su2_is_identity(SU2Element g);
double    su2_trace(SU2Element g);
double    su2_distance(SU2Element g, SU2Element h);
void      su2_to_matrix(SU2Element g, double U[2][2]);
void      su2_from_euler_angles(double alpha, double beta, double gamma,
                                SU2Element *out);

/* SU(2) Lie algebra su(2) */
SU2Algebra su2_alg_make(double x, double y, double z);
SU2Element su2_exp(SU2Algebra X);
SU2Algebra su2_log(SU2Element g);
SU2Algebra su2_alg_add(SU2Algebra X, SU2Algebra Y);
SU2Algebra su2_alg_scale(double s, SU2Algebra X);
SU2Algebra su2_alg_commutator(SU2Algebra X, SU2Algebra Y);
double     su2_alg_norm(SU2Algebra X);
double     su2_alg_inner(SU2Algebra X, SU2Algebra Y);
void       su2_alg_to_matrix(SU2Algebra X, double m[2][2]);

/* ============================================================================
 * SO(3) Group Operations (L1-L2)
 * =========================================================================== */

/**
 * Construct SO(3) element.
 * Theorem source: SO(3) is the rotation group of R^3.
 * Complexity: O(1) for Rodrigues, O(27) for matrix mult.
 */
SO3Element so3_identity(void);
SO3Element so3_rodrigues(double axis_x, double axis_y, double axis_z,
                          double angle);
SO3Element so3_multiply(SO3Element R1, SO3Element R2);
SO3Element so3_inverse(SO3Element R);
SO3Element so3_from_euler_zyx(double roll, double pitch, double yaw);
SO3Element so3_from_axis_angle(double ax, double ay, double az, double angle);
void       so3_to_axis_angle(SO3Element R, double *ax, double *ay, double *az,
                             double *angle);
void       so3_apply(SO3Element R, double v[3], double out[3]);
double     so3_trace(SO3Element R);
double     so3_distance(SO3Element R1, SO3Element R2);

/* SO(3) Lie algebra so(3) */
SO3Algebra so3_alg_make(double x, double y, double z);
SO3Element so3_exp(SO3Algebra X);
SO3Algebra so3_log(SO3Element R);
SO3Algebra so3_alg_commutator(SO3Algebra X, SO3Algebra Y);
double     so3_alg_norm(SO3Algebra X);
void       so3_alg_to_matrix(SO3Algebra X, double m[3][3]);
SO3Algebra so3_alg_scale(double s, SO3Algebra X);

/* SU(2) <-> SO(3) double cover (L3) */
SO3Element su2_to_so3(SU2Element q);
SU2Element so3_to_su2(SO3Element R);
SO3Algebra so3_alg_from_vector(double v[3]);

/* ============================================================================
 * Generic Lie Group Element Operations (L2)
 * =========================================================================== */

LieGroupElement lie_identity(LieGroupType t);
LieGroupElement lie_multiply(LieGroupElement g, LieGroupElement h);
LieGroupElement lie_inverse(LieGroupElement g);
LieAlgebraElement lie_log(LieGroupElement g);
LieGroupElement lie_exp(LieAlgebraElement X);
LieAlgebraElement lie_alg_commutator(LieAlgebraElement X, LieAlgebraElement Y);
double lie_distance(LieGroupElement g, LieGroupElement h);
double lie_trace(LieGroupElement g);
int    lie_dimension(LieGroupType t);

#endif /* LIE_GROUP_H */
