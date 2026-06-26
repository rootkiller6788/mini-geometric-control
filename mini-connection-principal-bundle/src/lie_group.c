#include "lie_group.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * U(1) Group Operations (abelian)
 * U(1) is isomorphic to S^1 = {e^{i*theta} | theta in [0, 2*pi)}
 * Group law: e^{i*theta} * e^{i*phi} = e^{i*(theta+phi)}
 * Lie algebra: u(1) is isomorphic to R, exp(a) = e^{i*a}
 * Reference: Nakahara (2003), section 2.2
 * ================================================================== */

static double normalize_angle(double theta) {
    theta = fmod(theta, 2.0 * M_PI);
    if (theta < 0.0) theta += 2.0 * M_PI;
    return theta;
}

U1Element u1_make(double theta) {
    U1Element g;
    g.theta = normalize_angle(theta);
    return g;
}

U1Element u1_identity(void) {
    U1Element g;
    g.theta = 0.0;
    return g;
}

U1Element u1_multiply(U1Element g, U1Element h) {
    return u1_make(g.theta + h.theta);
}

U1Element u1_inverse(U1Element g) {
    return u1_make(-g.theta);
}

U1Element u1_power(U1Element g, int n) {
    return u1_make(g.theta * (double)n);
}

bool u1_is_identity(U1Element g) {
    double d = normalize_angle(g.theta);
    return (d < 1e-12 || d > 2.0 * M_PI - 1e-12);
}

double u1_distance(U1Element g, U1Element h) {
    double diff = normalize_angle(g.theta - h.theta);
    if (diff > M_PI) diff = 2.0 * M_PI - diff;
    return diff;
}

double u1_trace(U1Element g) {
    return 2.0 * cos(g.theta);
}

U1Algebra u1_alg_make(double a) {
    U1Algebra X;
    X.a = a;
    return X;
}

U1Element u1_exp(U1Algebra X) {
    return u1_make(X.a);
}

U1Algebra u1_log(U1Element g) {
    U1Algebra X;
    double theta = normalize_angle(g.theta);
    if (theta > M_PI) theta -= 2.0 * M_PI;
    X.a = theta;
    return X;
}

U1Algebra u1_alg_add(U1Algebra X, U1Algebra Y) {
    U1Algebra Z;
    Z.a = X.a + Y.a;
    return Z;
}

U1Algebra u1_alg_scale(double s, U1Algebra X) {
    U1Algebra Z;
    Z.a = s * X.a;
    return Z;
}

double u1_alg_norm(U1Algebra X) {
    return fabs(X.a);
}

/* ==================================================================
 * SU(2) Group Operations (non-abelian)
 * SU(2) = {a + bi + cj + dk | a^2+b^2+c^2+d^2=1} isomorphic to S^3
 * Lie algebra su(2) is isomorphic to R^3 via Pauli matrices
 * Reference: Hall (2015), section 1.2; Nakahara (2003), section 2.3
 * ================================================================== */

SU2Element su2_make(double a, double b, double c, double d) {
    SU2Element g;
    double norm = sqrt(a*a + b*b + c*c + d*d);
    if (norm < 1e-15) {
        g.a = 1.0; g.b = 0.0; g.c = 0.0; g.d = 0.0;
    } else {
        g.a = a / norm; g.b = b / norm; g.c = c / norm; g.d = d / norm;
    }
    return g;
}

SU2Element su2_identity(void) {
    SU2Element g;
    g.a = 1.0; g.b = 0.0; g.c = 0.0; g.d = 0.0;
    return g;
}

SU2Element su2_multiply(SU2Element g, SU2Element h) {
    SU2Element r;
    r.a = g.a*h.a - g.b*h.b - g.c*h.c - g.d*h.d;
    r.b = g.a*h.b + g.b*h.a + g.c*h.d - g.d*h.c;
    r.c = g.a*h.c - g.b*h.d + g.c*h.a + g.d*h.b;
    r.d = g.a*h.d + g.b*h.c - g.c*h.b + g.d*h.a;
    return r;
}

SU2Element su2_inverse(SU2Element g) {
    SU2Element r;
    r.a = g.a; r.b = -g.b; r.c = -g.c; r.d = -g.d;
    return r;
}

SU2Algebra su2_adjoint_action(SU2Element g, SU2Algebra X) {
    SU2Algebra Y;
    double a = g.a, b = g.b, c = g.c, d = g.d;
    double x = X.x, y = X.y, z = X.z;
    double p0 = -b*x - c*y - d*z;
    double p1 =  a*x + c*z - d*y;
    double p2 =  a*y - b*z + d*x;
    double p3 =  a*z + b*y - c*x;
    Y.x = p0*(-b) + p1*a + p2*d - p3*c;
    Y.y = p0*(-c) - p1*d + p2*a + p3*b;
    Y.z = p0*(-d) + p1*c - p2*b + p3*a;
    return Y;
}

bool su2_is_identity(SU2Element g) {
    return (fabs(g.a - 1.0) < 1e-12 && fabs(g.b) < 1e-12 &&
            fabs(g.c) < 1e-12 && fabs(g.d) < 1e-12);
}

double su2_trace(SU2Element g) {
    return 2.0 * g.a;
}

double su2_distance(SU2Element g, SU2Element h) {
    double dot = g.a*h.a + g.b*h.b + g.c*h.c + g.d*h.d;
    if (dot > 1.0) dot = 1.0;
    if (dot < -1.0) dot = -1.0;
    double d = acos(fabs(dot));
    return 2.0 * d;
}

void su2_to_matrix(SU2Element g, double U[2][2]) {
    U[0][0] = g.a;  U[0][1] = g.b;
    U[1][0] = g.c;  U[1][1] = g.d;
}

void su2_from_euler_angles(double alpha, double beta, double gamma,
                            SU2Element *out) {
    double hb = beta * 0.5;
    double hs = (alpha + gamma) * 0.5;
    double hd = (alpha - gamma) * 0.5;
    out->a = cos(hb) * cos(hs);
    out->b = -sin(hb) * cos(hd);
    out->c = -sin(hb) * sin(hd);
    out->d = -cos(hb) * sin(hs);
}

SU2Algebra su2_alg_make(double x, double y, double z) {
    SU2Algebra X;
    X.x = x; X.y = y; X.z = z;
    return X;
}

SU2Element su2_exp(SU2Algebra X) {
    double norm = sqrt(X.x*X.x + X.y*X.y + X.z*X.z);
    SU2Element g;
    if (norm < 1e-15) {
        g.a = 1.0; g.b = 0.0; g.c = 0.0; g.d = 0.0;
        return g;
    }
    double s = sin(norm) / norm;
    g.a = cos(norm);
    g.b = s * X.x; g.c = s * X.y; g.d = s * X.z;
    return g;
}

SU2Algebra su2_log(SU2Element g) {
    SU2Algebra X;
    double a = g.a;
    if (a > 1.0) a = 1.0;
    if (a < -1.0) a = -1.0;
    double theta = acos(a);
    double scale = (theta < 1e-15) ? 1.0 : theta / sin(theta);
    X.x = scale * g.b; X.y = scale * g.c; X.z = scale * g.d;
    if (fabs(a + 1.0) < 1e-12 && fabs(g.b) < 1e-12 &&
        fabs(g.c) < 1e-12 && fabs(g.d) < 1e-12) {
        X.x = M_PI; X.y = 0.0; X.z = 0.0;
    }
    return X;
}

SU2Algebra su2_alg_add(SU2Algebra X, SU2Algebra Y) {
    SU2Algebra Z;
    Z.x = X.x + Y.x; Z.y = X.y + Y.y; Z.z = X.z + Y.z;
    return Z;
}

SU2Algebra su2_alg_scale(double s, SU2Algebra X) {
    SU2Algebra Z;
    Z.x = s * X.x; Z.y = s * X.y; Z.z = s * X.z;
    return Z;
}

SU2Algebra su2_alg_commutator(SU2Algebra X, SU2Algebra Y) {
    SU2Algebra Z;
    Z.x = X.y * Y.z - X.z * Y.y;
    Z.y = X.z * Y.x - X.x * Y.z;
    Z.z = X.x * Y.y - X.y * Y.x;
    return Z;
}

double su2_alg_norm(SU2Algebra X) {
    return sqrt(X.x*X.x + X.y*X.y + X.z*X.z);
}

double su2_alg_inner(SU2Algebra X, SU2Algebra Y) {
    return X.x*Y.x + X.y*Y.y + X.z*Y.z;
}

void su2_alg_to_matrix(SU2Algebra X, double m[2][2]) {
    m[0][0] = X.z;   m[0][1] = X.x + X.y;
    m[1][0] = X.x - X.y;  m[1][1] = -X.z;
}

/* ==================================================================
 * SO(3) Group Operations
 * SO(3) = {R in M_{3x3} | R^T R = I, det(R) = 1}
 * Lie algebra so(3) = skew-symmetric 3x3 matrices, isomorphic to R^3
 * Rodrigues formula: exp(theta*n_hat) = I + sin(theta)*n_hat + (1-cos(theta))*n_hat^2
 * Reference: Marsden & Ratiu (1999), section 9.2
 * ================================================================== */

SO3Element so3_identity(void) {
    SO3Element R;
    R.m[0][0] = 1.0; R.m[0][1] = 0.0; R.m[0][2] = 0.0;
    R.m[1][0] = 0.0; R.m[1][1] = 1.0; R.m[1][2] = 0.0;
    R.m[2][0] = 0.0; R.m[2][1] = 0.0; R.m[2][2] = 1.0;
    return R;
}

SO3Element so3_rodrigues(double axis_x, double axis_y, double axis_z,
                          double angle) {
    double n = sqrt(axis_x*axis_x + axis_y*axis_y + axis_z*axis_z);
    if (n < 1e-15) return so3_identity();
    double kx = axis_x / n, ky = axis_y / n, kz = axis_z / n;
    double ct = cos(angle), st = sin(angle), vt = 1.0 - ct;
    SO3Element R;
    R.m[0][0] = ct + kx*kx*vt;
    R.m[0][1] = kx*ky*vt - kz*st;
    R.m[0][2] = kx*kz*vt + ky*st;
    R.m[1][0] = ky*kx*vt + kz*st;
    R.m[1][1] = ct + ky*ky*vt;
    R.m[1][2] = ky*kz*vt - kx*st;
    R.m[2][0] = kz*kx*vt - ky*st;
    R.m[2][1] = kz*ky*vt + kx*st;
    R.m[2][2] = ct + kz*kz*vt;
    return R;
}

SO3Element so3_multiply(SO3Element R1, SO3Element R2) {
    SO3Element R;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) sum += R1.m[i][k] * R2.m[k][j];
            R.m[i][j] = sum;
        }
    }
    return R;
}

SO3Element so3_inverse(SO3Element R) {
    SO3Element Rt;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Rt.m[i][j] = R.m[j][i];
    return Rt;
}

SO3Element so3_from_euler_zyx(double roll, double pitch, double yaw) {
    double cr = cos(roll),  sr = sin(roll);
    double cp = cos(pitch), sp = sin(pitch);
    double cy = cos(yaw),   sy = sin(yaw);
    SO3Element R;
    R.m[0][0] = cy * cp;
    R.m[0][1] = cy * sp * sr - sy * cr;
    R.m[0][2] = cy * sp * cr + sy * sr;
    R.m[1][0] = sy * cp;
    R.m[1][1] = sy * sp * sr + cy * cr;
    R.m[1][2] = sy * sp * cr - cy * sr;
    R.m[2][0] = -sp;
    R.m[2][1] = cp * sr;
    R.m[2][2] = cp * cr;
    return R;
}

SO3Element so3_from_axis_angle(double ax, double ay, double az, double angle) {
    return so3_rodrigues(ax, ay, az, angle);
}

void so3_to_axis_angle(SO3Element R, double *ax, double *ay, double *az,
                       double *angle) {
    double trace = R.m[0][0] + R.m[1][1] + R.m[2][2];
    double ct = (trace - 1.0) * 0.5;
    if (ct > 1.0) ct = 1.0;
    if (ct < -1.0) ct = -1.0;
    *angle = acos(ct);
    if (*angle < 1e-12) {
        *ax = 1.0; *ay = 0.0; *az = 0.0; *angle = 0.0;
    } else if (fabs(*angle - M_PI) < 1e-12) {
        double x = sqrt(fmax(0.0, (R.m[0][0] + 1.0) * 0.5));
        double y = sqrt(fmax(0.0, (R.m[1][1] + 1.0) * 0.5));
        double z = sqrt(fmax(0.0, (R.m[2][2] + 1.0) * 0.5));
        *ax = (R.m[2][1] > 0) ? x : -x;
        *ay = (R.m[0][2] > 0) ? y : -y;
        *az = (R.m[1][0] > 0) ? z : -z;
        double nrm = sqrt((*ax)*(*ax) + (*ay)*(*ay) + (*az)*(*az));
        if (nrm > 1e-12) { *ax /= nrm; *ay /= nrm; *az /= nrm; }
    } else {
        double s = 0.5 / sin(*angle);
        *ax = (R.m[2][1] - R.m[1][2]) * s;
        *ay = (R.m[0][2] - R.m[2][0]) * s;
        *az = (R.m[1][0] - R.m[0][1]) * s;
    }
}

void so3_apply(SO3Element R, double v[3], double out[3]) {
    out[0] = R.m[0][0]*v[0] + R.m[0][1]*v[1] + R.m[0][2]*v[2];
    out[1] = R.m[1][0]*v[0] + R.m[1][1]*v[1] + R.m[1][2]*v[2];
    out[2] = R.m[2][0]*v[0] + R.m[2][1]*v[1] + R.m[2][2]*v[2];
}

double so3_trace(SO3Element R) {
    return R.m[0][0] + R.m[1][1] + R.m[2][2];
}

double so3_distance(SO3Element R1, SO3Element R2) {
    SO3Element Rdiff = so3_multiply(so3_inverse(R1), R2);
    double ax, ay, az, angle;
    so3_to_axis_angle(Rdiff, &ax, &ay, &az, &angle);
    return angle;
}

SO3Algebra so3_alg_make(double x, double y, double z) {
    SO3Algebra X;
    X.x = x; X.y = y; X.z = z;
    return X;
}

SO3Element so3_exp(SO3Algebra X) {
    double norm = sqrt(X.x*X.x + X.y*X.y + X.z*X.z);
    if (norm < 1e-15) return so3_identity();
    return so3_rodrigues(X.x, X.y, X.z, norm);
}

SO3Algebra so3_log(SO3Element R) {
    SO3Algebra X;
    double ax, ay, az, angle;
    so3_to_axis_angle(R, &ax, &ay, &az, &angle);
    X.x = angle * ax; X.y = angle * ay; X.z = angle * az;
    return X;
}

SO3Algebra so3_alg_commutator(SO3Algebra X, SO3Algebra Y) {
    SO3Algebra Z;
    Z.x = X.y * Y.z - X.z * Y.y;
    Z.y = X.z * Y.x - X.x * Y.z;
    Z.z = X.x * Y.y - X.y * Y.x;
    return Z;
}

SO3Algebra so3_alg_scale(double s, SO3Algebra X) {
    SO3Algebra Z;
    Z.x = s * X.x; Z.y = s * X.y; Z.z = s * X.z;
    return Z;
}

double so3_alg_norm(SO3Algebra X) {
    return sqrt(X.x*X.x + X.y*X.y + X.z*X.z);
}

void so3_alg_to_matrix(SO3Algebra X, double m[3][3]) {
    m[0][0] = 0.0;   m[0][1] = -X.z; m[0][2] = X.y;
    m[1][0] = X.z;   m[1][1] = 0.0;  m[1][2] = -X.x;
    m[2][0] = -X.y;  m[2][1] = X.x;  m[2][2] = 0.0;
}

SO3Algebra so3_alg_from_vector(double v[3]) {
    SO3Algebra X;
    X.x = v[0]; X.y = v[1]; X.z = v[2];
    return X;
}

/* ==================================================================
 * SU(2) <-> SO(3) Double Cover
 *
 * The map Phi: SU(2) -> SO(3) is a 2:1 covering homomorphism.
 * For unit quaternion q = (w, x, y, z):
 *   R = [[1-2(y^2+z^2), 2(xy-wz),    2(xz+wy)    ],
 *        [2(xy+wz),     1-2(x^2+z^2), 2(yz-wx)    ],
 *        [2(xz-wy),     2(yz+wx),     1-2(x^2+y^2)]]
 * Reference: Hall (2015), Proposition 1.11
 * ================================================================== */

SO3Element su2_to_so3(SU2Element q) {
    double w = q.a, x = q.b, y = q.c, z = q.d;
    SO3Element R;
    R.m[0][0] = 1.0 - 2.0*(y*y + z*z);
    R.m[0][1] = 2.0*(x*y - w*z);
    R.m[0][2] = 2.0*(x*z + w*y);
    R.m[1][0] = 2.0*(x*y + w*z);
    R.m[1][1] = 1.0 - 2.0*(x*x + z*z);
    R.m[1][2] = 2.0*(y*z - w*x);
    R.m[2][0] = 2.0*(x*z - w*y);
    R.m[2][1] = 2.0*(y*z + w*x);
    R.m[2][2] = 1.0 - 2.0*(x*x + y*y);
    return R;
}

SU2Element so3_to_su2(SO3Element R) {
    double trace = R.m[0][0] + R.m[1][1] + R.m[2][2];
    SU2Element q;
    if (trace > 0.0) {
        double s = sqrt(trace + 1.0) * 2.0;
        q.a = s * 0.25;
        q.b = (R.m[2][1] - R.m[1][2]) / s;
        q.c = (R.m[0][2] - R.m[2][0]) / s;
        q.d = (R.m[1][0] - R.m[0][1]) / s;
    } else if (R.m[0][0] > R.m[1][1] && R.m[0][0] > R.m[2][2]) {
        double s = sqrt(1.0 + R.m[0][0] - R.m[1][1] - R.m[2][2]) * 2.0;
        q.a = (R.m[2][1] - R.m[1][2]) / s;
        q.b = s * 0.25;
        q.c = (R.m[0][1] + R.m[1][0]) / s;
        q.d = (R.m[0][2] + R.m[2][0]) / s;
    } else if (R.m[1][1] > R.m[2][2]) {
        double s = sqrt(1.0 + R.m[1][1] - R.m[0][0] - R.m[2][2]) * 2.0;
        q.a = (R.m[0][2] - R.m[2][0]) / s;
        q.b = (R.m[0][1] + R.m[1][0]) / s;
        q.c = s * 0.25;
        q.d = (R.m[1][2] + R.m[2][1]) / s;
    } else {
        double s = sqrt(1.0 + R.m[2][2] - R.m[0][0] - R.m[1][1]) * 2.0;
        q.a = (R.m[1][0] - R.m[0][1]) / s;
        q.b = (R.m[0][2] + R.m[2][0]) / s;
        q.c = (R.m[1][2] + R.m[2][1]) / s;
        q.d = s * 0.25;
    }
    if (q.a < 0.0) {
        q.a = -q.a; q.b = -q.b; q.c = -q.c; q.d = -q.d;
    }
    return q;
}

/* ==================================================================
 * Generic Lie Group Element Operations (type-erased dispatch)
 * ================================================================== */

LieGroupElement lie_identity(LieGroupType t) {
    LieGroupElement e;
    e.type = t;
    switch (t) {
        case LIE_U1:  e.elem.u1  = u1_identity();  break;
        case LIE_SU2: e.elem.su2 = su2_identity(); break;
        case LIE_SO3: e.elem.so3 = so3_identity(); break;
        default:
            e.type = LIE_U1;
            e.elem.u1 = u1_identity();
            break;
    }
    return e;
}

LieGroupElement lie_multiply(LieGroupElement g, LieGroupElement h) {
    if (g.type != h.type) return g;
    LieGroupElement r;
    r.type = g.type;
    switch (g.type) {
        case LIE_U1:
            r.elem.u1 = u1_multiply(g.elem.u1, h.elem.u1);
            break;
        case LIE_SU2:
            r.elem.su2 = su2_multiply(g.elem.su2, h.elem.su2);
            break;
        case LIE_SO3:
            r.elem.so3 = so3_multiply(g.elem.so3, h.elem.so3);
            break;
        default:
            r = g;
            break;
    }
    return r;
}

LieGroupElement lie_inverse(LieGroupElement g) {
    LieGroupElement r;
    r.type = g.type;
    switch (g.type) {
        case LIE_U1:
            r.elem.u1 = u1_inverse(g.elem.u1);
            break;
        case LIE_SU2:
            r.elem.su2 = su2_inverse(g.elem.su2);
            break;
        case LIE_SO3:
            r.elem.so3 = so3_inverse(g.elem.so3);
            break;
        default:
            r = g;
            break;
    }
    return r;
}

LieAlgebraElement lie_log(LieGroupElement g) {
    LieAlgebraElement X;
    X.type = g.type;
    switch (g.type) {
        case LIE_U1:
            X.alg.u1 = u1_log(g.elem.u1);
            break;
        case LIE_SU2:
            X.alg.su2 = su2_log(g.elem.su2);
            break;
        case LIE_SO3:
            X.alg.so3 = so3_log(g.elem.so3);
            break;
        default:
            X.type = LIE_U1;
            X.alg.u1.a = 0.0;
            break;
    }
    return X;
}

LieGroupElement lie_exp(LieAlgebraElement X) {
    LieGroupElement g;
    g.type = X.type;
    switch (X.type) {
        case LIE_U1:
            g.elem.u1 = u1_exp(X.alg.u1);
            break;
        case LIE_SU2:
            g.elem.su2 = su2_exp(X.alg.su2);
            break;
        case LIE_SO3:
            g.elem.so3 = so3_exp(X.alg.so3);
            break;
        default:
            g.type = LIE_U1;
            g.elem.u1 = u1_identity();
            break;
    }
    return g;
}

LieAlgebraElement lie_alg_commutator(LieAlgebraElement X, LieAlgebraElement Y) {
    if (X.type != Y.type) return X;
    LieAlgebraElement Z;
    Z.type = X.type;
    switch (X.type) {
        case LIE_U1:
            Z.alg.u1.a = 0.0;  /* U(1) is abelian */
            break;
        case LIE_SU2:
            Z.alg.su2 = su2_alg_commutator(X.alg.su2, Y.alg.su2);
            break;
        case LIE_SO3:
            Z.alg.so3 = so3_alg_commutator(X.alg.so3, Y.alg.so3);
            break;
        default:
            Z.alg.u1.a = 0.0;
            break;
    }
    return Z;
}

double lie_distance(LieGroupElement g, LieGroupElement h) {
    if (g.type != h.type) return INFINITY;
    switch (g.type) {
        case LIE_U1:  return u1_distance(g.elem.u1, h.elem.u1);
        case LIE_SU2: return su2_distance(g.elem.su2, h.elem.su2);
        case LIE_SO3: return so3_distance(g.elem.so3, h.elem.so3);
        default: return 0.0;
    }
}

double lie_trace(LieGroupElement g) {
    switch (g.type) {
        case LIE_U1:  return u1_trace(g.elem.u1);
        case LIE_SU2: return su2_trace(g.elem.su2);
        case LIE_SO3: return so3_trace(g.elem.so3);
        default: return 0.0;
    }
}

int lie_dimension(LieGroupType t) {
    switch (t) {
        case LIE_U1:  return 1;
        case LIE_SU2: return 3;
        case LIE_SO3: return 3;
        case LIE_R:   return 1;
        case LIE_RN:  return -1;
        default: return 0;
    }
}
