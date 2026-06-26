#include "feedback_linearization.h"
#include "frobenius.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

#define DEFAULT_TOL 1e-8

RelativeDegreeResult compute_relative_degree(const NonlinearSystem *sys,
                                              const Vector *x,
                                              int max_order) {
    RelativeDegreeResult result;
    memset(&result, 0, sizeof(result));
    result.rho = -1;
    result.is_well_defined = false;
    if (!sys || !x || max_order <= 0) return result;
    if (sys->m != 1 || sys->p != 1) return result;
    VectorField *g = sys->g[0];
    VectorField *f = sys->f;
    ScalarField *h = sys->h[0];
    if (!f || !g || !h) return result;
    for (int k = 0; k < max_order; k++) {
        double LgLfk_h = lie_derivative_mixed(g, f, h, x, k);
        if (fabs(LgLfk_h) > DEFAULT_TOL) {
            result.rho = k + 1;
            result.is_well_defined = true;
            result.Lg_Lf_power = LgLfk_h;
            break;
        }
    }
    if (result.is_well_defined && result.rho > 0) {
        result.Lf_powers_h = compute_all_lie_derivatives(f, h, x, result.rho);
    }
    return result;
}


IOLinearizationResult input_output_linearize(const NonlinearSystem *sys,
                                              const Vector *x) {
    IOLinearizationResult result;
    memset(&result, 0, sizeof(result));
    result.is_feedback_linearizable = false;
    if (!sys || !x) return result;
    RelativeDegreeResult rd = compute_relative_degree(sys, x, sys->n);
    if (!rd.is_well_defined) return result;
    int rho = rd.rho;
    double Lf_rho_h = lie_derivative_kth(sys->f, sys->h[0], x, rho);
    double Lg_Lf_rho1_h = rd.Lg_Lf_power;
    if (fabs(Lg_Lf_rho1_h) < 1e-12) {
        if (rd.Lf_powers_h) free(rd.Lf_powers_h);
        return result;
    }
    result.is_feedback_linearizable = true;
    result.alpha = -Lf_rho_h / Lg_Lf_rho1_h;
    result.beta = 1.0 / Lg_Lf_rho1_h;
    if (rd.Lf_powers_h) free(rd.Lf_powers_h);
    return result;
}


ISLinearizationResult input_state_linearize(const NonlinearSystem *sys) {
    ISLinearizationResult result;
    memset(&result, 0, sizeof(result));
    result.is_exactly_linearizable = false;
    if (!sys || sys->m != 1) return result;
    int n = sys->n;
    VectorField *f = sys->f;
    VectorField *g = sys->g[0];
    if (!f || !g) return result;
    Vector *x0 = sys->equilibrium;
    Vector *temp_x = NULL;
    if (!x0) { temp_x = vector_create(n); vector_set(temp_x, 0.0); x0 = temp_x; }
    int rank = controllability_rank(f, g, x0, DEFAULT_TOL);
    if (temp_x) vector_free(temp_x);
    if (rank < n) { result.largest_controllable_index = rank - 1; return result; }
    x0 = sys->equilibrium; temp_x = NULL;
    if (!x0) { temp_x = vector_create(n); vector_set(temp_x, 0.0); x0 = temp_x; }
    int max_k = compute_max_involutive_subdist(f, g, x0);
    if (temp_x) vector_free(temp_x);
    result.largest_controllable_index = max_k;
    if (max_k < n - 1) return result;
    result.is_exactly_linearizable = true;
    result.phi = linearizing_coordinates(sys);
    result.f_transformed = vector_field_create(n, NULL, NULL);
    result.g_transformed = vector_field_create(n, NULL, NULL);
    return result;
}

ZeroDynamicsResult analyze_zero_dynamics(const NonlinearSystem *sys,
                                          const Vector *x) {
    ZeroDynamicsResult result;
    memset(&result, 0, sizeof(result));
    result.stability = ZERO_DYN_NONE;
    if (!sys || !x) return result;
    RelativeDegreeResult rd = compute_relative_degree(sys, x, sys->n);
    if (!rd.is_well_defined) return result;
    int rho = rd.rho;
    int n = sys->n;
    result.dim_zero_dynamics = n - rho;
    if (rd.Lf_powers_h) free(rd.Lf_powers_h);
    if (result.dim_zero_dynamics == 0) {
        result.stability = ZERO_DYN_NONE;
        return result;
    }
    NormalFormResult nf = compute_normal_form(sys, x);
    if (nf.eta && result.dim_zero_dynamics > 0) {
        result.zero_eigenvalues = (double*)malloc(
            (size_t)result.dim_zero_dynamics * sizeof(double));
        result.lyapunov_exponent = 0.0;
        if (result.zero_eigenvalues) {
            Matrix *J_q = matrix_create(
                result.dim_zero_dynamics, result.dim_zero_dynamics);
            double trace_jq = matrix_trace(J_q);
            double det_jq = matrix_determinant(J_q);
            if (result.dim_zero_dynamics == 1) {
                result.zero_eigenvalues[0] = trace_jq;
                result.lyapunov_exponent = trace_jq;
            } else if (result.dim_zero_dynamics == 2) {
                double disc = trace_jq * trace_jq - 4.0 * det_jq;
                double sqrt_disc = (disc >= 0) ? sqrt(disc) : 0.0;
                result.zero_eigenvalues[0] = (trace_jq + sqrt_disc) / 2.0;
                result.zero_eigenvalues[1] = (trace_jq - sqrt_disc) / 2.0;
                result.lyapunov_exponent = fmax(
                    result.zero_eigenvalues[0], result.zero_eigenvalues[1]);
            } else {
                double max_real = -DBL_MAX;
                for (int i = 0; i < result.dim_zero_dynamics; i++) {
                    double center = J_q->data[
                        i * result.dim_zero_dynamics + i];
                    double radius = 0.0;
                    for (int j = 0; j < result.dim_zero_dynamics; j++)
                        if (j != i) radius += fabs(J_q->data[
                            i * result.dim_zero_dynamics + j]);
                    if (center + radius > max_real)
                        max_real = center + radius;
                }
                result.lyapunov_exponent = max_real;
            }
            if (result.lyapunov_exponent < -DEFAULT_TOL)
                result.stability = ZERO_DYN_STABLE;
            else if (result.lyapunov_exponent > DEFAULT_TOL)
                result.stability = ZERO_DYN_UNSTABLE;
            else
                result.stability = ZERO_DYN_MARGINAL;
            matrix_free(J_q);
        }
    }
    result.zero_manifold = vector_create(rho);
    if (nf.xi) vector_free(nf.xi);
    if (nf.eta) vector_free(nf.eta);
    return result;
}


NormalFormResult compute_normal_form(const NonlinearSystem *sys,
                                      const Vector *x) {
    NormalFormResult result;
    memset(&result, 0, sizeof(result));
    if (!sys || !x) return result;
    RelativeDegreeResult rd = compute_relative_degree(sys, x, sys->n);
    if (!rd.is_well_defined) {
        if (rd.Lf_powers_h) free(rd.Lf_powers_h);
        return result;
    }
    int rho = rd.rho;
    int n = sys->n;
    result.xi = vector_create(rho);
    for (int k = 0; k < rho; k++)
        result.xi->data[k] = lie_derivative_kth(sys->f, sys->h[0], x, k);
    int dim_eta = n - rho;
    if (dim_eta > 0) {
        result.eta = vector_create(dim_eta);
        Vector *gx = vector_create(n);
        vector_field_eval(sys->g[0], x, gx);
        double g_norm_sq = vector_dot(gx, gx);
        if (g_norm_sq > 1e-12) {
            Matrix *W = matrix_create(n, n);
            int filled = 0;
            bool *used = (bool*)calloc((size_t)n, sizeof(bool));
            for (int i = 0; i < n && filled < dim_eta; i++) {
                Vector *w = vector_create(n);
                w->data[i] = 1.0;
                double dot_gw = vector_dot(gx, w);
                for (int j = 0; j < n; j++)
                    w->data[j] -= (dot_gw / g_norm_sq) * gx->data[j];
                double nrm = vector_norm(w);
                if (nrm > 1e-12 && !used[i]) {
                    vector_scale(w, 1.0 / nrm, w);
                    used[i] = true;
                    for (int j = 0; j < n; j++)
                        W->data[j * n + filled] = w->data[j];
                    filled++;
                }
                vector_free(w);
            }
            for (int i = 0; i < dim_eta; i++) {
                double eta_i = 0.0;
                for (int j = 0; j < n; j++)
                    eta_i += W->data[j * n + i] * x->data[j];
                result.eta->data[i] = eta_i;
            }
            matrix_free(W); free(used);
        } else {
            for (int i = 0; i < dim_eta && i < n; i++)
                result.eta->data[i] = x->data[i];
        }
        vector_free(gx);
    }
    result.h_normal = NULL;
    if (rd.Lf_powers_h) free(rd.Lf_powers_h);
    return result;
}

MIMOLinearizationResult mimo_linearize(const NonlinearSystem *sys,
                                        const Vector *x) {
    MIMOLinearizationResult result;
    memset(&result, 0, sizeof(result));
    if (!sys || !x) return result;
    int m = sys->m;
    result.num_inputs = m;
    result.num_outputs = sys->p;
    if (m != sys->p) {
        result.has_vector_relative_degree = false;
        return result;
    }
    result.relative_degrees = (int*)malloc((size_t)m * sizeof(int));
    if (!result.relative_degrees) return result;
    bool all_defined = true;
    int total_relative_degree = 0;
    for (int i = 0; i < m; i++) {
        int rho_i = -1;
        for (int r = 1; r <= sys->n; r++) {
            bool found = false;
            for (int j = 0; j < m; j++) {
                double Lgj_Lf_hi = lie_derivative_mixed(
                    sys->g[j], sys->f, sys->h[i], x, r - 1);
                if (fabs(Lgj_Lf_hi) > DEFAULT_TOL) {
                    found = true; break;
                }
            }
            if (found) { rho_i = r; break; }
        }
        if (rho_i < 0) { all_defined = false; break; }
        result.relative_degrees[i] = rho_i;
        total_relative_degree += rho_i;
    }
    result.has_vector_relative_degree = all_defined;
    if (all_defined) {
        result.decoupling_matrix = matrix_create(m, m);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                result.decoupling_matrix->data[i * m + j] =
                    lie_derivative_mixed(sys->g[j], sys->f, sys->h[i], x,
                                          result.relative_degrees[i] - 1);
        double det_E = matrix_determinant(result.decoupling_matrix);
        result.decoupling_matrix_invertible = (fabs(det_E) > DEFAULT_TOL);
        if (result.decoupling_matrix_invertible) {
            result.A = matrix_create(total_relative_degree,
                                      total_relative_degree);
            result.B = matrix_create(total_relative_degree, m);
            int offset = 0;
            for (int i = 0; i < m; i++) {
                int rho_i = result.relative_degrees[i];
                for (int k = 0; k < rho_i - 1; k++)
                    result.A->data[(offset + k) * total_relative_degree
                                    + (offset + k + 1)] = 1.0;
                result.B->data[(offset + rho_i - 1) * m + i] = 1.0;
                offset += rho_i;
            }
        }
    }
    return result;
}


double feedback_law_evaluate(const IOLinearizationResult *io_result,
                              double v) {
    if (!io_result || !io_result->is_feedback_linearizable) return 0.0;
    return io_result->alpha + io_result->beta * v;
}

bool check_minimum_phase(const ZeroDynamicsResult *zd_result) {
    if (!zd_result) return false;
    return (zd_result->stability == ZERO_DYN_STABLE) ||
           (zd_result->stability == ZERO_DYN_NONE);
}

bool full_state_feedback_linearizable(const NonlinearSystem *sys) {
    if (!sys || sys->m != 1) return false;
    ISLinearizationResult result = input_state_linearize(sys);
    bool linearizable = result.is_exactly_linearizable;
    if (result.phi) diffeomorphism_free(result.phi);
    if (result.f_transformed) vector_field_free(result.f_transformed);
    if (result.g_transformed) vector_field_free(result.g_transformed);
    return linearizable;
}

Matrix *compute_decoupling_matrix(const NonlinearSystem *sys,
                                   const Vector *x,
                                   const int *relative_degrees) {
    if (!sys || !x || !relative_degrees) return NULL;
    int m = sys->m;
    if (m != sys->p) return NULL;
    Matrix *E = matrix_create(m, m);
    if (!E) return NULL;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            E->data[i * m + j] = lie_derivative_mixed(
                sys->g[j], sys->f, sys->h[i], x,
                relative_degrees[i] - 1);
    return E;
}

Diffeomorphism *linearizing_coordinates(const NonlinearSystem *sys) {
    if (!sys) return NULL;
    int n = sys->n;
    Diffeomorphism *phi = identity_diffeomorphism(n);
    return phi;
}

IOLinearizationResult free_output_feedback_linearize(
        const NonlinearSystem *sys, int max_extend) {
    IOLinearizationResult result;
    memset(&result, 0, sizeof(result));
    result.is_feedback_linearizable = false;
    if (!sys || max_extend <= 0) return result;
    Vector *x_default = vector_create(sys->n);
    if (sys->equilibrium)
        memcpy(x_default->data, sys->equilibrium->data,
               (size_t)sys->n * sizeof(double));
    result = input_output_linearize(sys, x_default);
    vector_free(x_default);
    return result;
}
