/**
 * @file flatness_core.h
 * @brief Core definitions for differential flatness theory
 *
 * Reference: Fliess, Levine, Martin, Rouchon (1995)
 *   "Flatness and defect of non-linear systems"
 *   International Journal of Control, 61(6), 1327-1361.
 *
 * L1: differentially flat system, flat outputs, relative degree
 * L2: exact linearization, trajectory generation without ODE
 * L3: vector fields, state-space, polynomials
 * L4: Brunovsky normal form, Fliess equivalence theorem
 * L6: n-trailer, quadrotor, crane, space robot
 */

#ifndef FLATNESS_CORE_H
#define FLATNESS_CORE_H

#include <stddef.h>
#include <stdbool.h>

#define FLAT_MAX_STATE_DIM      20
#define FLAT_MAX_INPUT_DIM      10
#define FLAT_MAX_OUTPUT_DIM     10
#define FLAT_MAX_DERIV_ORDER    8
#define FLAT_MAX_VECTOR_FIELDS  12
#define FLAT_MAX_TERMS          64
#define FLAT_MAX_TRAJ_POINTS    1024

/** Monomial: c * x1^p1 * ... * xn^pn */
typedef struct {
    double coefficient;
    int    powers[FLAT_MAX_STATE_DIM];
} FlatMonomial;

/** Multivariate polynomial */
typedef struct {
    FlatMonomial terms[FLAT_MAX_TERMS];
    int          num_terms;
    int          num_vars;
} FlatPolynomial;

/** Smooth vector field f: R^n -> R^n */
typedef struct {
    FlatPolynomial components[FLAT_MAX_STATE_DIM];
    int            dimension;
    char           name[32];
} FlatVectorField;

/** Control-affine system: xdot = f(x) + sum g_i(x)*u_i, y = h(x) */
typedef struct {
    FlatVectorField  drift;
    FlatVectorField  control_fields[FLAT_MAX_INPUT_DIM];
    FlatPolynomial   output_map[FLAT_MAX_OUTPUT_DIM];
    int              state_dim;
    int              input_dim;
    int              output_dim;
} FlatControlAffineSystem;

/** B-spline for smooth trajectory parameterization */
typedef struct {
    double *coefficients;
    double *knot_vector;
    int     num_coeffs;
    int     order;
    int     num_derivatives;
} FlatSpline;

/** Flat output trajectory over [t0, tf] */
typedef struct {
    FlatSpline    *y_splines;
    double         t0;
    double         tf;
    int            num_components;
    int            deriv_order;
} FlatTrajectory;

/** Computed state-input trajectory from flat parameterization */
typedef struct {
    double *state_trajectory;
    double *input_trajectory;
    double *time_points;
    int     num_points;
    int     state_dim;
    int     input_dim;
} FlatStateInputTrajectory;

/** Endogenous dynamic feedback compensator */
typedef struct {
    FlatVectorField alpha;
    FlatPolynomial  beta[FLAT_MAX_INPUT_DIM];
    int             comp_order;
    int             num_virtual_inputs;
} FlatDynamicCompensator;

/* ---- Core API: Flatness verification ---- */
bool flat_system_is_flat(const FlatControlAffineSystem *sys);
int  flat_compute_relative_degree(const FlatControlAffineSystem *sys, int *rel_deg);
int  flat_compute_defect(const FlatControlAffineSystem *sys);
bool flat_verify_flat_outputs(const FlatControlAffineSystem *sys,
                               const FlatPolynomial *candidate_y,
                               int num_candidates);
int  flat_design_endo_feedback(const FlatControlAffineSystem *sys,
                                const FlatPolynomial *flat_outputs,
                                FlatDynamicCompensator *compensator);
int  flat_to_brunovsky(const FlatControlAffineSystem *sys,
                       const FlatPolynomial *flat_outputs,
                       int *bruno_dims, double *transform_mat);

/* ---- Canonical flat system constructors ---- */
int flat_make_n_trailer_system(int n, FlatControlAffineSystem *sys);
int flat_make_quadrotor_system(FlatControlAffineSystem *sys);
int flat_make_crane_system(FlatControlAffineSystem *sys);
int flat_make_space_robot_system(int n, FlatControlAffineSystem *sys);

#endif /* FLATNESS_CORE_H */

/* ================================================================
 * L1-L6 Extended API: Additional flatness operations
 * ================================================================ */

/* ---- Vector field operations ---- */
void flat_vf_init(FlatVectorField*vf,int dim,const char*name);
int flat_vf_evaluate(const FlatVectorField*vf,const double*x,double*out);
int flat_vf_jacobian(const FlatVectorField*vf,const double*x,double*J);
int flat_vf_add(const FlatVectorField*a,const FlatVectorField*b,FlatVectorField*out);
int flat_vf_scale(const FlatVectorField*a,double c,FlatVectorField*out);
int flat_vf_subtract(const FlatVectorField*a,const FlatVectorField*b,FlatVectorField*out);
int flat_vf_copy(const FlatVectorField*src,FlatVectorField*dst);

/* ---- Polynomial operations ---- */
void flat_poly_init(FlatPolynomial*p,int n_vars);
void flat_poly_clear(FlatPolynomial*p);
int flat_poly_add_term(FlatPolynomial*p,double coeff,const int*powers);
double flat_poly_evaluate(const FlatPolynomial*p,const double*x);
double flat_poly_evaluate_deriv(const FlatPolynomial*p,const double*x,int var);
int flat_poly_add(const FlatPolynomial*a,const FlatPolynomial*b,FlatPolynomial*out);
int flat_poly_subtract(const FlatPolynomial*a,const FlatPolynomial*b,FlatPolynomial*out);
int flat_poly_multiply(const FlatPolynomial*a,const FlatPolynomial*b,FlatPolynomial*out);
int flat_poly_scale(const FlatPolynomial*a,double c,FlatPolynomial*out);
int flat_poly_partial_deriv(const FlatPolynomial*p,int var_idx,FlatPolynomial*out);
int flat_poly_gradient(const FlatPolynomial*p,FlatPolynomial*grad);
int flat_poly_degree(const FlatPolynomial*p);
int flat_poly_var_degree(const FlatPolynomial*p,int var_idx);
int flat_poly_copy(const FlatPolynomial*src,FlatPolynomial*dst);
bool flat_poly_is_zero(const FlatPolynomial*p);
bool flat_poly_is_constant(const FlatPolynomial*p);
bool flat_poly_is_linear(const FlatPolynomial*p);
void flat_poly_print(const FlatPolynomial*p,char*buf,int bufsize);

/* ---- System analysis ---- */
int flat_system_evaluate_dynamics(const FlatControlAffineSystem*sys,const double*x,const double*u,double*xdot);
int flat_system_jacobian(const FlatControlAffineSystem*sys,const double*x,const double*u,double*df_dx,double*df_du);
int flat_controllability_matrix(const double*A,const double*B,int n,int m,double*C);
int flat_matrix_rank(const double*M,int rows,int cols,double tol);
int flat_controllability_matrix_rank(const FlatControlAffineSystem*sys,int*rank);
bool flat_system_is_driftless(const FlatControlAffineSystem*sys);
bool flat_system_is_control_affine(const FlatControlAffineSystem*sys);
bool flat_system_is_square(const FlatControlAffineSystem*sys);
int flat_system_get_state_dim(const FlatControlAffineSystem*sys);
int flat_system_get_input_dim(const FlatControlAffineSystem*sys);
int flat_system_get_output_dim(const FlatControlAffineSystem*sys);
int flat_kalman_decomposition(double*A,double*B,int n,int m,int*nc);
int flat_pbh_test(const double*A,const double*B,int n,int m);
bool flat_linear_controllability_check(const FlatControlAffineSystem*sys);
int flat_poly_substitute(const FlatPolynomial*p,const FlatPolynomial*subs,FlatPolynomial*out);
int flat_poly_integrate(const FlatPolynomial*p,int var,FlatPolynomial*out);
double flat_poly_inner_product(const FlatPolynomial*a,const FlatPolynomial*b,const double*domain);
int flat_poly_reduce(FlatPolynomial*p,FlatPolynomial*basis,int n_basis);
bool flat_poly_is_sos(const FlatPolynomial*p);
int flat_poly_rk4_step(const FlatPolynomial*f,int n,const double*x,double dt,double*x_next);
int flat_find_equilibria(const FlatControlAffineSystem*sys,double*equilibria,int max_eq);
int flat_bifurcation_detect(const FlatControlAffineSystem*sys,double param,double*param_crit);

/* ---- L5: Observers for flat systems ---- */
int flat_luenberger_observer(const double*A,const double*B,const double*C,int n,int m,int p,double*L);
int flat_extended_kalman_filter(void(*f)(const double*,const double*,double*,void*),void(*h)(const double*,double*,void*),void*ctx,int n,int m,int p,double*x,double*P,const double*u,const double*z,double Q,double R,double dt);

/* ---- L4-L5: Additional matrix decompositions ---- */
int flat_qr_decomposition(double*A,int m,int n,double*Q,double*R);
int flat_svd_decomposition(const double*A,int m,int n,double*U,double*S,double*Vt);
double flat_condition_number(const double*A,int n);
int flat_linear_solve(const double*A,const double*b,int n,double*x);
int flat_eigenvalues_2x2(const double*A,double*lambda_re,double*lambda_im);
int flat_cholesky_decomposition(const double*A,int n,double*L);
int flat_ldlt_decomposition(const double*A,int n,double*L,double*D);
int flat_forward_substitution(const double*L,const double*b,int n,double*y);
int flat_backward_substitution(const double*U,const double*y,int n,double*x);
int flat_newton_solve(void(*F)(const double*,double*,void*),void(*JF)(const double*,double*,void*),void*ctx,int n,double*x,int max_iter,double tol);
int flat_cubic_spline_interp(const double*x,const double*y,int n,double*a,double*b,double*c,double*d);
int flat_pseudo_inverse(const double*A,int m,int n,double*A_pinv);
double flat_golden_section_search(double(*f)(double,void*),void*ctx,double a,double b,double tol);
int flat_bfgs_optimize(double(*f)(const double*,void*),void(*grad)(const double*,double*,void*),void*ctx,int n,double*x,double tol,int max_iter);

/* ---- L6: Additional canonical system models ---- */
int flat_make_inverted_pendulum(FlatControlAffineSystem*sys);
int flat_make_pvtol_system(FlatControlAffineSystem*sys);
int flat_make_ball_beam_system(FlatControlAffineSystem*sys);
int flat_make_induction_motor_system(FlatControlAffineSystem*sys);
int flat_make_maglev_system(FlatControlAffineSystem*sys);
int flat_make_boost_converter(FlatControlAffineSystem*sys);
int flat_make_auv_system(FlatControlAffineSystem*sys);
int flat_make_satellite_attitude(FlatControlAffineSystem*sys);
int flat_make_wind_turbine(FlatControlAffineSystem*sys);
int flat_make_flexible_joint(FlatControlAffineSystem*sys);

/* ---- L4-L8: Advanced analysis functions ---- */
int flat_check_absolute_flatness(const FlatControlAffineSystem*sys);
int flat_verify_jacobi_identity(const FlatVectorField*f,const FlatVectorField*g,const FlatVectorField*h,double*tol);

int flat_count_required_extensions(const FlatControlAffineSystem*sys);
int flat_compute_characteristic_indices(const FlatControlAffineSystem*sys,int*indices,int max_indices);
int flat_output_degree(const FlatControlAffineSystem*sys,int output_idx,int*degree);
int flat_is_statically_linearizable(const FlatControlAffineSystem*sys);
double flat_decoupling_matrix_det(const FlatControlAffineSystem*sys,const double*x);
int flat_verify_trajectory(const FlatControlAffineSystem*sys,const double*x,const double*u,int n_steps,double dt,double*tol);
double flat_trajectory_power(const FlatControlAffineSystem*sys,const FlatStateInputTrajectory*traj);
double flat_tracking_error(const double*ref,const double*actual,int n,int n_steps);
int flat_monte_carlo_flatness_check(const FlatControlAffineSystem*sys,int n_samples,double*prob);
double flat_energy_measure(const FlatControlAffineSystem*sys,const double*x);
int flat_differential_weight(const FlatControlAffineSystem*sys,double*weights);
int flat_structural_distance(const FlatControlAffineSystem*sys,double*distance);
int flat_controllability_gramian(const double*A,const double*B,int n,int m,double T,double*Wc);
int flat_observability_gramian(const double*A,const double*C,int n,int p,double T,double*Wo);
int flat_lyapunov_solve(const double*A,const double*Q,int n,double*X);
int flat_jordan_2x2(const double*A,double*P,double*J);
int flat_discrete_time_flatness(const double*Ad,const double*Bd,int n,int m);
int flat_zoh_discretize(const double*A,const double*B,int n,int m,double Ts,double*Ad,double*Bd);
int flat_bisimulation_check(const FlatControlAffineSystem*sys1,const FlatControlAffineSystem*sys2,double*sim);
int flat_morphology_class(const FlatControlAffineSystem*sys,int*class_id);
int flat_input_to_state_linearization(const FlatControlAffineSystem*sys,const double*x,double*T);

/* ---- L4: Lie algebra analysis ---- */
int flat_lie_algebra_dimension(const FlatVectorField*fields,int n_fields,int max_dim,int*dim);
int flat_lie_algebra_nilpotency_check(const FlatVectorField*fields,int n_fields,int max_depth);
int flat_lie_algebra_solvability_check(const FlatVectorField*fields,int n_fields);
int flat_generate_field_basis(int n,int idx,FlatVectorField*vf);
int flat_structure_constants(const FlatVectorField*fields,int n_fields,double*constants);
int flat_is_symmetry(const FlatVectorField*vf,const FlatVectorField*f);

/* ---- L5: Control design functions ---- */
int flat_lqr_design(const double*A,const double*B,const double*Q,const double*R,int n,int m,double*K);
int flat_acker_pole_placement(const double*A,const double*B,int n,const double*poles,double*K);
int flat_mpc_control(const FlatControlAffineSystem*sys,const double*x_current,const FlatTrajectory*ref,int horizon,double*U_plan);
int flat_feedforward_control(const FlatControlAffineSystem*sys,const FlatTrajectory*traj,double t,double*U_ff);
int flat_sliding_mode_control(const FlatControlAffineSystem*sys,const double*x,const double*xd,double lambda,double eta,double*U);
int flat_backstepping_control(const FlatControlAffineSystem*sys,const double*x,const double*xd,double*k,double*U);
int flat_adaptive_control(const FlatControlAffineSystem*sys,const double*x,const double*xd,double*theta_hat,double gamma,double*U);
int flat_gain_scheduling(const FlatControlAffineSystem*sys,const double*scheduling_var,int n_schedules,double*K_schedule);
int flat_hinf_control(const FlatControlAffineSystem*sys,const double*x,double gamma,double*U);
int flat_state_transform_linearization(const FlatControlAffineSystem*sys,const double*z,double*T,double*alpha,double*beta);
int flat_approx_feedback_linearization(const FlatControlAffineSystem*sys,const double*x_op,double*K_approx);
int flat_partial_feedback_linearization(const FlatControlAffineSystem*sys,int*linearizable_dims);
int flat_virtual_constraints(const FlatControlAffineSystem*sys,const FlatTrajectory*ref,double*constraint_params);
int flat_time_scale_separation(const FlatControlAffineSystem*sys,double epsilon,double*slow_manifold,double*fast_dynamics);
int flat_receding_horizon_control(const FlatControlAffineSystem*sys,const double*x0,const FlatTrajectory*ref,int horizon,double*U_seq);
int flat_invariant_set(const FlatControlAffineSystem*sys,const double*gains,int kappa,double*set_radius);
int flat_lpv_embedding(const FlatControlAffineSystem*sys,double*scheduling_params,int n_params,double*A_lpv,double*B_lpv);

/* ---- L4-L8: Additional advanced analysis functions ---- */
int flat_poly_resultant(const FlatPolynomial*a,const FlatPolynomial*b,double*result);
double flat_poly_discriminant(const FlatPolynomial*p);
int flat_poly_sturm_sequence(const FlatPolynomial*p,FlatPolynomial*seq,int*seq_len);
int flat_poly_sign_changes(const FlatPolynomial*p,int*num_changes);
int flat_poly_budan_fourier(const FlatPolynomial*p,double a,double b,int*num_roots);
int flat_poly_companion_matrix(const FlatPolynomial*p,double*C);
int flat_routh_hurwitz(const double*coeffs,int degree,int*is_stable);
int flat_kharitonov(const double*coeffs_lo,const double*coeffs_hi,int degree,int*is_stable);
int flat_bilinear_transform(const double*num,const double*den,int order,double Ts,double*num_z,double*den_z);
int flat_pade_approximation(double delay,int order,double*num,double*den);
int flat_youla_parameterization(const double*A,const double*B,const double*C,int n,int m,int p,double*Q_params);
double flat_gap_metric_robustness(const double*A_nom,const double*B_nom,const double*C_nom,const double*A_pert,const double*B_pert,const double*C_pert,int n,double*stability_margin);
int flat_iqc_analysis(const double*A,const double*B,const double*C,const double*D,int n,int m,int p,double*multiplier);
int flat_sos_decomposition(const FlatPolynomial*p,FlatPolynomial*sos_terms,int max_terms,int*num_terms);

/* ---- L5: Additional trajectory planning algorithms ---- */
int flat_lp_trajectory(const double*A,const double*b,const double*c,int m,int n,double*x_opt);
int flat_qp_solve(const double*Q,const double*c,const double*A,const double*b,int n,int m,double*x);
int flat_dp_trajectory(const double*cost_matrix,int n_states,int n_steps,double*cost_to_go,int*optimal_path);
int flat_astar_path(const double*grid,int width,int height,const double*start,const double*goal,double*path,int*path_len);
int flat_rrt_plan(const FlatControlAffineSystem*sys,const double*start,const double*goal,int max_nodes,double step_size,double*tree,int*num_nodes);
int flat_prm_plan(const FlatControlAffineSystem*sys,int n_samples,int k_neighbors,double*roadmap,int*num_edges);
int flat_smooth_trajectory(FlatTrajectory*traj,double alpha,int iterations);
int flat_time_optimal(const FlatControlAffineSystem*sys,const double*start,const double*goal,double vmax,double amax,FlatTrajectory*traj);
int flat_tube_trajectory(const FlatControlAffineSystem*sys,const FlatTrajectory*nominal,double tube_radius,FlatTrajectory*tube);
int flat_stochastic_trajectory(const FlatControlAffineSystem*sys,const double*start,const double*goal,int n_samples,double noise_std,FlatTrajectory*traj);


int flat_agriculture_waypoints(const double*field_boundary,int n_vertices,double row_spacing,double**coverage_waypoints,int*num_waypoints);
int flat_ship_dynamic_positioning(const double*surge_sway_yaw_setpoint,double env_force,double env_direction,FlatTrajectory*traj);
int flat_wind_turbine_pitch_control(double wind_speed,double rated_power,double rated_speed,double*pitch_trajectory,int n_steps);
int flat_lane_change_trajectory(double lane_width,double vehicle_speed,double maneuver_time,FlatTrajectory*traj);
int flat_rocket_landing_trajectory(const double*initial_pos,const double*initial_vel,const double*landing_pad,double g,double max_thrust,FlatTrajectory*traj);
int flat_debris_removal_path(const double*debris_field,int n_debris,const double*robot_base,int priority_order,FlatTrajectory*traj);

/* ---- L5: Bezier, Hermite, and S-curve trajectories ---- */
int flat_bezier_curve(const double*control_points,int n_cp,int dim,double*curve,int n_samples);
int flat_hermite_spline(const double*t,const double*p,const double*v,int n,double*coeffs);
int flat_minimum_time_trajectory(const double*start,const double*goal,double vmax,double amax,double*traj_times);
int flat_s_curve_trajectory(double t,double t1,double t2,double t3,double vmax,double amax,double jmax,double*pos,double*vel,double*acc);
int flat_formation_trajectory(int n_vehicles,const double*formation_offsets,const FlatTrajectory*leader_traj,FlatTrajectory*follower_trajs);
int flat_funnel_trajectory(const FlatControlAffineSystem*sys,const double*nominal_path,int n_wp,double funnel_radius,double*funnel_bounds);

/* ---- L4: Additional matrix and system analysis ---- */
int flat_controllability_gramian(const double*A,const double*B,int n,int m,double T,double*Wc);
int flat_observability_gramian(const double*A,const double*C,int n,int p,double T,double*Wo);
int flat_lyapunov_solve(const double*A,const double*Q,int n,double*X);
int flat_jordan_2x2(const double*A,double*P,double*J);
int flat_discrete_time_flatness(const double*Ad,const double*Bd,int n,int m);
int flat_zoh_discretize(const double*A,const double*B,int n,int m,double Ts,double*Ad,double*Bd);
double flat_h2_norm(const double*A,const double*B,const double*C,int n,int m,int p);
double flat_hinf_norm(const double*A,const double*B,const double*C,int n,int m,int p);
double flat_gap_metric(const double*A1,const double*B1,const double*C1,const double*A2,const double*B2,const double*C2,int n1,int n2);
double flat_nu_gap_metric(const double*A1,const double*B1,const double*C1,const double*A2,const double*B2,const double*C2,int n);
int flat_duality_map(const double*A,const double*B,const double*C,int n,int m,int p,double*A_dual,double*B_dual,double*C_dual);
int flat_controllability_staircase(double*A,double*B,int n,int m,double*Ac,double*Bc,int*ci);
int flat_minimal_realization(const double*A,const double*B,const double*C,int n,int m,int p,int*n_min);
int flat_transmission_zeros(const double*A,const double*B,const double*C,const double*D,int n,int m,int p,double*freq,int max_z,int*nz);
int flat_balanced_truncation(const double*A,const double*B,const double*C,int n,int m,int p,int r,double*Ar,double*Br,double*Cr);
