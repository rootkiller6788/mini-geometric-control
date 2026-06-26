#ifndef TRAJECTORY_GENERATION_H
#define TRAJECTORY_GENERATION_H
#include "flatness_core.h"

/* L5: B-spline basis */
double flat_bspline_basis(int i,int p,const double*knots,double t);
int flat_bspline_basis_all(double t,int p,const double*knots,int n,double*basis_vals,int*first_nonzero);
double flat_bspline_eval(const FlatSpline*spline,double t,int deriv);
int flat_bspline_derivatives(const FlatSpline*spline,double t,int max_deriv,double*vals);

/* L5: Polynomial trajectory */
double flat_poly_eval_traj(const double*coeffs,int degree,double t);
int flat_poly_derivatives_traj(const double*coeffs,int degree,double t,int max_deriv,double*vals);
int flat_poly_fit_waypoints(const double*waypoints,int n_waypoints,const double*times,int n_segments,int continuity_order,double**coeffs_out);

/* L2: Flat output to state/input mapping */
int flat_map_to_state_input(const FlatControlAffineSystem*sys,const FlatTrajectory*flat_y,int n_samples,const double*times,FlatStateInputTrajectory*result);

/* L5: Trajectory optimization */
typedef enum{FLAT_COST_MINIMUM_SNAP,FLAT_COST_MINIMUM_JERK,FLAT_COST_MINIMUM_ACCEL,FLAT_COST_MINIMUM_ENERGY,FLAT_COST_MINIMUM_TIME}FlatCostType;
typedef struct{FlatCostType cost_type;double*waypoints;int n_waypoints;double*waypoint_times;double t0;double tf;int num_flat_outputs;double*state_bounds_lo;double*state_bounds_hi;double*input_bounds_lo;double*input_bounds_hi;int poly_degree;}FlatTrajectoryProblem;
int flat_trajectory_optimize(const FlatTrajectoryProblem*problem,FlatTrajectory*traj);

/* L7: Collision-free planning */
typedef struct{double center[3];double radius;}FlatObstacle;
int flat_collision_free_trajectory(const FlatControlAffineSystem*sys,const double*start,const double*goal,const FlatObstacle*obstacles,int n_obs,FlatTrajectory*traj);

/* Memory */
void flat_trajectory_free(FlatTrajectory*traj);
void flat_state_input_traj_free(FlatStateInputTrajectory*traj);
#endif
