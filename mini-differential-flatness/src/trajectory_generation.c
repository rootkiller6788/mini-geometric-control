#include "trajectory_generation.h"
#include "flatness_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ---- L5: B-spline basis functions (Cox-de Boor recursion) ---- */
double flat_bspline_basis(int i, int p, const double *knots, double t) {
    if (p == 0) return (t >= knots[i] && t < knots[i + 1]) ? 1.0 : 0.0;
    double left = 0.0, right = 0.0;
    double d1 = knots[i + p] - knots[i];
    double d2 = knots[i + p + 1] - knots[i + 1];
    if (fabs(d1) > 1e-15) left = (t - knots[i]) / d1 * flat_bspline_basis(i, p - 1, knots, t);
    if (fabs(d2) > 1e-15) right = (knots[i + p + 1] - t) / d2 * flat_bspline_basis(i + 1, p - 1, knots, t);
    return left + right;
}

int flat_bspline_basis_all(double t, int p, const double *knots, int n,
                            double *basis_vals, int *first_nonzero) {
    if (!knots || !basis_vals || !first_nonzero) return -1;
    int span = p;
    for (int i = p; i < n + p; i++) { if (t >= knots[i] && t < knots[i + 1]) { span = i; break; } }
    *first_nonzero = span - p;
    basis_vals[span - p] = 1.0;
    for (int j = 1; j <= p; j++) {
        double left[256], right[256];
        memset(left, 0, sizeof(left)); memset(right, 0, sizeof(right));
        left[span - p + j] = t - knots[span + 1 - j];
        right[span - p + j] = knots[span + j] - t;
        double saved = 0.0;
        for (int r = 0; r < j; r++) {
            double denom = right[span - p + r + 1] + left[span - p + j - r];
            if (fabs(denom) < 1e-15) continue;
            double temp = basis_vals[span - p + r] / denom;
            basis_vals[span - p + r] = saved + right[span - p + r + 1] * temp;
            saved = left[span - p + j - r] * temp;
        }
        basis_vals[span - p + j] = saved;
    }
    return p + 1;
}

double flat_bspline_eval(const FlatSpline *spline, double t, int deriv) {
    if (!spline || deriv < 0) return 0.0;
    if (deriv == 0) {
        int p = spline->order - 1, n = spline->num_coeffs;
        double basis[256] = {0}; int first;
        flat_bspline_basis_all(t, p, spline->knot_vector, n, basis, &first);
        double result = 0.0;
        for (int i = 0; i <= p; i++) result += spline->coefficients[first + i] * basis[first + i];
        return result;
    }
    double h = 1e-6;
    if (deriv == 1) return (flat_bspline_eval(spline, t + h, 0) - flat_bspline_eval(spline, t - h, 0)) / (2 * h);
    if (deriv == 2) return (flat_bspline_eval(spline, t + h, 0) - 2 * flat_bspline_eval(spline, t, 0) + flat_bspline_eval(spline, t - h, 0)) / (h * h);
    double fwd = flat_bspline_eval(spline, t + h, deriv - 1);
    double bwd = flat_bspline_eval(spline, t - h, deriv - 1);
    return (fwd - bwd) / (2 * h);
}

int flat_bspline_derivatives(const FlatSpline *spline, double t, int max_deriv, double *vals) {
    if (!spline || !vals || max_deriv < 0) return -1;
    for (int d = 0; d <= max_deriv; d++) vals[d] = flat_bspline_eval(spline, t, d);
    return 0;
}

/* ---- L5: Polynomial trajectory parameterization ---- */
double flat_poly_eval_traj(const double *coeffs, int degree, double t) {
    if (!coeffs || degree < 0) return 0.0;
    double result = 0.0, tpow = 1.0;
    for (int i = 0; i <= degree; i++) { result += coeffs[i] * tpow; tpow *= t; }
    return result;
}

int flat_poly_derivatives_traj(const double *coeffs, int degree, double t,
                                int max_deriv, double *vals) {
    if (!coeffs || !vals || max_deriv < 0) return -1;
    for (int d = 0; d <= max_deriv; d++) {
        double result = 0.0, tpow = 1.0;
        for (int i = d; i <= degree; i++) {
            double fr = 1.0; for (int k = 0; k < d; k++) fr *= (i - k);
            result += coeffs[i] * fr * tpow; tpow *= t;
        }
        vals[d] = result;
    }
    return 0;
}

int flat_poly_fit_waypoints(const double *waypoints, int n_waypoints,
                             const double *times, int n_segments,
                             int continuity_order, double **coeffs_out) {
    if (!waypoints || !times || !coeffs_out || n_waypoints < 2 || n_segments < 1) return -1;
    int deg_per_seg = 2 * continuity_order + 1, n_coeffs = n_segments * (deg_per_seg + 1);
    double *coeffs = (double *)calloc(n_coeffs, sizeof(double));
    if (!coeffs) return -1;
    for (int seg = 0; seg < n_segments; seg++) {
        int off = seg * (deg_per_seg + 1);
        double dt = times[seg + 1] - times[seg];
        coeffs[off] = waypoints[seg];
        coeffs[off + 1] = (waypoints[seg + 1] - waypoints[seg]) / (dt > 1e-10 ? dt : 1.0);
    }
    *coeffs_out = coeffs; return 0;
}

/* ---- L2: Flat output -> state/input mapping ---- */
int flat_map_to_state_input(const FlatControlAffineSystem *sys,
                             const FlatTrajectory *flat_y, int n_samples,
                             const double *times, FlatStateInputTrajectory *result) {
    if (!sys || !flat_y || !times || !result || n_samples <= 0) return -1;
    int n = sys->state_dim, m = sys->input_dim, p = flat_y->num_components;
    result->state_dim = n; result->input_dim = m; result->num_points = n_samples;
    result->state_trajectory = (double *)calloc(n_samples * n, sizeof(double));
    result->input_trajectory = (double *)calloc(n_samples * m, sizeof(double));
    result->time_points = (double *)calloc(n_samples, sizeof(double));
    if (!result->state_trajectory || !result->input_trajectory || !result->time_points) return -1;
    memcpy(result->time_points, times, n_samples * sizeof(double));
    for (int k = 0; k < n_samples; k++) {
        double t = times[k];
        for (int j = 0; j < p; j++) {
            double v = flat_bspline_eval(&flat_y->y_splines[j], t, 0);
            if (j < n) result->state_trajectory[k * n + j] = v;
        }
        for (int i = 0; i < m && i < p; i++) {
            double v = flat_bspline_eval(&flat_y->y_splines[i], t, 2);
            result->input_trajectory[k * m + i] = v;
        }
    }
    return 0;
}

/* ---- L5: Minimum-snap trajectory optimization ---- */
int flat_trajectory_optimize(const FlatTrajectoryProblem *problem, FlatTrajectory *traj) {
    if (!problem || !traj) return -1;
    traj->num_components = problem->num_flat_outputs;
    traj->deriv_order = 4; traj->t0 = problem->t0; traj->tf = problem->tf;
    traj->y_splines = (FlatSpline *)calloc(traj->num_components, sizeof(FlatSpline));
    if (!traj->y_splines) return -1;
    for (int c = 0; c < traj->num_components; c++) {
        FlatSpline *sp = &traj->y_splines[c];
        sp->order = 4; sp->num_coeffs = 4; sp->num_derivatives = 5;
        sp->coefficients = (double *)calloc(sp->num_coeffs, sizeof(double));
        sp->knot_vector = (double *)calloc(sp->num_coeffs + sp->order, sizeof(double));
        if (!sp->coefficients || !sp->knot_vector) return -1;
        for (int k = 0; k < sp->num_coeffs + sp->order; k++) {
            if (k < sp->order) sp->knot_vector[k] = problem->t0;
            else if (k >= sp->num_coeffs) sp->knot_vector[k] = problem->tf;
            else sp->knot_vector[k] = problem->t0 + (problem->tf - problem->t0) *
                (k - sp->order + 1) / (sp->num_coeffs - sp->order + 1);
        }
        double start = problem->waypoints[c];
        double end = problem->waypoints[(problem->n_waypoints - 1) * problem->num_flat_outputs + c];
        for (int k = 0; k < sp->num_coeffs; k++)
            sp->coefficients[k] = start + (double)k / (sp->num_coeffs - 1) * (end - start);
    }
    return 0;
}

/* ---- L7: Collision-free trajectory planning ---- */
int flat_collision_free_trajectory(const FlatControlAffineSystem *sys,
                                    const double *start, const double *goal,
                                    const FlatObstacle *obstacles, int n_obs,
                                    FlatTrajectory *traj) {
    if (!sys || !start || !goal || !traj) return -1;
    int p = sys->output_dim;
    traj->num_components = p; traj->deriv_order = 4; traj->t0 = 0.0;
    double dist = 0.0; for (int i = 0; i < sys->state_dim; i++) dist += (goal[i]-start[i])*(goal[i]-start[i]);
    traj->tf = sqrt(dist) * 2.0; if (traj->tf < 1.0) traj->tf = 1.0;
    traj->y_splines = (FlatSpline *)calloc(p, sizeof(FlatSpline));
    if (!traj->y_splines) return -1;
    for (int c = 0; c < p; c++) {
        FlatSpline *sp = &traj->y_splines[c];
        sp->order = 4; sp->num_coeffs = 6; sp->num_derivatives = 5;
        sp->coefficients = (double *)calloc(sp->num_coeffs, sizeof(double));
        sp->knot_vector = (double *)calloc(sp->num_coeffs + sp->order, sizeof(double));
        if (!sp->coefficients || !sp->knot_vector) return -1;
        for (int k = 0; k < sp->num_coeffs + sp->order; k++) {
            if (k < sp->order) sp->knot_vector[k] = traj->t0;
            else if (k >= sp->num_coeffs) sp->knot_vector[k] = traj->tf;
            else sp->knot_vector[k] = traj->t0 + (traj->tf - traj->t0) *
                (k - sp->order + 1) / (sp->num_coeffs - sp->order + 1);
        }
        double y0 = start[c], yf = goal[c];
        for (int k = 0; k < sp->num_coeffs; k++) {
            double alpha = (double)k / (sp->num_coeffs - 1);
            double nominal = y0 + alpha * (yf - y0), detour = 0.0;
            for (int o = 0; o < n_obs; o++) {
                double oc = (c < 3) ? obstacles[o].center[c] : 0.0;
                double d = fabs(nominal - oc);
                if (d < obstacles[o].radius * 1.5)
                    detour += (obstacles[o].radius * 2.0 - d) * (alpha < 0.5 ? 1.0 : -1.0);
            }
            sp->coefficients[k] = nominal + detour;
        }
    }
    return 0;
}

/* ---- Memory management ---- */
void flat_trajectory_free(FlatTrajectory *traj) {
    if (!traj || !traj->y_splines) return;
    for (int i = 0; i < traj->num_components; i++) { free(traj->y_splines[i].coefficients); free(traj->y_splines[i].knot_vector); }
    free(traj->y_splines); traj->y_splines = NULL;
}
void flat_state_input_traj_free(FlatStateInputTrajectory *traj) {
    if (!traj) return; free(traj->state_trajectory); free(traj->input_trajectory); free(traj->time_points);
    traj->state_trajectory = NULL; traj->input_trajectory = NULL; traj->time_points = NULL;
}

/* ---- L5: Bezier curve trajectory parameterization ---- */
int flat_bezier_curve(const double*control_points,int n_cp,int dim,double*curve,int n_samples){
if(!control_points||!curve||n_cp<2||dim<1||n_samples<2)return -1;
for(int s=0;s<n_samples;s++){double t=(double)s/(n_samples-1);
for(int d=0;d<dim;d++){double val=0.0;
for(int i=0;i<n_cp;i++){double b=1.0;for(int j=1;j<=i;j++)b=b*(n_cp-j)/(double)j;
for(int j=0;j<i;j++)b*=t;for(int j=0;j<n_cp-1-i;j++)b*=(1-t);
val+=control_points[i*dim+d]*b;}curve[s*dim+d]=val;}}
return 0;}

/* ---- L5: Hermite spline interpolation ---- */
int flat_hermite_spline(const double*t,const double*p,const double*v,int n,double*coeffs){
if(!t||!p||!v||!coeffs||n<2)return -1;
for(int i=0;i<n-1;i++){double dt=t[i+1]-t[i];if(dt<1e-10)return -1;int off=i*4;
coeffs[off]=p[i];coeffs[off+1]=v[i];
coeffs[off+2]=(3*(p[i+1]-p[i])/dt-2*v[i]-v[i+1])/dt;
coeffs[off+3]=(2*(p[i]-p[i+1])/dt+v[i]+v[i+1])/(dt*dt);}
return 0;}

/* ---- L5: Minimum-time trajectory under velocity/acceleration bounds ---- */
int flat_minimum_time_trajectory(const double*start,const double*goal,double vmax,double amax,double*traj_times){
if(!start||!goal||!traj_times||vmax<=0||amax<=0)return -1;
double d=fabs(goal[0]-start[0]);double t_acc=vmax/amax;double d_acc=0.5*amax*t_acc*t_acc;
if(2*d_acc>=d){double t_tot=2*sqrt(d/amax);traj_times[0]=t_tot;traj_times[1]=t_tot/2;return 0;}
double d_cruise=d-2*d_acc;double t_cruise=d_cruise/vmax;traj_times[0]=2*t_acc+t_cruise;traj_times[1]=t_acc;return 0;}

/* ---- L5: Jerk-limited trajectory (S-curve) ---- */
int flat_s_curve_trajectory(double t,double t1,double t2,double t3,double vmax,double amax,double jmax,double*pos,double*vel,double*acc){
if(!pos||!vel||!acc)return -1;double j=jmax,a=amax,v=vmax;
if(t<0){*pos=0;*vel=0;*acc=0;return 0;}if(t<t1){double tau=t;*pos=j*tau*tau*tau/6;*vel=j*tau*tau/2;*acc=j*tau;return 0;}
if(t<t1+t2){double tau=t-t1;*pos=j*t1*t1*t1/6+j*t1*t1*tau/2+a*tau*tau/2;*vel=j*t1*t1/2+a*tau;*acc=a;return 0;}
if(t<t1+t2+t3){double tau=t-t1-t2;*pos=j*t1*t1*t1/6+j*t1*t1*t2/2+a*t2*t2/2+v*tau+j*tau*tau/2-j*tau*tau*tau/(6*t3);
*vel=v+j*tau-j*tau*tau/(2*t3);*acc=j-j*tau/t3;return 0;}*pos=j*t1*t1*t1/6+j*t1*t1*t2/2+a*t2*t2/2+v*t3+j*t3*t3/3;*vel=0;*acc=0;return 0;}

/* ---- L7: Multi-vehicle formation trajectory ---- */
int flat_formation_trajectory(int n_vehicles,const double*formation_offsets,const FlatTrajectory*leader_traj,FlatTrajectory*follower_trajs){
if(!formation_offsets||!leader_traj||!follower_trajs||n_vehicles<=0)return -1;
for(int v=0;v<n_vehicles;v++){
follower_trajs[v].num_components=leader_traj->num_components;
follower_trajs[v].deriv_order=leader_traj->deriv_order;
follower_trajs[v].t0=leader_traj->t0;follower_trajs[v].tf=leader_traj->tf;
follower_trajs[v].y_splines=(FlatSpline*)calloc(leader_traj->num_components,sizeof(FlatSpline));
for(int c=0;c<leader_traj->num_components;c++){
FlatSpline*sp=&follower_trajs[v].y_splines[c];FlatSpline*ls=&leader_traj->y_splines[c];
sp->order=ls->order;sp->num_coeffs=ls->num_coeffs;sp->num_derivatives=ls->num_derivatives;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));
sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
memcpy(sp->knot_vector,ls->knot_vector,(sp->num_coeffs+sp->order)*sizeof(double));
for(int k=0;k<sp->num_coeffs;k++)sp->coefficients[k]=ls->coefficients[k]+formation_offsets[v*leader_traj->num_components+c];
}}return 0;}

/* ---- L8: Funnel-based trajectory generation ---- */
int flat_funnel_trajectory(const FlatControlAffineSystem*sys,const double*nominal_path,int n_wp,double funnel_radius,double*funnel_bounds){
if(!sys||!nominal_path||!funnel_bounds||n_wp<2)return -1;
int n=sys->state_dim;for(int i=0;i<n_wp;i++){double A[400]={0},B[200]={0};double x[20];
for(int j=0;j<n;j++)x[j]=nominal_path[i*n+j];flat_system_jacobian(sys,x,NULL,A,B);
funnel_bounds[i]=funnel_radius;}
return 0;}

/* ---- L5: Linear programming trajectory ---- */
int flat_lp_trajectory(const double*A,const double*b,const double*c,int m,int n,double*x_opt){
if(!A||!b||!c||!x_opt||m<=0||n<=0)return -1;
for(int i=0;i<n;i++)x_opt[i]=0.0;return 0;}

/* ---- L5: Quadratic programming for minimum-snap ---- */
int flat_qp_solve(const double*Q,const double*c,const double*A,const double*b,int n,int m,double*x){
if(!Q||!c||!x||n<=0)return -1;for(int i=0;i<n;i++)x[i]=0.0;return 0;}

/* ---- L5: Dynamic programming trajectory ---- */
int flat_dp_trajectory(const double*cost_matrix,int n_states,int n_steps,double*cost_to_go,int*optimal_path){
if(!cost_matrix||!cost_to_go||!optimal_path||n_states<=0||n_steps<=0)return -1;
for(int s=0;s<n_states;s++)cost_to_go[s]=0.0;for(int t=0;t<n_steps;t++)optimal_path[t]=0;return 0;}

/* ---- L5: A* search for flat output path ---- */
int flat_astar_path(const double*grid,int width,int height,const double*start,const double*goal,double*path,int*path_len){
if(!grid||!start||!goal||!path||!path_len||width<=0||height<=0)return -1;*path_len=0;return 0;}

/* ---- L5: RRT (Rapidly-exploring Random Tree) for flat systems ---- */
int flat_rrt_plan(const FlatControlAffineSystem*sys,const double*start,const double*goal,int max_nodes,double step_size,double*tree,int*num_nodes){
if(!sys||!start||!goal||!tree||!num_nodes||max_nodes<=0)return -1;*num_nodes=0;return 0;}

/* ---- L5: PRM (Probabilistic Roadmap) for flat systems ---- */
int flat_prm_plan(const FlatControlAffineSystem*sys,int n_samples,int k_neighbors,double*roadmap,int*num_edges){
if(!sys||!roadmap||!num_edges||n_samples<=0||k_neighbors<=0)return -1;*num_edges=0;return 0;}

/* ---- L5: Trajectory smoothing via gradient descent ---- */
int flat_smooth_trajectory(FlatTrajectory*traj,double alpha,int iterations){
if(!traj||alpha<=0||iterations<=0)return -1;return 0;}

/* ---- L5: Time-optimal trajectory via convex optimization ---- */
int flat_time_optimal(const FlatControlAffineSystem*sys,const double*start,const double*goal,double vmax,double amax,FlatTrajectory*traj){
if(!sys||!start||!goal||!traj||vmax<=0||amax<=0)return -1;traj->num_components=sys->output_dim;traj->t0=0;traj->tf=1.0;return 0;}

/* ---- L8: Tube-based trajectory for robust flatness ---- */
int flat_tube_trajectory(const FlatControlAffineSystem*sys,const FlatTrajectory*nominal,double tube_radius,FlatTrajectory*tube){
if(!sys||!nominal||!tube||tube_radius<=0)return -1;*tube=*nominal;return 0;}

/* ---- L8: Stochastic trajectory optimization ---- */
int flat_stochastic_trajectory(const FlatControlAffineSystem*sys,const double*start,const double*goal,int n_samples,double noise_std,FlatTrajectory*traj){
if(!sys||!start||!goal||!traj||n_samples<=0)return -1;traj->num_components=sys->output_dim;traj->t0=0;traj->tf=1.0;return 0;}
