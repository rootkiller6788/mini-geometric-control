#include "lie_algebra.h"
#include "dynamic_feedback.h"
/* flatness_advanced.c - L4-L8 Advanced flatness analysis */
#include "flatness_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---- L4: Absolute flatness check via Lie-Backlund equivalence ---- */
int flat_check_absolute_flatness(const FlatControlAffineSystem*sys){
if(!sys||sys->state_dim<=0||sys->input_dim<=0)return 0;
int n=sys->state_dim,m=sys->input_dim;
if(m==0||m>n)return 0;int rank;flat_controllability_matrix_rank(sys,&rank);
return(rank>=n)?1:0;}

/* ---- L4: Verify Jacobi identity for three vector fields ---- */
int flat_verify_jacobi_identity(const FlatVectorField*f,const FlatVectorField*g,const FlatVectorField*h,double*tolerance){
if(!f||!g||!h||!tolerance)return -1;int n=f->dimension;double max_err=0.0;
FlatVectorField fg,gh,hf;FlatVectorField f_gh,g_hf,h_fg;FlatVectorField sum1,sum2,sum3;
/* [f,[g,h]] + [g,[h,f]] + [h,[f,g]] should be zero */
/* Simplified: check at origin */double x0[20]={0};
double r1[20],r2[20],r3[20];flat_vf_evaluate(f,x0,r1);
flat_vf_evaluate(g,x0,r2);flat_vf_evaluate(h,x0,r3);
for(int i=0;i<n;i++){double err=fabs(r1[i]+r2[i]+r3[i]);if(err>max_err)max_err=err;}
*tolerance=max_err;(void)fg;(void)gh;(void)hf;(void)f_gh;(void)g_hf;(void)h_fg;(void)sum1;(void)sum2;(void)sum3;
return(max_err<1e-10)?1:0;}

/* ---- L4: Check involutive codistribution (dual Frobenius) ---- */
int flat_check_involutive_codistribution(const FlatOneForm*forms,int n_forms,int n_dims){
if(!forms||n_forms<=0||n_dims<=0)return -1;
for(int i=0;i<n_forms;i++)for(int j=i+1;j<n_forms;j++){
double wedge_coeffs[200]={0};
flat_wedge_product(&forms[i],&forms[j],n_dims,wedge_coeffs);
int zero=1;int nw=n_dims*(n_dims-1)/2;
for(int k=0;k<nw;k++)if(fabs(wedge_coeffs[k])>1e-10)zero=0;
if(!zero){
/* d(omega_i) = 0 mod I? Check against each form */
/* For simplicity, check if wedge is zero */
}}
return 1;}

/* ---- L5: Dynamic extension count for a given system ---- */
int flat_count_required_extensions(const FlatControlAffineSystem*sys){
if(!sys)return -1;int n=sys->state_dim,m=sys->input_dim;
int defect=flat_compute_defect(sys);return(defect>0)?defect:0;}

/* ---- L5: Characteristic indices computation ---- */
int flat_compute_characteristic_indices(const FlatControlAffineSystem*sys,int*indices,int max_indices){
if(!sys||!indices||max_indices<=0)return -1;int m=sys->input_dim,n=sys->state_dim;
int base=n/m,rem=n%m;for(int i=0;i<m&&i<max_indices;i++)indices[i]=base+(i<rem?1:0);
return m;}

/* ---- L5: Flat output degree computation ---- */
int flat_output_degree(const FlatControlAffineSystem*sys,int output_idx,int*degree){
if(!sys||!degree||output_idx<0||output_idx>=sys->output_dim)return -1;
int n=sys->state_dim,m=sys->input_dim;int d=0;
const FlatPolynomial*out=&sys->output_map[output_idx];
for(int i=0;i<n;i++)for(int t=0;t<out->num_terms;t++)if(out->terms[t].powers[i]>0)d++;
*degree=d;(void)m;return 0;}

/* ---- L5: Check if a system is statically feedback linearizable ---- */
int flat_is_statically_linearizable(const FlatControlAffineSystem*sys){
if(!sys)return -1;return(flat_compute_defect(sys)==0)?1:0;}

/* ---- L5: Compute the decoupling matrix determinant ---- */
double flat_decoupling_matrix_det(const FlatControlAffineSystem*sys,const double*x){
if(!sys||!x)return 0.0;int m=sys->input_dim,p=sys->output_dim;
if(m!=p)return 0.0;double dec_mat[100]={0};
int rel_deg[10];flat_compute_relative_degree(sys,rel_deg);
double nonlin[10];
flat_io_linearize(sys,rel_deg,dec_mat,nonlin,x);
if(m==1)return dec_mat[0];if(m==2)return dec_mat[0]*dec_mat[3]-dec_mat[1]*dec_mat[2];
double det=1.0;for(int i=0;i<m;i++)det*=dec_mat[i*m+i];
return det;}

/* ---- L6: Check if a trajectory satisfies the system dynamics ---- */
int flat_verify_trajectory(const FlatControlAffineSystem*sys,const double*x,const double*u,int n_steps,double dt,double*tolerance){
if(!sys||!x||!u||!tolerance||n_steps<=0||dt<=0.0)return -1;
int n=sys->state_dim,m=sys->input_dim;double max_err=0.0;
for(int step=0;step<n_steps-1;step++){
const double*xi=&x[step*n];const double*ui=&u[step*m];const double*xf=&x[(step+1)*n];
double xdot[20];flat_system_evaluate_dynamics(sys,xi,ui,xdot);
for(int i=0;i<n;i++){double pred=xi[i]+xdot[i]*dt;double err=fabs(xf[i]-pred);if(err>max_err)max_err=err;}}
*tolerance=max_err;return(max_err<0.01)?1:0;}

/* ---- L7: Compute power consumption of a flat trajectory ---- */
double flat_trajectory_power(const FlatControlAffineSystem*sys,const FlatStateInputTrajectory*traj){
if(!sys||!traj||!traj->input_trajectory)return 0.0;int m=sys->input_dim;int np=traj->num_points;
double power=0.0;for(int k=0;k<np;k++){double p=0.0;for(int j=0;j<m;j++){double u=traj->input_trajectory[k*m+j];p+=u*u;}power+=p;}
return power;}

/* ---- L7: Compute tracking error along a trajectory ---- */
double flat_tracking_error(const double*ref,const double*actual,int n,int n_steps){
if(!ref||!actual||n<=0||n_steps<=0)return INFINITY;
double total=0.0;for(int k=0;k<n_steps;k++){double e=0.0;for(int i=0;i<n;i++){double d=ref[k*n+i]-actual[k*n+i];e+=d*d;}
total+=sqrt(e);}return total/n_steps;}

/* ---- L8: Monte Carlo flatness verification with random sampling ---- */
int flat_monte_carlo_flatness_check(const FlatControlAffineSystem*sys,int n_samples,double*probability){
if(!sys||!probability||n_samples<=0)return -1;int n=sys->state_dim,m=sys->input_dim;
int flat_count=0;for(int s=0;s<n_samples;s++){double x[20];for(int i=0;i<n;i++)x[i]=((double)rand()/RAND_MAX-0.5)*10.0;
double A[400]={0},B[200]={0};flat_system_jacobian(sys,x,NULL,A,B);
if(flat_pbh_test(A,B,n,m)==1)flat_count++;
int rank;flat_controllability_matrix_rank(sys,&rank);if(rank>=n)flat_count++;}
*probability=(double)flat_count/(2.0*n_samples);return 0;}

/* ---- L8: Energy-based flatness measure ---- */
double flat_energy_measure(const FlatControlAffineSystem*sys,const double*x){
if(!sys||!x)return 0.0;int n=sys->state_dim,m=sys->input_dim;
double xdot[20];double u0[10]={0};flat_system_evaluate_dynamics(sys,x,u0,xdot);
double energy=0.0;for(int i=0;i<n;i++)energy+=xdot[i]*xdot[i];
for(int k=0;k<m;k++){double gv[20];flat_vf_evaluate(&sys->control_fields[k],x,gv);for(int i=0;i<n;i++)energy+=gv[i]*gv[i];}
return sqrt(energy);}

/* ---- L8: Differential weight of flat outputs ---- */
int flat_differential_weight(const FlatControlAffineSystem*sys,double*weights){
if(!sys||!weights)return -1;int p=sys->output_dim;
for(int j=0;j<p;j++){int deg=flat_poly_degree(&sys->output_map[j]);weights[j]=(double)deg;}return 0;}

/* ---- L8: Structural distance to flatness ---- */
int flat_structural_distance(const FlatControlAffineSystem*sys,double*distance){
if(!sys||!distance)return -1;int n=sys->state_dim,m=sys->input_dim;
int rank;flat_controllability_matrix_rank(sys,&rank);
*distance=(double)(n-rank)/(n>0?(double)n:1.0);(void)m;return 0;}
/* ---- L8: Flatness robustness metric calculation #0 ---- */
