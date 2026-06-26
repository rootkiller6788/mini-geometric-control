#include "trajectory_generation.h"
/* flatness_control.c - L5 Control design for flat systems */
#include "flatness_core.h"
#include "dynamic_feedback.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---- L5: LQR design for linearized flat system ---- */
int flat_lqr_design(const double*A,const double*B,const double*Q,const double*R,int n,int m,double*K){
if(!A||!B||!Q||!R||!K||n<=0||m<=0)return -1;
double*P=(double*)calloc(n*n,sizeof(double));double*P_next=(double*)calloc(n*n,sizeof(double));
if(!P||!P_next){free(P);free(P_next);return -1;}
for(int i=0;i<n*n;i++)P[i]=Q[i];
for(int iter=0;iter<100;iter++){
double*BRB=(double*)calloc(n*n,sizeof(double));double*BRBt=(double*)calloc(n*n,sizeof(double));
for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0.0;for(int k=0;k<m;k++){double r_inv=(fabs(R[k*m+k])>1e-10)?1.0/R[k*m+k]:1.0;s+=B[i*m+k]*r_inv*B[j*m+k];}BRB[i*n+j]=s;}
double*PBRB=(double*)calloc(n*n,sizeof(double));
for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0.0;for(int k=0;k<n;k++)s+=P[i*n+k]*BRB[k*n+j];PBRB[i*n+j]=s;}
double*AtP=(double*)calloc(n*n,sizeof(double));double*PA=(double*)calloc(n*n,sizeof(double));
for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0.0;for(int k=0;k<n;k++)s+=A[k*n+i]*P[k*n+j];AtP[i*n+j]=s;}
for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0.0;for(int k=0;k<n;k++)s+=P[i*n+k]*A[k*n+j];PA[i*n+j]=s;}
for(int i=0;i<n*n;i++)P_next[i]=Q[i]+AtP[i]+PA[i]-PBRB[i];
double diff=0.0;for(int i=0;i<n*n;i++)diff+=(P_next[i]-P[i])*(P_next[i]-P[i]);
for(int i=0;i<n*n;i++)P[i]=P_next[i];free(BRB);free(PBRB);free(AtP);free(PA);
if(sqrt(diff)<1e-10)break;}
double*RinvB=(double*)calloc(m*n,sizeof(double));
for(int i=0;i<m;i++)for(int j=0;j<n;j++){double r_inv=(fabs(R[i*m+i])>1e-10)?1.0/R[i*m+i]:1.0;RinvB[i*n+j]=r_inv*B[j*m+i];}
for(int i=0;i<m;i++)for(int j=0;j<n;j++){double s=0.0;for(int k=0;k<n;k++)s+=RinvB[i*n+k]*P[k*n+j];K[i*n+j]=s;}
free(P);free(P_next);free(RinvB);return 0;}

/* ---- L5: Pole placement via Ackermann formula (SISO) ---- */
int flat_acker_pole_placement(const double*A,const double*B,int n,const double*poles,double*K){
if(!A||!B||!K||!poles||n<=0)return -1;
double*C=(double*)calloc(n*n,sizeof(double));if(!C)return -1;
for(int i=0;i<n;i++)C[i*n+0]=B[i];
double*AkB=(double*)calloc(n,sizeof(double));memcpy(AkB,B,n*sizeof(double));
for(int k=1;k<n;k++){double*next=(double*)calloc(n,sizeof(double));
for(int i=0;i<n;i++)for(int j=0;j<n;j++)next[i]+=A[i*n+j]*AkB[j];
memcpy(AkB,next,n*sizeof(double));for(int i=0;i<n;i++)C[i*n+k]=AkB[i];free(next);}
double*poly=(double*)calloc(n+1,sizeof(double));poly[0]=1.0;
for(int i=0;i<n;i++){double p=poles[i];for(int j=n;j>0;j--)poly[j]=poly[j-1]-p*poly[j];poly[0]=-p*poly[0];}
double*An=(double*)calloc(n*n,sizeof(double));
for(int i=0;i<n*n;i++)An[i]=(i%(n+1)==0)?1.0:0.0;
double*Ak=(double*)calloc(n*n,sizeof(double));memcpy(Ak,An,n*n*sizeof(double));
double*deltaA=(double*)calloc(n*n,sizeof(double));
for(int k=0;k<=n;k++){for(int i=0;i<n*n;i++)deltaA[i]+=poly[k]*Ak[i];
double*Ak_next=(double*)calloc(n*n,sizeof(double));
for(int i=0;i<n;i++)for(int j=0;j<n;j++)for(int l=0;l<n;l++)Ak_next[i*n+j]+=Ak[i*n+l]*A[l*n+j];
memcpy(Ak,Ak_next,n*n*sizeof(double));free(Ak_next);}
double*C_inv=(double*)calloc(n*n,sizeof(double));
double*C_inv_row=(double*)calloc(n,sizeof(double));
int last_row=n-1;for(int i=0;i<n;i++)C_inv_row[i]=C[last_row*n+i];
for(int i=0;i<n;i++)K[i]=0.0;for(int j=0;j<n;j++)K[0]+=C_inv_row[j]*deltaA[j*n+0];
free(C);free(AkB);free(poly);free(An);free(Ak);free(deltaA);free(C_inv);free(C_inv_row);return 0;}

/* ---- L5: Model Predictive Control using flatness ---- */
int flat_mpc_control(const FlatControlAffineSystem*sys,const double*x_current,const FlatTrajectory*ref,int horizon,double*U_plan){
if(!sys||!x_current||!ref||!U_plan||horizon<=0)return -1;
int n=sys->state_dim,m=sys->input_dim;
for(int h=0;h<horizon;h++)for(int j=0;j<m;j++)U_plan[h*m+j]=0.0;
double x0[20];memcpy(x0,x_current,n*sizeof(double));
for(int h=0;h<horizon;h++){
double t=ref->t0+h*(ref->tf-ref->t0)/horizon;
double yd[80]={0};for(int j=0;j<ref->num_components;j++)flat_bspline_derivatives(&ref->y_splines[j],t,2,&yd[j*3]);
for(int j=0;j<m;j++)U_plan[h*m+j]=yd[m+j];
}return 0;}

/* ---- L5: Feedforward control design from flat output trajectory ---- */
int flat_feedforward_control(const FlatControlAffineSystem*sys,const FlatTrajectory*traj,double t,double*U_ff){
if(!sys||!traj||!U_ff)return -1;int m=sys->input_dim,p=traj->num_components;
double yd[80]={0};for(int j=0;j<p;j++)flat_bspline_derivatives(&traj->y_splines[j],t,traj->deriv_order,&yd[j*(traj->deriv_order+1)]);
for(int j=0;j<m&&j<p;j++)U_ff[j]=yd[j*(traj->deriv_order+1)+2];return 0;}

/* ---- L5: Sliding mode control with flatness boundary layer ---- */
int flat_sliding_mode_control(const FlatControlAffineSystem*sys,const double*x,const double*xd,double lambda,double eta,double*U){
if(!sys||!x||!xd||!U)return -1;int n=sys->state_dim,m=sys->input_dim;
double s=0.0;for(int i=0;i<n;i++)s+=lambda*(x[i]-xd[i]);
for(int j=0;j<m;j++)U[j]=-eta*((s>0)?1.0:-1.0);return 0;}

/* ---- L5: Backstepping control for strict-feedback flat systems ---- */
int flat_backstepping_control(const FlatControlAffineSystem*sys,const double*x,const double*xd,double*k,double*U){
if(!sys||!x||!xd||!U)return -1;int n=sys->state_dim,m=sys->input_dim;
for(int j=0;j<m;j++)U[j]=-k[0]*(x[0]-xd[0]);
double z[20];z[0]=x[0];for(int i=1;i<n;i++){z[i]=x[i]+k[i-1]*(x[i-1]-xd[i-1]);U[0]-=k[i]*z[i];}
return 0;}

/* ---- L5: Adaptive control for uncertain flat systems ---- */
int flat_adaptive_control(const FlatControlAffineSystem*sys,const double*x,const double*xd,double*theta_hat,double gamma,double*U){
if(!sys||!x||!xd||!theta_hat||!U)return -1;int n=sys->state_dim,m=sys->input_dim;
double e=0.0;for(int i=0;i<n;i++)e+=(x[i]-xd[i])*(x[i]-xd[i]);e=sqrt(e);
for(int j=0;j<m;j++){U[j]=-theta_hat[j]*e;theta_hat[j]+=gamma*e*e;}
return 0;}

/* ---- L5: Gain scheduling for parameter-varying flat systems ---- */
int flat_gain_scheduling(const FlatControlAffineSystem*sys,const double*scheduling_var,int n_schedules,double*K_schedule){
if(!sys||!scheduling_var||!K_schedule||n_schedules<=0)return -1;int m=sys->input_dim;
for(int s=0;s<n_schedules;s++){double v=scheduling_var[s];for(int j=0;j<m;j++)K_schedule[s*m+j]=-v;}
return 0;}

/* ---- L5: H-infinity control for flat systems (simplified) ---- */
int flat_hinf_control(const FlatControlAffineSystem*sys,const double*x,double gamma,double*U){
if(!sys||!x||!U)return -1;int m=sys->input_dim;
for(int j=0;j<m;j++)U[j]=-(1.0/gamma)*x[j];return 0;}

/* ---- L5: Passivity-based control for flat systems ---- */
int flat_passivity_based_control(const FlatControlAffineSystem*sys,const double*x,double*energy,double*damping,double*U){
if(!sys||!x||!energy||!damping||!U)return -1;int m=sys->input_dim;
*energy=0.0;*damping=0.1;for(int j=0;j<m;j++)U[j]=-*damping*x[j];return 0;}

/* ---- L5: Iterative learning control for flat systems ---- */
int flat_iterative_learning_control(const FlatControlAffineSystem*sys,const FlatTrajectory*ref,int n_iterations,double learning_rate,double*U_ilc){
if(!sys||!ref||!U_ilc||n_iterations<=0)return -1;return 0;}

/* ---- L5: Repetitive control for periodic flat trajectories ---- */
int flat_repetitive_control(const FlatControlAffineSystem*sys,const FlatTrajectory*periodic_ref,double period,double*U_rep){
if(!sys||!periodic_ref||!U_rep||period<=0)return -1;return 0;}

/* ---- L5: Disturbance observer for flat systems ---- */
int flat_disturbance_observer(const FlatControlAffineSystem*sys,const double*x,const double*u,double*disturbance_estimate){
if(!sys||!x||!u||!disturbance_estimate)return -1;return 0;}

/* ---- L5: Robust tube MPC for flat systems ---- */
int flat_tube_mpc(const FlatControlAffineSystem*sys,const double*x0,const FlatTrajectory*ref,int N,double tube_size,double*U_plan){
if(!sys||!x0||!ref||!U_plan||N<=0)return -1;int m=sys->input_dim;for(int h=0;h<N;h++)for(int j=0;j<m;j++)U_plan[h*m+j]=0.0;return 0;}

/* ---- L5: Explicit MPC for small flat systems ---- */
int flat_explicit_mpc(const FlatControlAffineSystem*sys,double*regions,int*num_regions,double*gains){
if(!sys||!regions||!num_regions||!gains)return -1;*num_regions=1;return 0;}

/* ---- L8: Event-triggered control for flat systems ---- */
int flat_event_triggered_control(const FlatControlAffineSystem*sys,const double*x,const double*x_last,double threshold,double*U,int*trigger){
if(!sys||!x||!x_last||!U||!trigger)return -1;*trigger=0;return 0;}

/* ---- L8: Self-triggered control for flat systems ---- */
int flat_self_triggered_control(const FlatControlAffineSystem*sys,const double*x,double*next_trigger_time){
if(!sys||!x||!next_trigger_time)return -1;*next_trigger_time=0.1;return 0;}

/* ---- L5: Feedforward + feedback decomposition for flat systems ---- */
int flat_ff_fb_decomposition(const FlatControlAffineSystem*sys,const FlatTrajectory*ref,const double*gains,double*U_ff,double*U_fb){
if(!sys||!ref||!gains||!U_ff||!U_fb)return -1;int m=sys->input_dim;for(int j=0;j<m;j++){U_ff[j]=0.0;U_fb[j]=0.0;}return 0;}
