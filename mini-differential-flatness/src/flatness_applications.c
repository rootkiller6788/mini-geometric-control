#include "trajectory_generation.h"
#include "dynamic_feedback.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "flatness_applications.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int flat_quadrotor_min_snap(const QuadrotorParams*params,const double*wp,int n,const double*st,FlatTrajectory*t){
if(!params||!wp||!t||n<2)return -1;t->num_components=4;t->deriv_order=4;t->t0=0.0;t->tf=0.0;
for(int i=0;i<n-1;i++)t->tf+=(st?st[i]:1.0);t->y_splines=(FlatSpline*)calloc(4,sizeof(FlatSpline));if(!t->y_splines)return -1;
for(int c=0;c<4;c++){FlatSpline*sp=&t->y_splines[c];sp->order=4;sp->num_coeffs=n+2;sp->num_derivatives=5;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
if(!sp->coefficients||!sp->knot_vector)return -1;
for(int k=0;k<sp->num_coeffs+sp->order;k++){if(k<sp->order)sp->knot_vector[k]=0.0;else if(k>=sp->num_coeffs)sp->knot_vector[k]=t->tf;
else sp->knot_vector[k]=t->tf*(k-sp->order+1)/(sp->num_coeffs-sp->order+1);}
for(int k=0;k<sp->num_coeffs;k++){int wi=(k<n)?k:n-1;sp->coefficients[k]=wp[wi*4+c];}}return 0;}

int flat_quadrotor_rotor_speeds(const QuadrotorParams*params,const double*u,double*omega){
if(!params||!u||!omega)return -1;double F=u[0],tx=u[1],ty=u[2],tz=u[3];
double l=params->arm_length,kf=params->kf,km=params->km;
double w1=F/(4*kf)+tx/(2*kf*l)+ty/(2*kf*l)+tz/(4*km);omega[0]=(w1>0)?sqrt(w1):0;
double w2=F/(4*kf)-tx/(2*kf*l)+ty/(2*kf*l)-tz/(4*km);omega[1]=(w2>0)?sqrt(w2):0;
double w3=F/(4*kf)-tx/(2*kf*l)-ty/(2*kf*l)+tz/(4*km);omega[2]=(w3>0)?sqrt(w3):0;
double w4=F/(4*kf)+tx/(2*kf*l)-ty/(2*kf*l)-tz/(4*km);omega[3]=(w4>0)?sqrt(w4):0;return 0;}

int flat_trailer_parking(const TrailerParams*params,const double*start,const double*goal,FlatTrajectory*traj){
if(!params||!start||!goal||!traj)return -1;int nt=params->n_trailers;double d=0;for(int i=0;i<2;i++){double dd=goal[i]-start[i];d+=dd*dd;}
double T=sqrt(d)/params->max_velocity;if(T<2.0)T=2.0;traj->num_components=2;traj->deriv_order=3;traj->t0=0.0;traj->tf=T;
traj->y_splines=(FlatSpline*)calloc(2,sizeof(FlatSpline));if(!traj->y_splines)return -1;
for(int c=0;c<2;c++){FlatSpline*sp=&traj->y_splines[c];sp->order=4;sp->num_coeffs=6;sp->num_derivatives=4;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
if(!sp->coefficients||!sp->knot_vector)return -1;
for(int k=0;k<sp->num_coeffs+sp->order;k++){if(k<sp->order)sp->knot_vector[k]=0.0;else if(k>=sp->num_coeffs)sp->knot_vector[k]=traj->tf;
else sp->knot_vector[k]=traj->tf*(k-sp->order+1)/(sp->num_coeffs-sp->order+1);}
double y0=start[c],yf=goal[c];for(int k=0;k<sp->num_coeffs;k++){double a=(double)k/(sp->num_coeffs-1);sp->coefficients[k]=y0*(1.0-a)+yf*a;}}
(void)nt;return 0;}

int flat_trailer_hitch_angles(const TrailerParams*params,const FlatTrajectory*traj,double t,double*ha){
if(!params||!traj||!ha)return -1;int nt=params->n_trailers;
double yd[8]={0};for(int j=0;j<2;j++)flat_bspline_derivatives(&traj->y_splines[j],t,3,&yd[j*4]);
double xd=yd[1],ydv=yd[5],xdd=yd[2],ydd=yd[6];double vsq=xd*xd+ydv*ydv;
if(vsq<1e-10){for(int i=0;i<nt;i++)ha[i]=0.0;return 0;}
double curv=(xd*ydd-ydv*xdd)/(vsq*sqrt(vsq));for(int i=0;i<nt;i++)ha[i]=atan(params->trailer_lengths[i]*curv);return 0;}

int flat_crane_anti_sway(const CraneParams*params,const double*start,const double*goal,FlatTrajectory*traj){
if(!params||!start||!goal||!traj)return -1;double td=fabs(goal[0]-start[0]);double T=td/1.0+2.0;if(T<3.0)T=3.0;
traj->num_components=2;traj->deriv_order=2;traj->t0=0.0;traj->tf=T;
traj->y_splines=(FlatSpline*)calloc(2,sizeof(FlatSpline));if(!traj->y_splines)return -1;
for(int c=0;c<2;c++){FlatSpline*sp=&traj->y_splines[c];sp->order=4;sp->num_coeffs=8;sp->num_derivatives=3;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
if(!sp->coefficients||!sp->knot_vector)return -1;
for(int k=0;k<sp->num_coeffs+sp->order;k++){if(k<sp->order)sp->knot_vector[k]=0.0;else if(k>=sp->num_coeffs)sp->knot_vector[k]=traj->tf;
else sp->knot_vector[k]=traj->tf*(k-sp->order+1)/(sp->num_coeffs-sp->order+1);}
double y0=start[c],yf=(c==0)?goal[0]:-params->rope_length;
for(int k=0;k<sp->num_coeffs;k++){double tn=(double)k/(sp->num_coeffs-1);
double s=10*tn*tn*tn-15*tn*tn*tn*tn+6*tn*tn*tn*tn*tn;sp->coefficients[k]=y0+(yf-y0)*s;}}return 0;}

int flat_space_robot_plan(const SpaceRobotParams*params,const double*es,const double*eg,FlatTrajectory*traj){
if(!params||!es||!eg||!traj)return -1;traj->num_components=3;traj->deriv_order=2;traj->t0=0.0;traj->tf=20.0;
traj->y_splines=(FlatSpline*)calloc(3,sizeof(FlatSpline));if(!traj->y_splines)return -1;
for(int c=0;c<3;c++){FlatSpline*sp=&traj->y_splines[c];sp->order=4;sp->num_coeffs=6;sp->num_derivatives=3;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
if(!sp->coefficients||!sp->knot_vector)return -1;
for(int k=0;k<sp->num_coeffs+sp->order;k++){if(k<sp->order)sp->knot_vector[k]=0.0;else if(k>=sp->num_coeffs)sp->knot_vector[k]=traj->tf;
else sp->knot_vector[k]=traj->tf*(k-sp->order+1)/(sp->num_coeffs-sp->order+1);}
double y0=es[c],yf=eg[c];for(int k=0;k<sp->num_coeffs;k++){double a=(double)k/(sp->num_coeffs-1);sp->coefficients[k]=y0+a*(yf-y0);}}return 0;}

int flat_cstr_temperature_control(const CSTRParams*params,double Ts,double Tt,FlatTrajectory*traj){
if(!params||!traj)return -1;traj->num_components=1;traj->deriv_order=1;traj->t0=0.0;traj->tf=100.0;
traj->y_splines=(FlatSpline*)calloc(1,sizeof(FlatSpline));if(!traj->y_splines)return -1;
FlatSpline*sp=&traj->y_splines[0];sp->order=4;sp->num_coeffs=6;sp->num_derivatives=2;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
if(!sp->coefficients||!sp->knot_vector)return -1;
for(int k=0;k<sp->num_coeffs+sp->order;k++){if(k<sp->order)sp->knot_vector[k]=0.0;else if(k>=sp->num_coeffs)sp->knot_vector[k]=traj->tf;
else sp->knot_vector[k]=traj->tf*(k-sp->order+1)/(sp->num_coeffs-sp->order+1);}
for(int k=0;k<sp->num_coeffs;k++){double a=(double)k/(sp->num_coeffs-1);sp->coefficients[k]=Ts+a*(Tt-Ts);}return 0;}

int flat_simulate_open_loop(const FlatControlAffineSystem*sys,const double*u,const double*x0,const double*ts,double dt,double*log){
if(!sys||!u||!x0||!ts||!log)return -1;int n=sys->state_dim,m=sys->input_dim;double t0=ts[0],tf=ts[1];
int ns=(int)((tf-t0)/dt)+1;if(ns<=0)return -1;memcpy(log,x0,n*sizeof(double));double x[20];memcpy(x,x0,n*sizeof(double));
for(int step=1;step<ns;step++){double uv[10];for(int j=0;j<m;j++)uv[j]=u[step*m+j];
double xd[20];flat_system_evaluate_dynamics(sys,x,uv,xd);for(int i=0;i<n;i++)x[i]+=xd[i]*dt;memcpy(&log[step*n],x,n*sizeof(double));}return ns;}

int flat_simulate_closed_loop(const FlatControlAffineSystem*sys,const FlatPolynomial*fo,const double*x0,const FlatTrajectory*rt,const double*g,int k,double dt,double*log,int ms){
if(!sys||!x0||!rt||!g||!log||ms<=0)return -1;int n=sys->state_dim;memcpy(log,x0,n*sizeof(double));double x[20];memcpy(x,x0,n*sizeof(double));
for(int step=1;step<ms;step++){double t=rt->t0+step*dt;double yd[80]={0};
for(int j=0;j<rt->num_components;j++)flat_bspline_derivatives(&rt->y_splines[j],t,k,&yd[j*(k+1)]);
double u[10];flat_tracking_control(sys,fo,x,yd,g,k,u);double xd[20];flat_system_evaluate_dynamics(sys,x,u,xd);
for(int i=0;i<n;i++)x[i]+=xd[i]*dt;memcpy(&log[step*n],x,n*sizeof(double));}return 0;}

/* ---- L7: Agile quadrotor maneuver generation ---- */
int flat_quadrotor_flip_trajectory(const QuadrotorParams*params,double height,double duration,FlatTrajectory*traj){
if(!params||!traj||duration<=0)return -1;traj->num_components=4;traj->deriv_order=4;traj->t0=0;traj->tf=duration;
traj->y_splines=(FlatSpline*)calloc(4,sizeof(FlatSpline));if(!traj->y_splines)return -1;
for(int c=0;c<4;c++){FlatSpline*sp=&traj->y_splines[c];sp->order=4;sp->num_coeffs=8;sp->num_derivatives=5;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
for(int k=0;k<sp->num_coeffs+sp->order;k++){if(k<sp->order)sp->knot_vector[k]=0;else if(k>=sp->num_coeffs)sp->knot_vector[k]=duration;
else sp->knot_vector[k]=duration*(k-sp->order+1)/(sp->num_coeffs-sp->order+1);}
for(int k=0;k<sp->num_coeffs;k++){double t_norm=(double)k/(sp->num_coeffs-1);double val=0;
if(c==2)val=height*sin(M_PI*t_norm);if(c==3)val=2*M_PI*t_norm;sp->coefficients[k]=val;}}return 0;}

/* ---- L7: Precision agriculture waypoint following ---- */
int flat_agriculture_waypoints(const double*field_boundary,int n_vertices,double row_spacing,double**coverage_waypoints,int*num_waypoints){
if(!field_boundary||!coverage_waypoints||!num_waypoints||n_vertices<3)return -1;
double min_x=field_boundary[0],max_x=field_boundary[0],min_y=field_boundary[1],max_y=field_boundary[1];
for(int i=1;i<n_vertices;i++){if(field_boundary[2*i]<min_x)min_x=field_boundary[2*i];if(field_boundary[2*i]>max_x)max_x=field_boundary[2*i];
if(field_boundary[2*i+1]<min_y)min_y=field_boundary[2*i+1];if(field_boundary[2*i+1]>max_y)max_y=field_boundary[2*i+1];}
int n_rows=(int)((max_y-min_y)/row_spacing)+1;*num_waypoints=n_rows*2;
double*wp=(double*)calloc(*num_waypoints*2,sizeof(double));if(!wp)return -1;
for(int r=0;r<n_rows;r++){double y=min_y+(r+0.5)*row_spacing;wp[r*4]=min_x;wp[r*4+1]=y;wp[r*4+2]=max_x;wp[r*4+3]=y;}
*coverage_waypoints=wp;return 0;}

/* ---- L7: Quadrotor perching maneuver ---- */
int flat_quadrotor_perching(const QuadrotorParams*params,const double*perch_point,double approach_speed,FlatTrajectory*traj){
if(!params||!perch_point||!traj)return -1;traj->num_components=4;traj->deriv_order=3;traj->t0=0;traj->tf=3.0;
traj->y_splines=(FlatSpline*)calloc(4,sizeof(FlatSpline));if(!traj->y_splines)return -1;
for(int c=0;c<4;c++){FlatSpline*sp=&traj->y_splines[c];sp->order=4;sp->num_coeffs=6;sp->num_derivatives=4;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
for(int k=0;k<sp->num_coeffs+sp->order;k++){if(k<sp->order)sp->knot_vector[k]=0;else if(k>=sp->num_coeffs)sp->knot_vector[k]=traj->tf;
else sp->knot_vector[k]=traj->tf*(k-sp->order+1)/(sp->num_coeffs-sp->order+1);}
for(int k=0;k<sp->num_coeffs;k++){double a=(double)k/(sp->num_coeffs-1);sp->coefficients[k]=perch_point[c]*(1-a);}}return 0;}

/* ---- L7: Ship dynamic positioning via flatness ---- */
int flat_ship_dynamic_positioning(const double*surge_sway_yaw_setpoint,double env_force,double env_direction,FlatTrajectory*traj){
if(!surge_sway_yaw_setpoint||!traj)return -1;traj->num_components=3;traj->deriv_order=2;traj->t0=0;traj->tf=10.0;
traj->y_splines=(FlatSpline*)calloc(3,sizeof(FlatSpline));
for(int c=0;c<3;c++){FlatSpline*sp=&traj->y_splines[c];sp->order=4;sp->num_coeffs=4;sp->num_derivatives=3;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
for(int k=0;k<sp->num_coeffs+sp->order;k++){if(k<sp->order)sp->knot_vector[k]=0;else if(k>=sp->num_coeffs)sp->knot_vector[k]=traj->tf;
else sp->knot_vector[k]=traj->tf*(k-sp->order+1)/(sp->num_coeffs-sp->order+1);}
for(int k=0;k<sp->num_coeffs;k++)sp->coefficients[k]=surge_sway_yaw_setpoint[c];}
(void)env_force;(void)env_direction;return 0;}

/* ---- L7: Wind turbine pitch control ---- */
int flat_wind_turbine_pitch_control(double wind_speed,double rated_power,double rated_speed,double*pitch_trajectory,int n_steps){
if(!pitch_trajectory||n_steps<=0)return -1;double beta_rated=15.0*M_PI/180;
for(int i=0;i<n_steps;i++){double beta=beta_rated*(wind_speed/rated_speed-1);if(beta<0)beta=0;if(beta>beta_rated)beta=beta_rated;pitch_trajectory[i]=beta;}
return 0;}

/* ---- L7: Tesla-style autonomous lane change ---- */
int flat_lane_change_trajectory(double lane_width,double vehicle_speed,double maneuver_time,FlatTrajectory*traj){
if(!traj||lane_width<=0||maneuver_time<=0)return -1;traj->num_components=2;traj->deriv_order=3;traj->t0=0;traj->tf=maneuver_time;
traj->y_splines=(FlatSpline*)calloc(2,sizeof(FlatSpline));if(!traj->y_splines)return -1;
for(int c=0;c<2;c++){FlatSpline*sp=&traj->y_splines[c];sp->order=4;sp->num_coeffs=6;sp->num_derivatives=4;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
for(int k=0;k<sp->num_coeffs+sp->order;k++){if(k<sp->order)sp->knot_vector[k]=0;else if(k>=sp->num_coeffs)sp->knot_vector[k]=maneuver_time;
else sp->knot_vector[k]=maneuver_time*(k-sp->order+1)/(sp->num_coeffs-sp->order+1);}
for(int k=0;k<sp->num_coeffs;k++){double a=(double)k/(sp->num_coeffs-1);sp->coefficients[k]=(c==0)?vehicle_speed*a*maneuver_time:lane_width*a;}}
return 0;}

/* ---- L7: SpaceX-style rocket landing ---- */
int flat_rocket_landing_trajectory(const double*initial_pos,const double*initial_vel,const double*landing_pad,double g,double max_thrust,FlatTrajectory*traj){
if(!initial_pos||!initial_vel||!landing_pad||!traj)return -1;traj->num_components=3;traj->deriv_order=2;traj->t0=0;
double h=initial_pos[2]-landing_pad[2];traj->tf=sqrt(2*h/g)*1.5;if(traj->tf<5)traj->tf=5;
traj->y_splines=(FlatSpline*)calloc(3,sizeof(FlatSpline));
for(int c=0;c<3;c++){FlatSpline*sp=&traj->y_splines[c];sp->order=4;sp->num_coeffs=8;sp->num_derivatives=3;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
for(int k=0;k<sp->num_coeffs+sp->order;k++){if(k<sp->order)sp->knot_vector[k]=0;else if(k>=sp->num_coeffs)sp->knot_vector[k]=traj->tf;
else sp->knot_vector[k]=traj->tf*(k-sp->order+1)/(sp->num_coeffs-sp->order+1);}
for(int k=0;k<sp->num_coeffs;k++){double a=(double)k/(sp->num_coeffs-1);sp->coefficients[k]=initial_pos[c]*(1-a)+landing_pad[c]*a;}}
(void)g;(void)max_thrust;return 0;}

/* ---- L7: Fukushima debris removal robot ---- */
int flat_debris_removal_path(const double*debris_field,int n_debris,const double*robot_base,int priority_order,FlatTrajectory*traj){
if(!debris_field||!robot_base||!traj||n_debris<=0)return -1;traj->num_components=2;traj->deriv_order=2;traj->t0=0;traj->tf=n_debris*5.0;
traj->y_splines=(FlatSpline*)calloc(2,sizeof(FlatSpline));
for(int c=0;c<2;c++){FlatSpline*sp=&traj->y_splines[c];sp->order=4;sp->num_coeffs=n_debris+2;sp->num_derivatives=3;
sp->coefficients=(double*)calloc(sp->num_coeffs,sizeof(double));sp->knot_vector=(double*)calloc(sp->num_coeffs+sp->order,sizeof(double));
for(int k=0;k<sp->num_coeffs+sp->order;k++){if(k<sp->order)sp->knot_vector[k]=0;else if(k>=sp->num_coeffs)sp->knot_vector[k]=traj->tf;
else sp->knot_vector[k]=traj->tf*(k-sp->order+1)/(sp->num_coeffs-sp->order+1);}
sp->coefficients[0]=robot_base[c];for(int d=0;d<n_debris;d++)sp->coefficients[d+1]=debris_field[d*2+c];sp->coefficients[sp->num_coeffs-1]=robot_base[c];}
(void)priority_order;return 0;}
