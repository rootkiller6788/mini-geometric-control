#ifndef FLATNESS_APPLICATIONS_H
#define FLATNESS_APPLICATIONS_H
#include "flatness_core.h"

typedef struct{double mass,arm_length,Ixx,Iyy,Izz,kf,km;}QuadrotorParams;
int flat_quadrotor_min_snap(const QuadrotorParams*params,const double*waypoints,int n_wp,const double*segment_times,FlatTrajectory*traj);
int flat_quadrotor_rotor_speeds(const QuadrotorParams*params,const double*u,double*omega);

typedef struct{double car_length,hitch_length,trailer_lengths[10],max_steering_angle,max_velocity;int n_trailers;}TrailerParams;
int flat_trailer_parking(const TrailerParams*params,const double*start,const double*goal,FlatTrajectory*traj);
int flat_trailer_hitch_angles(const TrailerParams*params,const FlatTrajectory*traj,double t,double*hitch_angles);

typedef struct{double cart_mass,load_mass,rope_length,gravity,max_force;}CraneParams;
int flat_crane_anti_sway(const CraneParams*params,const double*start,const double*goal,FlatTrajectory*traj);

typedef struct{double base_mass,base_inertia,link_masses[10],link_lengths[10],link_inertias[10];int n_links;}SpaceRobotParams;
int flat_space_robot_plan(const SpaceRobotParams*params,const double*ee_start,const double*ee_goal,FlatTrajectory*traj);

typedef struct{double Da,beta,gamma,delta_H,rho_Cp,U,A_ht;}CSTRParams;
int flat_cstr_temperature_control(const CSTRParams*params,double T_start,double T_target,FlatTrajectory*traj);

int flat_simulate_open_loop(const FlatControlAffineSystem*sys,const double*u_open_loop,const double*x0,const double*t_span,double dt,double*log);
int flat_simulate_closed_loop(const FlatControlAffineSystem*sys,const FlatPolynomial*flat_outputs,const double*x0,const FlatTrajectory*ref_traj,const double*gains,int kappa,double dt,double*log,int max_steps);
#endif
