/* flatness_models.c - L6 Additional canonical flat system models */
#include "flatness_core.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* ---- Inverted pendulum on a cart (4 states, 1 input) ---- */
int flat_make_inverted_pendulum(FlatControlAffineSystem*sys){
if(!sys)return -1;int n=4,m=1;sys->state_dim=n;sys->input_dim=m;sys->output_dim=m;
int p0[20]={0};flat_poly_add_term(&sys->drift.components[3],-9.81,p0);
flat_poly_add_term(&sys->control_fields[0].components[2],1.0,p0);
int px[20]={0};px[0]=1;flat_poly_add_term(&sys->output_map[0],1.0,px);
return 0;}

/* ---- PVTOL aircraft (6 states, 2 inputs) ---- */
int flat_make_pvtol_system(FlatControlAffineSystem*sys){
if(!sys)return -1;int n=6,m=2;sys->state_dim=n;sys->input_dim=m;sys->output_dim=m;
int p0[20]={0};flat_poly_add_term(&sys->drift.components[5],9.81,p0);
flat_poly_add_term(&sys->control_fields[0].components[3],1.0,p0);
flat_poly_add_term(&sys->control_fields[0].components[4],1.0,p0);
flat_poly_add_term(&sys->control_fields[1].components[5],1.0,p0);
int pv[20]={0};pv[0]=1;flat_poly_add_term(&sys->output_map[0],1.0,pv);
pv[0]=0;pv[1]=1;flat_poly_add_term(&sys->output_map[1],1.0,pv);
return 0;}

/* ---- Ball and beam system (4 states, 1 input) ---- */
int flat_make_ball_beam_system(FlatControlAffineSystem*sys){
if(!sys)return -1;int n=4,m=1;sys->state_dim=n;sys->input_dim=m;sys->output_dim=m;
int p0[20]={0};flat_poly_add_term(&sys->drift.components[2],-4.905,p0);
flat_poly_add_term(&sys->control_fields[0].components[3],1.0,p0);
flat_poly_add_term(&sys->output_map[0],1.0,p0);return 0;}

/* ---- Induction motor (5 states, 2 inputs) ---- */
int flat_make_induction_motor_system(FlatControlAffineSystem*sys){
if(!sys)return -1;int n=5,m=2;sys->state_dim=n;sys->input_dim=m;sys->output_dim=m;
int p0[20]={0};flat_poly_add_term(&sys->control_fields[0].components[3],1.0,p0);
flat_poly_add_term(&sys->control_fields[1].components[4],1.0,p0);
int po[20]={0};flat_poly_add_term(&sys->output_map[0],1.0,po);
po[4]=1;po[0]=0;flat_poly_add_term(&sys->output_map[1],1.0,po);
return 0;}

/* ---- Magnetic levitation system (3 states, 1 input) ---- */
int flat_make_maglev_system(FlatControlAffineSystem*sys){
if(!sys)return -1;int n=3,m=1;sys->state_dim=n;sys->input_dim=m;sys->output_dim=m;
int p0[20]={0};flat_poly_add_term(&sys->drift.components[2],9.81,p0);
flat_poly_add_term(&sys->control_fields[0].components[2],1.0,p0);
int po[20]={0};po[0]=1;flat_poly_add_term(&sys->output_map[0],1.0,po);
return 0;}

/* ---- DC-DC boost converter (2 states, 1 input) ---- */
int flat_make_boost_converter(FlatControlAffineSystem*sys){
if(!sys)return -1;int n=2,m=1;sys->state_dim=n;sys->input_dim=m;sys->output_dim=m;
int p0[20]={0};flat_poly_add_term(&sys->control_fields[0].components[0],1.0,p0);
flat_poly_add_term(&sys->output_map[0],1.0,p0);return 0;}

/* ---- Underwater vehicle (6 DOF, 6 inputs) ---- */
int flat_make_auv_system(FlatControlAffineSystem*sys){
if(!sys)return -1;int n=12,m=6;sys->state_dim=n;sys->input_dim=m;sys->output_dim=m;
int p0[20]={0};
for(int i=0;i<6;i++){flat_poly_add_term(&sys->drift.components[i],1.0,p0);
flat_poly_add_term(&sys->control_fields[i].components[6+i],1.0,p0);
int po[20]={0};po[i]=1;flat_poly_add_term(&sys->output_map[i],1.0,po);}
return 0;}

/* ---- Satellite attitude dynamics (6 states, 3 inputs) ---- */
int flat_make_satellite_attitude(FlatControlAffineSystem*sys){
if(!sys)return -1;int n=6,m=3;sys->state_dim=n;sys->input_dim=m;sys->output_dim=m;
int p0[20]={0};
for(int i=0;i<3;i++){flat_poly_add_term(&sys->drift.components[i],1.0,p0);
flat_poly_add_term(&sys->control_fields[i].components[3+i],1.0,p0);
int po[20]={0};po[i]=1;flat_poly_add_term(&sys->output_map[i],1.0,po);}
return 0;}

/* ---- Wind turbine (3 states, 2 inputs) ---- */
int flat_make_wind_turbine(FlatControlAffineSystem*sys){
if(!sys)return -1;int n=3,m=2;sys->state_dim=n;sys->input_dim=m;sys->output_dim=m;
int p0[20]={0};flat_poly_add_term(&sys->control_fields[0].components[1],1.0,p0);
flat_poly_add_term(&sys->control_fields[1].components[2],1.0,p0);
int po[20]={0};po[1]=1;flat_poly_add_term(&sys->output_map[0],1.0,po);
po[1]=0;po[2]=1;flat_poly_add_term(&sys->output_map[1],1.0,po);
return 0;}

/* ---- Flexible joint robot (4 states, 1 input) ---- */
int flat_make_flexible_joint(FlatControlAffineSystem*sys){
if(!sys)return -1;int n=4,m=1;sys->state_dim=n;sys->input_dim=m;sys->output_dim=m;
int p0[20]={0};flat_poly_add_term(&sys->drift.components[2],-1.0,p0);
flat_poly_add_term(&sys->control_fields[0].components[3],1.0,p0);
int po[20]={0};po[0]=1;flat_poly_add_term(&sys->output_map[0],1.0,po);
return 0;}