#include "flatness_core.h"
#include "trajectory_generation.h"
#include "flatness_applications.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void) {
  printf("=== mini-differential-flatness Examples ===\n\n");

  printf("Example 1: Quadrotor Minimum-Snap Trajectory\n");
  QuadrotorParams qp={1,0.25,0.01,0.01,0.02,1,0.1};
  double wp[]={0,0,0,0, 1,0,0,0, 1,1,0,0, 0,1,0,0, 0,0,1,0};
  double st[]={2,2,2,2}; FlatTrajectory qt;
  if(flat_quadrotor_min_snap(&qp,wp,5,st,&qt)==0){
    printf("  Generated %d-component trajectory [%.1f, %.1f]s\n",qt.num_components,qt.t0,qt.tf);
    printf("  Flat outputs: x, y, z, yaw\n");
    flat_trajectory_free(&qt); }

  printf("\nExample 2: 3-Trailer Vehicle Parking\n");
  TrailerParams tp={1,0.5,{1,1,1},0.5,2,3};
  double s[]={0,0},g[]={5,3}; FlatTrajectory tt;
  if(flat_trailer_parking(&tp,s,g,&tt)==0){
    printf("  Parking trajectory [%.1f, %.1f]s\n",tt.t0,tt.tf);
    flat_trajectory_free(&tt); }

  printf("\nExample 3: Overhead Crane Anti-Sway Transport\n");
  CraneParams cp={5,2,3,9.81,100};
  double cs[]={0,-3},cg[]={10,-3}; FlatTrajectory ct;
  if(flat_crane_anti_sway(&cp,cs,cg,&ct)==0){
    printf("  Anti-sway trajectory [%.1f, %.1f]s\n",ct.t0,ct.tf);
    flat_trajectory_free(&ct); }

  printf("\nExample 4: Space Robot End-Effector Planning\n");
  SpaceRobotParams srp={500,100,{10,10,10},{1,1,1},{1,1,1},3};
  double es[]={0,0,0},eg[]={2,1,0.5}; FlatTrajectory srt;
  if(flat_space_robot_plan(&srp,es,eg,&srt)==0){
    printf("  EE trajectory [%.1f, %.1f]s\n",srt.t0,srt.tf);
    flat_trajectory_free(&srt); }

  printf("\nExample 5: Flatness Analysis of Quadrotor System\n");
  FlatControlAffineSystem sys; flat_make_quadrotor_system(&sys);
  printf("  State dim: %d, Input dim: %d\n",sys.state_dim,sys.input_dim);
  printf("  Is flat: %s\n",flat_system_is_flat(&sys)?"yes":"no");
  int d=flat_compute_defect(&sys);
  printf("  Defect: %d (0 = statically linearizable)\n",d);
  int rank; flat_controllability_matrix_rank(&sys,&rank);
  printf("  Controllability rank: %d / %d\n",rank,sys.state_dim);

  printf("\n=== All examples completed ===\n");
  return 0;
}
