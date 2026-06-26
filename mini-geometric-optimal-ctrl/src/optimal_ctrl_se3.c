#include "geo_optctrl_core.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double mass;
    double I[3];
    double R_init[3][3];
    double pos_init[3];
    double vel_init[3];
    double Omega_init[3];
    double R_target[3][3];
    double pos_target[3];
    double T;
    int n_steps;
} SE3OptimalProblem;

double geo_se3_pose_error(const double R1[3][3], const double p1[3], const double R2[3][3], const double p2[3]) {
    extern double geo_attitude_error(const double R1[3][3], const double R2[3][3]);
    double re = geo_attitude_error(R1, R2);
    double pe = sqrt((p2[0]-p1[0])*(p2[0]-p1[0])+(p2[1]-p1[1])*(p2[1]-p1[1])+(p2[2]-p1[2])*(p2[2]-p1[2]));
    return re + 0.1*pe;
}

void geo_quadrotor_motor_mixing(const double R[3][3], const double tau[3], double ft, double L, double kt, double w[4]) {
    double tr=tau[0], tp=tau[1], ty=tau[2];
    w[0]=ft/4.0-tp/(2.0*L)-ty/(4.0*kt);
    w[1]=ft/4.0-tr/(2.0*L)+ty/(4.0*kt);
    w[2]=ft/4.0+tp/(2.0*L)-ty/(4.0*kt);
    w[3]=ft/4.0+tr/(2.0*L)+ty/(4.0*kt);
    for (int i=0;i<4;i++) if(w[i]<0.0)w[i]=0.0;
}
