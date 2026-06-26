#define _USE_MATH_DEFINES
#include "geo_optctrl_core.h"
#include "lie_group_methods.h"
#include "symplectic_geometry.h"
#include "pontryagin_maximum.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ATT_MAX_STEPS 1000

typedef struct {
    double I[3];
    double R_init[3][3];
    double Omega_init[3];
    double R_target[3][3];
    double Omega_target[3];
    double tau_max[3];
    double T;
    int    n_steps;
} SO3AttitudeProblem;

double geo_attitude_error(const double R1[3][3], const double R2[3][3]) {
    double Re[3][3], RT[3][3];
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) RT[i][j] = R1[j][i];
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
        Re[i][j] = 0.0;
        for (int k=0;k<3;k++) Re[i][j] += RT[i][k] * R2[k][j];
    }
    double omega[3];
    if (geo_so3_log(Re, omega) != 0) return M_PI;
    return sqrt(omega[0]*omega[0] + omega[1]*omega[1] + omega[2]*omega[2]);
}

int geo_attitude_energy_optimal(const SO3AttitudeProblem *prob,
                                 double *R_traj, double *Omega_traj,
                                 double *tau_traj) {
    if (!prob || !R_traj || !Omega_traj || !tau_traj) return -1;
    if (prob->n_steps <= 0 || prob->n_steps > ATT_MAX_STEPS) return -1;
    int N = prob->n_steps;
    double dt = prob->T / (double)N;
    double I[3]; memcpy(I, prob->I, 3*sizeof(double));
    double R_curr[3][3], Omega_curr[3];
    memcpy(R_curr, prob->R_init, 9*sizeof(double));
    memcpy(Omega_curr, prob->Omega_init, 3*sizeof(double));
    for (int k = 0; k < N; k++) {
        memcpy(&R_traj[k*9], R_curr, 9*sizeof(double));
        memcpy(&Omega_traj[k*3], Omega_curr, 3*sizeof(double));
        double R_err[3][3], omega_err[3], RT[3][3];
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) RT[i][j] = R_curr[j][i];
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
            R_err[i][j]=0.0;
            for (int kk=0;kk<3;kk++) R_err[i][j]+=RT[i][kk]*prob->R_target[kk][j];
        }
        geo_so3_log(R_err, omega_err);
        double Kp=10.0, Kd=5.0, tau[3], IO[3], gyro[3];
        IO[0]=I[0]*Omega_curr[0]; IO[1]=I[1]*Omega_curr[1]; IO[2]=I[2]*Omega_curr[2];
        gyro[0]=Omega_curr[1]*IO[2]-Omega_curr[2]*IO[1];
        gyro[1]=Omega_curr[2]*IO[0]-Omega_curr[0]*IO[2];
        gyro[2]=Omega_curr[0]*IO[1]-Omega_curr[1]*IO[0];
        for (int i=0;i<3;i++) {
            tau[i] = -Kp*omega_err[i] - Kd*Omega_curr[i] + gyro[i];
            if (tau[i]>prob->tau_max[i]) tau[i]=prob->tau_max[i];
            if (tau[i]<-prob->tau_max[i]) tau[i]=-prob->tau_max[i];
        }
        memcpy(&tau_traj[k*3], tau, 3*sizeof(double));
        double dOmega[3];
        geo_euler_poincare_rigid_body(I, Omega_curr, tau, dOmega);
        for (int i=0;i<3;i++) Omega_curr[i] += dt*dOmega[i];
        double dR[3][3];
        geo_so3_kinematics(R_curr, Omega_curr, dR);
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) R_curr[i][j]+=dt*dR[i][j];
        double r0[3]={R_curr[0][0],R_curr[1][0],R_curr[2][0]};
        double r1[3]={R_curr[0][1],R_curr[1][1],R_curr[2][1]};
        double r2[3]={R_curr[0][2],R_curr[1][2],R_curr[2][2]};
        double nr0=sqrt(r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2]);
        for (int i=0;i<3;i++) r0[i]/=nr0;
        double dot01=r0[0]*r1[0]+r0[1]*r1[1]+r0[2]*r1[2];
        for (int i=0;i<3;i++) r1[i]-=dot01*r0[i];
        double nr1=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]);
        for (int i=0;i<3;i++) r1[i]/=nr1;
        r2[0]=r0[1]*r1[2]-r0[2]*r1[1];
        r2[1]=r0[2]*r1[0]-r0[0]*r1[2];
        r2[2]=r0[0]*r1[1]-r0[1]*r1[0];
        for (int i=0;i<3;i++){R_curr[i][0]=r0[i];R_curr[i][1]=r1[i];R_curr[i][2]=r2[i];}
    }
    memcpy(&R_traj[N*9], R_curr, 9*sizeof(double));
    memcpy(&Omega_traj[N*3], Omega_curr, 3*sizeof(double));
    return (int)(geo_attitude_error(R_curr, prob->R_target)*1000.0);
}

double geo_attitude_min_time(const SO3AttitudeProblem *prob,
                              double *R_traj, double *Omega_traj,
                              double *tau_traj, double tol) {
    if (!prob) return -1.0;
    double T_low=0.1, T_high=100.0;
    SO3AttitudeProblem p = *prob;
    for (int iter=0;iter<30;iter++) {
        double T_mid=0.5*(T_low+T_high);
        p.T=T_mid; p.n_steps=(int)(T_mid*100);
        if (p.n_steps<10) p.n_steps=10;
        if (p.n_steps>ATT_MAX_STEPS) p.n_steps=ATT_MAX_STEPS;
        int err = geo_attitude_energy_optimal(&p,R_traj,Omega_traj,tau_traj);
        if (err<100) T_high=T_mid; else T_low=T_mid;
        if (T_high-T_low<tol) break;
    }
    return T_high;
}

void geo_apollo_attitude_maneuver(double roll_deg, double pitch_deg, double yaw_deg,
                                  const double R_init[3][3], double R_final[3][3]) {
    (void)R_init;
    double cr=cos(roll_deg*M_PI/180.0), sr=sin(roll_deg*M_PI/180.0);
    double cp=cos(pitch_deg*M_PI/180.0), sp=sin(pitch_deg*M_PI/180.0);
    double cy=cos(yaw_deg*M_PI/180.0), sy=sin(yaw_deg*M_PI/180.0);
    R_final[0][0]=cy*cp; R_final[0][1]=cy*sp*sr-sy*cr; R_final[0][2]=cy*sp*cr+sy*sr;
    R_final[1][0]=sy*cp; R_final[1][1]=sy*sp*sr+cy*cr; R_final[1][2]=sy*sp*cr-cy*sr;
    R_final[2][0]=-sp;   R_final[2][1]=cp*sr;           R_final[2][2]=cp*cr;
}

void geo_reaction_wheel_desaturation(const double h_rw[3],
                                      const double b_field[3],
                                      double m_dipole[3]) {
    double b2=b_field[0]*b_field[0]+b_field[1]*b_field[1]+b_field[2]*b_field[2];
    if (b2<1e-12) { m_dipole[0]=m_dipole[1]=m_dipole[2]=0.0; return; }
    m_dipole[0]=(b_field[1]*h_rw[2]-b_field[2]*h_rw[1])/b2;
    m_dipole[1]=(b_field[2]*h_rw[0]-b_field[0]*h_rw[2])/b2;
    m_dipole[2]=(b_field[0]*h_rw[1]-b_field[1]*h_rw[0])/b2;
}