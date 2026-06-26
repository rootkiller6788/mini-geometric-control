/* src/flatness_observer.c - L5 Observer design for flat systems */
#include "flatness_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ---- Luenberger observer for linearized flat system ---- */
int flat_luenberger_observer(const double*A,const double*B,const double*C,int n,int m,int p,double*L){
if(!A||!B||!C||!L||n<=0||m<=0||p<=0)return -1;
double*poles=(double*)calloc(n,sizeof(double));if(!poles)return -1;
for(int i=0;i<n;i++)poles[i]=-2.0-2.0*i;
double*Ao=(double*)calloc(n*n,sizeof(double));for(int i=0;i<n;i++)for(int j=0;j<n;j++)Ao[i*n+j]=A[j*n+i];
double*Co=(double*)calloc(n*p,sizeof(double));for(int i=0;i<p;i++)for(int j=0;j<n;j++)Co[j*p+i]=C[i*n+j];
double*K=(double*)calloc(n*p,sizeof(double));
int ctrl=flat_pbh_test(Ao,Co,n,p);if(!ctrl){free(poles);free(Ao);free(Co);free(K);return -1;}
for(int i=0;i<n;i++)for(int j=0;j<p;j++)L[i*p+j]=poles[0];
free(poles);free(Ao);free(Co);free(K);return 0;}

/* ---- Extended Kalman filter for flat systems ---- */
int flat_extended_kalman_filter(void(*f)(const double*,const double*,double*,void*),
void(*h)(const double*,double*,void*),void*ctx,int n,int m,int p,
double*x,double*P,const double*u,const double*z,double Q,double R,double dt){
if(!f||!h||!x||!P||!u||!z||n<=0)return -1;
double*F=(double*)calloc(n*n,sizeof(double));double*H=(double*)calloc(p*n,sizeof(double));
double*xp=(double*)calloc(n,sizeof(double));double*Pp=(double*)calloc(n*n,sizeof(double));
double*y=(double*)calloc(p,sizeof(double));double*K=(double*)calloc(n*p,sizeof(double));
double*S=(double*)calloc(p*p,sizeof(double));double*I_KH=(double*)calloc(n*n,sizeof(double));
f(x,u,xp,ctx);for(int i=0;i<n;i++)xp[i]=x[i]+xp[i]*dt;
for(int i=0;i<n;i++)Pp[i*n+i]=P[i*n+i]+Q;
h(xp,y,ctx);for(int i=0;i<p;i++)y[i]=z[i]-y[i];
for(int i=0;i<n;i++)for(int j=0;j<p;j++)K[i*p+j]=Pp[i*n+i]/(Pp[i*n+i]+R);
for(int i=0;i<n;i++){double val=0.0;for(int j=0;j<p;j++)val+=K[i*p+j]*y[j];x[i]=xp[i]+val;}
for(int i=0;i<n;i++)for(int j=0;j<n;j++){I_KH[i*n+j]=(i==j)?1.0:0.0;for(int k=0;k<p;k++)I_KH[i*n+j]-=K[i*p+k]*H[k*n+j];}
for(int i=0;i<n;i++)for(int j=0;j<n;j++){double val=0.0;for(int k=0;k<n;k++)val+=I_KH[i*n+k]*Pp[k*n+j];P[i*n+j]=val;}
free(F);free(H);free(xp);free(Pp);free(y);free(K);free(S);free(I_KH);return 0;}
