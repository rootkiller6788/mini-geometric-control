#include "flatness_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* flatness_numerical.c - L5 Numerical methods for flatness analysis */

/* QR decomposition via Householder reflections */
int flat_qr_decomposition(double*A,int m,int n,double*Q,double*R){
if(!A||!Q||!R||m<=0||n<=0)return -1;int i,j,k;for(i=0;i<m*n;i++)R[i]=A[i];
for(i=0;i<m*m;i++)Q[i]=(i%(m+1)==0)?1.0:0.0;
for(k=0;k<n&&k<m;k++){double norm=0.0;for(i=k;i<m;i++)norm+=R[i*n+k]*R[i*n+k];norm=sqrt(norm);if(norm<1e-15)continue;
double beta=(R[k*n+k]>0)?-norm:norm,sigma=R[k*n+k]-beta;double v[100]={0};for(i=k;i<m;i++)v[i]=R[i*n+k]/sigma;v[k]=1.0;
double tau=-sigma/beta;for(j=k;j<n;j++){double s=0.0;for(i=k;i<m;i++)s+=v[i]*R[i*n+j];s*=tau;for(i=k;i<m;i++)R[i*n+j]-=s*v[i];}
for(j=0;j<m;j++){double s=0.0;for(i=k;i<m;i++)s+=v[i]*Q[j*m+i];s*=tau;for(i=0;i<m;i++)Q[j*m+i]-=s*v[i];}}return 0;}

/* SVD via Golub-Kahan */
int flat_svd_decomposition(const double*A,int m,int n,double*U,double*S,double*Vt){
if(!A||!U||!S||!Vt||m<=0||n<=0)return -1;int md=(m<n)?m:n,i,j;
for(i=0;i<m*m;i++)U[i]=(i%(m+1)==0)?1.0:0.0;for(i=0;i<n*n;i++)Vt[i]=(i%(n+1)==0)?1.0:0.0;
double*B=(double*)calloc(m*n,sizeof(double));if(!B)return -1;memcpy(B,A,m*n*sizeof(double));
for(int iter=0;iter<md*10;iter++){
for(i=0;i<md-1;i++){double a=B[i*n+i],b=B[(i+1)*n+i],r=sqrt(a*a+b*b);if(r<1e-15)continue;
double c=a/r,s=b/r;for(j=0;j<n;j++){double t1=B[i*n+j],t2=B[(i+1)*n+j];B[i*n+j]=c*t1+s*t2;B[(i+1)*n+j]=-s*t1+c*t2;}
for(j=0;j<m;j++){double t1=U[j*m+i],t2=U[j*m+i+1];U[j*m+i]=c*t1+s*t2;U[j*m+i+1]=-s*t1+c*t2;}}
for(i=0;i<md-1;i++){double a=B[i*n+i],b=B[i*n+i+1],r=sqrt(a*a+b*b);if(r<1e-15)continue;
double c=a/r,s=b/r;for(j=0;j<m;j++){double t1=B[j*n+i],t2=B[j*n+i+1];B[j*n+i]=c*t1+s*t2;B[j*n+i+1]=-s*t1+c*t2;}
for(j=0;j<n;j++){double t1=Vt[i*n+j],t2=Vt[(i+1)*n+j];Vt[i*n+j]=c*t1+s*t2;Vt[(i+1)*n+j]=-s*t1+c*t2;}}
}for(i=0;i<md;i++)S[i]=fabs(B[i*n+i]);free(B);return 0;}

/* Condition number */
double flat_condition_number(const double*A,int n){
if(!A||n<=0)return INFINITY;double S[100]={0},U[400]={0},Vt[400]={0};
flat_svd_decomposition(A,n,n,U,S,Vt);double sm=S[0],sM=S[0];
for(int i=1;i<n;i++){if(S[i]>sM)sM=S[i];if(S[i]<sm&&S[i]>1e-15)sm=S[i];}return(sm>1e-15)?sM/sm:INFINITY;}

/* Linear solve Ax=b via Gaussian elimination */
int flat_linear_solve(const double*A,const double*b,int n,double*x){
if(!A||!b||!x||n<=0)return -1;double*M=(double*)calloc(n*(n+1),sizeof(double));if(!M)return -1;
int i,j;for(i=0;i<n;i++){for(j=0;j<n;j++)M[i*(n+1)+j]=A[i*n+j];M[i*(n+1)+n]=b[i];}
for(int col=0;col<n;col++){int pv=col;for(i=col+1;i<n;i++)if(fabs(M[i*(n+1)+col])>fabs(M[pv*(n+1)+col]))pv=i;
if(fabs(M[pv*(n+1)+col])<1e-15){free(M);return -1;}
if(pv!=col)for(j=0;j<=n;j++){double t=M[col*(n+1)+j];M[col*(n+1)+j]=M[pv*(n+1)+j];M[pv*(n+1)+j]=t;}
double p=M[col*(n+1)+col];for(j=col;j<=n;j++)M[col*(n+1)+j]/=p;
for(i=0;i<n;i++){if(i==col)continue;double f=M[i*(n+1)+col];for(j=col;j<=n;j++)M[i*(n+1)+j]-=f*M[col*(n+1)+j];}}
for(i=0;i<n;i++)x[i]=M[i*(n+1)+n];free(M);return 0;}

/* Eigenvalues of 2x2 matrix */
int flat_eigenvalues_2x2(const double*A,double*lr,double*li){
if(!A||!lr||!li)return -1;double tr=A[0]+A[3],det=A[0]*A[3]-A[1]*A[2],disc=tr*tr-4*det;
if(disc>=0){lr[0]=(tr+sqrt(disc))/2;li[0]=0;lr[1]=(tr-sqrt(disc))/2;li[1]=0;}
else{lr[0]=tr/2;li[0]=sqrt(-disc)/2;lr[1]=tr/2;li[1]=-sqrt(-disc)/2;}return 0;}

/* Cholesky decomposition A = L*L^T for symmetric positive definite A */
int flat_cholesky_decomposition(const double*A,int n,double*L){
if(!A||!L||n<=0)return -1;memset(L,0,n*n*sizeof(double));
for(int i=0;i<n;i++){for(int j=0;j<=i;j++){double s=0.0;for(int k=0;k<j;k++)s+=L[i*n+k]*L[j*n+k];
if(i==j){double v=A[i*n+i]-s;if(v<=0)return -1;L[i*n+i]=sqrt(v);}else L[i*n+j]=(A[i*n+j]-s)/L[j*n+j];}}return 0;}

/* LDL^T decomposition for symmetric matrices */
int flat_ldlt_decomposition(const double*A,int n,double*L,double*D){
if(!A||!L||!D||n<=0)return -1;memset(L,0,n*n*sizeof(double));memset(D,0,n*sizeof(double));
for(int i=0;i<n;i++){double s=0.0;for(int k=0;k<i;k++)s+=L[i*n+k]*L[i*n+k]*D[k];D[i]=A[i*n+i]-s;
if(fabs(D[i])<1e-15)return -1;L[i*n+i]=1.0;
for(int j=i+1;j<n;j++){s=0.0;for(int k=0;k<i;k++)s+=L[j*n+k]*L[i*n+k]*D[k];L[j*n+i]=(A[j*n+i]-s)/D[i];}}return 0;}

/* Forward substitution: Ly = b */
int flat_forward_substitution(const double*L,const double*b,int n,double*y){
if(!L||!b||!y||n<=0)return -1;for(int i=0;i<n;i++){double s=0.0;for(int j=0;j<i;j++)s+=L[i*n+j]*y[j];y[i]=(b[i]-s)/L[i*n+i];}return 0;}

/* Backward substitution: Ux = y */
int flat_backward_substitution(const double*U,const double*y,int n,double*x){
if(!U||!y||!x||n<=0)return -1;for(int i=n-1;i>=0;i--){double s=0.0;for(int j=i+1;j<n;j++)s+=U[i*n+j]*x[j];x[i]=(y[i]-s)/U[i*n+i];}return 0;}

/* Newton method for nonlinear equations */
int flat_newton_solve(void(*F)(const double*,double*,void*),void(*JF)(const double*,double*,void*),
void*ctx,int n,double*x,int max_iter,double tol){
if(!F||!JF||!x||n<=0)return -1;double*dx=(double*)calloc(n,sizeof(double));double*Fv=(double*)calloc(n,sizeof(double));
double*J=(double*)calloc(n*n,sizeof(double));if(!dx||!Fv||!J){free(dx);free(Fv);free(J);return -1;}
for(int iter=0;iter<max_iter;iter++){F(x,Fv,ctx);double err=0.0;for(int i=0;i<n;i++)err+=Fv[i]*Fv[i];
if(sqrt(err)<tol){free(dx);free(Fv);free(J);return iter+1;}
JF(x,J,ctx);if(flat_linear_solve(J,Fv,n,dx)!=0){free(dx);free(Fv);free(J);return -1;}
for(int i=0;i<n;i++)x[i]-=dx[i];}free(dx);free(Fv);free(J);return -1;}

/* Cubic spline interpolation */
int flat_cubic_spline_interp(const double*x,const double*y,int n,double*a,double*b,double*c,double*d){
if(!x||!y||!a||!b||!c||!d||n<2)return -1;
double*h=(double*)calloc(n-1,sizeof(double));double*alpha=(double*)calloc(n-1,sizeof(double));
double*l=(double*)calloc(n,sizeof(double));double*mu=(double*)calloc(n,sizeof(double));double*z=(double*)calloc(n,sizeof(double));
for(int i=0;i<n-1;i++)h[i]=x[i+1]-x[i];
for(int i=1;i<n-1;i++)alpha[i]=(3/h[i])*(y[i+1]-y[i])-(3/h[i-1])*(y[i]-y[i-1]);
l[0]=1;mu[0]=0;z[0]=0;for(int i=1;i<n-1;i++){l[i]=2*(x[i+1]-x[i-1])-h[i-1]*mu[i-1];mu[i]=h[i]/l[i];z[i]=(alpha[i]-h[i-1]*z[i-1])/l[i];}
l[n-1]=1;z[n-1]=0;c[n-1]=0;
for(int j=n-2;j>=0;j--){c[j]=z[j]-mu[j]*c[j+1];b[j]=(y[j+1]-y[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;d[j]=(c[j+1]-c[j])/(3*h[j]);}
for(int i=0;i<n;i++)a[i]=y[i];free(h);free(alpha);free(l);free(mu);free(z);return 0;}

/* Pseudoinverse via SVD */
int flat_pseudo_inverse(const double*A,int m,int n,double*A_pinv){
if(!A||!A_pinv||m<=0||n<=0)return -1;double*U=(double*)calloc(m*m,sizeof(double));
double*S=(double*)calloc((m<n?m:n),sizeof(double));double*Vt=(double*)calloc(n*n,sizeof(double));
if(!U||!S||!Vt){free(U);free(S);free(Vt);return -1;}flat_svd_decomposition(A,m,n,U,S,Vt);
double*Si=(double*)calloc(n*m,sizeof(double));int md=(m<n)?m:n;
for(int i=0;i<md;i++)if(S[i]>1e-10)Si[i*n+i]=1.0/S[i];
double*V=(double*)calloc(n*n,sizeof(double));for(int i=0;i<n;i++)for(int j=0;j<n;j++)V[i*n+j]=Vt[j*n+i];
double*VS=(double*)calloc(n*m,sizeof(double));for(int i=0;i<n;i++)for(int j=0;j<m;j++){double s=0.0;for(int k=0;k<n;k++)s+=V[i*n+k]*Si[k*m+j];VS[i*m+j]=s;}
for(int i=0;i<n;i++)for(int j=0;j<m;j++){double s=0.0;for(int k=0;k<m;k++)s+=VS[i*m+k]*U[j*m+k];A_pinv[i*m+j]=s;}
free(U);free(S);free(Vt);free(Si);free(V);free(VS);return 0;}

/* Golden section search */
double flat_golden_section_search(double(*f)(double,void*),void*ctx,double a,double b,double tol){
double phi=(sqrt(5)-1)/2;double c=b-phi*(b-a),d=a+phi*(b-a);
while(fabs(b-a)>tol){if(f(c,ctx)<f(d,ctx)){b=d;d=c;c=b-phi*(b-a);}else{a=c;c=d;d=a+phi*(b-a);}}return(a+b)/2;}

/* BFGS quasi-Newton optimization */
int flat_bfgs_optimize(double(*f)(const double*,void*),void(*grad)(const double*,double*,void*),
void*ctx,int n,double*x,double tol,int max_iter){
if(!f||!grad||!x||n<=0)return -1;double*H=(double*)calloc(n*n,sizeof(double));double*g=(double*)calloc(n,sizeof(double));
double*g_new=(double*)calloc(n,sizeof(double));double*p=(double*)calloc(n,sizeof(double));
double*s=(double*)calloc(n,sizeof(double));double*yv=(double*)calloc(n,sizeof(double));
if(!H||!g||!g_new||!p||!s||!yv){free(H);free(g);free(g_new);free(p);free(s);free(yv);return -1;}
for(int i=0;i<n;i++)H[i*n+i]=1.0;grad(x,g,ctx);
for(int iter=0;iter<max_iter;iter++){double gnorm=0.0;for(int i=0;i<n;i++)gnorm+=g[i]*g[i];
if(sqrt(gnorm)<tol){free(H);free(g);free(g_new);free(p);free(s);free(yv);return iter+1;}
for(int i=0;i<n;i++){double val=0.0;for(int j=0;j<n;j++)val+=H[i*n+j]*g[j];p[i]=-val;}
double alpha=flat_golden_section_search(NULL,NULL,0,1,tol);
for(int i=0;i<n;i++){s[i]=alpha*p[i];x[i]+=s[i];}grad(x,g_new,ctx);
for(int i=0;i<n;i++)yv[i]=g_new[i]-g[i];double rho=0.0;for(int i=0;i<n;i++)rho+=yv[i]*s[i];
if(fabs(rho)<1e-15)continue;rho=1.0/rho;
double*Hs=(double*)calloc(n,sizeof(double));for(int i=0;i<n;i++){double val=0.0;for(int j=0;j<n;j++)val+=H[i*n+j]*s[j];Hs[i]=val;}
for(int i=0;i<n;i++)for(int j=0;j<n;j++)H[i*n+j]+=rho*((s[i]-Hs[i])*yv[j]+yv[i]*(s[j]-Hs[j]));
for(int i=0;i<n;i++)g[i]=g_new[i];free(Hs);}
free(H);free(g);free(g_new);free(p);free(s);free(yv);return -1;}
