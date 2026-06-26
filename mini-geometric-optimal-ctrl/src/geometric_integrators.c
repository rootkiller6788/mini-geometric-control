#include "geometric_integrators.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

void geo_gauss_legendre_tableau(GeoButcherTableau *t, int stages) {
    if (!t || stages < 1 || stages > 3) return;
    t->stages = stages;
    memset(t->a, 0, sizeof(t->a));
    memset(t->b, 0, sizeof(t->b));
    memset(t->c, 0, sizeof(t->c));
    if (stages == 1) {
        t->c[0] = 0.5; t->a[0][0] = 0.5; t->b[0] = 1.0;
    } else if (stages == 2) {
        double s3 = sqrt(3.0);
        t->c[0] = 0.5 - s3/6.0;
        t->c[1] = 0.5 + s3/6.0;
        t->a[0][0] = 0.25;
        t->a[0][1] = 0.25 - s3/6.0;
        t->a[1][0] = 0.25 + s3/6.0;
        t->a[1][1] = 0.25;
        t->b[0] = 0.5; t->b[1] = 0.5;
    } else {
        double s15 = sqrt(15.0);
        t->c[0] = 0.5 - s15/10.0;
        t->c[1] = 0.5;
        t->c[2] = 0.5 + s15/10.0;
        t->a[0][0]=5.0/36.0; t->a[0][1]=2.0/9.0-s15/15.0; t->a[0][2]=5.0/36.0-s15/30.0;
        t->a[1][0]=5.0/36.0+s15/24.0; t->a[1][1]=2.0/9.0; t->a[1][2]=5.0/36.0-s15/24.0;
        t->a[2][0]=5.0/36.0+s15/30.0; t->a[2][1]=2.0/9.0+s15/15.0; t->a[2][2]=5.0/36.0;
        t->b[0]=5.0/18.0; t->b[1]=4.0/9.0; t->b[2]=5.0/18.0;
    }
}double geo_discrete_lagrangian_midpoint(GeoScalarFunction L, void *ctx, const double *qk, const double *qk1, double h, int n) {
    if (!L||!qk||!qk1||n<=0) return 0.0;
    double *qm=(double*)malloc(n*sizeof(double));
    for (int i=0;i<n;i++) qm[i]=0.5*(qk[i]+qk1[i]);
    double r=h*L(qm,n,ctx); free(qm); return r;
}

int geo_del_step(GeoDiscreteLagrangian L_d, void *ctx, const double *qp, const double *qc, double *qn, double h, int n, double tol, int mi) {
    if (!L_d||!qp||!qc||!qn||n<=0) return -1;
    if (mi<=0) mi=20;
    for (int i=0;i<n;i++) qn[i]=2.0*qc[i]-qp[i];
    double *res=(double*)malloc(n*sizeof(double));
    for (int it=0;it<mi;it++) {
        double rn=0.0, hf=1e-6;
        for (int i=0;i<n;i++) {
            double *x=(double*)malloc(n*sizeof(double));
            memcpy(x,qc,n*sizeof(double));
            x[i]=qc[i]+hf; double Lp=L_d(qp,x,h,n,ctx);
            x[i]=qc[i]-hf; double Lm=L_d(qp,x,h,n,ctx);
            free(x); double d2=(Lp-Lm)/(2.0*hf);
            x=(double*)malloc(n*sizeof(double));
            memcpy(x,qc,n*sizeof(double));
            x[i]=qc[i]+hf; Lp=L_d(x,qn,h,n,ctx);
            x[i]=qc[i]-hf; Lm=L_d(x,qn,h,n,ctx);
            free(x); double d1=(Lp-Lm)/(2.0*hf);
            res[i]=d2+d1; rn+=res[i]*res[i];
        }
        rn=sqrt(rn); if (rn<tol) { free(res); return 0; }
        for (int i=0;i<n;i++) qn[i]-=0.5*res[i];
    }
    free(res); return -1;
}

void geo_discrete_legendre(const GeoDiscreteLagrangian L_d, void *ctx, const double *qk, const double *qk1, double *pk, double *pk1, double h, int n) {
    if (!L_d||!qk||!qk1||!pk||!pk1||n<=0) return;
    double hf=1e-6;
    for (int i=0;i<n;i++) {
        double *x=(double*)malloc(n*sizeof(double));
        memcpy(x,qk,n*sizeof(double));
        x[i]=qk[i]+hf; double Lp=L_d(x,qk1,h,n,ctx);
        x[i]=qk[i]-hf; double Lm=L_d(x,qk1,h,n,ctx);
        pk[i]=-(Lp-Lm)/(2.0*hf); free(x);
        x=(double*)malloc(n*sizeof(double));
        memcpy(x,qk1,n*sizeof(double));
        x[i]=qk1[i]+hf; Lp=L_d(qk,x,h,n,ctx);
        x[i]=qk1[i]-hf; Lm=L_d(qk,x,h,n,ctx);
        pk1[i]=(Lp-Lm)/(2.0*hf); free(x);
    }
}

int geo_avf_step(GeoHamiltonian H, void *ctx, double *x, double h, int n, int nq, double tol, int mi) {
    if (!H||!x||n<=0) return -1;
    if (nq<=0) nq=3; if (mi<=0) mi=20;
    int N=2*n;
    double *xn=(double*)malloc(N*sizeof(double));
    memcpy(xn,x,N*sizeof(double));
    for (int it=0;it<mi;it++) {
        double *avf=(double*)calloc(N,sizeof(double));
        for (int q=0;q<nq;q++) {
            double tau=(q+0.5)/(double)nq;
            double *xt=(double*)malloc(N*sizeof(double));
            for (int i=0;i<N;i++) xt[i]=(1.0-tau)*x[i]+tau*xn[i];
            double *Xv=(double*)malloc(N*sizeof(double));
            geo_hamiltonian_vector_field(H,xt,Xv,n,1e-6,ctx);
            for (int i=0;i<N;i++) avf[i]+=Xv[i]/nq;
            free(xt); free(Xv);
        }
        double *xc=(double*)malloc(N*sizeof(double));
        for (int i=0;i<N;i++) xc[i]=x[i]+h*avf[i];
        double diff=0.0;
        for (int i=0;i<N;i++) diff+=(xc[i]-xn[i])*(xc[i]-xn[i]);
        memcpy(xn,xc,N*sizeof(double));
        free(avf); free(xc);
        if (sqrt(diff)<tol) break;
    }
    memcpy(x,xn,N*sizeof(double)); free(xn);
    return 0;
}

int geo_discrete_gradient_step(GeoHamiltonian H, void *ctx, double *x, double h, int n, double tol, int mi) {
    if (!H||!x||n<=0) return -1;
    if (mi<=0) mi=20; int N=2*n;
    double *xn=(double*)malloc(N*sizeof(double));
    double *xb=(double*)malloc(N*sizeof(double));
    memcpy(xn,x,N*sizeof(double));
    for (int it=0;it<mi;it++) {
        memcpy(xb,x,N*sizeof(double));
        for (int i=0;i<N;i++) {
            xb[i]=xn[i]; double Hb=H(xb,n,ctx);
            xb[i]=x[i];  double Ho=H(xb,n,ctx);
            double den=xn[i]-x[i];
            double dg=(fabs(den)>1e-14)?(Hb-Ho)/den:0.0;
            if (i<n) xn[i]=x[i]+h*dg; else xn[i]=x[i]-h*dg;
        }
        double diff=0.0;
        for (int i=0;i<N;i++) diff+=(xn[i]-x[i])*(xn[i]-x[i]);
        if (sqrt(diff)<tol) break;
    }
    memcpy(x,xn,N*sizeof(double)); free(xn); free(xb);
    return 0;
}

void geo_yoshida_weights(int order, double *w, int *n_w) {
    if (!w||!n_w) return;
    if (order==4) {
        double w0=-cbrt(2.0)/(2.0-cbrt(2.0));
        double w1=1.0/(2.0-cbrt(2.0));
        w[0]=w1;w[1]=w0;w[2]=w1; *n_w=3;
    } else if (order==6) {
        double w0=-pow(2.0,1.0/7.0)/(2.0-pow(2.0,1.0/7.0));
        double w1=1.0/(2.0-pow(2.0,1.0/7.0));
        w[0]=w1;w[1]=w0;w[2]=w1;w[3]=w0;w[4]=w1;w[5]=w0;w[6]=w1; *n_w=7;
    } else { w[0]=1.0; *n_w=1; }
}

void geo_symplectic_euler(char v, GeoScalarFunction T, void *tc, GeoScalarFunction V, void *vc, double *q, double *p, double h, int n) {
    if (!T||!V||!q||!p||n<=0) return;
    double *g=(double*)malloc(n*sizeof(double));
    if (v==65||v==97) {
        double *x=(double*)malloc(n*sizeof(double));
        memcpy(x,q,n*sizeof(double));
        for (int i=0;i<n;i++) {
            x[i]=q[i]+1e-6; double Vp=V(x,n,vc);
            x[i]=q[i]-1e-6; double Vm=V(x,n,vc);
            g[i]=(Vp-Vm)/(2.0*1e-6); x[i]=q[i];
        }
        free(x); for (int i=0;i<n;i++) p[i]-=h*g[i];
        x=(double*)malloc(n*sizeof(double));
        memcpy(x,p,n*sizeof(double));
        for (int i=0;i<n;i++) {
            x[i]=p[i]+1e-6; double Tp=T(x,n,tc);
            x[i]=p[i]-1e-6; double Tm=T(x,n,tc);
            g[i]=(Tp-Tm)/(2.0*1e-6); x[i]=p[i];
        }
        free(x); for (int i=0;i<n;i++) q[i]+=h*g[i];
    } else {
        double *x=(double*)malloc(n*sizeof(double));
        memcpy(x,p,n*sizeof(double));
        for (int i=0;i<n;i++) {
            x[i]=p[i]+1e-6; double Tp=T(x,n,tc);
            x[i]=p[i]-1e-6; double Tm=T(x,n,tc);
            g[i]=(Tp-Tm)/(2.0*1e-6); x[i]=p[i];
        }
        free(x); for (int i=0;i<n;i++) q[i]+=h*g[i];
        x=(double*)malloc(n*sizeof(double));
        memcpy(x,q,n*sizeof(double));
        for (int i=0;i<n;i++) {
            x[i]=q[i]+1e-6; double Vp=V(x,n,vc);
            x[i]=q[i]-1e-6; double Vm=V(x,n,vc);
            g[i]=(Vp-Vm)/(2.0*1e-6); x[i]=q[i];
        }
        free(x); for (int i=0;i<n;i++) p[i]-=h*g[i];
    }
    free(g);
}

int geo_dmoc_solve(GeoDiscreteLagrangian Ld, void *Lc, GeoDiscreteDynamics fd, void *fc, const double *q0, const double *qT, int N, double *ug, double *qt, double *ut, int n, int m, double tol, int mi) {
    if (!Ld||!fd||!q0||!qT||!ug||!qt||!ut) return -1;
    if (N<=0||n<=0) return -1; (void)mi;
    memcpy(qt,q0,n*sizeof(double));
    memcpy(ut,ug,N*m*sizeof(double));
    for (int k=0;k<N;k++) fd(&qt[k*n],&ut[k*m],&qt[(k+1)*n],n,m,fc);
    double err=0.0;
    for (int i=0;i<n;i++) { double d=qt[N*n+i]-qT[i]; err+=d*d; }
    if (sqrt(err)>tol) return -1;
    return 0;
}
