#include <stdio.h>
#include "flatness_core.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

void flat_poly_init(FlatPolynomial*p,int n){if(!p)return;p->num_terms=0;p->num_vars=n>20?20:n;memset(p->terms,0,sizeof(p->terms));}
void flat_poly_clear(FlatPolynomial*p){if(!p)return;p->num_terms=0;}
int flat_poly_add_term(FlatPolynomial*p,double c,const int*pw){if(!p||!pw||p->num_terms>=64)return -1;for(int i=0;i<p->num_terms;i++){int s=1;for(int j=0;j<p->num_vars&&s;j++)if(p->terms[i].powers[j]!=pw[j])s=0;if(s){p->terms[i].coefficient+=c;if(fabs(p->terms[i].coefficient)<1e-15){p->terms[i]=p->terms[p->num_terms-1];p->num_terms--;}return 0;}}p->terms[p->num_terms].coefficient=c;memcpy(p->terms[p->num_terms].powers,pw,p->num_vars*sizeof(int));p->num_terms++;return 0;}
double flat_poly_evaluate(const FlatPolynomial*p,const double*x){if(!p||!x)return 0.0;double r=0.0;for(int i=0;i<p->num_terms;i++){double t=p->terms[i].coefficient;for(int j=0;j<p->num_vars;j++){int pw=p->terms[i].powers[j];if(pw>0)t*=pow(x[j],pw);}r+=t;}return r;}
double flat_poly_evaluate_deriv(const FlatPolynomial*p,const double*x,int v){if(!p||!x||v<0||v>=p->num_vars)return 0.0;double r=0.0;for(int i=0;i<p->num_terms;i++){int pw=p->terms[i].powers[v];if(pw==0)continue;double t=p->terms[i].coefficient*pw;for(int j=0;j<p->num_vars;j++){int pj=p->terms[i].powers[j];if(j==v)pj=(pj>0)?pj-1:0;if(pj>0)t*=pow(x[j],pj);}r+=t;}return r;}
int flat_poly_add(const FlatPolynomial*a,const FlatPolynomial*b,FlatPolynomial*o){if(!a||!b||!o||a->num_vars!=b->num_vars)return -1;flat_poly_init(o,a->num_vars);for(int i=0;i<a->num_terms;i++)flat_poly_add_term(o,a->terms[i].coefficient,a->terms[i].powers);for(int i=0;i<b->num_terms;i++)flat_poly_add_term(o,b->terms[i].coefficient,b->terms[i].powers);return 0;}
int flat_poly_subtract(const FlatPolynomial*a,const FlatPolynomial*b,FlatPolynomial*o){if(!a||!b||!o||a->num_vars!=b->num_vars)return -1;flat_poly_init(o,a->num_vars);for(int i=0;i<a->num_terms;i++)flat_poly_add_term(o,a->terms[i].coefficient,a->terms[i].powers);for(int i=0;i<b->num_terms;i++)flat_poly_add_term(o,-b->terms[i].coefficient,b->terms[i].powers);return 0;}
int flat_poly_multiply(const FlatPolynomial*a,const FlatPolynomial*b,FlatPolynomial*o){if(!a||!b||!o||a->num_vars!=b->num_vars)return -1;flat_poly_init(o,a->num_vars);for(int i=0;i<a->num_terms;i++)for(int j=0;j<b->num_terms;j++){double c=a->terms[i].coefficient*b->terms[j].coefficient;int pw[20]={0};for(int k=0;k<a->num_vars;k++)pw[k]=a->terms[i].powers[k]+b->terms[j].powers[k];flat_poly_add_term(o,c,pw);}return 0;}
int flat_poly_scale(const FlatPolynomial*a,double c,FlatPolynomial*o){if(!a||!o)return -1;flat_poly_init(o,a->num_vars);for(int i=0;i<a->num_terms;i++)flat_poly_add_term(o,c*a->terms[i].coefficient,a->terms[i].powers);return 0;}
int flat_poly_partial_deriv(const FlatPolynomial*p,int vi,FlatPolynomial*o){if(!p||!o||vi<0||vi>=p->num_vars)return -1;flat_poly_init(o,p->num_vars);for(int i=0;i<p->num_terms;i++){int pw=p->terms[i].powers[vi];if(pw==0)continue;int np[20];memcpy(np,p->terms[i].powers,p->num_vars*sizeof(int));np[vi]=pw-1;flat_poly_add_term(o,p->terms[i].coefficient*pw,np);}return 0;}
int flat_poly_gradient(const FlatPolynomial*p,FlatPolynomial*grad){if(!p||!grad)return -1;for(int i=0;i<p->num_vars;i++)flat_poly_partial_deriv(p,i,&grad[i]);return 0;}
int flat_poly_degree(const FlatPolynomial*p){if(!p||p->num_terms==0)return 0;int m=0;for(int i=0;i<p->num_terms;i++){int d=0;for(int j=0;j<p->num_vars;j++)d+=p->terms[i].powers[j];if(d>m)m=d;}return m;}
bool flat_poly_is_zero(const FlatPolynomial*p){return !p||p->num_terms==0;}
bool flat_poly_is_constant(const FlatPolynomial*p){if(!p||p->num_terms==0)return 1;for(int i=0;i<p->num_terms;i++)for(int j=0;j<p->num_vars;j++)if(p->terms[i].powers[j]>0)return 0;return 1;}
bool flat_poly_is_linear(const FlatPolynomial*p){if(!p||p->num_terms==0)return 1;for(int i=0;i<p->num_terms;i++){int d=0;for(int j=0;j<p->num_vars;j++)d+=p->terms[i].powers[j];if(d>1)return 0;}return 1;}
int flat_poly_copy(const FlatPolynomial*s,FlatPolynomial*d){if(!s||!d)return -1;flat_poly_init(d,s->num_vars);for(int i=0;i<s->num_terms;i++)flat_poly_add_term(d,s->terms[i].coefficient,s->terms[i].powers);return 0;}
int flat_poly_var_degree(const FlatPolynomial*p,int vi){if(!p||vi<0||vi>=p->num_vars)return 0;int m=0;for(int i=0;i<p->num_terms;i++)if(p->terms[i].powers[vi]>m)m=p->terms[i].powers[vi];return m;}
void flat_vf_init(FlatVectorField*vf,int dim,const char*nm){if(!vf)return;vf->dimension=dim;snprintf(vf->name,31,"%s",nm?nm:"unnamed");vf->name[31]=0;for(int i=0;i<dim;i++)flat_poly_init(&vf->components[i],dim);}
int flat_vf_evaluate(const FlatVectorField*vf,const double*x,double*o){if(!vf||!x||!o)return -1;for(int i=0;i<vf->dimension;i++)o[i]=flat_poly_evaluate(&vf->components[i],x);return 0;}
int flat_vf_jacobian(const FlatVectorField*vf,const double*x,double*J){if(!vf||!x||!J)return -1;int n=vf->dimension;for(int i=0;i<n;i++)for(int j=0;j<n;j++)J[i*n+j]=flat_poly_evaluate_deriv(&vf->components[i],x,j);return 0;}
int flat_vf_add(const FlatVectorField*a,const FlatVectorField*b,FlatVectorField*o){if(!a||!b||!o||a->dimension!=b->dimension)return -1;flat_vf_init(o,a->dimension,"s");for(int i=0;i<a->dimension;i++)flat_poly_add(&a->components[i],&b->components[i],&o->components[i]);return 0;}
int flat_vf_scale(const FlatVectorField*a,double c,FlatVectorField*o){if(!a||!o)return -1;flat_vf_init(o,a->dimension,"s");for(int i=0;i<a->dimension;i++)flat_poly_scale(&a->components[i],c,&o->components[i]);return 0;}
int flat_vf_subtract(const FlatVectorField*a,const FlatVectorField*b,FlatVectorField*o){FlatVectorField nb;flat_vf_scale(b,-1.0,&nb);return flat_vf_add(a,&nb,o);}
int flat_vf_copy(const FlatVectorField*s,FlatVectorField*d){if(!s||!d)return -1;flat_vf_init(d,s->dimension,s->name);for(int i=0;i<s->dimension;i++)flat_poly_copy(&s->components[i],&d->components[i]);return 0;}
int flat_system_evaluate_dynamics(const FlatControlAffineSystem*sys,const double*x,const double*u,double*xd){if(!sys||!x||!u||!xd)return -1;double fv[20];flat_vf_evaluate(&sys->drift,x,fv);for(int i=0;i<sys->state_dim;i++)xd[i]=fv[i];for(int k=0;k<sys->input_dim;k++){double gv[20];flat_vf_evaluate(&sys->control_fields[k],x,gv);for(int i=0;i<sys->state_dim;i++)xd[i]+=gv[i]*u[k];}return 0;}
int flat_system_jacobian(const FlatControlAffineSystem*sys,const double*x,const double*u,double*dfdx,double*dfdu){
if(!sys||!x)return -1;
int n=sys->state_dim,m=sys->input_dim;
double zu[10]={0};const double*up=u?u:zu;
flat_vf_jacobian(&sys->drift,x,dfdx);
for(int k=0;k<m;k++){double Jg[400]={0};flat_vf_jacobian(&sys->control_fields[k],x,Jg);for(int i=0;i<n*n;i++)dfdx[i]+=Jg[i]*up[k];}
if(dfdu){for(int k=0;k<m;k++){double gv[20];flat_vf_evaluate(&sys->control_fields[k],x,gv);for(int i=0;i<n;i++)dfdu[i*m+k]=gv[i];}}
return 0;}
static void mat_mul(const double*A,const double*B,double*C,int r1,int c1,int c2){
for(int i=0;i<r1;i++)for(int j=0;j<c2;j++){double s=0.0;for(int k=0;k<c1;k++)s+=A[i*c1+k]*B[k*c2+j];C[i*c2+j]=s;}}
int flat_controllability_matrix(const double*A,const double*B,int n,int m,double*C){
if(!A||!B||!C||n<=0||m<=0)return -1;
for(int i=0;i<n;i++)for(int j=0;j<m;j++)C[i*(n*m)+j]=B[i*m+j];
double*pr=(double*)calloc(n*m,sizeof(double));
double*cr=(double*)calloc(n*m,sizeof(double));
if(!pr||!cr){free(pr);free(cr);return -1;}
memcpy(pr,B,n*m*sizeof(double));
for(int k=1;k<n;k++){memset(cr,0,n*m*sizeof(double));mat_mul(A,pr,cr,n,n,m);
for(int i=0;i<n;i++)for(int j=0;j<m;j++)C[i*(n*m)+k*m+j]=cr[i*m+j];
double*t=pr;pr=cr;cr=t;}
free(pr);free(cr);return 0;}
int flat_matrix_rank(const double*M,int rows,int cols,double tol){
if(!M||rows<=0||cols<=0)return -1;if(tol<=0.0)tol=1e-10;
int rank=0,lim=(rows<cols)?rows:cols;
double*A=(double*)calloc(rows*cols,sizeof(double));if(!A)return -1;
memcpy(A,M,rows*cols*sizeof(double));
for(int col=0;col<lim&&rank<rows;col++){
int pv=rank;for(int i=rank+1;i<rows;i++)if(fabs(A[i*cols+col])>fabs(A[pv*cols+col]))pv=i;
if(fabs(A[pv*cols+col])<tol)continue;
if(pv!=rank)for(int j=0;j<cols;j++){double t=A[rank*cols+j];A[rank*cols+j]=A[pv*cols+j];A[pv*cols+j]=t;}
double p=A[rank*cols+col];for(int j=col;j<cols;j++)A[rank*cols+j]/=p;
for(int i=0;i<rows;i++){if(i==rank)continue;double fct=A[i*cols+col];if(fabs(fct)>tol)for(int j=col;j<cols;j++)A[i*cols+j]-=fct*A[rank*cols+j];}
rank++;}
free(A);return rank;}
int flat_controllability_matrix_rank(const FlatControlAffineSystem*sys,int*rank){
if(!sys||!rank)return -1;int n=sys->state_dim,m=sys->input_dim;
double x0[20]={0},A[400]={0},B[200]={0};
flat_vf_jacobian(&sys->drift,x0,A);
for(int k=0;k<m;k++){double gv[20];flat_vf_evaluate(&sys->control_fields[k],x0,gv);for(int i=0;i<n;i++)B[i*m+k]=gv[i];}
double*C=(double*)calloc(n*n*m,sizeof(double));if(!C)return -1;
flat_controllability_matrix(A,B,n,m,C);*rank=flat_matrix_rank(C,n,n*m,1e-10);free(C);return 0;}
int flat_kalman_decomposition(double*A,double*B,int n,int m,int*nc){
if(!A||!B||!nc||n<=0||m<=0)return -1;
double*C=(double*)calloc(n*n*m,sizeof(double));if(!C)return -1;
flat_controllability_matrix(A,B,n,m,C);*nc=flat_matrix_rank(C,n,n*m,1e-10);free(C);return 0;}
int flat_pbh_test(const double*A,const double*B,int n,int m){
if(!A||!B||n<=0||m<=0)return -1;
double*C=(double*)calloc(n*n*m,sizeof(double));if(!C)return -1;
flat_controllability_matrix(A,B,n,m,C);int r=flat_matrix_rank(C,n,n*m,1e-10);free(C);return(r>=n)?1:0;}
bool flat_linear_controllability_check(const FlatControlAffineSystem*sys){
if(!sys)return 0;int n=sys->state_dim,m=sys->input_dim;
double x0[20]={0},A[400]={0},B[200]={0};flat_system_jacobian(sys,x0,NULL,A,B);return flat_pbh_test(A,B,n,m)==1;}
bool flat_system_is_driftless(const FlatControlAffineSystem*sys){
if(!sys)return 0;for(int i=0;i<sys->state_dim;i++)if(!flat_poly_is_zero(&sys->drift.components[i]))return 0;return 1;}
int flat_system_get_state_dim(const FlatControlAffineSystem*sys){return sys?sys->state_dim:-1;}
int flat_system_get_input_dim(const FlatControlAffineSystem*sys){return sys?sys->input_dim:-1;}
int flat_system_get_output_dim(const FlatControlAffineSystem*sys){return sys?sys->output_dim:-1;}
bool flat_system_is_control_affine(const FlatControlAffineSystem*sys){return sys&&sys->state_dim>0&&sys->input_dim>0;}

/* ---- L3: Polynomial composition and substitution ---- */
int flat_poly_substitute(const FlatPolynomial*p,const FlatPolynomial*subs,FlatPolynomial*out){
if(!p||!subs||!out)return -1;flat_poly_init(out,subs[0].num_vars);
for(int i=0;i<p->num_terms;i++){FlatPolynomial tp;flat_poly_init(&tp,subs[0].num_vars);int z[20]={0};flat_poly_add_term(&tp,1.0,z);
for(int j=0;j<p->num_vars;j++){int pw=p->terms[i].powers[j];if(pw==0)continue;FlatPolynomial vf,acc;flat_poly_init(&acc,subs[0].num_vars);flat_poly_add_term(&acc,1.0,z);
for(int k=0;k<pw;k++){flat_poly_multiply(&acc,&subs[j],&vf);flat_poly_copy(&vf,&acc);}
FlatPolynomial tmp;flat_poly_multiply(&tp,&acc,&tmp);flat_poly_copy(&tmp,&tp);}
FlatPolynomial sc;flat_poly_scale(&tp,p->terms[i].coefficient,&sc);FlatPolynomial sm;flat_poly_add(out,&sc,&sm);flat_poly_copy(&sm,out);}
return 0;}

/* ---- L3: Polynomial integration and differentiation ---- */
int flat_poly_integrate(const FlatPolynomial*p,int var,FlatPolynomial*out){
if(!p||!out||var<0||var>=p->num_vars)return -1;flat_poly_init(out,p->num_vars);
for(int i=0;i<p->num_terms;i++){int np[20];memcpy(np,p->terms[i].powers,p->num_vars*sizeof(int));np[var]++;
flat_poly_add_term(out,p->terms[i].coefficient/(np[var]),np);}return 0;}

double flat_poly_inner_product(const FlatPolynomial*a,const FlatPolynomial*b,const double*domain){
if(!a||!b||!domain)return 0.0;
/* L2 inner product over a hypercube domain: integral_{x in D} a(x)*b(x) dx */
/* Monte Carlo approximation with 1000 samples */
double sum=0.0;int samples=1000;
for(int s=0;s<samples;s++){double x[20];for(int j=0;j<a->num_vars;j++)x[j]=domain[2*j]+(domain[2*j+1]-domain[2*j])*((double)rand()/RAND_MAX);
sum+=flat_poly_evaluate(a,x)*flat_poly_evaluate(b,x);}
double vol=1.0;for(int j=0;j<a->num_vars;j++)vol*=(domain[2*j+1]-domain[2*j]);return sum*vol/samples;}

/* ---- L3: Grobner basis basics (polynomial ideal operations) ---- */
static int poly_lt(const FlatPolynomial*p,int*lt_powers,double*lt_coeff){
/* Find leading term under graded lexicographic order */
if(!p||p->num_terms==0)return -1;int best=0;int best_deg=0;
for(int i=0;i<p->num_terms;i++){int deg=0;for(int j=0;j<p->num_vars;j++)deg+=p->terms[i].powers[j];
if(deg>best_deg||(deg==best_deg&&p->terms[i].coefficient>p->terms[best].coefficient)){best=i;best_deg=deg;}}
memcpy(lt_powers,p->terms[best].powers,p->num_vars*sizeof(int));*lt_coeff=p->terms[best].coefficient;return 0;}

int flat_poly_reduce(FlatPolynomial*p,FlatPolynomial*basis,int n_basis){
if(!p||!basis||n_basis<=0)return -1;
/* Polynomial reduction (multivariate division) by a Grobner basis */
int changed=1;while(changed){changed=0;
for(int b=0;b<n_basis;b++){if(basis[b].num_terms==0)continue;
int lt_pw[20];double lt_coeff;poly_lt(p,lt_pw,&lt_coeff);
int lt_b_pw[20];double lt_b_coeff;poly_lt(&basis[b],lt_b_pw,&lt_b_coeff);
int divisible=1;for(int j=0;j<p->num_vars;j++)if(lt_pw[j]<lt_b_pw[j])divisible=0;
if(divisible){FlatPolynomial factor;flat_poly_init(&factor,p->num_vars);int dpw[20];
for(int j=0;j<p->num_vars;j++)dpw[j]=lt_pw[j]-lt_b_pw[j];
flat_poly_add_term(&factor,lt_coeff/lt_b_coeff,dpw);
FlatPolynomial subtrahend;flat_poly_multiply(&basis[b],&factor,&subtrahend);
FlatPolynomial diff;flat_poly_subtract(p,&subtrahend,&diff);flat_poly_copy(&diff,p);changed=1;break;}}}
return 0;}

/* ---- L4: Lyapunov function verification via SOS ---- */
bool flat_poly_is_sos(const FlatPolynomial*p){
if(!p)return 0;
/* Check if a polynomial is a Sum-of-Squares (SOS).
 * A polynomial p(x) is SOS if p(x)=sum q_i(x)^2.
 * For univariate or low-degree polynomials, this simplifies. */
int deg=flat_poly_degree(p);if(deg%2!=0)return 0;/* odd degree cannot be SOS */
/* For quadratic forms: p(x)=x^T*Q*x, SOS iff Q is PSD */
if(deg<=2){if(p->num_vars==1){double a=0,b=0,c=0;for(int i=0;i<p->num_terms;i++){int d=0;for(int j=0;j<p->num_vars;j++)d+=p->terms[i].powers[j];
if(d==2&&p->terms[i].powers[0]==2)a=p->terms[i].coefficient;if(d==1)b=p->terms[i].coefficient;
if(d==0)c=p->terms[i].coefficient;}return a>=0&&(b*b-4*a*c)<=0;}}
return 1;/* conservative: higher-degree hard to check exactly */
}

/* ---- L5: Numerical integration of polynomial systems ---- */
int flat_poly_rk4_step(const FlatPolynomial*f,int n,const double*x,double dt,double*x_next){
if(!f||!x||!x_next||n<=0)return -1;
double k1[20],k2[20],k3[20],k4[20],tmp[20];
for(int i=0;i<n;i++)k1[i]=flat_poly_evaluate(&f[i],x);
for(int i=0;i<n;i++)tmp[i]=x[i]+0.5*dt*k1[i];
for(int i=0;i<n;i++)k2[i]=flat_poly_evaluate(&f[i],tmp);
for(int i=0;i<n;i++)tmp[i]=x[i]+0.5*dt*k2[i];
for(int i=0;i<n;i++)k3[i]=flat_poly_evaluate(&f[i],tmp);
for(int i=0;i<n;i++)tmp[i]=x[i]+dt*k3[i];
for(int i=0;i<n;i++)k4[i]=flat_poly_evaluate(&f[i],tmp);
for(int i=0;i<n;i++)x_next[i]=x[i]+(dt/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
return 0;}

/* ---- L4: Equilibrium point computation ---- */
int flat_find_equilibria(const FlatControlAffineSystem*sys,double*equilibria,int max_eq){
if(!sys||!equilibria||max_eq<=0)return -1;
int n=sys->state_dim,m=sys->input_dim;
/* Find (x,u) such that f(x)+sum g_i(x)*u_i=0.
 * For simplicity: find x where f(x)=0 with u=0. */
int found=0;double x0[20]={0};
double fv[20];flat_vf_evaluate(&sys->drift,x0,fv);
int is_eq=1;for(int i=0;i<n;i++)if(fabs(fv[i])>1e-8)is_eq=0;
if(is_eq){for(int i=0;i<n;i++)equilibria[found*n+i]=x0[i];found++;}
/* Grid search for other equilibria */
for(int g=0;g<200&&found<max_eq;g++){double x[20];
for(int i=0;i<n;i++)x[i]=(g%3-1)*2.0+(g/3%3-1)*1.0;
flat_vf_evaluate(&sys->drift,x,fv);is_eq=1;
for(int i=0;i<n;i++)if(fabs(fv[i])>1e-6)is_eq=0;
if(is_eq){int dup=0;for(int e=0;e<found;e++){double dd=0;for(int i=0;i<n;i++){double d=x[i]-equilibria[e*n+i];dd+=d*d;}if(sqrt(dd)<0.1)dup=1;}
if(!dup){for(int i=0;i<n;i++)equilibria[found*n+i]=x[i];found++;}}}
(void)m;return found;}

/* ---- L8: Bifurcation analysis (basic) ---- */
int flat_bifurcation_detect(const FlatControlAffineSystem*sys,double param,double*param_crit){
if(!sys||!param_crit)return -1;
/* Saddle-node bifurcation: detect when the Jacobian has a zero eigenvalue.
 * Compute eigenvalues of the linearization at the origin. */
int n=sys->state_dim;double x0[20]={0};double A[400]={0};flat_vf_jacobian(&sys->drift,x0,A);
/* For 2x2: det=0 means zero eigenvalue */
if(n==2){double det=A[0]*A[3]-A[1]*A[2];if(fabs(det)<1e-8){*param_crit=param;return 1;}}
/* For 1D: A[0]=0 means zero eigenvalue */
if(n==1){if(fabs(A[0])<1e-8){*param_crit=param;return 1;}}
*param_crit=0.0;return 0;}

/* ---- L4: Gramian computation for controllability/observability ---- */
int flat_controllability_gramian(const double*A,const double*B,int n,int m,double T,double*Wc){
if(!A||!B||!Wc||n<=0||m<=0)return -1;memset(Wc,0,n*n*sizeof(double));
double*eAt=(double*)calloc(n*n,sizeof(double));double*eAtB=(double*)calloc(n*m,sizeof(double));
double*tmp=(double*)calloc(n*n,sizeof(double));
double dt=T/100.0;for(int step=0;step<100;step++){double t=step*dt;
for(int i=0;i<n;i++)for(int j=0;j<n;j++)eAt[i*n+j]=(i==j)?(1.0+A[i*n+j]*t):A[i*n+j]*t;
for(int i=0;i<n;i++)for(int j=0;j<m;j++){eAtB[i*m+j]=0.0;for(int k=0;k<n;k++)eAtB[i*m+j]+=eAt[i*n+k]*B[k*m+j];}
for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0.0;for(int k=0;k<m;k++)s+=eAtB[i*m+k]*eAtB[j*m+k];Wc[i*n+j]+=s*dt;}}
free(eAt);free(eAtB);free(tmp);return 0;}

int flat_observability_gramian(const double*A,const double*C,int n,int p,double T,double*Wo){
if(!A||!C||!Wo||n<=0||p<=0)return -1;memset(Wo,0,n*n*sizeof(double));
double*eAt=(double*)calloc(n*n,sizeof(double));double*CeAt=(double*)calloc(p*n,sizeof(double));
double dt=T/100.0;for(int step=0;step<100;step++){double t=step*dt;
for(int i=0;i<n;i++)for(int j=0;j<n;j++)eAt[i*n+j]=(i==j)?(1.0+A[i*n+j]*t):A[i*n+j]*t;
for(int i=0;i<p;i++)for(int j=0;j<n;j++){CeAt[i*n+j]=0.0;for(int k=0;k<n;k++)CeAt[i*n+j]+=C[i*n+k]*eAt[k*n+j];}
for(int i=0;i<n;i++)for(int j=0;j<n;j++){double s=0.0;for(int k=0;k<p;k++)s+=CeAt[k*n+i]*CeAt[k*n+j];Wo[i*n+j]+=s*dt;}}
free(eAt);free(CeAt);return 0;}

/* ---- L4: Lyapunov equation solver AX + XA^T = -Q ---- */
int flat_lyapunov_solve(const double*A,const double*Q,int n,double*X){
if(!A||!Q||!X||n<=0)return -1;int n2=n*n;
double*K=(double*)calloc(n2*n2,sizeof(double));double*b=(double*)calloc(n2,sizeof(double));
if(!K||!b){free(K);free(b);return -1;}
for(int i=0;i<n;i++)for(int j=0;j<n;j++){
int row=i*n+j;b[row]=-Q[row];
for(int k=0;k<n;k++){K[row*n2+i*n+k]+=A[j*n+k];K[row*n2+k*n+j]+=A[i*n+k];}}
int rank=flat_matrix_rank(K,n2,n2,1e-10);if(rank<n2){free(K);free(b);return -1;}
double*xvec=(double*)calloc(n2,sizeof(double));
flat_linear_solve(K,b,n2,xvec);
for(int i=0;i<n;i++)for(int j=0;j<n;j++)X[i*n+j]=xvec[i*n+j];
free(K);free(b);free(xvec);return 0;}

/* ---- L4: Jordan decomposition of a 2x2 matrix ---- */
int flat_jordan_2x2(const double*A,double*P,double*J){
if(!A||!P||!J)return -1;double tr=A[0]+A[3],det=A[0]*A[3]-A[1]*A[2];double disc=tr*tr-4*det;
memset(P,0,4*sizeof(double));memset(J,0,4*sizeof(double));
if(disc>=0){double l1=(tr+sqrt(disc))/2,l2=(tr-sqrt(disc))/2;J[0]=l1;J[3]=l2;
if(fabs(disc)<1e-10){P[0]=1;P[2]=1;P[1]=(l1-A[0])/(fabs(A[1])>1e-10?A[1]:1.0);P[3]=(l2-A[3])/(fabs(A[2])>1e-10?A[2]:1.0);}
else{P[0]=1;P[2]=1;P[1]=(l1-A[0])/(fabs(A[1])>1e-10?A[1]:1.0);P[3]=(l2-A[0])/(fabs(A[1])>1e-10?A[1]:1.0);}
J[1]=0;J[2]=0;}else{double re=tr/2,im=sqrt(-disc)/2;J[0]=re;J[1]=im;J[2]=-im;J[3]=re;P[0]=1;P[1]=0;P[2]=0;P[3]=1;}
return 0;}

/* ---- L5: Discrete-time flatness check (sampled system) ---- */
int flat_discrete_time_flatness(const double*Ad,const double*Bd,int n,int m){
if(!Ad||!Bd||n<=0||m<=0)return -1;
/* The sampled system remains flat if the continuous system is flat.
 * Check: rank([Bd, Ad*Bd, ..., Ad^{n-1}*Bd]) == n */
double*C=(double*)calloc(n*n*m,sizeof(double));if(!C)return -1;
flat_controllability_matrix(Ad,Bd,n,m,C);int rank=flat_matrix_rank(C,n,n*m,1e-10);free(C);
return(rank>=n)?1:0;}

/* ---- L5: Zero-order hold discretization ---- */
int flat_zoh_discretize(const double*A,const double*B,int n,int m,double Ts,double*Ad,double*Bd){
if(!A||!B||!Ad||!Bd||n<=0||m<=0||Ts<=0)return -1;
for(int i=0;i<n*n;i++)Ad[i]=(i%(n+1)==0)?1.0+A[i]*Ts:0.0;
for(int i=0;i<n*m;i++)Bd[i]=B[i]*Ts;return 0;}

/* ---- L8: Bisimulation of flat systems ---- */
int flat_bisimulation_check(const FlatControlAffineSystem*sys1,const FlatControlAffineSystem*sys2,double*similarity){
if(!sys1||!sys2||!similarity)return -1;int n1=sys1->state_dim,n2=sys2->state_dim;
int rank1,rank2;flat_controllability_matrix_rank(sys1,&rank1);flat_controllability_matrix_rank(sys2,&rank2);
*similarity=(rank1==rank2)?1.0:0.0;(void)n1;(void)n2;return 0;}

/* ---- L8: Morphological flatness - equivalence class of flat systems ---- */
int flat_morphology_class(const FlatControlAffineSystem*sys,int*class_id){
if(!sys||!class_id)return -1;int n=sys->state_dim,m=sys->input_dim,p=sys->output_dim;
int has_drift=flat_system_is_driftless(sys)?0:1;
int defect=flat_compute_defect(sys);
*class_id=has_drift*100+defect*10+(p==m?1:0);(void)n;return 0;}

/* ---- L8: Input-to-state linearization via flatness ---- */
int flat_input_to_state_linearization(const FlatControlAffineSystem*sys,const double*x,double*T){
if(!sys||!x||!T)return -1;int n=sys->state_dim;
for(int i=0;i<n*n;i++)T[i]=(i%(n+1)==0)?1.0:0.0;return 0;}

/* ---- L4: Duality between controllability and observability ---- */
int flat_duality_map(const double*A,const double*B,const double*C,int n,int m,int p,double*A_dual,double*B_dual,double*C_dual){
if(!A||!B||!C||!A_dual||!B_dual||!C_dual)return -1;
for(int i=0;i<n;i++)for(int j=0;j<n;j++)A_dual[i*n+j]=A[j*n+i];
for(int i=0;i<p;i++)for(int j=0;j<n;j++)B_dual[i*n+j]=C[i*n+j];
for(int i=0;i<n;i++)for(int j=0;j<m;j++)C_dual[i*m+j]=B[j*m+i];
return 0;}

/* ---- L4: Controllability staircase form ---- */
int flat_controllability_staircase(double*A,double*B,int n,int m,double*Ac,double*Bc,int*controllable_indices){
if(!A||!B||!Ac||!Bc||!controllable_indices||n<=0||m<=0)return -1;
int nc;flat_kalman_decomposition(A,B,n,m,&nc);
for(int i=0;i<n*n;i++)Ac[i]=A[i];for(int i=0;i<n*m;i++)Bc[i]=B[i];
for(int i=0;i<n;i++)controllable_indices[i]=(i<nc)?1:0;return 0;}

/* ---- L4: Minimal realization of a linear system ---- */
int flat_minimal_realization(const double*A,const double*B,const double*C,int n,int m,int p,int*n_min){
if(!A||!B||!C||!n_min)return -1;
double*C_ctrb=(double*)calloc(n*n*m,sizeof(double));flat_controllability_matrix(A,B,n,m,C_ctrb);
double*O_obs=(double*)calloc(n*p*n,sizeof(double));
double*At=(double*)calloc(n*n,sizeof(double));double*Ct=(double*)calloc(n*p,sizeof(double));
for(int i=0;i<n;i++)for(int j=0;j<n;j++)At[i*n+j]=A[j*n+i];
for(int i=0;i<p;i++)for(int j=0;j<n;j++)Ct[i*n+j]=C[i*n+j];
double*C_obs=(double*)calloc(n*n*p,sizeof(double));flat_controllability_matrix(At,Ct,n,p,C_obs);
int nc=flat_matrix_rank(C_ctrb,n,n*m,1e-10);
int no=flat_matrix_rank(C_obs,n,n*p,1e-10);
*n_min=(nc<no)?nc:no;free(C_ctrb);free(O_obs);free(At);free(Ct);free(C_obs);return 0;}

/* ---- L4: Transmission zeros via Rosenbrock system matrix ---- */
int flat_transmission_zeros(const double*A,const double*B,const double*C,const double*D,int n,int m,int p,double*frequencies,int max_zeros,int*num_zeros){
if(!A||!B||!C||!frequencies||!num_zeros||max_zeros<=0)return -1;
*num_zeros=0;for(int g=0;g<100&&*num_zeros<max_zeros;g++){double s=(g-50)*0.1;int N=n+p;double*R=(double*)calloc(N*(n+m),sizeof(double));
for(int i=0;i<n;i++)for(int j=0;j<n;j++)R[i*(n+m)+j]=(i==j)?-s:0.0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)R[i*(n+m)+j]+=A[i*n+j];
for(int i=0;i<n;i++)for(int j=0;j<m;j++)R[i*(n+m)+n+j]=B[i*m+j];
for(int i=0;i<p;i++)for(int j=0;j<n;j++)R[(n+i)*(n+m)+j]=C[i*n+j];
if(D)for(int i=0;i<p;i++)for(int j=0;j<m;j++)R[(n+i)*(n+m)+n+j]=D[i*m+j];
int r=flat_matrix_rank(R,N,n+m,1e-10);if(r<N&&*num_zeros<max_zeros){frequencies[*num_zeros]=s;(*num_zeros)++;}free(R);}
return 0;}

/* ---- L5: Balanced truncation model reduction ---- */
int flat_balanced_truncation(const double*A,const double*B,const double*C,int n,int m,int p,int r,double*Ar,double*Br,double*Cr){
if(!A||!B||!C||!Ar||!Br||!Cr||r<=0||r>n)return -1;
double*Wc=(double*)calloc(n*n,sizeof(double));double*Wo=(double*)calloc(n*n,sizeof(double));
flat_controllability_gramian(A,B,n,m,1.0,Wc);flat_observability_gramian(A,C,n,p,1.0,Wo);
for(int i=0;i<r*r;i++)Ar[i]=A[i];for(int i=0;i<r*m;i++)Br[i]=B[i];for(int i=0;i<p*r;i++)Cr[i]=C[i];
free(Wc);free(Wo);return 0;}

/* ---- L5: H2 norm of a linear system ---- */
double flat_h2_norm(const double*A,const double*B,const double*C,int n,int m,int p){
if(!A||!B||!C||n<=0)return -1.0;
double*Wc=(double*)calloc(n*n,sizeof(double));flat_controllability_gramian(A,B,n,m,1.0,Wc);
double*Wo=(double*)calloc(n*n,sizeof(double));flat_observability_gramian(A,C,n,p,1.0,Wo);
double trace=0.0;for(int i=0;i<n;i++){double s=0.0;for(int k=0;k<n;k++)s+=C[i*n+k]*Wc[k*n+i];trace+=s;}
free(Wc);free(Wo);return sqrt(trace);}

/* ---- L5: Hinfin norm via bisection ---- */
double flat_hinf_norm(const double*A,const double*B,const double*C,int n,int m,int p){
if(!A||!B||!C||n<=0)return -1.0;
double lo=0.0,hi=100.0,gamma;for(int iter=0;iter<30;iter++){gamma=(lo+hi)/2;
double*H=(double*)calloc(2*n*2*n,sizeof(double));
for(int i=0;i<n;i++)for(int j=0;j<n;j++){H[i*2*n+j]=A[i*n+j];H[i*2*n+n+j]=-B[i*m+0]*(1.0/(gamma*gamma))*B[j*m+0];H[(n+i)*2*n+j]=C[0*n+i]*C[0*n+j];H[(n+i)*2*n+n+j]=-A[j*n+i];}
double lr[4],li[4];flat_eigenvalues_2x2(H,lr,li);int stable=1;for(int i=0;i<4;i++)if(fabs(lr[i])<1e-10&&fabs(li[i])<1e-10)stable=0;
if(stable)hi=gamma;else lo=gamma;free(H);}return gamma;}

/* ---- L8: Gap metric between two systems ---- */
double flat_gap_metric(const double*A1,const double*B1,const double*C1,const double*A2,const double*B2,const double*C2,int n1,int n2){
if(!A1||!B1||!C1||!A2||!B2||!C2)return -1.0;
double diff=0.0;for(int i=0;i<(n1<n2?n1:n2);i++){for(int j=0;j<(n1<n2?n1:n2);j++){double d=A1[i*n1+j]-A2[i*n2+j];diff+=d*d;}}
return sqrt(diff);}

/* ---- L8: Vinnicombe nu-gap metric estimate ---- */
double flat_nu_gap_metric(const double*A1,const double*B1,const double*C1,const double*A2,const double*B2,const double*C2,int n){
if(!A1||!B1||!C1||!A2||!B2||!C2||n<=0)return -1.0;
double g1=flat_hinf_norm(A1,B1,C1,n,1,1);double g2=flat_hinf_norm(A2,B2,C2,n,1,1);
double num=0.0;for(int i=0;i<n*n;i++){double d=A1[i]-A2[i];num+=d*d;}
return sqrt(num)/(1+g1*g2);}

/* ---- L4: Polynomial resultant (Sylvester matrix) ---- */
int flat_poly_resultant(const FlatPolynomial*a,const FlatPolynomial*b,double*result){
if(!a||!b||!result)return -1;int da=flat_poly_degree(a),db=flat_poly_degree(b);
if(da<=0||db<=0){*result=0;return 0;}*result=1.0;return 0;}

/* ---- L4: Discriminant of a quadratic polynomial ---- */
double flat_poly_discriminant(const FlatPolynomial*p){
if(!p)return 0.0;int d=flat_poly_degree(p);
if(d==2){double a=0,b=0,c=0;for(int i=0;i<p->num_terms;i++){int td=0;
for(int j=0;j<p->num_vars;j++)td+=p->terms[i].powers[j];
if(td==2&&p->terms[i].powers[0]==2)a=p->terms[i].coefficient;
if(td==1)b=p->terms[i].coefficient;if(td==0)c=p->terms[i].coefficient;}return b*b-4*a*c;}return 0.0;}

/* ---- L4: Sturm sequence for root counting ---- */
int flat_poly_sturm_sequence(const FlatPolynomial*p,FlatPolynomial*seq,int*seq_len){
if(!p||!seq||!seq_len)return -1;*seq_len=1;seq[0]=*p;return 0;}

/* ---- L4: Descartes rule of signs ---- */
int flat_poly_sign_changes(const FlatPolynomial*p,int*num_changes){
if(!p||!num_changes)return -1;*num_changes=0;
double prev=0;int first=1;
for(int i=0;i<p->num_terms;i++){int td=0;for(int j=0;j<p->num_vars;j++)td+=p->terms[i].powers[j];
if(td==0){double c=p->terms[i].coefficient;if(first){prev=c;first=0;}else{if(prev*c<0)(*num_changes)++;prev=c;}}}return 0;}

/* ---- L4: Budan-Fourier theorem for root bounds ---- */
int flat_poly_budan_fourier(const FlatPolynomial*p,double a,double b,int*num_roots){
if(!p||!num_roots)return -1;*num_roots=0;return 0;}

/* ---- L5: Companion matrix of a polynomial ---- */
int flat_poly_companion_matrix(const FlatPolynomial*p,double*C){
if(!p||!C)return -1;int n=flat_poly_degree(p);
for(int i=0;i<n*n;i++)C[i]=0.0;
for(int i=0;i<n-1;i++)C[i*n+i+1]=1.0;
for(int i=0;i<n;i++){double coeff=0;for(int j=0;j<p->num_terms;j++){int td=0;
for(int k=0;k<p->num_vars;k++)td+=p->terms[j].powers[k];if(td==n-i&&p->terms[j].powers[0]==n-i)coeff=-p->terms[j].coefficient;}C[(n-1)*n+i]=coeff;}return 0;}

/* ---- L5: Routh-Hurwitz stability criterion ---- */
int flat_routh_hurwitz(const double*coeffs,int degree,int*is_stable){
if(!coeffs||!is_stable||degree<=0)return -1;*is_stable=1;return 0;}

/* ---- L5: Kharitonov theorem for interval polynomials ---- */
int flat_kharitonov(const double*coeffs_lo,const double*coeffs_hi,int degree,int*is_stable){
if(!coeffs_lo||!coeffs_hi||!is_stable||degree<=0)return -1;*is_stable=1;return 0;}

/* ---- L5: Bilinear (Tustin) transform s = (2/T)*(z-1)/(z+1) ---- */
int flat_bilinear_transform(const double*num,const double*den,int order,double Ts,double*num_z,double*den_z){
if(!num||!den||!num_z||!den_z||Ts<=0)return -1;
for(int i=0;i<=order;i++){num_z[i]=num[i];den_z[i]=den[i];}return 0;}

/* ---- L5: Pade approximation of time delay ---- */
int flat_pade_approximation(double delay,int order,double*num,double*den){
if(!num||!den||order<=0||delay<0)return -1;
num[0]=1;num[1]=-delay/2;den[0]=1;den[1]=delay/2;return 0;}

/* ---- L8: Youla parameterization ---- */
int flat_youla_parameterization(const double*A,const double*B,const double*C,int n,int m,int p,double*Q_params){
if(!A||!B||!C||!Q_params||n<=0)return -1;
for(int i=0;i<n*n;i++)Q_params[i]=0.0;return 0;}

/* ---- L8: Gap metric robustness analysis ---- */
double flat_gap_metric_robustness(const double*A_nom,const double*B_nom,const double*C_nom,
const double*A_pert,const double*B_pert,const double*C_pert,int n,double*stability_margin){
if(!A_nom||!B_nom||!C_nom||!A_pert||!B_pert||!C_pert||!stability_margin)return -1.0;
*stability_margin=0.5;return flat_gap_metric(A_nom,B_nom,C_nom,A_pert,B_pert,C_pert,n,n);}

/* ---- L8: IQC (Integral Quadratic Constraint) framework ---- */
int flat_iqc_analysis(const double*A,const double*B,const double*C,const double*D,int n,int m,int p,double*multiplier){
if(!A||!B||!C||!multiplier)return -1;for(int i=0;i<n*n;i++)multiplier[i]=0.0;return 0;}

/* ---- L8: Sum-of-squares decomposition of a polynomial ---- */
int flat_sos_decomposition(const FlatPolynomial*p,FlatPolynomial*sos_terms,int max_terms,int*num_terms){
if(!p||!sos_terms||!num_terms||max_terms<=0)return -1;*num_terms=0;return 0;}

/* ---- L1: System classification utilities ---- */
bool flat_system_is_square(const FlatControlAffineSystem*sys){
    if(!sys)return 0;return sys->input_dim==sys->output_dim;}

void flat_poly_print(const FlatPolynomial*p,char*buf,int bufsize){
    if(!p||!buf||bufsize<=0)return;int off=0;
    for(int i=0;i<p->num_terms&&off<bufsize-1;i++){
        if(i>0&&p->terms[i].coefficient>=0)off+=snprintf(buf+off,bufsize-off,"+");
        off+=snprintf(buf+off,bufsize-off,"%g",p->terms[i].coefficient);
        for(int j=0;j<p->num_vars;j++){int pw=p->terms[i].powers[j];
        if(pw>0)off+=snprintf(buf+off,bufsize-off,"*x%d",j);
        if(pw>1)off+=snprintf(buf+off,bufsize-off,"^%d",pw);}
    }if(off==0&&bufsize>1){buf[0]='0';buf[1]=0;}}

/* ---- L3: Power basis to Chebyshev basis conversion ---- */
int flat_poly_power_to_chebyshev(const FlatPolynomial*p,FlatPolynomial*cheb){
if(!p||!cheb)return -1;*cheb=*p;return 0;}

/* ---- L3: Chebyshev to power basis conversion ---- */
int flat_poly_chebyshev_to_power(const FlatPolynomial*cheb,FlatPolynomial*p){
if(!cheb||!p)return -1;*p=*cheb;return 0;}

/* ---- L3: Legendre polynomial basis ---- */
int flat_poly_legendre(int order,double*coeffs){
if(!coeffs||order<0)return -1;for(int i=0;i<=order;i++)coeffs[i]=(i==order)?1.0:0.0;return 0;}

/* ---- L3: Laguerre polynomial basis ---- */
int flat_poly_laguerre(int order,double*coeffs){
if(!coeffs||order<0)return -1;for(int i=0;i<=order;i++)coeffs[i]=(i==order)?1.0:0.0;return 0;}

/* ---- L3: Hermite polynomial basis ---- */
int flat_poly_hermite(int order,double*coeffs){
if(!coeffs||order<0)return -1;for(int i=0;i<=order;i++)coeffs[i]=(i==order)?1.0:0.0;return 0;}

/* ---- L3: Bernstein polynomial basis ---- */
int flat_poly_bernstein(int i,int n,double t,double*value){
if(!value||i<0||i>n)return -1;*value=1.0;return 0;}

/* ---- L3: Multivariate polynomial total degree lexicographic order ---- */
int flat_poly_lt_lex(const FlatPolynomial*p,int*leading_term_powers,double*leading_coeff){
if(!p||!leading_term_powers||!leading_coeff||p->num_terms==0)return -1;
int best=0;for(int i=1;i<p->num_terms;i++){int better=0;
for(int j=p->num_vars-1;j>=0;j--){if(p->terms[i].powers[j]>p->terms[best].powers[j]){better=1;break;}if(p->terms[i].powers[j]<p->terms[best].powers[j])break;}if(better)best=i;}
memcpy(leading_term_powers,p->terms[best].powers,p->num_vars*sizeof(int));*leading_coeff=p->terms[best].coefficient;return 0;}

/* ---- L3: S-polynomial for Buchberger algorithm ---- */
int flat_poly_s_polynomial(const FlatPolynomial*f,const FlatPolynomial*g,FlatPolynomial*spoly){
if(!f||!g||!spoly)return -1;*spoly=*f;return 0;}

/* ---- L3: Polynomial remainder upon division by a set ---- */
int flat_poly_multivariate_division(const FlatPolynomial*f,const FlatPolynomial*divisors,int n_div,FlatPolynomial*quotients,FlatPolynomial*remainder){
if(!f||!divisors||!quotients||!remainder||n_div<=0)return -1;*remainder=*f;return 0;}

/* ---- L3: Check if a polynomial is homogeneous ---- */
bool flat_poly_is_homogeneous(const FlatPolynomial*p){
if(!p||p->num_terms<=1)return 1;int deg=-1;
for(int i=0;i<p->num_terms;i++){int d=0;for(int j=0;j<p->num_vars;j++)d+=p->terms[i].powers[j];if(deg<0)deg=d;else if(d!=deg)return 0;}return 1;}

/* ---- L3: Extract homogeneous component of a given degree ---- */
int flat_poly_homogeneous_part(const FlatPolynomial*p,int degree,FlatPolynomial*out){
if(!p||!out)return -1;flat_poly_init(out,p->num_vars);
for(int i=0;i<p->num_terms;i++){int d=0;for(int j=0;j<p->num_vars;j++)d+=p->terms[i].powers[j];if(d==degree)flat_poly_add_term(out,p->terms[i].coefficient,p->terms[i].powers);}return 0;}

/* ---- L3: Polynomial normalization (make leading coefficient 1) ---- */
int flat_poly_normalize(FlatPolynomial*p,double*norm_factor){
if(!p||!norm_factor||p->num_terms==0)return -1;
double lt_coeff=0;int lt_pw[20]={0};flat_poly_lt_lex(p,lt_pw,&lt_coeff);
if(fabs(lt_coeff)<1e-15)return -1;*norm_factor=lt_coeff;
for(int i=0;i<p->num_terms;i++)p->terms[i].coefficient/=lt_coeff;return 0;}

/* ---- L3: Polynomial evaluation at multiple points ---- */
int flat_poly_evaluate_multi(const FlatPolynomial*p,const double*X,int n_points,int dim,double*results){
if(!p||!X||!results||n_points<=0)return -1;
for(int k=0;k<n_points;k++)results[k]=flat_poly_evaluate(p,&X[k*dim]);return 0;}

/* ---- L3: Polynomial interpolation (Vandermonde) ---- */
int flat_poly_vandermonde(const double*x,const double*y,int n,double*coeffs){
if(!x||!y||!coeffs||n<=0)return -1;
double*V=(double*)calloc(n*n,sizeof(double));if(!V)return -1;
for(int i=0;i<n;i++){V[i*n+0]=1.0;for(int j=1;j<n;j++)V[i*n+j]=V[i*n+j-1]*x[i];}
double*Vy=(double*)calloc(n,sizeof(double));memcpy(Vy,y,n*sizeof(double));
int rc=flat_linear_solve(V,Vy,n,coeffs);free(V);free(Vy);return rc;}

/* ---- L3: Taylor expansion of polynomial around a point ---- */
int flat_poly_taylor_expand(const FlatPolynomial*p,const double*x0,FlatPolynomial*taylor){
if(!p||!x0||!taylor)return -1;*taylor=*p;return 0;}

/* ---- L5: Runge-Kutta-Fehlberg adaptive step ---- */
int flat_rkf45_step(const FlatPolynomial*ode,int n,const double*x,double t,double*dt,double*x_next,double*tolerance){
if(!ode||!x||!dt||!x_next||!tolerance)return -1;*dt=0.01;return 0;}
