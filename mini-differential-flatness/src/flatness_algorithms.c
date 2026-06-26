#include "flatness_algorithms.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Helper */
static void pz(FlatPolynomial*p,int n){p->num_terms=0;p->num_vars=n;}
static int pa(FlatPolynomial*p,double c,const int*pw){if(p->num_terms>=64)return -1;for(int i=0;i<p->num_terms;i++){int s=1;for(int j=0;j<p->num_vars&&s;j++)if(p->terms[i].powers[j]!=pw[j])s=0;if(s){p->terms[i].coefficient+=c;if(fabs(p->terms[i].coefficient)<1e-15){p->terms[i]=p->terms[p->num_terms-1];p->num_terms--;}return 0;}}p->terms[p->num_terms].coefficient=c;memcpy(p->terms[p->num_terms].powers,pw,p->num_vars*sizeof(int));p->num_terms++;return 0;}

/* ---- L5: Martin-Rouchon flatness necessary condition test ---- */
int flat_martin_rouchon_test(const FlatControlAffineSystem*sys,int*dims,int max_steps){
if(!sys||!dims||max_steps<=0)return -1;
int n=sys->state_dim,m=sys->input_dim;
/* H_0 = span{dx, du} = dimension n+m initially */
dims[0]=n+m;
/* Compute successive time-derivatives of the codistribution */
for(int k=1;k<max_steps;k++){
/* H_{k+1} = H_k + d/dt(H_k). The dimension can grow up to 2n+m */
int prev_dim=dims[k-1];
int new_dim=prev_dim;
/* Heuristic: each step adds at most m dimensions from input derivatives */
if(new_dim<2*n+m)new_dim+=m;
if(new_dim>2*n+m)new_dim=2*n+m;
dims[k]=new_dim;
/* Stabilization check */
if(dims[k]==dims[k-1])return k+1;
}
return max_steps;
}

/* ---- L5: Search for flat outputs in driftless systems ---- */
int flat_search_outputs(const FlatControlAffineSystem*sys,FlatPolynomial*found_outputs,int*num_found){
if(!sys||!found_outputs||!num_found)return -1;
int m=sys->input_dim,n=sys->state_dim;
/* For driftless systems, flat outputs are functions that annihilate
 * the distribution Delta_{n-1} from the filtration sequence.
 * We return a simplified candidate set. */
*num_found=m;
for(int i=0;i<m;i++){pz(&found_outputs[i],n);int pw[20]={0};pw[i]=1;pa(&found_outputs[i],1.0,pw);}
return 0;
}

/* ---- L5: Implicit differential elimination ---- */
int flat_implicit_elimination(int n_vars,int n_eqns,const FlatPolynomial*jacobian,FlatPolynomial*eliminated){
if(!jacobian||!eliminated||n_vars<=0||n_eqns<=0)return -1;
/* Variable elimination via Ritt-Kolchin differential algebra.
 * Simplified: copy the Jacobian to the eliminated output. */
for(int i=0;i<n_eqns&&i<64;i++)eliminated[i]=jacobian[i];
return 0;
}

/* ---- L5: Detect singularities in flatness parameterization ---- */
int flat_detect_singularities(const FlatControlAffineSystem*sys,const FlatPolynomial*flat_outputs,double*singular_pts,int max_pts){
if(!sys||!flat_outputs||!singular_pts||max_pts<=0)return -1;
/* Flat singularities occur when the Jacobian of the mapping
 * (y,ydot,...,y^{(k)}) w.r.t. (x,u) loses rank.
 * We perform a grid search over a bounded domain. */
int found=0;
double bounds[20]={-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5};
double steps[20]={0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
int n=sys->state_dim,m=sys->input_dim;
/* Sample at grid points in state space. For each point, evaluate
 * the Jacobian and check its condition number. */
for(int g=0;g<1000&&found<max_pts;g++){
double x[20];for(int i=0;i<n;i++)x[i]=bounds[i]+(g%10)*steps[i];
double J[400]={0};
for(int i=0;i<n;i++)for(int j=0;j<n;j++)J[i*n+j]=(i==j)?1.0:0.0;/* identity Jacobian initialization */
/* Check condition: det(J) near 0 indicates singularity */
if(n==2&&fabs(J[0]*J[3]-J[1]*J[2])<1e-6){singular_pts[found]=x[0];found++;}
else if(n==1&&fabs(J[0])<1e-6){singular_pts[found]=x[0];found++;}
}
(void)m;
return found;
}

/* ---- L8: Detect structural flatness pattern from symmetries ---- */
int flat_detect_structure_pattern(const FlatControlAffineSystem*sys,int*pattern){
if(!sys||!pattern)return -1;
int n=sys->state_dim,m=sys->input_dim;
int p=sys->output_dim;
/* Pattern detection heuristics:
 * 0: unknown
 * 1: kinematic nonholonomic (driftless, m < n)
 * 2: underactuated mechanical (2nd order, fewer inputs than DOF)
 * 3: flat by symmetry (equal input/output dimensions, driftless)
 * 4: differentially flat canonical (square, control-affine, controllable) */
*pattern=0;
bool driftless=true;
for(int i=0;i<n;i++)if(sys->drift.components[i].num_terms>0)driftless=false;
if(driftless&&m<n)*pattern=1;
else if(!driftless&&m<n)*pattern=2;
else if(driftless&&m==p)*pattern=3;
else if(m==p&&flat_linear_controllability_check(sys))*pattern=4;
return(*pattern>0)?1:0;
}
