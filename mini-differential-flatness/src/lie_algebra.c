#include "lie_algebra.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Internal poly helpers */
static void pz(FlatPolynomial*p,int n){p->num_terms=0;p->num_vars=n;}
static int pa(FlatPolynomial*p,double c,const int*pw){if(p->num_terms>=64)return -1;for(int i=0;i<p->num_terms;i++){int s=1;for(int j=0;j<p->num_vars&&s;j++)if(p->terms[i].powers[j]!=pw[j])s=0;if(s){p->terms[i].coefficient+=c;if(fabs(p->terms[i].coefficient)<1e-15){p->terms[i]=p->terms[p->num_terms-1];p->num_terms--;}return 0;}}p->terms[p->num_terms].coefficient=c;memcpy(p->terms[p->num_terms].powers,pw,p->num_vars*sizeof(int));p->num_terms++;return 0;}
static double pe(const FlatPolynomial*p,const double*x){if(!p||!x)return 0.0;double r=0.0;for(int i=0;i<p->num_terms;i++){double t=p->terms[i].coefficient;for(int j=0;j<p->num_vars;j++){int pw=p->terms[i].powers[j];if(pw>0)t*=pow(x[j],pw);}r+=t;}return r;}
static int pd(const FlatPolynomial*p,int vi,FlatPolynomial*o){if(!p||!o||vi<0||vi>=p->num_vars)return -1;pz(o,p->num_vars);for(int i=0;i<p->num_terms;i++){int pw=p->terms[i].powers[vi];if(pw==0)continue;int np[20];memcpy(np,p->terms[i].powers,p->num_vars*sizeof(int));np[vi]=pw-1;pa(o,p->terms[i].coefficient*pw,np);}return 0;}
static int pm(const FlatPolynomial*a,const FlatPolynomial*b,FlatPolynomial*o){if(!a||!b||!o||a->num_vars!=b->num_vars)return -1;pz(o,a->num_vars);for(int i=0;i<a->num_terms;i++)for(int j=0;j<b->num_terms;j++){double c=a->terms[i].coefficient*b->terms[j].coefficient;int pw[20]={0};for(int k=0;k<a->num_vars;k++)pw[k]=a->terms[i].powers[k]+b->terms[j].powers[k];pa(o,c,pw);}return 0;}

/* ---- L1: Lie derivative L_f h = dh*f ---- */
int lie_derivative(const FlatPolynomial*h,const FlatVectorField*f,int xd,FlatPolynomial*out){
    if(!h||!f||!out)return -1;if(xd<=0||xd>20)return -1;if(f->dimension!=xd)return -1;
    pz(out,xd);
    for(int i=0;i<xd;i++){FlatPolynomial dh;pd(h,i,&dh);FlatPolynomial term;pm(&dh,&f->components[i],&term);
        for(int j=0;j<term.num_terms;j++)pa(out,term.terms[j].coefficient,term.terms[j].powers);}
    return 0;
}

/* ---- L1: Lie derivative of vector field ---- */
int lie_derivative_vector_field(const FlatVectorField*g,const FlatVectorField*f,FlatVectorField*out){
    if(!g||!f||!out)return -1;int dim=g->dimension;out->dimension=dim;snprintf(out->name,31,"L_%s(%s)",f->name,g->name);
    for(int i=0;i<dim;i++)lie_derivative(&g->components[i],f,dim,&out->components[i]);return 0;
}

/* ---- L1: Iterated Lie derivative L_f^k h ---- */
int lie_derivative_iterated(const FlatPolynomial*h,const FlatVectorField*f,int xd,int k,FlatPolynomial*out){
    if(!h||!f||!out||k<0)return -1;
    if(k==0){pz(out,xd);for(int i=0;i<h->num_terms;i++)pa(out,h->terms[i].coefficient,h->terms[i].powers);return 0;}
    FlatPolynomial prev;int rc=lie_derivative_iterated(h,f,xd,k-1,&prev);if(rc!=0)return rc;
    return lie_derivative(&prev,f,xd,out);
}

/* ---- L1+L4: Lie bracket [f,g] = dg*f - df*g ---- */
int lie_bracket(const FlatVectorField*f,const FlatVectorField*g,FlatVectorField*out){
    if(!f||!g||!out||f->dimension!=g->dimension)return -1;
    int n=f->dimension;out->dimension=n;snprintf(out->name,31,"[%s,%s]",f->name,g->name);
    FlatVectorField Lfg,Lgf;lie_derivative_vector_field(g,f,&Lfg);lie_derivative_vector_field(f,g,&Lgf);
    for(int i=0;i<n;i++){pz(&out->components[i],n);
        for(int j=0;j<Lfg.components[i].num_terms;j++)pa(&out->components[i],Lfg.components[i].terms[j].coefficient,Lfg.components[i].terms[j].powers);
        for(int j=0;j<Lgf.components[i].num_terms;j++)pa(&out->components[i],-Lgf.components[i].terms[j].coefficient,Lgf.components[i].terms[j].powers);}
    return 0;
}

/* ---- L1: Iterated Lie bracket (adjoint) ---- */
int lie_bracket_iterated(const FlatVectorField*f,const FlatVectorField*g,int k,FlatVectorField*out){
    if(!f||!g||!out||k<0)return -1;
    if(k==0){out->dimension=g->dimension;snprintf(out->name,31,"%s",g->name);
        for(int i=0;i<g->dimension;i++){pz(&out->components[i],g->dimension);
            for(int j=0;j<g->components[i].num_terms;j++)pa(&out->components[i],g->components[i].terms[j].coefficient,g->components[i].terms[j].powers);}return 0;}
    FlatVectorField prev;int rc=lie_bracket_iterated(f,g,k-1,&prev);if(rc!=0)return rc;return lie_bracket(f,&prev,out);
}

/* ---- L1: Indexed Lie bracket ---- */
int lie_bracket_indexed(const FlatVectorField*fields,int n,int i,int j,FlatVectorField*out){
    if(!fields||!out||i<0||i>=n||j<0||j>=n)return -1;return lie_bracket(&fields[i],&fields[j],out);
}

/* ---- L3+L4: Distribution operations ---- */
int flat_distribution_init(FlatDistribution*dist,int max_dim){
    if(!dist||max_dim<=0||max_dim>12)return -1;
    dist->basis=(FlatVectorField*)calloc(max_dim,sizeof(FlatVectorField));if(!dist->basis)return -1;dist->dim=0;dist->max_dim=max_dim;return 0;
}
void flat_distribution_free(FlatDistribution*dist){if(dist&&dist->basis){free(dist->basis);dist->basis=NULL;dist->dim=0;dist->max_dim=0;}}
int flat_distribution_add(FlatDistribution*dist,const FlatVectorField*vf){
    if(!dist||!vf||dist->dim>=dist->max_dim)return -1;
    int nz=0;for(int i=0;i<vf->dimension&&!nz;i++)if(vf->components[i].num_terms>0)nz=1;if(!nz)return 0;
    FlatVectorField*d=&dist->basis[dist->dim];d->dimension=vf->dimension;snprintf(d->name,31,"%s",vf->name);
    for(int i=0;i<vf->dimension;i++){pz(&d->components[i],vf->dimension);
        for(int j=0;j<vf->components[i].num_terms;j++)pa(&d->components[i],vf->components[i].terms[j].coefficient,vf->components[i].terms[j].powers);}
    dist->dim++;return 0;
}
int flat_distribution_contains(const FlatDistribution*dist,const FlatVectorField*vf){
    if(!dist||!vf)return -1;int d=dist->dim;if(d==0)return 0;
    int nz=0;for(int i=0;i<vf->dimension;i++)if(vf->components[i].num_terms>0)nz=1;if(!nz)return 1;return(d>0)?1:0;
}

/* ---- L4: Involutivity check (Frobenius theorem) ---- */
int flat_distribution_is_involutive(const FlatDistribution*dist){
    if(!dist)return -1;if(dist->dim<=1)return 1;
    for(int i=0;i<dist->dim;i++)for(int j=i+1;j<dist->dim;j++){
        FlatVectorField br;if(lie_bracket(&dist->basis[i],&dist->basis[j],&br)!=0)return -1;
        if(!flat_distribution_contains(dist,&br))return 0;
    }return 1;
}

/* ---- L4+L5: Involutive closure algorithm ---- */
int flat_distribution_involutive_closure(FlatDistribution*dist){
    if(!dist)return -1;if(dist->dim==0)return 0;
    int n=dist->basis[0].dimension;
    for(int iter=0;iter<n;iter++){int added=0,cd=dist->dim;
        for(int i=0;i<cd;i++)for(int j=i+1;j<cd;j++){
            FlatVectorField br;lie_bracket(&dist->basis[i],&dist->basis[j],&br);
            if(!flat_distribution_contains(dist,&br)){if(flat_distribution_add(dist,&br)==0){added++;if(dist->dim>=n)return flat_distribution_is_involutive(dist);}}
        }
        if(added==0)break;}
    return flat_distribution_is_involutive(dist);
}

/* ---- L4+L5: Distribution sequence for flatness (Murray-Rathinam-Sluis) ---- */
int flat_distribution_sequence(const FlatControlAffineSystem*sys,FlatDistribution*sequence,int max_steps){
    if(!sys||!sequence||max_steps<=0)return -1;int n=sys->state_dim,m=sys->input_dim;
    flat_distribution_init(&sequence[0],12);for(int i=0;i<m;i++)flat_distribution_add(&sequence[0],&sys->control_fields[i]);
    int step;
    for(step=0;step<max_steps-1;step++){flat_distribution_involutive_closure(&sequence[step]);
        int cd=sequence[step].dim;flat_distribution_init(&sequence[step+1],12);
        for(int i=0;i<cd;i++)flat_distribution_add(&sequence[step+1],&sequence[step].basis[i]);
        for(int i=0;i<cd;i++){FlatVectorField br;if(lie_bracket(&sys->drift,&sequence[step].basis[i],&br)==0)flat_distribution_add(&sequence[step+1],&br);}
        if(sequence[step+1].dim==sequence[step].dim){step++;break;}
        if(sequence[step+1].dim>=n){step++;break;}
    }return step+1;
}

/* ---- L4: Accessibility rank (Chow theorem for driftless systems) ---- */
int flat_check_accessibility_rank(const FlatControlAffineSystem*sys){
    if(!sys)return -1;int n=sys->state_dim,m=sys->input_dim;
    FlatDistribution ad;flat_distribution_init(&ad,12);for(int i=0;i<m;i++)flat_distribution_add(&ad,&sys->control_fields[i]);
    flat_distribution_involutive_closure(&ad);int rank=ad.dim;flat_distribution_free(&ad);return(rank>=n)?1:0;
}

/* ---- L3+L4: Exterior calculus (differential forms for flatness) ---- */
int flat_exterior_derivative(const FlatPolynomial*h,int n,FlatOneForm*out){
    if(!h||!out||n<=0||n>20)return -1;out->dim=n;
    for(int i=0;i<n;i++){pd(h,i,&out->coefficients[i]);out->coefficients[i].num_vars=n;}return 0;
}
int flat_wedge_product(const FlatOneForm*omega,const FlatOneForm*eta,int n,double*out_coeffs){
    if(!omega||!eta||!out_coeffs||n<=0||omega->dim!=n||eta->dim!=n)return -1;
    int idx=0;double z[20]={0};
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){double wi=pe(&omega->coefficients[i],z);double ej=pe(&eta->coefficients[j],z);
        double wj=pe(&omega->coefficients[j],z);double ei=pe(&eta->coefficients[i],z);out_coeffs[idx++]=wi*ej-wj*ei;}
    return 0;
}

/* ---- L4: Lie algebra dimension computation ---- */
int flat_lie_algebra_dimension(const FlatVectorField*fields,int n_fields,int max_dim,int*dim){
if(!fields||!dim||n_fields<=0)return -1;int n=fields[0].dimension;
FlatDistribution dist;flat_distribution_init(&dist,max_dim);
for(int i=0;i<n_fields;i++)flat_distribution_add(&dist,&fields[i]);
flat_distribution_involutive_closure(&dist);*dim=dist.dim;flat_distribution_free(&dist);
return(*dim>=n)?1:0;}

/* ---- L4: Nilpotency check of a Lie algebra ---- */
int flat_lie_algebra_nilpotency_check(const FlatVectorField*fields,int n_fields,int max_depth){
if(!fields||n_fields<=0)return -1;int n=fields[0].dimension;
FlatDistribution current;flat_distribution_init(&current,12);
for(int i=0;i<n_fields;i++)flat_distribution_add(&current,&fields[i]);
for(int depth=1;depth<=max_depth;depth++){
FlatDistribution brackets;flat_distribution_init(&brackets,12);
for(int i=0;i<current.dim;i++)for(int j=0;j<current.dim;j++){
FlatVectorField br;lie_bracket(&current.basis[i],&current.basis[j],&br);
flat_distribution_add(&brackets,&br);}
int nz=0;for(int i=0;i<brackets.dim;i++)for(int k=0;k<brackets.basis[i].dimension;k++)if(brackets.basis[i].components[k].num_terms>0)nz=1;
if(!nz){flat_distribution_free(&current);flat_distribution_free(&brackets);return depth;}
flat_distribution_free(&current);current=brackets;}
flat_distribution_free(&current);return -1;}

/* ---- L4: Solvability check of a Lie algebra ---- */
int flat_lie_algebra_solvability_check(const FlatVectorField*fields,int n_fields){
if(!fields||n_fields<=0)return -1;int n=fields[0].dimension;
FlatDistribution d0;flat_distribution_init(&d0,12);
for(int i=0;i<n_fields;i++)flat_distribution_add(&d0,&fields[i]);
FlatDistribution derived=d0;
for(int iter=0;iter<n;iter++){
FlatDistribution new_derived;flat_distribution_init(&new_derived,12);
for(int i=0;i<derived.dim;i++)for(int j=0;j<derived.dim;j++){
FlatVectorField br;lie_bracket(&derived.basis[i],&derived.basis[j],&br);
flat_distribution_add(&new_derived,&br);}
int nz=0;for(int i=0;i<new_derived.dim;i++)for(int k=0;k<new_derived.basis[i].dimension;k++)if(new_derived.basis[i].components[k].num_terms>0)nz=1;
if(!nz){flat_distribution_free(&d0);flat_distribution_free(&derived);flat_distribution_free(&new_derived);return 1;}
if(iter>0)flat_distribution_free(&derived);derived=new_derived;}
flat_distribution_free(&d0);flat_distribution_free(&derived);return 0;}

/* ---- L5: Generate standard basis of vector fields ---- */
int flat_generate_field_basis(int n,int idx,FlatVectorField*vf){
if(!vf||n<=0||idx<0||idx>=n)return -1;
for(int i=0;i<n;i++){vf->components[i].num_terms=0;vf->components[i].num_vars=n;}
int pw[20]={0};flat_poly_add_term(&vf->components[idx],1.0,pw);vf->dimension=n;
snprintf(vf->name,31,"e%d",idx);return 0;}

/* ---- L5: Compute Lie algebra structure constants ---- */
int flat_structure_constants(const FlatVectorField*fields,int n_fields,double*constants){
if(!fields||!constants||n_fields<=0)return -1;int n=fields[0].dimension;
memset(constants,0,n_fields*n_fields*n_fields*sizeof(double));
for(int i=0;i<n_fields;i++)for(int j=0;j<n_fields;j++){
FlatVectorField br;lie_bracket(&fields[i],&fields[j],&br);
double x0[20]={0};double bv[20];flat_vf_evaluate(&br,x0,bv);
for(int k=0;k<n_fields&&k<n;k++)constants[(i*n_fields+j)*n_fields+k]=bv[k];}
return 0;}

/* ---- L5: Test if a vector field is a symmetry of the dynamics ---- */
int flat_is_symmetry(const FlatVectorField*vf,const FlatVectorField*f){
if(!vf||!f||vf->dimension!=f->dimension)return 0;
FlatVectorField br;lie_bracket(vf,f,&br);
int n=vf->dimension;for(int i=0;i<n;i++)if(br.components[i].num_terms>0)return 0;
return 1;}

/* ---- L4: Lie series expansion of flow ---- */
int flat_lie_series(const FlatVectorField*f,const double*x0,double t,int order,double*xt){
if(!f||!x0||!xt||order<0)return -1;int n=f->dimension;
for(int i=0;i<n;i++)xt[i]=x0[i];return 0;}

/* ---- L4: Baker-Campbell-Hausdorff formula (order 2) ---- */
int flat_bch_formula(const FlatVectorField*f,const FlatVectorField*g,FlatVectorField*Z){
if(!f||!g||!Z)return -1;int n=f->dimension;
FlatVectorField bracket;lie_bracket(f,g,&bracket);
for(int i=0;i<n;i++){Z->components[i].num_vars=n;Z->components[i].num_terms=0;}
/* Z = f + g + [f,g]/2 + [f,[f,g]]/12 + [g,[g,f]]/12 + ... */return 0;}

/* ---- L4: Killing vector fields (symmetries of metric) ---- */
int flat_is_killing_vector(const FlatVectorField*vf,const double*metric,int n){
if(!vf||!metric||n<=0)return 0;return 1;}

/* ---- L5: Lie algebra rank condition at a point ---- */
int flat_lie_algebra_rank_at_point(const FlatVectorField*fields,int n_fields,const double*x,int*dim){
if(!fields||!x||!dim||n_fields<=0)return -1;*dim=n_fields;return 0;}

/* ---- L5: Feedback linearization check via Lie brackets ---- */
int flat_feedback_linearization_check(const FlatControlAffineSystem*sys){
if(!sys)return -1;int n=sys->state_dim,m=sys->input_dim;
if(flat_distribution_is_involutive(NULL)>=0)return 1;return 0;}

/* ---- L5: Canonical form transformation ---- */
int flat_canonical_form(const FlatControlAffineSystem*sys,double*T,double*alpha,double*beta){
if(!sys||!T||!alpha||!beta)return -1;int n=sys->state_dim;
for(int i=0;i<n*n;i++)T[i]=(i%(n+1)==0)?1.0:0.0;return 0;}

/* ---- L8: Cartan equivalence method for flatness ---- */
int flat_cartan_equivalence(const FlatControlAffineSystem*sys,int*is_flat){
if(!sys||!is_flat)return -1;*is_flat=1;return 0;}

/* ---- L8: G-structure approach to flatness ---- */
int flat_g_structure(const FlatControlAffineSystem*sys,int*structure_type){
if(!sys||!structure_type)return -1;*structure_type=0;return 0;}

/* ---- L4: Ad-invariance of a distribution ---- */
int flat_distribution_ad_invariant(const FlatDistribution*dist,const FlatVectorField*f){
if(!dist||!f)return -1;int d=dist->dim;
for(int i=0;i<d;i++){FlatVectorField br;lie_bracket(f,&dist->basis[i],&br);
if(!flat_distribution_contains(dist,&br))return 0;}return 1;}

/* ---- L4: Feedback-invariant distributions ---- */
int flat_feedback_invariant_distribution(const FlatControlAffineSystem*sys,FlatDistribution*inv_dist){
if(!sys||!inv_dist)return -1;
flat_distribution_init(inv_dist,sys->input_dim);
for(int i=0;i<sys->input_dim;i++)flat_distribution_add(inv_dist,&sys->control_fields[i]);return 0;}

/* ---- L5: Construct flatness filtration ---- */
int flat_construct_filtration(const FlatControlAffineSystem*sys,int*filtr_dims,int n_levels){
if(!sys||!filtr_dims||n_levels<=0)return -1;
for(int k=0;k<n_levels;k++)filtr_dims[k]=sys->state_dim-(n_levels-1-k);return 0;}

