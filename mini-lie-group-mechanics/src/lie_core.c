/* ==============================================================
 * lie_core.c -- Lie Group and Lie Algebra Core Implementation
 *
 * Implements: matrix ops, vector ops, Lie group/algebra element API,
 * bracket, Jacobi identity, hat/vee maps, exponential/log maps
 * (Rodrigues, scaling+squaring), adjoint/coadjoint, Killing form,
 * BCH formula to 3rd order.
 *
 * Reference: Marsden & Ratiu (1999), Murray et al. (1994),
 *            Hall (2015), Higham (2009)
 * ============================================================== */

#include "lie_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define LIE_EPS 1e-15

static double lie_abs(double x) { return (x < 0) ? -x : x; }

/* ---- Matrix create/destroy ---- */
LieMatrix* lie_matrix_create(int rows, int cols) {
    assert(rows > 0 && cols > 0);
    LieMatrix *m = (LieMatrix*)malloc(sizeof(LieMatrix));
    assert(m);
    m->data = (double*)calloc((size_t)(rows * cols), sizeof(double));
    assert(m->data);
    m->rows = rows; m->cols = cols;
    return m;
}
void lie_matrix_free(LieMatrix *m) { if(m){free(m->data);free(m);} }
LieMatrix* lie_matrix_copy(const LieMatrix *src) {
    assert(src);
    LieMatrix *d = lie_matrix_create(src->rows, src->cols);
    memcpy(d->data, src->data, (size_t)(src->rows*src->cols)*sizeof(double));
    return d;
}
void lie_matrix_set(LieMatrix *m, int row, int col, double val) {
    assert(m && row>=0 && row<m->rows && col>=0 && col<m->cols);
    m->data[row*m->cols+col] = val;
}
double lie_matrix_get(const LieMatrix *m, int row, int col) {
    assert(m && row>=0 && row<m->rows && col>=0 && col<m->cols);
    return m->data[row*m->cols+col];
}
LieMatrix* lie_matrix_identity(int n) {
    LieMatrix *m = lie_matrix_create(n,n);
    for(int i=0;i<n;i++) m->data[i*n+i]=1.0;
    return m;
}
LieMatrix* lie_matrix_zeros(int rows, int cols) { return lie_matrix_create(rows,cols); }
LieMatrix* lie_matrix_diag(int n, const double *diag) {
    LieMatrix *m = lie_matrix_create(n,n);
    for(int i=0;i<n;i++) m->data[i*n+i]=diag[i];
    return m;
}
void lie_matrix_fill(LieMatrix *m, double val) {
    assert(m);
    int sz=m->rows*m->cols;
    for(int i=0;i<sz;i++) m->data[i]=val;
}
LieMatrix* lie_matrix_add(const LieMatrix *a, const LieMatrix *b) {
    assert(a&&b&&a->rows==b->rows&&a->cols==b->cols);
    LieMatrix *c=lie_matrix_create(a->rows,a->cols);
    int sz=a->rows*a->cols;
    for(int i=0;i<sz;i++) c->data[i]=a->data[i]+b->data[i];
    return c;
}
LieMatrix* lie_matrix_sub(const LieMatrix *a, const LieMatrix *b) {
    assert(a&&b&&a->rows==b->rows&&a->cols==b->cols);
    LieMatrix *c=lie_matrix_create(a->rows,a->cols);
    int sz=a->rows*a->cols;
    for(int i=0;i<sz;i++) c->data[i]=a->data[i]-b->data[i];
    return c;
}
LieMatrix* lie_matrix_mul(const LieMatrix *a, const LieMatrix *b) {
    assert(a&&b&&a->cols==b->rows);
    int n=a->rows,m=b->cols,inner=a->cols;
    LieMatrix *c=lie_matrix_create(n,m);
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++){
            double s=0.0;
            for(int k=0;k<inner;k++)
                s+=a->data[i*a->cols+k]*b->data[k*b->cols+j];
            c->data[i*m+j]=s;
        }
    return c;
}
LieMatrix* lie_matrix_scale(const LieMatrix *m, double s) {
    assert(m);
    LieMatrix *c=lie_matrix_create(m->rows,m->cols);
    int sz=m->rows*m->cols;
    for(int i=0;i<sz;i++) c->data[i]=m->data[i]*s;
    return c;
}
LieMatrix* lie_matrix_transpose(const LieMatrix *m) {
    assert(m);
    LieMatrix *t=lie_matrix_create(m->cols,m->rows);
    for(int i=0;i<m->rows;i++)
        for(int j=0;j<m->cols;j++)
            t->data[j*m->rows+i]=m->data[i*m->cols+j];
    return t;
}
double lie_matrix_trace(const LieMatrix *m) {
    assert(m&&m->rows==m->cols);
    double tr=0.0;
    for(int i=0;i<m->rows;i++) tr+=m->data[i*m->cols+i];
    return tr;
}
double lie_matrix_det_2x2(const LieMatrix *m) {
    assert(m&&m->rows==2&&m->cols==2);
    return m->data[0]*m->data[3]-m->data[1]*m->data[2];
}
double lie_matrix_det_3x3(const LieMatrix *m) {
    assert(m&&m->rows==3&&m->cols==3);
    double *d=m->data;
    return d[0]*(d[4]*d[8]-d[5]*d[7])-d[1]*(d[3]*d[8]-d[5]*d[6])+d[2]*(d[3]*d[7]-d[4]*d[6]);
}
double lie_matrix_det_4x4(const LieMatrix *m) {
    assert(m&&m->rows==4&&m->cols==4);
    double *a=m->data;
    double m00=a[5]*(a[10]*a[15]-a[11]*a[14])-a[6]*(a[9]*a[15]-a[11]*a[13])+a[7]*(a[9]*a[14]-a[10]*a[13]);
    double m01=a[4]*(a[10]*a[15]-a[11]*a[14])-a[6]*(a[8]*a[15]-a[11]*a[12])+a[7]*(a[8]*a[14]-a[10]*a[12]);
    double m02=a[4]*(a[9]*a[15]-a[11]*a[13])-a[5]*(a[8]*a[15]-a[11]*a[12])+a[7]*(a[8]*a[13]-a[9]*a[12]);
    double m03=a[4]*(a[9]*a[14]-a[10]*a[13])-a[5]*(a[8]*a[14]-a[10]*a[12])+a[6]*(a[8]*a[13]-a[9]*a[12]);
    return a[0]*m00-a[1]*m01+a[2]*m02-a[3]*m03;
}
double lie_matrix_det(const LieMatrix *m) {
    assert(m&&m->rows==m->cols);
    if(m->rows==2)return lie_matrix_det_2x2(m);
    if(m->rows==3)return lie_matrix_det_3x3(m);
    if(m->rows==4)return lie_matrix_det_4x4(m);
    int n=m->rows;
    LieMatrix *LU=lie_matrix_copy(m);
    double det=1.0;
    for(int k=0;k<n-1;k++){
        double mx=lie_abs(LU->data[k*n+k]);int mr=k;
        for(int i=k+1;i<n;i++){double v=lie_abs(LU->data[i*n+k]);if(v>mx){mx=v;mr=i;}}
        if(mx<LIE_EPS){lie_matrix_free(LU);return 0.0;}
        if(mr!=k){
            det=-det;
            for(int j=0;j<n;j++){double t=LU->data[k*n+j];LU->data[k*n+j]=LU->data[mr*n+j];LU->data[mr*n+j]=t;}
        }
        for(int i=k+1;i<n;i++){
            LU->data[i*n+k]/=LU->data[k*n+k];
            for(int j=k+1;j<n;j++)LU->data[i*n+j]-=LU->data[i*n+k]*LU->data[k*n+j];
        }
    }
    for(int i=0;i<n;i++)det*=LU->data[i*n+i];
    lie_matrix_free(LU);
    return det;
}
double lie_matrix_norm_frobenius(const LieMatrix *m) {
    assert(m);
    double s=0.0;int sz=m->rows*m->cols;
    for(int i=0;i<sz;i++)s+=m->data[i]*m->data[i];
    return sqrt(s);
}
double lie_matrix_norm_infinity(const LieMatrix *m) {
    assert(m);
    double mr=0.0;
    for(int i=0;i<m->rows;i++){
        double rs=0.0;
        for(int j=0;j<m->cols;j++)rs+=lie_abs(m->data[i*m->cols+j]);
        if(rs>mr)mr=rs;
    }
    return mr;
}
bool lie_matrix_is_identity(const LieMatrix *m, double tol) {
    if(!m||m->rows!=m->cols)return false;
    int n=m->rows;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){
        double e=(i==j)?1.0:0.0;
        if(lie_abs(m->data[i*n+j]-e)>tol)return false;
    }
    return true;
}
bool lie_matrix_is_skew_symmetric(const LieMatrix *m, double tol) {
    if(!m||m->rows!=m->cols)return false;
    int n=m->rows;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)
        if(lie_abs(m->data[i*n+j]+m->data[j*n+i])>tol)return false;
    return true;
}
bool lie_matrix_is_orthogonal(const LieMatrix *m, double tol) {
    if(!m||m->rows!=m->cols)return false;
    int n=m->rows;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){
        double dot=0.0;
        for(int k=0;k<n;k++)dot+=m->data[i*n+k]*m->data[j*n+k];
        double e=(i==j)?1.0:0.0;
        if(lie_abs(dot-e)>tol)return false;
    }
    return true;
}
bool lie_matrix_is_special_orthogonal(const LieMatrix *m, double tol) {
    return lie_matrix_is_orthogonal(m,tol)&&lie_abs(lie_matrix_det(m)-1.0)<tol;
}
bool lie_matrix_approx_equal(const LieMatrix *a, const LieMatrix *b, double tol) {
    if(!a||!b||a->rows!=b->rows||a->cols!=b->cols)return false;
    int sz=a->rows*a->cols;
    for(int i=0;i<sz;i++)if(lie_abs(a->data[i]-b->data[i])>tol)return false;
    return true;
}
LieMatrix* lie_matrix_inverse_2x2(const LieMatrix *m) {
    double det=lie_matrix_det_2x2(m);assert(lie_abs(det)>LIE_EPS);
    LieMatrix *inv=lie_matrix_create(2,2);
    double id=1.0/det;
    inv->data[0]=m->data[3]*id;inv->data[1]=-m->data[1]*id;
    inv->data[2]=-m->data[2]*id;inv->data[3]=m->data[0]*id;
    return inv;
}
LieMatrix* lie_matrix_inverse_3x3(const LieMatrix *m) {
    double det=lie_matrix_det_3x3(m);assert(lie_abs(det)>LIE_EPS);
    LieMatrix *inv=lie_matrix_create(3,3);
    double *a=m->data,*d=inv->data,id=1.0/det;
    d[0]=(a[4]*a[8]-a[5]*a[7])*id;d[1]=(a[2]*a[7]-a[1]*a[8])*id;d[2]=(a[1]*a[5]-a[2]*a[4])*id;
    d[3]=(a[5]*a[6]-a[3]*a[8])*id;d[4]=(a[0]*a[8]-a[2]*a[6])*id;d[5]=(a[2]*a[3]-a[0]*a[5])*id;
    d[6]=(a[3]*a[7]-a[4]*a[6])*id;d[7]=(a[1]*a[6]-a[0]*a[7])*id;d[8]=(a[0]*a[4]-a[1]*a[3])*id;
    return inv;
}
LieMatrix* lie_matrix_inverse_4x4(const LieMatrix *m) {
    double det=lie_matrix_det_4x4(m);assert(lie_abs(det)>LIE_EPS);
    int n=4;
    LieMatrix *aug=lie_matrix_create(n,2*n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)aug->data[i*2*n+j]=m->data[i*n+j];
        aug->data[i*2*n+n+i]=1.0;
    }
    for(int k=0;k<n;k++){
        double pv=aug->data[k*2*n+k];
        if(lie_abs(pv)<LIE_EPS){
            int sw=-1;
            for(int i=k+1;i<n;i++)if(lie_abs(aug->data[i*2*n+k])>LIE_EPS){sw=i;break;}
            assert(sw>=0);
            for(int j=0;j<2*n;j++){
                double t=aug->data[k*2*n+j];
                aug->data[k*2*n+j]=aug->data[sw*2*n+j];
                aug->data[sw*2*n+j]=t;
            }
            pv=aug->data[k*2*n+k];
        }
        for(int j=0;j<2*n;j++)aug->data[k*2*n+j]/=pv;
        for(int i=0;i<n;i++){
            if(i==k)continue;
            double f=aug->data[i*2*n+k];
            for(int j=0;j<2*n;j++)aug->data[i*2*n+j]-=f*aug->data[k*2*n+j];
        }
    }
    LieMatrix *inv=lie_matrix_create(n,n);
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)inv->data[i*n+j]=aug->data[i*2*n+n+j];
    lie_matrix_free(aug);
    return inv;
}
LieMatrix* lie_matrix_inverse(const LieMatrix *m) {
    if(m->rows==2)return lie_matrix_inverse_2x2(m);
    if(m->rows==3)return lie_matrix_inverse_3x3(m);
    if(m->rows==4)return lie_matrix_inverse_4x4(m);
    int n=m->rows;
    assert(lie_abs(lie_matrix_det(m))>LIE_EPS);
    LieMatrix *LU=lie_matrix_copy(m);
    LieMatrix *inv=lie_matrix_identity(n);
    int *piv=(int*)malloc((size_t)n*sizeof(int));
    assert(piv);
    for(int i=0;i<n;i++)piv[i]=i;
    for(int k=0;k<n-1;k++){
        double mx=lie_abs(LU->data[k*n+k]);int mr=k;
        for(int i=k+1;i<n;i++){double v=lie_abs(LU->data[i*n+k]);if(v>mx){mx=v;mr=i;}}
        if(mr!=k){
            int t=piv[k];piv[k]=piv[mr];piv[mr]=t;
            for(int j=0;j<n;j++){
                double tv=LU->data[k*n+j];LU->data[k*n+j]=LU->data[mr*n+j];LU->data[mr*n+j]=tv;
            }
        }
        for(int i=k+1;i<n;i++){
            LU->data[i*n+k]/=LU->data[k*n+k];
            for(int j=k+1;j<n;j++)LU->data[i*n+j]-=LU->data[i*n+k]*LU->data[k*n+j];
        }
    }
    double *b=(double*)malloc((size_t)n*sizeof(double));
    double *x=(double*)malloc((size_t)n*sizeof(double));
    for(int col=0;col<n;col++){
        for(int i=0;i<n;i++)b[i]=(i==col)?1.0:0.0;
        for(int i=0;i<n;i++){x[i]=b[piv[i]];for(int j=0;j<i;j++)x[i]-=LU->data[i*n+j]*x[j];}
        for(int i=n-1;i>=0;i--){for(int j=i+1;j<n;j++)x[i]-=LU->data[i*n+j]*x[j];x[i]/=LU->data[i*n+i];}
        for(int i=0;i<n;i++)inv->data[i*n+col]=x[i];
    }
    free(b);free(x);free(piv);lie_matrix_free(LU);
    return inv;
}
void lie_matrix_solve_linear(LieMatrix *A, double *b, int n, double *x) {
    assert(A&&A->rows==n&&A->cols==n);
    double *Ab=(double*)malloc((size_t)(n*(n+1))*sizeof(double));
    assert(Ab);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)Ab[i*(n+1)+j]=A->data[i*n+j];
        Ab[i*(n+1)+n]=b[i];
    }
    for(int k=0;k<n;k++){
        double mx=lie_abs(Ab[k*(n+1)+k]);int mr=k;
        for(int i=k+1;i<n;i++){double v=lie_abs(Ab[i*(n+1)+k]);if(v>mx){mx=v;mr=i;}}
        if(mr!=k)for(int j=0;j<=n;j++){
            double t=Ab[k*(n+1)+j];Ab[k*(n+1)+j]=Ab[mr*(n+1)+j];Ab[mr*(n+1)+j]=t;
        }
        for(int i=k+1;i<n;i++){
            double f=Ab[i*(n+1)+k]/Ab[k*(n+1)+k];
            for(int j=k;j<=n;j++)Ab[i*(n+1)+j]-=f*Ab[k*(n+1)+j];
        }
    }
    for(int i=n-1;i>=0;i--){
        x[i]=Ab[i*(n+1)+n];
        for(int j=i+1;j<n;j++)x[i]-=Ab[i*(n+1)+j]*x[j];
        x[i]/=Ab[i*(n+1)+i];
    }
    free(Ab);
}

/* ---- Vector operations ---- */
LieVector* lie_vector_create(int size) {
    assert(size>0);
    LieVector *v=(LieVector*)malloc(sizeof(LieVector));
    assert(v);
    v->data=(double*)calloc((size_t)size,sizeof(double));
    assert(v->data);
    v->size=size;
    return v;
}
void lie_vector_free(LieVector *v){if(v){free(v->data);free(v);}}
void lie_vector_set(LieVector *v,int idx,double val){assert(v&&idx>=0&&idx<v->size);v->data[idx]=val;}
double lie_vector_get(const LieVector *v,int idx){assert(v&&idx>=0&&idx<v->size);return v->data[idx];}
LieVector* lie_vector_copy(const LieVector *src){
    assert(src);
    LieVector *d=lie_vector_create(src->size);
    memcpy(d->data,src->data,(size_t)src->size*sizeof(double));
    return d;
}
double lie_vector_norm(const LieVector *v){return sqrt(lie_vector_norm_sq(v));}
double lie_vector_norm_sq(const LieVector *v){
    assert(v);
    double s=0.0;
    for(int i=0;i<v->size;i++)s+=v->data[i]*v->data[i];
    return s;
}
double lie_vector_dot(const LieVector *a,const LieVector *b){
    assert(a&&b&&a->size==b->size);
    double s=0.0;
    for(int i=0;i<a->size;i++)s+=a->data[i]*b->data[i];
    return s;
}
LieVector* lie_vector_cross3(const LieVector *a,const LieVector *b){
    assert(a&&b&&a->size==3&&b->size==3);
    LieVector *c=lie_vector_create(3);
    c->data[0]=a->data[1]*b->data[2]-a->data[2]*b->data[1];
    c->data[1]=a->data[2]*b->data[0]-a->data[0]*b->data[2];
    c->data[2]=a->data[0]*b->data[1]-a->data[1]*b->data[0];
    return c;
}
LieVector* lie_vector_add(const LieVector *a,const LieVector *b){
    assert(a&&b&&a->size==b->size);
    LieVector *c=lie_vector_create(a->size);
    for(int i=0;i<a->size;i++)c->data[i]=a->data[i]+b->data[i];
    return c;
}
LieVector* lie_vector_sub(const LieVector *a,const LieVector *b){
    assert(a&&b&&a->size==b->size);
    LieVector *c=lie_vector_create(a->size);
    for(int i=0;i<a->size;i++)c->data[i]=a->data[i]-b->data[i];
    return c;
}
LieVector* lie_vector_scale(const LieVector *v,double s){
    assert(v);
    LieVector *c=lie_vector_create(v->size);
    for(int i=0;i<v->size;i++)c->data[i]=v->data[i]*s;
    return c;
}
void lie_vector_normalize(LieVector *v){
    double n=lie_vector_norm(v);
    assert(n>LIE_EPS);
    for(int i=0;i<v->size;i++)v->data[i]/=n;
}

/* ---- Lie Group Element ---- */
LieGroupElement* lie_group_create(LieGroupType gtype, const char *name) {
    LieGroupElement *g=(LieGroupElement*)malloc(sizeof(LieGroupElement));
    assert(g);
    g->group_type=gtype;
    g->name=name?strdup(name):strdup("unnamed");
    g->dim = 0; g->algebra_dim = 0; g->mat_size = 0;
    switch(gtype){
        case LIE_GROUP_SO2:g->dim=1;g->algebra_dim=1;g->mat_size=2;break;
        case LIE_GROUP_SO3:g->dim=3;g->algebra_dim=3;g->mat_size=3;break;
        case LIE_GROUP_SE2:g->dim=3;g->algebra_dim=3;g->mat_size=3;break;
        case LIE_GROUP_SE3:g->dim=6;g->algebra_dim=6;g->mat_size=4;break;
        case LIE_GROUP_SU2:g->dim=3;g->algebra_dim=3;g->mat_size=2;break;
        case LIE_GROUP_GL3:g->dim=9;g->algebra_dim=9;g->mat_size=3;break;
        case LIE_GROUP_SL3:g->dim=8;g->algebra_dim=8;g->mat_size=3;break;
        default:break;
    }
    g->g=lie_matrix_identity(g->mat_size);
    g->params=lie_vector_create(g->algebra_dim);
    return g;
}
void lie_group_free(LieGroupElement *g){if(g){free(g->name);lie_matrix_free(g->g);lie_vector_free(g->params);free(g);}}
void lie_group_set_identity(LieGroupElement *g){
    assert(g);
    for(int i=0;i<g->mat_size;i++)for(int j=0;j<g->mat_size;j++)g->g->data[i*g->mat_size+j]=(i==j)?1.0:0.0;
    for(int i=0;i<g->algebra_dim;i++)g->params->data[i]=0.0;
}
LieGroupElement* lie_group_mul(const LieGroupElement *g1, const LieGroupElement *g2){
    assert(g1&&g2&&g1->mat_size==g2->mat_size);
    LieGroupElement *g=lie_group_create(g1->group_type,"product");
    lie_matrix_free(g->g);
    g->g=lie_matrix_mul(g1->g,g2->g);
    return g;
}
LieGroupElement* lie_group_inverse(const LieGroupElement *g){
    assert(g);
    LieGroupElement *inv=lie_group_create(g->group_type,"inverse");
    lie_matrix_free(inv->g);
    if(g->group_type==LIE_GROUP_SO2||g->group_type==LIE_GROUP_SO3)
        inv->g=lie_matrix_transpose(g->g);
    else
        inv->g=lie_matrix_inverse(g->g);
    return inv;
}
void lie_group_invert_inplace(LieGroupElement *g){
    LieMatrix *inv=lie_group_inverse(g)->g;
    lie_matrix_free(g->g);g->g=inv;
}
LieVector* lie_group_act_vector(const LieGroupElement *g, const LieVector *x){
    assert(g&&x&&x->size==g->g->cols);
    LieVector *y=lie_vector_create(g->g->rows);
    for(int i=0;i<g->g->rows;i++){
        double s=0.0;
        for(int j=0;j<g->g->cols;j++)s+=g->g->data[i*g->g->cols+j]*x->data[j];
        y->data[i]=s;
    }
    return y;
}

/* ---- Lie Algebra Element ---- */
LieAlgebraElement* lie_algebra_create(LieAlgebraType atype, const char *name) {
    LieAlgebraElement *xi=(LieAlgebraElement*)malloc(sizeof(LieAlgebraElement));
    assert(xi);
    xi->algebra_type=atype;
    xi->name=name?strdup(name):strdup("unnamed");
    xi->dim = 0; xi->mat = NULL;
    switch(atype){
        case LIE_ALG_SO2:xi->dim=1;xi->mat=lie_matrix_zeros(2,2);break;
        case LIE_ALG_SO3:xi->dim=3;xi->mat=lie_matrix_zeros(3,3);break;
        case LIE_ALG_SE2:xi->dim=3;xi->mat=lie_matrix_zeros(3,3);break;
        case LIE_ALG_SE3:xi->dim=6;xi->mat=lie_matrix_zeros(4,4);break;
        case LIE_ALG_SU2:xi->dim=3;xi->mat=lie_matrix_zeros(2,2);break;
        case LIE_ALG_GL3:xi->dim=9;xi->mat=lie_matrix_zeros(3,3);break;
        case LIE_ALG_SL3:xi->dim=8;xi->mat=lie_matrix_zeros(3,3);break;
        default:break;
    }
    xi->coords=lie_vector_create(xi->dim);
    return xi;
}
void lie_algebra_free(LieAlgebraElement *xi){if(xi){free(xi->name);lie_matrix_free(xi->mat);lie_vector_free(xi->coords);free(xi);}}
LieAlgebraElement* lie_algebra_copy(const LieAlgebraElement *src){
    assert(src);
    LieAlgebraElement *d=lie_algebra_create(src->algebra_type,src->name);
    memcpy(d->coords->data,src->coords->data,(size_t)src->dim*sizeof(double));
    if(src->mat&&d->mat&&src->mat->rows==d->mat->rows){
        int sz=src->mat->rows*src->mat->cols;
        memcpy(d->mat->data,src->mat->data,(size_t)sz*sizeof(double));
    }
    return d;
}
void lie_algebra_set_zero(LieAlgebraElement *xi){
    assert(xi);
    if(xi->mat)lie_matrix_fill(xi->mat,0.0);
    for(int i=0;i<xi->dim;i++)xi->coords->data[i]=0.0;
}

/* ---- Lie Bracket [xi, eta] ---- */
LieAlgebraElement* lie_bracket(const LieAlgebraElement *xi, const LieAlgebraElement *eta) {
    assert(xi&&eta&&xi->algebra_type==eta->algebra_type);
    LieAlgebraElement *r=lie_algebra_create(xi->algebra_type,"bracket");
    if(xi->algebra_type==LIE_ALG_SO3){
        LieVector *cross=lie_vector_cross3(xi->coords,eta->coords);
        for(int i=0;i<3;i++)r->coords->data[i]=cross->data[i];
        lie_vector_free(cross);
        LieMatrix *hat=lie_so3_hat(r->coords);
        memcpy(r->mat->data,hat->data,9*sizeof(double));
        lie_matrix_free(hat);
    }else if(xi->algebra_type==LIE_ALG_SE3){
        LieVector w1={.data=xi->coords->data+3,.size=3},w2={.data=eta->coords->data+3,.size=3};
        LieVector v1={.data=xi->coords->data,.size=3},v2={.data=eta->coords->data,.size=3};
        LieVector *w1xv2=lie_vector_cross3(&w1,&v2),*w2xv1=lie_vector_cross3(&w2,&v1);
        LieVector *w1xw2=lie_vector_cross3(&w1,&w2);
        for(int i=0;i<3;i++){
            r->coords->data[i]=w1xv2->data[i]-w2xv1->data[i];
            r->coords->data[i+3]=w1xw2->data[i];
        }
        lie_vector_free(w1xv2);lie_vector_free(w2xv1);lie_vector_free(w1xw2);
        LieMatrix *hat=lie_se3_hat(r->coords);
        memcpy(r->mat->data,hat->data,16*sizeof(double));
        lie_matrix_free(hat);
    }else{
        if(xi->mat&&eta->mat&&xi->mat->rows==xi->mat->cols){
            LieMatrix *AB=lie_matrix_mul(xi->mat,eta->mat);
            LieMatrix *BA=lie_matrix_mul(eta->mat,xi->mat);
            lie_matrix_free(r->mat);
            r->mat=lie_matrix_sub(AB,BA);
            lie_matrix_free(AB);lie_matrix_free(BA);
        }
    }
    return r;
}
bool lie_jacobi_identity_check(const LieAlgebraElement *xi, const LieAlgebraElement *eta,
    const LieAlgebraElement *zeta, double tol){
    LieAlgebraElement *b1=lie_bracket(eta,zeta),*b2=lie_bracket(zeta,xi),*b3=lie_bracket(xi,eta);
    LieAlgebraElement *j1=lie_bracket(xi,b1),*j2=lie_bracket(eta,b2),*j3=lie_bracket(zeta,b3);
    double n2=0.0;
    for(int i=0;i<xi->dim;i++){
        double s=j1->coords->data[i]+j2->coords->data[i]+j3->coords->data[i];
        n2+=s*s;
    }
    lie_algebra_free(b1);lie_algebra_free(b2);lie_algebra_free(b3);
    lie_algebra_free(j1);lie_algebra_free(j2);lie_algebra_free(j3);
    return sqrt(n2)<tol;
}

/* ---- Hat and Vee Maps ---- */
LieMatrix* lie_so3_hat(const LieVector *omega){
    assert(omega&&omega->size==3);
    LieMatrix *hat=lie_matrix_zeros(3,3);
    hat->data[1]=-omega->data[2];hat->data[2]=omega->data[1];
    hat->data[3]=omega->data[2];hat->data[5]=-omega->data[0];
    hat->data[6]=-omega->data[1];hat->data[7]=omega->data[0];
    return hat;
}
LieVector* lie_so3_vee(const LieMatrix *omega_hat){
    assert(omega_hat&&omega_hat->rows==3&&omega_hat->cols==3);
    LieVector *omega=lie_vector_create(3);
    omega->data[0]=omega_hat->data[7];omega->data[1]=omega_hat->data[2];omega->data[2]=omega_hat->data[3];
    return omega;
}
LieMatrix* lie_se3_hat(const LieVector *xi){
    assert(xi&&xi->size==6);
    LieMatrix *hat=lie_matrix_zeros(4,4);
    hat->data[1]=-xi->data[5];hat->data[2]=xi->data[4];hat->data[3]=xi->data[0];
    hat->data[4]=xi->data[5];hat->data[6]=-xi->data[3];hat->data[7]=xi->data[1];
    hat->data[8]=-xi->data[4];hat->data[9]=xi->data[3];hat->data[11]=xi->data[2];
    return hat;
}
LieVector* lie_se3_vee(const LieMatrix *xi_hat){
    assert(xi_hat&&xi_hat->rows==4&&xi_hat->cols==4);
    LieVector *xi=lie_vector_create(6);
    xi->data[0]=xi_hat->data[3];xi->data[1]=xi_hat->data[7];xi->data[2]=xi_hat->data[11];
    xi->data[3]=xi_hat->data[9];xi->data[4]=xi_hat->data[2];xi->data[5]=xi_hat->data[4];
    return xi;
}
LieMatrix* lie_su2_from_vector(const LieVector *v){
    assert(v&&v->size==3);
    LieMatrix *m=lie_matrix_zeros(2,2);
    m->data[1]=v->data[1];m->data[2]=-v->data[1];
    return m;
}
LieVector* lie_su2_to_vector(const LieMatrix *u){
    assert(u&&u->rows==2&&u->cols==2);
    LieVector *v=lie_vector_create(3);
    v->data[1]=u->data[1];
    return v;
}

/* ---- Exponential Map ---- */
LieMatrix* lie_matrix_exponential(const LieMatrix *A){
    assert(A&&A->rows==A->cols);
    int n=A->rows;
    double nA=lie_matrix_norm_infinity(A);
    int s=(nA>1.0)?(int)ceil(log2(nA)):0;
    if (s < 0) s = 0;
    if (s > 20) s = 20;
    double scale=1.0/(double)(1<<s);
    LieMatrix *B=lie_matrix_scale(A,scale);
    LieMatrix *expB=lie_matrix_identity(n);
    LieMatrix *term=lie_matrix_identity(n);
    double fact=1.0;
    for(int k=1;k<=6;k++){
        LieMatrix *nt=lie_matrix_mul(term,B);
        lie_matrix_free(term);term=nt;
        fact*=(double)k;
        for(int i=0;i<n*n;i++)expB->data[i]+=term->data[i]/fact;
    }
    lie_matrix_free(term);
    LieMatrix *result=lie_matrix_copy(expB);
    for(int k=0;k<s;k++){
        LieMatrix *sq=lie_matrix_mul(result,result);
        lie_matrix_free(result);result=sq;
    }
    lie_matrix_free(B);lie_matrix_free(expB);
    return result;
}
LieGroupElement* lie_exp_so3(const LieVector *omega){
    assert(omega&&omega->size==3);
    double theta=lie_vector_norm(omega);
    LieGroupElement *R=lie_group_create(LIE_GROUP_SO3,"exp_so3");
    if(theta<LIE_EPS){lie_group_set_identity(R);return R;}
    LieMatrix *wh=lie_so3_hat(omega);
    LieMatrix *wh2=lie_matrix_mul(wh,wh);
    double c=cos(theta),s=sin(theta);
    double a=s/theta,b2=(1.0-c)/(theta*theta);
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)
        R->g->data[i*3+j]=(i==j?1.0:0.0)+a*wh->data[i*3+j]+b2*wh2->data[i*3+j];
    memcpy(R->params->data,omega->data,3*sizeof(double));
    lie_matrix_free(wh);lie_matrix_free(wh2);
    return R;
}
LieGroupElement* lie_exp_se3(const LieVector *xi){
    assert(xi&&xi->size==6);
    double theta=sqrt(xi->data[3]*xi->data[3]+xi->data[4]*xi->data[4]+xi->data[5]*xi->data[5]);
    LieGroupElement *g=lie_group_create(LIE_GROUP_SE3,"exp_se3");
    if(theta<LIE_EPS){
        for(int i=0;i<4;i++)for(int j=0;j<4;j++){
            double val=(i==j)?1.0:0.0;
            if(i<3&&j==3)val=xi->data[i];
            g->g->data[i*4+j]=val;
        }
        return g;
    }
    LieMatrix *wh=lie_so3_hat(&(LieVector){.data=xi->data+3,.size=3});
    LieMatrix *wh2=lie_matrix_mul(wh,wh);
    double ct=cos(theta),st=sin(theta),t2=theta*theta,t3=t2*theta;
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)
        g->g->data[i*4+j]=(i==j?1.0:0.0)+(st/theta)*wh->data[i*3+j]+((1.0-ct)/t2)*wh2->data[i*3+j];
    for(int i=0;i<3;i++){
        double tr=xi->data[i];
        for(int j=0;j<3;j++){
            tr+=((1.0-ct)/t2)*wh->data[i*3+j]*xi->data[j];
            tr+=((theta-st)/t3)*wh2->data[i*3+j]*xi->data[j];
        }
        g->g->data[i*4+3]=tr;
    }
    g->g->data[15]=1.0;
    memcpy(g->params->data,xi->data,6*sizeof(double));
    lie_matrix_free(wh);lie_matrix_free(wh2);
    return g;
}
LieGroupElement* lie_group_exp(const LieAlgebraElement *xi){
    assert(xi);
    if(xi->algebra_type==LIE_ALG_SO3)return lie_exp_so3(xi->coords);
    if(xi->algebra_type==LIE_ALG_SE3)return lie_exp_se3(xi->coords);
    if(xi->mat){
        LieMatrix *em=lie_matrix_exponential(xi->mat);
        LieGroupElement *g=lie_group_create((LieGroupType)(xi->algebra_type),xi->name);
        lie_matrix_free(g->g);g->g=em;
        return g;
    }
    return NULL;
}

/* ---- Logarithm Map ---- */
LieVector* lie_log_so3(const LieGroupElement *R){
    assert(R&&R->group_type==LIE_GROUP_SO3);
    double tr=lie_matrix_trace(R->g);
    double cth=(tr-1.0)/2.0;
    if (cth > 1.0) cth = 1.0;
    if (cth < -1.0) cth = -1.0;
    double theta=acos(cth);
    LieVector *omega=lie_vector_create(3);
    if(theta<LIE_EPS)return omega;
    double factor=theta/(2.0*sin(theta));
    omega->data[0]=factor*(R->g->data[7]-R->g->data[5]);
    omega->data[1]=factor*(R->g->data[2]-R->g->data[6]);
    omega->data[2]=factor*(R->g->data[3]-R->g->data[1]);
    return omega;
}
LieVector* lie_log_se3(const LieGroupElement *g){
    assert(g&&g->group_type==LIE_GROUP_SE3);
    LieMatrix *Rm=lie_matrix_create(3,3);
    LieVector *p=lie_vector_create(3);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++)Rm->data[i*3+j]=g->g->data[i*4+j];
        p->data[i]=g->g->data[i*4+3];
    }
    LieGroupElement Rt={.group_type=LIE_GROUP_SO3,.g=Rm,.mat_size=3};
    LieVector *omega=lie_log_so3(&Rt);
    double theta=lie_vector_norm(omega);
    LieVector *xi=lie_vector_create(6);
    if(theta<LIE_EPS){
        memcpy(xi->data,p->data,3*sizeof(double));
    }else{
        LieMatrix *wh=lie_so3_hat(omega),*wh2=lie_matrix_mul(wh,wh);
        for(int i=0;i<3;i++){
            double val=p->data[i];
            for(int j=0;j<3;j++)val-=0.5*wh->data[i*3+j]*p->data[j];
            double alpha=1.0/(theta*theta)-(1.0+cos(theta))/(2.0*theta*sin(theta));
            for(int j=0;j<3;j++)val+=alpha*wh2->data[i*3+j]*p->data[j];
            xi->data[i]=val;xi->data[i+3]=omega->data[i];
        }
        lie_matrix_free(wh);lie_matrix_free(wh2);
    }
    lie_matrix_free(Rm);lie_vector_free(p);lie_vector_free(omega);
    return xi;
}
LieMatrix* lie_matrix_logarithm(const LieMatrix *G){
    assert(G&&G->rows==G->cols);
    int n=G->rows;
    LieMatrix *Gk=lie_matrix_copy(G);
    LieMatrix *I=lie_matrix_identity(n);
    int s=0;
    while(s<15){
        LieMatrix *diff=lie_matrix_sub(Gk,I);
        double nd=lie_matrix_norm_frobenius(diff);lie_matrix_free(diff);
        if(nd<=0.5)break;
        LieMatrix *Y=lie_matrix_copy(Gk);
        for(int iter=0;iter<8;iter++){
            LieMatrix *Yinv=lie_matrix_inverse(Y),*sum=lie_matrix_add(Y,Yinv);
            LieMatrix *nY=lie_matrix_scale(sum,0.5);
            lie_matrix_free(Y);Y=nY;
            lie_matrix_free(Yinv);lie_matrix_free(sum);
        }
        lie_matrix_free(Gk);Gk=Y;
        s++;
    }
    LieMatrix *X=lie_matrix_sub(Gk,I);
    LieMatrix *logG=lie_matrix_zeros(n,n);
    LieMatrix *term=lie_matrix_copy(X);
    double sign=1.0;
    for(int k=1;k<=15;k++){
        for(int i=0;i<n*n;i++)logG->data[i]+=sign*term->data[i]/(double)k;
        sign=-sign;
        LieMatrix *nt=lie_matrix_mul(term,X);
        lie_matrix_free(term);term=nt;
    }
    double factor=(double)(1<<s);
    LieMatrix *result=lie_matrix_scale(logG,factor);
    lie_matrix_free(Gk);lie_matrix_free(I);lie_matrix_free(X);lie_matrix_free(term);lie_matrix_free(logG);
    return result;
}

/* ---- Adjoint and Coadjoint Representations ---- */
LieAlgebraElement* lie_Ad(const LieGroupElement *g, const LieAlgebraElement *eta){
    assert(g&&eta);
    LieMatrix *g_eta=lie_matrix_mul(g->g,eta->mat);
    LieGroupElement *gi=lie_group_inverse(g);
    LieMatrix *Ad_eta=lie_matrix_mul(g_eta,gi->g);
    LieAlgebraElement *r=lie_algebra_create(eta->algebra_type,"Ad");
    lie_matrix_free(r->mat);r->mat=Ad_eta;
    lie_matrix_free(g_eta);lie_group_free(gi);
    return r;
}
LieMatrix* lie_Ad_matrix_so3(const LieMatrix *R){
    assert(R&&R->rows==3&&R->cols==3);
    return lie_matrix_copy(R);
}
LieMatrix* lie_Ad_matrix_se3(const LieMatrix *R, const LieVector *p){
    assert(R&&R->rows==3&&R->cols==3);
    LieMatrix *Ad=lie_matrix_zeros(6,6);
    for(int i=0;i<3;i++)for(int j=0;j<3;j++){
        Ad->data[i*6+j]=R->data[i*3+j];
        Ad->data[(i+3)*6+(j+3)]=R->data[i*3+j];
    }
    LieMatrix *ph=lie_so3_hat(p),*phR=lie_matrix_mul(ph,R);
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)Ad->data[(i+3)*6+j]=phR->data[i*3+j];
    lie_matrix_free(ph);lie_matrix_free(phR);
    return Ad;
}
LieAlgebraElement* lie_ad(const LieAlgebraElement *xi, const LieAlgebraElement *eta){
    return lie_bracket(xi,eta);
}
LieMatrix* lie_ad_matrix_so3(const LieVector *omega){
    return lie_so3_hat(omega);
}
LieMatrix* lie_ad_matrix_se3(const LieVector *xi){
    assert(xi&&xi->size==6);
    LieMatrix *ad=lie_matrix_zeros(6,6);
    LieVector w={.data=xi->data+3,.size=3},v={.data=xi->data,.size=3};
    LieMatrix *wh=lie_so3_hat(&w),*vh=lie_so3_hat(&v);
    for(int i=0;i<3;i++)for(int j=0;j<3;j++){
        ad->data[i*6+j]=wh->data[i*3+j];
        ad->data[i*6+j+3]=vh->data[i*3+j];
        ad->data[(i+3)*6+j+3]=wh->data[i*3+j];
    }
    lie_matrix_free(wh);lie_matrix_free(vh);
    return ad;
}
LieVector* lie_Ad_star_so3(const LieMatrix *R, const LieVector *mu){
    assert(R&&mu&&mu->size==3);
    LieVector *r=lie_vector_create(3);
    for(int i=0;i<3;i++){double s=0.0;for(int j=0;j<3;j++)s+=R->data[j*3+i]*mu->data[j];r->data[i]=s;}
    return r;
}
LieVector* lie_ad_star_so3(const LieVector *omega, const LieVector *mu){
    return lie_vector_cross3(mu,omega);
}

/* ---- Killing Form ---- */
double lie_killing_form(const LieAlgebraElement *xi, const LieAlgebraElement *eta){
    assert(xi&&eta);
    if(xi->algebra_type==LIE_ALG_SO3)return -2.0*lie_vector_dot(xi->coords,eta->coords);
    LieMatrix *adxi=lie_ad_matrix_so3(xi->coords),*adeta=lie_ad_matrix_so3(eta->coords);
    LieMatrix *prod=lie_matrix_mul(adxi,adeta);
    double k=lie_matrix_trace(prod);
    lie_matrix_free(adxi);lie_matrix_free(adeta);lie_matrix_free(prod);
    return k;
}

/* ---- BCH Formula ---- */
LieAlgebraElement* lie_bch_2nd_order(const LieAlgebraElement *A, const LieAlgebraElement *B){
    assert(A&&B&&A->algebra_type==B->algebra_type);
    LieAlgebraElement *br=lie_bracket(A,B);
    LieAlgebraElement *r=lie_algebra_create(A->algebra_type,"bch2");
    for(int i=0;i<A->dim;i++)r->coords->data[i]=A->coords->data[i]+B->coords->data[i]+0.5*br->coords->data[i];
    lie_algebra_free(br);
    return r;
}
LieAlgebraElement* lie_bch_3rd_order(const LieAlgebraElement *A, const LieAlgebraElement *B){
    assert(A&&B&&A->algebra_type==B->algebra_type);
    LieAlgebraElement *bAB=lie_bracket(A,B),*bAAB=lie_bracket(A,bAB);
    LieAlgebraElement *bBA=lie_bracket(B,A),*bBBA=lie_bracket(B,bBA);
    LieAlgebraElement *r=lie_algebra_create(A->algebra_type,"bch3");
    for(int i=0;i<A->dim;i++)
        r->coords->data[i]=A->coords->data[i]+B->coords->data[i]+0.5*bAB->coords->data[i]
            +(1.0/12.0)*(bAAB->coords->data[i]+bBBA->coords->data[i]);
    lie_algebra_free(bAB);lie_algebra_free(bAAB);lie_algebra_free(bBA);lie_algebra_free(bBBA);
    return r;
}

/* ---- Display ---- */
void lie_matrix_print(const LieMatrix *m, const char *label){
    if(!m){printf("%s: NULL\n",label?label:"matrix");return;}
    printf("%s (%dx%d):\n",label?label:"matrix",m->rows,m->cols);
    for(int i=0;i<m->rows;i++){
        printf("  ");
        for(int j=0;j<m->cols;j++)printf("%12.6f ",m->data[i*m->cols+j]);
        printf("\n");
    }
}
void lie_vector_print(const LieVector *v, const char *label){
    if(!v){printf("%s: NULL\n",label?label:"vector");return;}
    printf("%s (dim=%d): [",label?label:"vector",v->size);
    for(int i=0;i<v->size;i++)printf("%12.6f%s",v->data[i],(i<v->size-1)?", ":"");
    printf("]\n");
}
void lie_group_print(const LieGroupElement *g){
    if(!g){printf("Group: NULL\n");return;}
    printf("Group: %s (type=%d,dim=%d)\n",g->name,g->group_type,g->dim);
    lie_matrix_print(g->g,"  matrix");
}
void lie_algebra_print(const LieAlgebraElement *xi){
    if(!xi){printf("Algebra: NULL\n");return;}
    printf("Algebra: %s (type=%d,dim=%d)\n",xi->name,xi->algebra_type,xi->dim);
    lie_vector_print(xi->coords,"  coords");
}
