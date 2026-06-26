#ifndef FLATNESS_ANALYSIS_H
#define FLATNESS_ANALYSIS_H
#include "flatness_core.h"

/* L8: Orbital flatness */
int flat_orbital_check(const FlatControlAffineSystem*sys,const double*x_ref,int n_ref);
int flat_resolve_singularity(const FlatControlAffineSystem*sys,const double*singular_pt,double*regularized_output);

/* L8: Implicit systems */
typedef struct{double(*F)(const double*x,const double*xdot,const double*u,int idx);int n_states,n_inputs,n_equations;}FlatImplicitSystem;
int flat_implicit_to_explicit(const FlatImplicitSystem*imp,FlatControlAffineSystem*exp);

/* L8: Numerical flatness */
int flat_numerical_flatness_check(const FlatControlAffineSystem*sys,double tol,int max_iter,double*confidence);

/* L8: Time-varying flatness */
typedef struct{FlatVectorField drift;FlatVectorField control_fields[10];double t;}FlatTimeVaryingSystem;
int flat_time_varying_check(const FlatTimeVaryingSystem*sys,double t,double*result);

/* L8: Symmetry-based flat output discovery */
int flat_symmetry_invariant_outputs(const FlatControlAffineSystem*sys,int*symmetry_type);
int flat_construct_flat_outputs_dynamic(const FlatControlAffineSystem*sys,int max_deriv,FlatPolynomial*outputs,int*num_outputs);

/* L8: Density evolution under flatness */
int flat_frobenius_perron(const FlatControlAffineSystem*sys,const double*density,int n_grid,double*evolved_density);

/* L8: Reachable set approximation */
int flat_reachable_set(const FlatControlAffineSystem*sys,double t_horizon,int n_samples,double*boundary_points,int*num_boundary);

/* L8: Cascade flatness */
int flat_cascade_check(const FlatControlAffineSystem*sys1,const FlatControlAffineSystem*sys2,int*is_flat);

/* L4-L8: Additional polynomial tools */
int flat_poly_substitute(const FlatPolynomial*p,const FlatPolynomial*subs,FlatPolynomial*out);
int flat_poly_integrate(const FlatPolynomial*p,int var,FlatPolynomial*out);
double flat_poly_inner_product(const FlatPolynomial*a,const FlatPolynomial*b,const double*domain);
int flat_poly_reduce(FlatPolynomial*p,FlatPolynomial*basis,int n_basis);
bool flat_poly_is_sos(const FlatPolynomial*p);
int flat_poly_rk4_step(const FlatPolynomial*f,int n,const double*x,double dt,double*x_next);
int flat_find_equilibria(const FlatControlAffineSystem*sys,double*equilibria,int max_eq);
int flat_bifurcation_detect(const FlatControlAffineSystem*sys,double param,double*param_crit);

#endif
