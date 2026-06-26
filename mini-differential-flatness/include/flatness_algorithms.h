#ifndef FLATNESS_ALGORITHMS_H
#define FLATNESS_ALGORITHMS_H
#include "flatness_core.h"
#include "lie_algebra.h"

int flat_controllability_matrix_rank(const FlatControlAffineSystem*sys,int*rank);
bool flat_linear_controllability_check(const FlatControlAffineSystem*sys);
int flat_martin_rouchon_test(const FlatControlAffineSystem*sys,int*dims,int max_steps);
int flat_search_outputs(const FlatControlAffineSystem*sys,FlatPolynomial*found_outputs,int*num_found);
int flat_implicit_elimination(int n_vars,int n_eqns,const FlatPolynomial*jacobian,FlatPolynomial*eliminated);
int flat_pbh_test(const double*A,const double*B,int n,int m);
int flat_kalman_decomposition(double*A,double*B,int n,int m,int*nc);
int flat_detect_singularities(const FlatControlAffineSystem*sys,const FlatPolynomial*flat_outputs,double*singular_pts,int max_pts);
int flat_detect_structure_pattern(const FlatControlAffineSystem*sys,int*pattern);
#endif
