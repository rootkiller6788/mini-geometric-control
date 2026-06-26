#ifndef DYNAMIC_FEEDBACK_H
#define DYNAMIC_FEEDBACK_H
#include "flatness_core.h"

/* L5: Dynamic extension (add integrators to inputs) */
int flat_dynamic_extend(FlatControlAffineSystem*sys,int input_channel);
int flat_dynamic_extend_to_flat(FlatControlAffineSystem*sys,int*num_extensions);

/* L4: Input-output linearization (Isidori) */
int flat_io_linearize(const FlatControlAffineSystem*sys,const int*rel_deg,double*decoupling_mat,double*nonlin_vec,const double*x_current);

/* L4: Full state linearization via dynamic feedback */
int flat_full_state_linearize(const FlatControlAffineSystem*sys,const FlatPolynomial*flat_outputs,double*K,int*n_gains);

/* L5: Flatness-based tracking control */
int flat_tracking_control(const FlatControlAffineSystem*sys,const FlatPolynomial*flat_outputs,const double*x_current,const double*y_desired_derivs,const double*gains,int kappa,double*u_out);

/* L8: Quasi-static feedback (trajectory-dependent) */
typedef struct{double(*alpha)(const double*x,const double*y_ref,int d);double(*beta)(const double*x,const double*y_ref,int d);int deriv_order;}FlatQuasiStaticFeedback;
int flat_design_quasi_static(const FlatControlAffineSystem*sys,const FlatPolynomial*flat_outputs,FlatQuasiStaticFeedback*qsf);

/* L4: Error dynamics and stability */
int flat_error_dynamics_matrix(const double*gains,int kappa,double*A);
int flat_place_poles(const double*poles,int kappa,double*gains);
#endif
