#ifndef SYMPLECTIC_GEOMETRY_H
#define SYMPLECTIC_GEOMETRY_H

#include "hamiltonian_control.h"

/*=============================================================================
 * Symplectic geometry structures for control theory.
 * Reference: Arnold, Mathematical Methods of Classical Mechanics
 *            da Silva, Lectures on Symplectic Geometry
 *===========================================================================*/

/* L2: Symplectic form omega on R^{2n}.
 * In canonical coordinates: omega = sum_i dq_i AND dp_i
 * omega(v, w) = v^T J w where J is the symplectic matrix. */
double symplectic_form(const double *v, const double *w, int n_dof);

/* L2: Symplectic area (sum of projected areas onto (q_i, p_i) planes).
 * For 2D: just the oriented area in phase plane. */
double symplectic_area(const phase_point_t *points, int n_points);

/* L2: Check if a linear map M preserves symplectic structure: M^T J M = J.
 * Returns Frobenius norm of M^T J M - J. */
double check_symplectic_map(double **M, int dim2n);

/* L3: Darboux coordinates -- locally, every symplectic manifold is R^{2n}.
 * This function finds a local linear transformation to canonical form.
 * Given a bilinear form omega (as 2n x 2n matrix), find basis where
 * omega becomes the standard symplectic matrix J.
 * Returns 0 on success. */
int darboux_basis(double **omega, int dim2n, double **basis);

/* L3: Symplectic eigenvalues of a positive-definite matrix under J.
 * Williamson theorem: for SPD A of size 2n, there exists symplectic S
 * such that S^T A S = diag(D, D) where D = diag(d1,...,dn).
 * d_i are the symplectic eigenvalues.
 * Reference: Williamson (1936), de Gosson §2.3 */
int symplectic_eigenvalues(double **A, int n_dof, double *symplectic_evals);

/* L3: Symplectic capacity (Gromov width).
 * For an ellipsoid E = {z: z^T A z <= 1}, the symplectic capacity
 * is c(E) = pi / d_max where d_max is the largest symplectic eigenvalue.
 * Reference: Gromov (1985), Hofer-Zehnder §3.5 */
double symplectic_capacity_ellipsoid(double **A, int n_dof);

/* L4: Darboux theorem -- local existence of canonical coordinates.
 * Given a closed nondegenerate 2-form, construct local canonical (q,p).
 * This version works for linear (constant) symplectic forms. */
int darboux_local_coordinates(double **omega_constant, int dim2n,
                              double **transform_matrix);

/* L5: Symplectic Gram-Schmidt.
 * Given 2k vectors v_1,...,v_{2k}, produce a symplectic basis
 * e_1,...,e_k, f_1,...,f_k such that omega(e_i, f_j) = delta_ij
 * and all other symplectic pairings vanish.
 * Reference: McDuff-Salamon §2.1 */
int symplectic_gram_schmidt(double **vectors, int n_vectors, int n_dof,
                            double **symplectic_basis_out);

/* L5: Linear symplectic reduction (Marsden-Weinstein for linear spaces).
 * Given a constraint surface (kernel of constraint matrix C),
 * compute the reduced symplectic space.
 * Reference: Marsden-Weinstein (1974) */
int linear_symplectic_reduction(const symplectic_matrix_t *J,
                                double **constraint_matrix, int n_constraints,
                                int n_dof, double **reduced_J, int *reduced_dim);

/* L6: Lagrangian submanifolds.
 * A subspace L of R^{2n} is Lagrangian if dim L = n and omega|_L = 0.
 * Check if a set of n vectors spans a Lagrangian subspace. */
int is_lagrangian_subspace(double **basis_vectors, int n, int n_dof);

/* L8: Symplectic capacity with holomorphic curves (introductory).
 * For simple domains, compute an upper bound on symplectic capacity
 * via the displacement energy. */
double displacement_energy_bound(double r, int n_dof);

#endif /* SYMPLECTIC_GEOMETRY_H */
