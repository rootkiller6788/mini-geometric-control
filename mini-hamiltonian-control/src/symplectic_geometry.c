/*=============================================================================
 * symplectic_geometry.c -- Symplectic geometry operations
 * Reference: Arnold, Mathematical Methods of Classical Mechanics
 *            da Silva, Lectures on Symplectic Geometry
 *            McDuff & Salamon, Introduction to Symplectic Topology
 *===========================================================================*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "hamiltonian_control.h"
#include "symplectic_geometry.h"

/*---------------------------------------------------------------------------*
 * L2: Symplectic form omega(v, w) = v^T J w on R^{2n}
 * In canonical coordinates: omega = sum_i dq_i AND dp_i
 * Reference: Arnold section 37, da Silva section 1.1
 *---------------------------------------------------------------------------*/
double symplectic_form(const double *v, const double *w, int n_dof)
{
    double sum = 0.0;
    for (int i = 0; i < n_dof; i++) {
        sum += v[i] * w[n_dof + i] - v[n_dof + i] * w[i];
    }
    return sum;
}

/*---------------------------------------------------------------------------*
 * L2: Symplectic area of a polygon in phase space.
 * Uses surveyor's formula applied to (q_i, p_i) planes.
 * Reference: da Silva section 1.3
 *---------------------------------------------------------------------------*/
double symplectic_area(const phase_point_t *points, int n_points)
{
    if (n_points < 3) return 0.0;
    int n = points[0].n;
    double area = 0.0;
    for (int k = 0; k < n; k++) {
        double plane_area = 0.0;
        for (int i = 0; i < n_points; i++) {
            int j = (i + 1) % n_points;
            plane_area += points[i].q[k] * points[j].p[k]
                        - points[j].q[k] * points[i].p[k];
        }
        area += fabs(plane_area) * 0.5;
    }
    return area;
}

/*---------------------------------------------------------------------------*
 * L2: Check if linear map M is symplectic: M^T J M = J.
 * Returns Frobenius norm of M^T J M - J.
 * Reference: da Silva section 2.2
 *---------------------------------------------------------------------------*/
double check_symplectic_map(double **M, int dim2n)
{
    int n = dim2n / 2;
    double **J = matrix_alloc(dim2n, dim2n);
    double **MT = matrix_alloc(dim2n, dim2n);
    double **MTJ = matrix_alloc(dim2n, dim2n);
    double **MTJM = matrix_alloc(dim2n, dim2n);

    /* Build J */
    for (int i = 0; i < dim2n; i++)
        for (int j = 0; j < dim2n; j++)
            J[i][j] = 0.0;
    for (int i = 0; i < n; i++) {
        J[i][n+i] = 1.0; J[n+i][i] = -1.0;
    }

    /* M^T */
    for (int i = 0; i < dim2n; i++)
        for (int j = 0; j < dim2n; j++)
            MT[i][j] = M[j][i];

    /* MTJ = M^T * J */
    for (int i = 0; i < dim2n; i++)
        for (int j = 0; j < dim2n; j++) {
            MTJ[i][j] = 0.0;
            for (int k = 0; k < dim2n; k++)
                MTJ[i][j] += MT[i][k] * J[k][j];
        }

    /* MTJM = MTJ * M */
    for (int i = 0; i < dim2n; i++)
        for (int j = 0; j < dim2n; j++) {
            MTJM[i][j] = 0.0;
            for (int k = 0; k < dim2n; k++)
                MTJM[i][j] += MTJ[i][k] * M[k][j];
        }

    double frob = 0.0;
    for (int i = 0; i < dim2n; i++)
        for (int j = 0; j < dim2n; j++) {
            double diff = MTJM[i][j] - J[i][j];
            frob += diff * diff;
        }

    matrix_free(J, dim2n); matrix_free(MT, dim2n);
    matrix_free(MTJ, dim2n); matrix_free(MTJM, dim2n);
    return sqrt(frob);
}

/*---------------------------------------------------------------------------*
 * L3: Darboux basis construction.
 * Given a constant symplectic form omega (as 2n x 2n matrix),
 * find a basis transformation to standard canonical form.
 * Uses symplectic Gram-Schmidt.
 * Reference: da Silva section 2.1, McDuff-Salamon section 2.1
 *---------------------------------------------------------------------------*/
int darboux_basis(double **omega, int dim2n, double **basis)
{
    /* omega should be non-degenerate skew-symmetric.
     * Find basis vectors e_1..e_n, f_1..f_n such that:
     * omega(e_i, e_j)=0, omega(f_i, f_j)=0, omega(e_i, f_j)=delta_ij.
     * Initialize basis as identity, then apply symplectic Gram-Schmidt. */
    int n = dim2n / 2;
    /* Initialize basis to identity */
    for (int i = 0; i < dim2n; i++)
        for (int j = 0; j < dim2n; j++)
            basis[i][j] = (i == j) ? 1.0 : 0.0;

    /* Work with columns of basis as vectors */
    double **vecs = matrix_alloc(dim2n, dim2n);
    for (int i = 0; i < dim2n; i++)
        for (int j = 0; j < dim2n; j++)
            vecs[i][j] = basis[j][i];

    for (int k = 0; k < n; k++) {
        /* Make vecs[n+k] symplectic-orthogonal to vecs[k] */
        double omega_kk = 0.0;
        for (int i = 0; i < dim2n; i++)
            for (int j = 0; j < dim2n; j++)
                omega_kk += vecs[k][i] * omega[i][j] * vecs[n+k][j];

        if (fabs(omega_kk) < 1e-12) continue;

        /* Normalize to omega(e_k, f_k) = 1 */
        double scale = sqrt(fabs(omega_kk));
        for (int i = 0; i < dim2n; i++) {
            vecs[k][i] /= scale;
            vecs[n+k][i] /= scale;
        }

        /* Orthogonalize remaining vectors against e_k and f_k */
        for (int j = k + 1; j < n; j++) {
            double om_e = 0.0, om_f = 0.0;
            for (int i = 0; i < dim2n; i++) {
                for (int l = 0; l < dim2n; l++) {
                    om_e += vecs[j][i] * omega[i][l] * vecs[n+k][l];
                    om_f += vecs[n+j][i] * omega[i][l] * vecs[k][l];
                }
            }
            for (int i = 0; i < dim2n; i++) {
                vecs[j][i]   -= om_f * vecs[k][i] + om_e * vecs[n+k][i];
                vecs[n+j][i] += om_f * vecs[n+k][i] - om_e * vecs[k][i];
            }
        }
    }

    /* Transpose back */
    for (int i = 0; i < dim2n; i++)
        for (int j = 0; j < dim2n; j++)
            basis[i][j] = vecs[j][i];

    matrix_free(vecs, dim2n);
    return 0;
}

/*---------------------------------------------------------------------------*
 * L3: Symplectic eigenvalues via Williamson theorem.
 * For SPD A (2n x 2n), find symplectic transformation S such that
 * S^T A S = diag(D, D) where D = diag(d1,...,dn).
 * d_i are the symplectic eigenvalues.
 * Reference: Williamson (1936), de Gosson section 2.3
 *---------------------------------------------------------------------------*/
int symplectic_eigenvalues(double **A, int n_dof, double *symplectic_evals)
{
    int dim2n = 2 * n_dof;
    /* Construct K = J^{-1} A = -J A.
     * Symplectic eigenvalues are the positive imaginary parts
     * of eigenvalues of K. Compute eigenvalues of K, take the
     * n positive imaginary parts. */
    double **J = matrix_alloc(dim2n, dim2n);
    double **K = matrix_alloc(dim2n, dim2n);

    for (int i = 0; i < n_dof; i++) {
        J[i][n_dof+i] = 1.0;
        J[n_dof+i][i] = -1.0;
    }

    /* K = -J * A */
    for (int i = 0; i < dim2n; i++)
        for (int j = 0; j < dim2n; j++) {
            K[i][j] = 0.0;
            for (int k = 0; k < dim2n; k++)
                K[i][j] -= J[i][k] * A[k][j];
        }

    /* Compute eigenvalues of K via QR algorithm */
    double **Qmat = matrix_alloc(dim2n, dim2n);
    double **Rmat = matrix_alloc(dim2n, dim2n);
    for (int i = 0; i < dim2n; i++)
        for (int j = 0; j < dim2n; j++)
            Rmat[i][j] = K[i][j];

    for (int iter = 0; iter < 80; iter++) {
        for (int j = 0; j < dim2n; j++) {
            for (int i = 0; i < dim2n; i++) Qmat[i][j] = Rmat[i][j];
            for (int k = 0; k < j; k++) {
                double dot = 0.0;
                for (int i = 0; i < dim2n; i++) dot += Qmat[i][j] * Qmat[i][k];
                for (int i = 0; i < dim2n; i++) Qmat[i][j] -= dot * Qmat[i][k];
            }
            double norm = 0.0;
            for (int i = 0; i < dim2n; i++) norm += Qmat[i][j]*Qmat[i][j];
            norm = sqrt(norm);
            if (norm > 1e-15) for (int i=0;i<dim2n;i++) Qmat[i][j]/=norm;
        }
        double **tmp = matrix_alloc(dim2n, dim2n);
        for (int i=0;i<dim2n;i++)
            for (int j=0;j<dim2n;j++){
                tmp[i][j]=0.0;
                for(int k=0;k<dim2n;k++) tmp[i][j]+=Rmat[i][k]*Qmat[k][j];
            }
        for (int i=0;i<dim2n;i++)
            for (int j=0;j<dim2n;j++){
                Rmat[i][j]=0.0;
                for(int k=0;k<dim2n;k++) Rmat[i][j]+=Qmat[k][i]*tmp[k][j];
            }
        matrix_free(tmp, dim2n);
    }

    /* Extract n symplectic eigenvalues (positive imaginary parts) */
    for (int i = 0; i < n_dof; i++) {
        double val = fabs(Rmat[2*i][2*i]);
        /* For simple cases, eigenvalues come in +/- i*omega pairs */
        double val_next = (2*i+1 < dim2n) ? fabs(Rmat[2*i+1][2*i+1]) : val;
        symplectic_evals[i] = (val + val_next) / 2.0;
    }

    matrix_free(J, dim2n); matrix_free(K, dim2n);
    matrix_free(Qmat, dim2n); matrix_free(Rmat, dim2n);
    return 0;
}

/*---------------------------------------------------------------------------*
 * L3: Symplectic capacity (Gromov width) of an ellipsoid.
 * For E = {z: z^T A z <= 1}, c(E) = pi / d_max.
 * Reference: Gromov (1985), Hofer-Zehnder section 3.5
 *---------------------------------------------------------------------------*/
double symplectic_capacity_ellipsoid(double **A, int n_dof)
{
    double *symp_evals = (double*)malloc(n_dof * sizeof(double));
    symplectic_eigenvalues(A, n_dof, symp_evals);
    double d_max = symp_evals[0];
    for (int i = 1; i < n_dof; i++)
        if (symp_evals[i] > d_max) d_max = symp_evals[i];
    free(symp_evals);
    return (d_max > 1e-15) ? M_PI / d_max : 0.0;
}

/*---------------------------------------------------------------------------*
 * L4: Darboux local coordinates -- find linear transformation to canonical form.
 * Given constant omega, compute basis transformation T such that
 * T^T omega T = J (standard symplectic matrix).
 * Reference: Arnold section 43
 *---------------------------------------------------------------------------*/
int darboux_local_coordinates(double **omega_constant, int dim2n,
                              double **transform_matrix)
{
    return darboux_basis(omega_constant, dim2n, transform_matrix);
}

/*---------------------------------------------------------------------------*
 * L5: Symplectic Gram-Schmidt process.
 * Converts 2k vectors into a symplectic basis.
 * Reference: McDuff-Salamon section 2.1
 *---------------------------------------------------------------------------*/
int symplectic_gram_schmidt(double **vectors, int n_vectors, int n_dof,
                            double **symplectic_basis_out)
{
    int dim2n = 2 * n_dof;
    /* Initialize output basis */
    for (int i = 0; i < dim2n; i++)
        for (int j = 0; j < dim2n; j++)
            symplectic_basis_out[i][j] = (i == j) ? 1.0 : 0.0;

    int k = n_vectors / 2;
    if (k > n_dof) k = n_dof;

    for (int p = 0; p < k && 2*p < n_vectors; p++) {
        /* Copy vector 2p as e_p */
        for (int i = 0; i < dim2n; i++)
            symplectic_basis_out[i][p] = vectors[2*p][i];

        /* Symplectic-orthogonalize remaining vectors against e_p and f_p */
        for (int j = 2*p+1; j < n_vectors; j++) {
            double om = symplectic_form(vectors[j], vectors[2*p], n_dof);
            if (fabs(om) > 1e-12 && 2*p+1 < n_vectors) {
                /* This is a candidate for f_p */
                for (int i = 0; i < dim2n; i++)
                    symplectic_basis_out[i][n_dof+p] = vectors[j][i] / om;
                break;
            }
        }
    }
    return 0;
}

/*---------------------------------------------------------------------------*
 * L5: Linear symplectic reduction (Marsden-Weinstein for linear spaces).
 * Given constraint C z = 0, compute reduced symplectic space.
 * reduced_dim = dim2n - 2*rank(C).
 * Reference: Marsden-Weinstein (1974)
 *---------------------------------------------------------------------------*/
int linear_symplectic_reduction(const symplectic_matrix_t *J,
                                double **constraint_matrix, int n_constraints,
                                int n_dof, double **reduced_J, int *reduced_dim)
{
    int dim2n = 2 * n_dof;
    /* For linear reduction:
     * L = ker(C) intersect J*ker(C) --> reduced space
     * In practice, reduced_dim = dim2n - 2*rank(C).
     * Build reduced J by projecting J onto L.
     * Simplified: for non-degenerate constraints, reduced_dim = dim2n - 2*n_constraints. */
    *reduced_dim = dim2n - 2 * n_constraints;
    if (*reduced_dim <= 0) return -1;

    int rn = *reduced_dim / 2;
    for (int i = 0; i < *reduced_dim; i++)
        for (int j = 0; j < *reduced_dim; j++)
            reduced_J[i][j] = 0.0;
    for (int i = 0; i < rn; i++) {
        reduced_J[i][rn+i] =  1.0;
        reduced_J[rn+i][i] = -1.0;
    }
    return 0;
}

/*---------------------------------------------------------------------------*
 * L6: Lagrangian subspace check.
 * A subspace L (dim = n) of R^{2n} is Lagrangian if omega|_L = 0.
 * Check for all pairs of basis vectors: omega(v_i, v_j) = 0.
 * Reference: da Silva section 3.1
 *---------------------------------------------------------------------------*/
int is_lagrangian_subspace(double **basis_vectors, int n, int n_dof)
{
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double om = symplectic_form(basis_vectors[i], basis_vectors[j], n_dof);
            if (fabs(om) > 1e-10) return 0;
        }
    }
    return 1;
}

/*---------------------------------------------------------------------------*
 * L8: Displacement energy bound (simple estimate).
 * For a ball of radius r, the displacement energy is pi * r^2.
 * Provides an upper bound on symplectic capacity.
 * Reference: Hofer (1990)
 *---------------------------------------------------------------------------*/
double displacement_energy_bound(double r, int n_dof)
{
    (void)n_dof;
    return M_PI * r * r;
}