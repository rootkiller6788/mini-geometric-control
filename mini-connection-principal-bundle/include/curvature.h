#ifndef CURVATURE_H
#define CURVATURE_H

#include <stdbool.h>
#include "lie_group.h"
#include "principal_bundle.h"
#include "connection.h"

/* ============================================================================
 * Curvature Theory on Principal Bundles
 *
 * The curvature 2-form Ω of a connection ω is defined by the Cartan
 * structure equation:
 *   Ω = dω + (1/2)[ω, ω]
 *
 * Equivalently, Ω(X, Y) = dω(X, Y) + [ω(X), ω(Y)].
 *
 * In terms of the gauge potential A:
 *   F_{μν} = ∂_μ A_ν − ∂_ν A_μ + [A_μ, A_ν]
 *
 * The Bianchi identity: DΩ = 0, i.e. dΩ + [ω, Ω] = 0.
 * In components: D_μ F_{νρ} + D_ν F_{ρμ} + D_ρ F_{μν} = 0.
 *
 * On the lattice, curvature is measured by the plaquette (Wilson loop
 * around an elementary square):
 *   U_{μν}(x) = U_μ(x) U_ν(x+μ̂) U_μ^{-1}(x+ν̂) U_ν^{-1}(x)
 *
 * The continuum field strength is recovered as a → 0:
 *   U_{μν}(x) = exp(i a^2 F_{μν}) + O(a^3)  for U(1)
 *   U_{μν}(x) = exp(i a^2 F_{μν}) + O(a^3)  for SU(N)
 *
 * References:
 *   Cartan (1928), Leçons sur la géométrie des espaces de Riemann
 *   Kobayashi & Nomizu (1963) Vol I, Ch.II §5
 *   Wilson (1974), Phys. Rev. D10, 2445
 *   Nakahara (2003), §10.3
 * ============================================================================ */

/* --- Plaquette on the lattice --- */
typedef struct {
    int site;            /* base site x */
    int mu;              /* first direction */
    int nu;              /* second direction */
    LieGroupElement U;   /* plaquette value U_{μν}(x) */
} Plaquette;

/* ============================================================================
 * Plaquette / Curvature Computation (L3-L4)
 * =========================================================================== */

/**
 * Compute the plaquette U_{μν}(x) for the elementary square in the μ-ν plane.
 * U_{μν}(x) = U_μ(x) · U_ν(x+μ̂) · U_μ^{-1}(x+ν̂) · U_ν^{-1}(x)
 *
 * This is the lattice discretization of the field strength tensor F_{μν}.
 * For U(1): U_{μν} = exp(i a^2 F_{μν}), so F_{μν} ≈ -i log(U_{μν}) / a^2.
 *
 * Theorem (Wilson, 1974): The Wilson action S = β Σ(1 - Re Tr U_{μν})
 * converges to the Yang-Mills action in the continuum limit a→0.
 *
 * Complexity: O(1) per plaquette
 */
Plaquette pb_plaquette(const PrincipalBundle *pb, int site, int mu, int nu);

/**
 * Compute field strength F_{μν}(x) ∈ g from plaquette.
 * F = -log(U_{μν}) / a^2 (for small a).
 * Complexity: O(1)
 */
LieAlgebraElement pb_field_strength(const PrincipalBundle *pb, int site,
                                     int mu, int nu);

/**
 * Compute all plaquettes on the lattice.
 * Returns array of size n_sites * (dim choose 2).
 * Complexity: O(n_sites * dim^2)
 */
Plaquette *pb_all_plaquettes(const PrincipalBundle *pb, int *n_plaq);

/** Free array returned by pb_all_plaquettes. */
void pb_plaquettes_free(Plaquette *plaqs);

/* ============================================================================
 * Wilson Action (L4-L5)
 * =========================================================================== */

/**
 * Compute the Wilson gauge action for a single plaquette.
 * S_{μν}(x) = β · (1 - (1/N) Re Tr U_{μν}(x))
 * where β = 2N/g^2 (for SU(N)) or β = 1/g^2 (for U(1)).
 *
 * Theorem (Wilson, 1974): As a → 0,
 *   β Σ (1 - Re Tr U_{μν}/N) → (1/4g^2) ∫ d^4x Tr F_{μν}F^{μν}
 *
 * Complexity: O(1)
 */
double pb_wilson_action_plaq(const PrincipalBundle *pb, int site,
                              int mu, int nu, double beta);

/** Compute total Wilson action summed over all plaquettes. O(n_links) */
double pb_wilson_action_total(const PrincipalBundle *pb, double beta);

/** Compute the topological charge density q(x) from plaquettes. O(1) */
double pb_topological_charge_density(const PrincipalBundle *pb, int site);

/* ============================================================================
 * Bianchi Identity Verification (L4)
 * =========================================================================== */

/**
 * Verify the Bianchi identity on the lattice.
 * D_μ F_{νρ} + D_ν F_{ρμ} + D_ρ F_{μν} = 0
 *
 * On the lattice, this becomes the coclosure condition:
 * The product of six plaquettes around a 3-cube equals identity.
 *
 * Theorem (Bianchi, 1902): The exterior covariant derivative of
 * the curvature vanishes identically. DΩ = dΩ + [ω, Ω] = 0.
 *
 * Complexity: O(1) per 3-cube
 */
bool pb_check_bianchi(const PrincipalBundle *pb, int site,
                      int mu, int nu, int rho, double tolerance);

/**
 * Verify Bianchi identity at all sites and all triplets of directions.
 * Returns the maximum deviation from identity across all 3-cubes.
 * Complexity: O(n_sites * dim^3)
 */
double pb_bianchi_max_deviation(const PrincipalBundle *pb);

/* ============================================================================
 * Chern Classes / Topological Invariants (L8)
 * =========================================================================== */

/**
 * Compute the first Chern number (U(1) bundles in 2D).
 * c_1 = (1/2π) ∫ d^2x F_{12}(x)
 *
 * On the lattice: c_1 = (1/2π) Σ_{sites} arg(U_{12}(x))
 *
 * Theorem (Chern-Weil): The Chern number is an integer topological invariant,
 * independent of the connection.
 *
 * Complexity: O(n_sites)
 */
double pb_chern_number_1(const PrincipalBundle *pb);

/**
 * Compute the second Chern number / instanton number (SU(2) bundles in 4D).
 * c_2 = (1/8π^2) ∫ d^4x Tr(F ∧ F)
 *
 * On the lattice: c_2 = (1/16π^2) Σ ε_{μνρσ} Tr(U_{μν}U_{ρσ})
 *
 * Theorem (Belavin et al., 1975): The instanton number is an integer
 * topological charge, characterizing the bundle's topology.
 *
 * Complexity: O(n_sites)
 */
double pb_chern_number_2(const PrincipalBundle *pb);

#endif /* CURVATURE_H */
