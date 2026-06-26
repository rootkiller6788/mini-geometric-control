# Knowledge Graph — Symmetry Reduction

## L1: Definitions

| # | Item | C Implementation | Lean Formalization |
|---|------|-----------------|-------------------|
| 1 | Lie group G | `SRLieGroup` struct | `structure LieAlgebra` |
| 2 | Lie algebra g = T_e G | `SRLieAlgebra` struct | `structure LieAlgebra` |
| 3 | Group element | `SRGroupElement` struct | — |
| 4 | Lie algebra element (infinitesimal generator) | `SRLieAlgebraElement` struct | — |
| 5 | Structure constants c_{ij}^k | `SRLieAlgebra.constants` | `bracket` field |
| 6 | Adjoint action Ad_g : g → g | `SRAdjointAction` struct | — |
| 7 | Coadjoint action Ad*_{g^{-1}} : g* → g* | `SRCoadjointAction` struct | `coadjointAction` def |
| 8 | Momentum map J : P → g* | `SRMomentumMap` struct | `MomentumMap` structure |
| 9 | Symplectic manifold (P, ω) | `SRSymplecticManifold` struct | `SymplecticManifold` structure |
| 10 | Poisson manifold | `SRPoissonManifold` struct | `LiePoisson` structure |
| 11 | Lie-Poisson bracket | `SRLiePoissonBracket` struct | `isCasimir` def |
| 12 | Coadjoint orbit O_μ | — | `coadjointOrbit` def |
| 13 | Marsden-Weinstein reduced space P_μ | `SRMarsdenWeinsteinSpace` struct | — |
| 14 | Semidirect product G = S ⋉ V | `SRSemidirectProduct` struct | — |
| 15 | Invariant Hamiltonian H(g·z) = H(z) | `SRInvariantHamiltonian` struct | — |
| 16 | Reduced control system | `SRReducedControlSystem` struct | — |

## L2: Core Concepts

| # | Item | Implementation |
|---|------|---------------|
| 1 | Group action (free, proper, Hamiltonian) | `SRActionType` enum |
| 2 | Momentum map equivariance | `sr_momentum_verify_equivariance` |
| 3 | Symplectic reduction | `sr_mw_reduce`, `sr_mw_compute_reduced_omega` |
| 4 | Cotangent bundle reduction T*G/G ≅ g* | `sr_cotangent_reduce` |
| 5 | Euler-Poincare reduction | `sr_ep_create`, `sr_ep_rhs` |
| 6 | Lagrange-Poincare reduction | `sr_lp_create`, `sr_lp_rhs` |
| 7 | Semidirect product reduction | `sr_semidirect_create` |
| 8 | Poisson reduction | `sr_poisson_reduction` |
| 9 | Reduction by stages | `sr_reduce_by_stages` |
| 10 | Singular reduction (stratified spaces) | `sr_singular_stratification` |
| 11 | Lie algebra cohomology | `sr_ce_coboundary_1form`, `sr_ce_2cocycle_condition` |
| 12 | Dual pair | `SRDualPair`, `sr_dual_pair_verify` |

## L3: Mathematical Structures

| # | Item | Implementation |
|---|------|---------------|
| 1 | so(3) Lie algebra | `sr_algebra_create_so3` |
| 2 | so(2) Lie algebra | `sr_algebra_create_so2` |
| 3 | su(2) Lie algebra | `sr_algebra_create_su2` |
| 4 | se(3) Lie algebra | `sr_algebra_create_se3` |
| 5 | Heisenberg algebra | `sr_algebra_create_heisenberg` |
| 6 | Abelian Lie algebra | `sr_algebra_create_abelian` |
| 7 | Lie bracket computation | `sr_algebra_bracket` |
| 8 | Killing form | `sr_algebra_killing_form` |
| 9 | Jacobi identity verification | `sr_algebra_verify_jacobi` |
| 10 | Exponential map (Rodrigues) | `sr_exp_map_so3`, `sr_exp_map_se3` |
| 11 | BCH formula (2nd order) | `sr_bch_formula_2nd` |
| 12 | Canonical symplectic form | `sr_symplectic_set_canonical` |
| 13 | Poisson tensor computation | `sr_poisson_tensor_at_point` |

## L4: Fundamental Laws

| # | Theorem | C Verification | Lean Formalization |
|---|---------|---------------|-------------------|
| 1 | Noether's theorem | `sr_noether_verify`, `sr_noether_conserved_quantities` | `noether_conservation` |
| 2 | Marsden-Weinstein reduction | `sr_mw_reduce`, `sr_mw_compute_reduced_omega` | `marsden_weinstein_reduction` |
| 3 | Lie-Poisson reduction | `sr_cotangent_reduce` | `lie_poisson_reduction` |
| 4 | Cotangent bundle reduction | `sr_cotangent_reduce` | `cotangent_bundle_reduction` |
| 5 | Euler-Poincare reduction | `sr_ep_rhs` | `euler_poincare_reduction` |
| 6 | Semidirect product reduction | `sr_semidirect_lie_poisson_rhs` | `semidirect_product_reduction` |
| 7 | Jacobi identity (bracket) | `sr_algebra_verify_jacobi` | — |
| 8 | KKS symplectic form | `sr_kks_symplectic_form` | `kks_symplectic_form` |
| 9 | Reduction by stages | `sr_reduce_by_stages` | `reduction_by_stages` |
| 10 | Energy-momentum stability | `sr_energy_momentum_stability` | `energy_momentum_stability_theorem` |
| 11 | Poisson reduction | `sr_poisson_reduction` | `poisson_reduction_theorem` |

## L5: Algorithms

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | Lie bracket computation | `sr_algebra_bracket` |
| 2 | Exponential map (SO(3) Rodrigues) | `sr_exp_map_so3` |
| 3 | BCH formula | `sr_bch_formula_2nd` |
| 4 | RK4 integration on coadjoint orbit | `sr_integrate_lie_poisson_rk4` |
| 5 | Symplectic Euler integrator | `sr_integrate_symplectic_euler` |
| 6 | Stormer-Verlet integrator | `sr_integrate_stormer_verlet` |
| 7 | Lie group variational integrator | `sr_integrate_lie_group_var` |
| 8 | RATTLE constrained integrator | `sr_integrate_rattle` |
| 9 | Momentum map computation | `sr_momentum_compute` |
| 10 | Isotropy subalgebra computation | `sr_isotropy_subalgebra_basis` |
| 11 | Coadjoint orbit sampling (Monte Carlo) | `sr_orbit_monte_carlo` |
| 12 | Relative equilibrium finding | `sr_find_relative_equilibrium` |
| 13 | Linearized reduced stability | `sr_linearize_reduced_dynamics` |

## L6: Canonical Problems

| # | Problem | Implementation |
|---|---------|---------------|
| 1 | Free rigid body (SO(3) reduction) | `sr_ep_rhs`, `sr_rigid_body_poisson` |
| 2 | Heavy top (semidirect SE(3) reduction) | `sr_heavy_top_bracket`, `sr_heavy_top_casimir` |
| 3 | Spherical pendulum (S^1 reduction) | `sr_spherical_pendulum_reduced_rhs` |
| 4 | Double spherical pendulum | `sr_double_spherical_rhs` |
| 5 | Spacecraft attitude control | `sr_spacecraft_reduced_rhs` |
| 6 | Underwater vehicle (Kirchhoff equations) | `sr_underwater_vehicle_reduced_rhs` |
| 7 | Quadrotor UAV | `sr_quadrotor_reduced_rhs` |
| 8 | Robot arm (cyclic coordinate reduction) | `sr_robot_arm_reduced_rhs` |

## L7: Applications

| # | Application | Implementation |
|---|------------|---------------|
| 1 | Spacecraft attitude control (NASA) | `sr_spacecraft_reduced_rhs` |
| 2 | Quadrotor UAV control (Tesla/SpaceX) | `sr_quadrotor_reduced_rhs` |
| 3 | Underwater vehicle dynamics (AUV/REMUS) | `sr_underwater_vehicle_reduced_rhs` |
| 4 | Robot arm control | `sr_robot_arm_reduced_rhs` |

## L8: Advanced Topics

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Controlled Euler-Poincare equations | `sr_ctrl_ep_rhs` |
| 2 | Controlled Lagrangians method | `sr_controlled_lagrangian_modify`, `sr_ctrl_matching_conditions` |
| 3 | Energy shaping via Casimir functions | `sr_energy_shape_casimirs` |
| 4 | Nonholonomic reduced dynamics | `sr_nonholonomic_reduced_rhs` |
| 5 | Inverse optimal reduced control | `sr_inverse_optimal_reduced` |
| 6 | LQR on reduced space | `sr_lqr_reduced` |
| 7 | Geometric phase computation | `sr_geometric_phase` |
| 8 | Hamiltonian-Hopf bifurcation detection | `sr_detect_hamiltonian_hopf` |
| 9 | Arnold stability analysis | `sr_arnold_stability` |
| 10 | Poisson reduction | `sr_poisson_reduction` |
| 11 | Reduction by stages | `sr_reduce_by_stages` |
| 12 | Singular reduction stratification | `sr_singular_stratification` |

## L9: Research Frontiers

| # | Topic | Status |
|---|-------|--------|
| 1 | Singular symplectic reduction (Sjamaar-Lerman) | Documented, Lean theorem |
| 2 | Optimal transport on coadjoint orbits | Documented, Lean theorem |
| 3 | Symplectic integrators on Lie groups | Partial implementation |
| 4 | Meta-complexity of symmetry detection | Documented only |
| 5 | Geometric deep learning on homogeneous spaces | Documented only |
