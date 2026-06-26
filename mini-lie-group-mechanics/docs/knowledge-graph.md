# Knowledge Graph -- mini-lie-group-mechanics

## L1: Definitions (Complete)

| Concept | C Implementation | Lean Formalization |
|---------|-----------------|-------------------|
| Lie group (matrix group) | `LieGroupElement` struct, `LieGroupType` enum | N/A (algebraic structure) |
| Lie algebra | `LieAlgebraElement` struct, `LieAlgebraType` enum | `Vec3` structure |
| Lie bracket [.,.] | `lie_bracket()` | `LieBracket` |
| Exponential map exp: g → G | `lie_exp_so3()`, `lie_exp_se3()`, `lie_matrix_exponential()` | `rodriguesForm` |
| Logarithm map log: G → g | `lie_log_so3()`, `lie_log_se3()`, `lie_matrix_logarithm()` | N/A |
| Adjoint representation Ad | `lie_Ad()`, `lie_Ad_matrix_so3()`, `lie_Ad_matrix_se3()` | Implicit |
| Coadjoint representation ad* | `lie_ad_star_so3()` | `discreteMomentumUpdate` |
| Hat/vee isomorphism | `lie_so3_hat()`, `lie_so3_vee()`, `lie_se3_hat()`, `lie_se3_vee()` | `Mat3x3.hat` |
| Killing form | `lie_killing_form()` | Implicit |
| SO(3), SE(3) groups | `LIE_GROUP_SO3`, `LIE_GROUP_SE3` | N/A |

## L2: Core Concepts (Complete)

| Concept | Implementation |
|---------|---------------|
| Group action on manifold | `lie_group_act_vector()`, `lie_manifold_state_create()` |
| Left/right translations | `lie_left_translate()`, `lie_right_translate()` |
| Tangent map | `lie_tangent_left_translate()` |
| Maurer-Cartan form | `lie_maurer_cartan_form()` |
| Momentum map | `lie_momentum_map_so3_position()`, `lie_angular_momentum_compute()` |
| Coadjoint orbits | `lie_coadjoint_orbit_radius_so3()`, `sameCoadjointOrbit` |
| Locked inertia tensor | `LieInertiaTensor`, `lie_inertia_apply()` |
| Symmetry reduction | Euler-Poincare, Lie-Poisson equations |

## L3: Mathematical Structures (Complete)

| Structure | Implementation |
|-----------|---------------|
| Matrix operations | `LieMatrix` with create, add, mul, inv, det, transpose, trace |
| Vector operations | `LieVector` with create, cross, dot, norm, normalize |
| so(3) skew-symmetric matrices | `lie_so3_hat()`, `lie_matrix_is_skew_symmetric()` |
| se(3) 4×4 matrices | `lie_se3_hat()`, `lie_se3_vee()` |
| Rotation matrices | Verified via `lie_matrix_is_special_orthogonal()` |
| Inertia tensor | `LieInertiaTensor` (3x3 symmetric positive definite) |
| Rigid body state | `LieRigidBodyState` |
| Heavy top state | `LieHeavyTopState` |
| Quadrotor state | `LieQuadrotorState` (SE(3) × R^6) |

## L4: Fundamental Laws/Theorems (Complete)

| Theorem | Code Verification | Lean Formalization |
|---------|------------------|-------------------|
| Rodrigues' formula | `lie_exp_so3()` implements exp via closed form | `rodrigues_is_rotation` |
| Euler-Rodrigues (log) | `lie_log_so3()` implements log via (R-R^T)^vee | N/A |
| Jacobi identity | `lie_jacobi_identity_check()` (numerical verification) | `bracket_jacobi` |
| Euler-Poincare theorem | `lie_euler_poincare_rhs()` | N/A |
| Noether's theorem | `lie_noether_conserved_quantity()` | N/A |
| Euler rigid body equations | `lie_rigid_body_rhs()` | `eulerRHS`, `energy_conservation_example` |
| Heavy top equations | `lie_heavy_top_rhs()` | N/A |
| Lie-Poisson bracket | `lie_lie_poisson_bracket_so3()` | N/A |
| BCH formula | `lie_bch_2nd_order()`, `lie_bch_3rd_order()` | N/A |
| Adjoint/codajoint properties | `lie_Ad()`, `lie_ad_star_so3()` | `bracket_antisym_direct` |

## L5: Algorithms/Methods (Complete)

| Algorithm | Implementation |
|-----------|---------------|
| Matrix exponential (scaling+squaring) | `lie_matrix_exponential()` |
| SO(3) exponential (Rodrigues closed form) | `lie_exp_so3()` |
| SE(3) exponential (closed form) | `lie_exp_se3()` |
| SO(3) logarithm | `lie_log_so3()` |
| SE(3) logarithm | `lie_log_se3()` |
| Matrix logarithm (inverse scaling+squaring) | `lie_matrix_logarithm()` |
| RKMK-RK4 integrator on SO(3) | `lie_rkmk_step_so3()`, `lie_rkmk4_integrate_so3()` |
| Lie-Euler integrator | `lie_euler_step_so3()` |
| Magnus expansion (2nd order) | `lie_magnus_expansion_order2()` |
| Variational integrator | `lie_variational_rigid_body_step()` |
| Geodesic interpolation on SO(3) | `lie_so3_geodesic_interpolate()` |
| BCH formula computation | `lie_bch_2nd_order()`, `lie_bch_3rd_order()` |

## L6: Canonical Problems (Complete)

| Problem | Example |
|---------|---------|
| Rigid body rotation | `example_rigid_body.c` — Euler equations on SO(3) |
| Heavy top dynamics | `example_heavy_top.c` — Lagrange-Poisson equations |
| Quadrotor control on SE(3) | `example_quadrotor.c` — Geometric tracking controller |

## L7: Applications (Complete — 5 applications)

| Application | File | Key Parameter |
|------------|------|--------------|
| Apollo Lunar Module | `lie_applications.c` | Ixx=15000, Iyy=35000, Izz=25000 kg·m² |
| SpaceX Falcon 9 separation | `lie_applications.c` | Stage 2 on SE(3): 111,500 kg |
| Mars helicopter (NASA) | `lie_applications.c` | 1.8 kg, Mars g=3.721 m/s² |
| GPS satellite orientation | `lie_applications.c` | MEO 20,200 km, gravity gradient |
| Tesla IMU navigation | `lie_applications.c` | SO(3) gyro propagation |

## L8: Advanced Topics (Partial+ — 3 topics)

| Topic | Implementation |
|-------|---------------|
| Geometric integration (RKMK) | `lie_rkmk_step_so3()`, `lie_rkmk4_integrate_so3()` |
| Variational integrators | `lie_variational_rigid_body_step()` |
| Magnus expansion | `lie_magnus_expansion_order2()` |
| Stochastic Hamiltonian systems | Not implemented |
| Lie group methods for PDEs | Not implemented |

## L9: Research Frontiers (Partial)

| Topic | Status |
|-------|--------|
| Geometric complexity theory (GCT) & Lie groups | Documented |
| Quantum Lie group symmetries | Documented |
| Infinite-dimensional Lie groups | Documented |
| Meta-stable geometric control | Documented |

