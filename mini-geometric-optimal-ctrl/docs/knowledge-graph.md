# Knowledge Graph ¡ª mini-geometric-optimal-ctrl

## L1: Definitions (Complete)
- Smooth manifold, atlas, chart (GeoChart, GeoAtlas)
- Tangent vector, cotangent vector (GeoTangentVector, GeoCotangentVector)
- Vector field, Lie bracket (GeoVectorField)
- Control-affine system (GeoControlAffineSystem)
- Optimal control problem in Bolza form (GeoOptimalControlProblem)
- Symplectic vector space, canonical symplectic form (GeoSymplecticSpace)
- Hamiltonian function, Hamiltonian vector field
- Lie group SO(3), Lie algebra so(3)
- Lie group SE(3), Lie algebra se(3)
- Sub-Riemannian manifold, horizontal distribution (GeoSubRiemannianManifold)
- Discrete Lagrangian (GeoDiscreteLagrangian)
- Poisson bracket, Liouville 1-form

## L2: Core Concepts (Complete)
- Lie bracket as failure of commuting flows
- Flow of vector field as one-parameter group
- Hamiltonian mechanics: Hamilton equations
- Symplectic flat/sharp maps (Omega-flat, Omega-sharp)
- Adjoint representations Ad and ad
- Exponential map on Lie groups
- Controllability via Lie algebra rank condition (LARC)
- Pontryagin Maximum Principle (PMP)
- Symplectic integrators (preserving omega)
- Energy-preserving integrators (AVF, discrete gradient)
- Discrete Noether theorem (via DEL)
- Yoshida composition for higher-order methods

## L3: Mathematical Structures (Complete)
- Atlas with finite chart covering (GeoAtlas)
- Tangent bundle TM, cotangent bundle T*M
- Product manifold M1 x M2
- Canonical Darboux coordinates (q,p)
- Lie group SO(3) as matrix manifold
- Exponential map as geodesic on Lie group
- Euler-Poincare reduced equations on so(3)*
- Sub-Riemannian structure (M, Delta, g)
- Heisenberg group as step-2 nilpotent Lie group
- Butcher tableau for Runge-Kutta methods
- Discrete mechanics on Q x Q

## L4: Fundamental Laws/Theorems (Complete)
- Pontryagin Maximum Principle (state-costate dynamics)
- Chow-Rashevskii controllability theorem
- Hormander bracket-generating condition
- Transversality conditions (free/fixed final state/time)
- Liouville theorem (Hamiltonian flow preserves symplectic form)
- Noether theorem (symmetries imply conservation laws)
- Jacobi identity (Lie bracket, Poisson bracket)
- Rodrigues formula for SO(3) exponential
- Euler-Poincare theorem (reduction by symmetry)
- Legendre transform as involution

## L5: Algorithms/Methods (Complete)
- RK4 integration of flows
- Finite-difference Lie bracket computation
- Stormer-Verlet symplectic integrator
- Symplectic Runge-Kutta (Gauss-Legendre)
- Discrete Euler-Lagrange (DEL) integrator
- Average Vector Field (AVF) energy-preserving method
- Discrete gradient method (Itoh-Abe)
- Yoshida composition for higher-order symplectic
- RK-Munthe-Kaas Lie group integrator
- Crouch-Grossman commutator-free integrator
- Single shooting for PMP boundary value problems
- Matrix exponential (scaling-and-squaring)
- Discrete Legendre transform

## L6: Canonical Problems (Complete)
- Optimal attitude control on SO(3) (reorientation)
- Sub-Riemannian Heisenberg geodesic problem
- Dido isoperimetric problem on Heisenberg group
- Kinematic car nonholonomic motion planning
- Parallel parking maneuver
- Brockett integrator (canonical nonholonomic)
- Arnold cat map on torus (chaotic dynamics)
- Rigid body Euler-Poincare dynamics
- Reaction wheel desaturation (spacecraft)

## L7: Applications (Partial+)
- NASA Apollo-style attitude maneuver (historical)
- SpaceX Falcon 9 landing burn computation
- Tesla optimal battery thrust optimization
- Quadrotor motor mixing (attitude to thrust)
- Reaction wheel momentum management
- Spacecraft attitude control (SO(3))

## L8: Advanced Topics (Partial+)
- Sub-Riemannian geometry on Heisenberg group
- Symplectic RK methods (Gauss-Legendre collocation)
- Discrete Mechanics and Optimal Control (DMOC)
- Variational integrators (DEL + discrete Legendre)
- Energy-preserving methods (AVF, discrete gradient)
- Lie group integrators (RK-MK, Crouch-Grossman)
- Conjugate point analysis for optimality

## L9: Research Frontiers (Partial)
- Geometric optimal control on SO(3) and SE(3)
- Structure-preserving discretization for optimal control
- Sub-Riemannian geodesic computation
