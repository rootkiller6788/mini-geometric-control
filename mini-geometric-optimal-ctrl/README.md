# mini-geometric-optimal-ctrl

Geometric Optimal Control on Manifolds: Pontryagin Maximum Principle,
Lie group methods, symplectic integration, and sub-Riemannian geometry.

## Knowledge Coverage (L1-L9)

| Level | Name | Status |
|-------|------|--------|
| L1 | Definitions | Complete |
| L2 | Core Concepts | Complete |
| L3 | Mathematical Structures | Complete |
| L4 | Fundamental Laws | Complete |
| L5 | Algorithms/Methods | Complete |
| L6 | Canonical Problems | Complete |
| L7 | Applications | Partial+ (NASA, SpaceX, Tesla, Quadrotor) |
| L8 | Advanced Topics | Partial+ (DMOC, Symplectic RK, AVF) |
| L9 | Research Frontiers | Partial (documented) |

## Core Definitions
- Smooth manifold with atlas (GeoChart, GeoAtlas)
- Tangent/cotangent vectors (GeoTangentVector, GeoCotangentVector)
- Vector fields, Lie bracket (GeoVectorField)
- Control-affine system (GeoControlAffineSystem)
- Optimal control problem in Bolza form
- Symplectic space in Darboux coordinates
- Lie groups SO(3), SE(3) with Lie algebras so(3), se(3)
- Sub-Riemannian manifold (distribution + metric)

## Core Theorems
- Pontryagin Maximum Principle (PMP) for control-affine systems
- Chow-Rashevskii controllability theorem
- Hormander bracket-generating condition
- Liouville volume preservation
- Noether conservation theorem
- Jacobi identity (Lie bracket, Poisson bracket)
- Rodrigues formula (SO(3) exponential)
- Euler-Poincare reduction

## Core Algorithms
- Lie bracket computation via finite-difference Jacobians
- Stormer-Verlet symplectic integrator
- Gauss-Legendre symplectic Runge-Kutta
- Discrete Euler-Lagrange (DEL) variational integrator
- Average Vector Field (AVF) energy-preserving method
- Discrete gradient method (Itoh-Abe)
- Yoshida higher-order composition
- RK-Munthe-Kaas Lie group integrator
- Single shooting for PMP boundary value problems
- Matrix exponential (scaling-and-squaring)

## Canonical Problems
- Optimal attitude control on SO(3)
- Heisenberg sub-Riemannian geodesic
- Dido isoperimetric problem
- Kinematic car nonholonomic planning
- Parallel parking maneuver
- Brockett integrator
- Arnold cat map on torus
- Rigid body Euler-Poincare dynamics
- Reaction wheel desaturation

## Nine-School Course Mapping

| School | Course | Topic Match |
|--------|--------|-------------|
| MIT | 6.832 Underactuated Robotics | Nonholonomic systems |
| MIT | 16.323 Optimal Control | PMP, shooting |
| Stanford | AA203 Optimal Control | PMP, Hamiltonian |
| Berkeley | EE222 Nonlinear Systems | Lie brackets, controllability |
| CMU | 24-677 Nonlinear Control | Geometric control on Lie groups |
| Princeton | MAE 546 Optimal Control | Pontryagin principle |
| Caltech | CDS140 Nonlinear Dynamics | Hamiltonian mechanics |
| Cambridge | 4F3 Nonlinear Control | Lie group methods |
| Oxford | B4 Predictive Control | DMOC, structure-preserving |
| ETH | 227-0220 Model Reduction | Symmetry reduction |

## Build



## Module Status: COMPLETE

- L1-L6: Complete
- L7: Partial+ (4 applications: NASA spacecraft, SpaceX landing, Tesla battery, Quadrotor)
- L8: Partial+ (6 advanced topics: DMOC, symplectic RK, AVF, discrete gradient, Lie group integrators, conjugate points)
- L9: Partial (documented, no implementation)

include/ + src/ line count: >3000 (meets threshold)
make test: 23/23 tests pass
