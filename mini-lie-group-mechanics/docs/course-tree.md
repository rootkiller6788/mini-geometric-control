# Course Tree -- mini-lie-group-mechanics

## Prerequisite Dependency Graph

```
Linear Algebra (matrices, determinants, eigenvalues)
    │
    ├── Group Theory (groups, subgroups, homomorphisms)
    │       │
    │       ├── Lie Groups (smooth manifold + group structure)
    │       │       │
    │       │       ├── Matrix Lie groups SO(3), SE(3), SU(2)
    │       │       ├── Lie algebra (tangent space at identity)
    │       │       ├── Exponential map exp: g → G
    │       │       ├── Logarithm map log: G → g
    │       │       ├── Adjoint representation Ad, ad
    │       │       ├── Coadjoint representation Ad*, ad*
    │       │       └── Killing form
    │       │
    │       └── Group Actions on Manifolds
    │               │
    │               ├── Momentum maps (Noether's theorem)
    │               ├── Symmetry reduction
    │               └── Coadjoint orbits
    │
    ├── Differential Geometry
    │       │
    │       ├── Tangent bundle TG, cotangent bundle T*G
    │       ├── Left/right trivialization
    │       ├── Maurer-Cartan form
    │       └── Connection and curvature
    │
    └── Classical Mechanics
            │
            ├── Lagrangian mechanics (TG → R)
            ├── Hamiltonian mechanics (T*G → R)
            ├── Euler-Poincare reduction (TG/G → g)
            ├── Lie-Poisson reduction (T*G/G → g*)
            └── Euler rigid body equations
                    │
                    ├── Free rigid body (Euler-Poinsot)
                    ├── Heavy top (Lagrange-Poisson)
                    └── Multi-body systems

Geometric Integration
    │
    ├── RKMK (Runge-Kutta-Munthe-Kaas)
    ├── Crouch-Grossman (composition of exponentials)
    ├── Magnus expansion
    └── Variational integrators
            │
            └── Discrete Euler-Poincare

Applications
    │
    ├── Spacecraft attitude (SO(3))
    ├── Quadrotor UAV (SE(3))
    ├── Underwater vehicles (SE(3) with added mass)
    ├── Robotic arms (SE(3) kinematics)
    └── Inertial navigation (SO(3) IMU propagation)

Research Frontiers (L9)
    │
    ├── Geometric Complexity Theory (GCT)
    ├── Quantum control on Lie groups
    ├── Equivariant neural networks
    ├── Stochastic geometric mechanics
    └── Infinite-dimensional Lie groups (EPDiff, Camassa-Holm)
```

## Module Dependencies

This module depends on:
- **Linear system theory** (matrix operations, notation)
- **Nonlinear dynamics** (state-space, ODEs)
- **Stability theory** (Lyapunov functions on manifolds)

This module is a prerequisite for:
- **mini-connection-principal-bundle** (Connection theory on Lie groups)
- **mini-hamiltonian-control** (Hamiltonian mechanics on symplectic manifolds)
- **mini-geometric-optimal-ctrl** (Optimal control on Lie groups)
- **mini-nonholonomic-planning** (Nonholonomic systems on principal bundles)
- **mini-symmetry-reduction** (Advanced symmetry reduction)

