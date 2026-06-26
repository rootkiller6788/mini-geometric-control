# Course Tree — Nonholonomic Motion Planning

## Prerequisites

This module depends on the following foundational topics:

```
mini-nonholonomic-planning
├── mini-lie-group-mechanics (SE(2)/SE(3) geometry)
│   └── Lie groups, exponential map, adjoint action
├── mini-geometric-optimal-ctrl (optimal control on manifolds)
│   └── Shooting methods, gradient-based optimization
├── mini-feedback-linearization-geo
│   └── Chained form conversion, differential flatness
├── mini-differential-flatness
│   └── Flatness-based planning (for N ≤ 1 trailers)
├── mini-connection-principal-bundle
│   └── Geometric phase, holonomy, connection forms
├── mini-symmetry-reduction
│   └── Symmetries in nonholonomic systems
├── mini-hamiltonian-control
│   └── Hamiltonian formulation, Poisson brackets
│
├── External prerequisites:
│   ├── Differential geometry (manifolds, vector fields, Lie derivatives)
│   ├── Control theory (controllability, feedback, Lyapunov stability)
│   ├── Robotics (kinematic models, motion planning)
│   └── Numerical methods (ODE integration, linear algebra)
```

## Dependency Tree

```
Differential Geometry
    ├── Manifolds & Tangent Bundles
    │   └── Configuration space Q, TQ, vector fields
    ├── Lie Groups & Lie Algebras
    │   ├── SE(2), SE(3), so(3)
    │   ├── Exponential map, logarithm
    │   └── Lie bracket [·,·]
    │
Control Theory
    ├── Controllability
    │   ├── STLC (small-time local controllability)
    │   ├── LARC (Lie Algebra Rank Condition)
    │   └── Chow-Rashevskii theorem
    ├── Feedback Stabilization
    │   ├── Brockett's necessary condition
    │   └── Time-varying feedback
    └── Optimal Control
        ├── Minimum-time problems
        └── Sub-Riemannian geodesics
            │
            └── Nonholonomic Motion Planning ★ (this module)
                ├── Sampling-based: RRT, PRM, RRT-Connect
                ├── Geometric: sinusoidal, nilpotent steering
                ├── Search-based: lattice A*, Dijkstra on PRM
                ├── Optimization: shortcut, gradient descent
                └── Models: unicycle, car, trailer, snakeboard, ball
```

## Learning Path

1. **Start here**: Unicycle model → Lie bracket → Geometric phase
2. **Controllability**: Frobenius test → Chow-Rashevskii → LARC
3. **Planning**: RRT → PRM → Steering methods
4. **Advanced**: Trailer systems → Nilpotent approximation → Optimal planning

## Difficulty Progression

- **Beginner**: Unicycle model, forward simulation, Lie bracket computation
- **Intermediate**: RRT planning, PRM construction, sinusoidal steering
- **Advanced**: Nilpotent approximation, trailer systems, geometric integration
- **Expert**: Optimal nonholonomic planning, sub-Riemannian geometry
