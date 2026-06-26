# Course Tree — Symmetry Reduction

## Prerequisite Knowledge Dependencies

```
Symmetry Reduction
├── Lie Group Theory
│   ├── Group actions (free, proper, transitive)
│   ├── Lie algebra, structure constants
│   ├── Adjoint/coadjoint representations
│   ├── Exponential map, BCH formula
│   └── Coadjoint orbits, KKS symplectic form
├── Symplectic Geometry
│   ├── Symplectic manifolds, Darboux theorem
│   ├── Hamiltonian vector fields
│   ├── Poisson brackets, Poisson manifolds
│   └── symplectic leaves (coadjoint orbits)
├── Hamiltonian Mechanics
│   ├── Phase space, canonical coordinates
│   ├── Noether's theorem
│   └── Integrability, action-angle variables
├── Geometric Control Theory
│   ├── Control-affine systems
│   ├── Controllability (Chow-Rashevskii)
│   ├── Optimal control (PMP)
│   └── Feedback stabilization
└── Numerical Analysis
    ├── RK4 integration
    ├── Geometric integrators (symplectic, variational)
    └── Constrained dynamics (RATTLE, SHAKE)
```

## Module Progression

1. **L1-L2**: Learn Lie groups/algebras, structure constants → `sr_algebra_create_so3`, `sr_algebra_bracket`
2. **L3**: Build symplectic/Poisson structures → `sr_symplectic_set_canonical`, `sr_poisson_set_lie_poisson`
3. **L4**: Prove reduction theorems → `sr_mw_reduce`, `sr_ep_rhs`, `sr_noether_verify`
4. **L5**: Implement algorithms → `sr_integrate_lie_poisson_rk4`, `sr_exp_map_so3`
5. **L6**: Solve canonical problems → `sr_rigid_body_poisson`, `sr_heavy_top_bracket`
6. **L7**: Apply to engineering → `sr_spacecraft_reduced_rhs`, `sr_quadrotor_reduced_rhs`
7. **L8**: Study advanced methods → `sr_controlled_lagrangian_modify`, `sr_geometric_phase`
8. **L9**: Explore frontiers → singular reduction, geometric deep learning

## Research Frontiers (L9)

- Singular reduction: stratification theory for non-regular momentum values (Sjamaar-Lerman 1991)
- Optimal transport on coadjoint orbits: Wasserstein geometry on g*
- Geometric deep learning: equivariant neural networks on homogeneous spaces
- Meta-complexity: detecting hidden symmetries in dynamical systems
- Symplectic integrators: higher-order variational integrators on Lie groups
