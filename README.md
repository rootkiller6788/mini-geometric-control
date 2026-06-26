# Mini Geometric Control

A collection of **from-scratch, zero-dependency C implementations** of Geometric Control Theory — the differential-geometric approach to nonlinear control systems pioneered by Brockett, Jurdjevic, Isidori, Nijmeijer, van der Schaft, and others. Each module translates foundational geometric-control formalisms into runnable C code, bridging modern control theory on manifolds with computational practice.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|--------|--------|-------------|
| [mini-connection-principal-bundle](mini-connection-principal-bundle/) | Ehresmann connections, connection 1-forms, curvature 2-forms, holonomy groups, parallel transport, principal G-bundles, lattice gauge theory (link variables, Wilson loops), gauge transformations | MIT 18.755 (Lie Groups & Lie Algebras), MIT 6.832 |
| [mini-differential-flatness](mini-differential-flatness/) | Differential flatness (Fliess–Levine–Martin–Rouchon), flat outputs, relative degree, trajectory generation without ODE integration, Brunovsky normal form, dynamic feedback linearization, multi-derivative chains | MIT 6.832 (Underactuated Robotics), Stanford AA 203 |
| [mini-feedback-linearization-geo](mini-feedback-linearization-geo/) | Input-output / full-state feedback linearization, Lie derivatives & Lie brackets, relative degree analysis, Frobenius theorem for involutive distributions, zero dynamics, nonlinear coordinate transforms, Isidori–Nijmeijer framework | MIT 6.241J (Dynamic Systems & Control), MIT 6.832, Stanford ENGR 205 |
| [mini-geometric-optimal-ctrl](mini-geometric-optimal-ctrl/) | Control-affine systems on smooth manifolds, tangent/cotangent bundles, Pontryagin Maximum Principle on manifolds, sub-Riemannian geometry, geometric integrators (Lie group methods), reachable sets & controllability (Chow's theorem) | MIT 6.832, Stanford AA 203, Caltech CDS 205 |
| [mini-hamiltonian-control](mini-hamiltonian-control/) | Hamiltonian systems on symplectic manifolds, Poisson geometry, Hamilton–Jacobi–Bellman (HJB) equation, energy shaping (PCHD systems), symplectic matrices, phase-space analysis, Pontryagin maximum principle (dual formulation) | MIT 2.151 (Advanced Dynamics), Stanford AA 242 |
| [mini-lie-group-mechanics](mini-lie-group-mechanics/) | Lie groups in mechanics (SO(3), SE(3), SU(2), SO(2), SE(2)), Lie algebra structure constants, group actions & adjoint/coadjoint orbits, Euler–Poincaré equations, Lie–Poisson reduction, geometric integration on Lie groups | MIT 2.032 (Dynamics), MIT 2.151, Caltech CDS 205 |
| [mini-nonholonomic-planning](mini-nonholonomic-planning/) | Nonholonomic constraints & distributions, involutivity & Frobenius theorem, chained-form conversion (Brockett), kinematic car / snake robot / ball-plate systems, motion planning under velocity constraints, Laumond–Murray framework | MIT 6.832, MIT 6.821 (Robotic Manipulation) |
| [mini-symmetry-reduction](mini-symmetry-reduction/) | Marsden–Weinstein symplectic reduction, momentum map, coadjoint orbits (Arnold), Poisson reduction, Euler–Poincaré reduction, free/proper group actions, Hamiltonian reduction, reduced phase spaces | MIT 18.755, Caltech CDS 205, MIT 2.151 |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module includes `docs/` with alignment to original geometric control source works (Brockett, Jurdjevic, Isidori, Marsden, Agrachev & Sachkov)
- **From-scratch implementations** — no existing libraries; each concept is built from geometric and control-theoretic primitives (Lie brackets, differential forms, symplectic structures)

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-lie-group-mechanics
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-geometric-control/
├── mini-connection-principal-bundle/ # Ehresmann connections, holonomy, curvature on principal bundles
├── mini-differential-flatness/       # Flat outputs, trajectory generation, Brunovsky normal form
├── mini-feedback-linearization-geo/  # Lie-derivative based feedback linearization, zero dynamics
├── mini-geometric-optimal-ctrl/      # PMP on manifolds, sub-Riemannian geometry, geometric integrators
├── mini-hamiltonian-control/         # Symplectic geometry, HJB equation, energy shaping
├── mini-lie-group-mechanics/         # SO(3)/SE(3) kinematics, Euler–Poincaré, Lie–Poisson reduction
├── mini-nonholonomic-planning/       # Distributions, chained forms, motion planning under constraints
└── mini-symmetry-reduction/          # Marsden–Weinstein reduction, momentum map, coadjoint orbits
```

## License

MIT
