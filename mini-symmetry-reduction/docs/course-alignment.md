# Course Alignment — Symmetry Reduction

## Nine-School Curriculum Mapping

| School | Course | Topic | Our Implementation |
|--------|--------|-------|-------------------|
| **MIT** | 6.832 Underactuated Robotics | Lie group symmetries, reduction | `sr_ep_create`, `sr_ctrl_ep_rhs` |
| **MIT** | 2.154/6.141 Nonlinear Control | Geometric nonlinear control | `sr_spacecraft_reduced_rhs` |
| **Stanford** | AA203 Optimal Control | Symmetry in optimal control | `sr_reduce_optimal_control`, `sr_lqr_reduced` |
| **Stanford** | AA274 Multi-agent Control | Group-theoretic methods | `sr_kinematic_rhs`, `sr_kinematic_controllable` |
| **Berkeley** | EE222 Nonlinear Systems | Reduction of mechanical systems | `sr_mw_reduce`, `sr_ep_rhs` |
| **CMU** | 24-677 Nonlinear Control | Geometric mechanics | `sr_lp_create`, `sr_lp_rhs` |
| **Princeton** | MAE 546 Optimal Control | Hamiltonian systems | `sr_integrate_symplectic_euler` |
| **Caltech** | CDS140 Nonlinear Dynamics | Symmetry and bifurcation | `sr_detect_hamiltonian_hopf` |
| **Cambridge** | Part II Geometry | Symplectic geometry | `sr_symplectic_create`, `sr_symplectic_set_canonical` |
| **Oxford** | C5.2 Geometric Mechanics | Lie-Poisson brackets | `sr_lpb_create`, `sr_rigid_body_poisson` |
| **ETH** | 151-0532 Geometric Mechanics | Poisson reduction | `sr_poisson_reduction` |

## Textbook Alignment

| Textbook | Chapter | Implementation |
|----------|---------|---------------|
| Marsden & Ratiu (1999) | Ch. 9-13 (Momentum Maps, Reduction) | `sr_momentum_*`, `sr_mw_*`, `sr_ep_*` |
| Arnold (1989) | Appendices 2, 5 (Lie groups, symplectic) | `sr_algebra_*`, `sr_symplectic_*` |
| Holm, Schmah & Stoica (2009) | Ch. 2-4 (EP reduction) | `sr_ep_*`, `sr_lp_*` |
| Abraham & Marsden (1978) | Ch. 4 (Hamiltonian systems with symmetry) | `sr_mw_*`, `sr_noether_*` |
| Bloch (2003) | Ch. 5-7 (Nonholonomic reduction) | `sr_nonholonomic_reduced_rhs` |
| Bullo & Lewis (2004) | Ch. 4-5 (Geometric control) | `sr_ctrl_*`, `sr_controlled_lagrangian_*` |
