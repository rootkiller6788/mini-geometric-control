# Course Alignment -- mini-lie-group-mechanics

## Nine-School Curriculum Mapping

| School | Course | Relevant Content |
|--------|--------|-----------------|
| **MIT** | 6.832 Underactuated Robotics | Lie groups for robot dynamics, SO(3)/SE(3) control |
| **MIT** | 16.323 Optimal Control | Geometric optimal control on Lie groups |
| **Stanford** | AA203 Optimal Control | Pontryagin on Lie groups, reduction |
| **Stanford** | CS237B Optimal Control | Geometric mechanics and control |
| **Berkeley** | EE222 Nonlinear Systems | Geometric nonlinear control |
| **Berkeley** | ME290 Geometric Control | Marsden's geometric control theory |
| **CMU** | 24-677 Nonlinear Control | Lie bracket analysis, controllability |
| **Princeton** | MAE 546 Optimal Control | Hamiltonian reduction, symplectic geometry |
| **Caltech** | CDS140 Nonlinear Dynamics | Euler-Poincare, rigid body, heavy top |
| **Caltech** | CDS270 Geometric Mechanics | Marsden & Ratiu textbook |
| **Cambridge** | 4F3 Nonlinear Control | Lie group methods in control |
| **Oxford** | C20 Adaptive Control | Geometric adaptive control |
| **ETH** | 227-0220 Model Reduction | Symmetry reduction, POD on Lie groups |

## Textbook Alignment

| Textbook | Chapters Covered |
|----------|-----------------|
| Marsden & Ratiu (1999) Intro to Mechanics and Symmetry | Ch.9 (Lie groups), Ch.11 (Momentum maps), Ch.13 (Euler-Poincare, Lie-Poisson), Ch.15 (Rigid bodies) |
| Murray, Li, Sastry (1994) Mathematical Intro to Robotic Manipulation | Appendix A (Lie groups), Ch.2 (Rigid body motion) |
| Hall (2015) Lie Groups, Lie Algebras, and Representations | Ch.3 (Matrix Lie groups), Ch.5 (Representations) |
| Hairer, Lubich, Wanner (2006) Geometric Numerical Integration | Ch.IV (Lie group methods), Ch.VIII (Variational integrators) |
| Holm, Schmah, Stoica (2009) Geometric Mechanics and Symmetry | Ch.2 (Lie groups), Ch.9 (Euler-Poincare), Ch.11 (Rigid bodies) |
| Bloch (2003) Nonholonomic Mechanics and Control | Ch.5 (Reduction), Ch.7 (Controllability) |
| Bullo & Lewis (2004) Geometric Control of Mechanical Systems | Ch.5 (Lie groups), Ch.7 (Controllability) |

## Keyword Mapping

| Concept | File Location |
|---------|--------------|
| Lie group | `include/lie_core.h` (LieGroupType enum) |
| Lie algebra | `include/lie_core.h` (LieAlgebraType enum) |
| Exponential map | `src/lie_core.c` (lie_exp_so3, lie_exp_se3) |
| Rodrigues formula | `src/lie_core.c` (lie_exp_so3) |
| Euler-Poincare | `src/lie_reduction.c` (lie_euler_poincare_rhs) |
| Rigid body | `src/lie_reduction.c` (lie_rigid_body_rhs) |
| Heavy top | `src/lie_reduction.c` (lie_heavy_top_rhs) |
| RKMK integrator | `src/lie_integration.c` (lie_rkmk_step_so3) |
| SE(3) quadrotor | `src/lie_mechanics.c` (lie_quadrotor_rhs) |
| Geometric controller | `src/lie_mechanics.c` (lie_quadrotor_geometric_controller) |

