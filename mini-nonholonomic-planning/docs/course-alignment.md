# Course Alignment — Nonholonomic Motion Planning

## Nine-School Curriculum Mapping

| University | Course | Topics Covered | Our Implementation |
|------------|--------|---------------|-------------------|
| **MIT** | 6.832 Underactuated Robotics | Nonholonomic constraints, Lie brackets, motion planning | `noplan_lie_bracket()`, RRT, sinusoidal steering |
| **MIT** | 2.12 Intro to Robotics | Wheeled mobile robots, car-like kinematics | `noplan_model_unicycle()`, `noplan_model_car()` |
| **Stanford** | AA274 Multi-agent Control | Nonholonomic vehicle models, formation control | Trailer models, lattice planning |
| **Berkeley** | EE222 Nonlinear Systems | Frobenius theorem, Chow-Rashevskii, controllability | `noplan_frobenius_theorem_test()`, LARC |
| **CMU** | 24-677 Nonlinear Control | Geometric control, Brockett's condition | `noplan_brockett_necessary()` |
| **Princeton** | MAE 546 Optimal Control | Trajectory optimization under constraints | `noplan_trajectory_optimize()` |
| **Caltech** | CDS 243 Geometric Mechanics | SE(2)/SE(3), Lie groups, geometric phase | `SE2`/`SE3` implementation |
| **Cambridge** | 4F3 Nonlinear & Predictive Ctrl | Nonholonomic path planning, RRT/PRM | Full RRT + PRM implementations |
| **Oxford** | B4 Predictive Control | Mobile robot planning, parking maneuvers | Parallel/garage parking benchmarks |
| **ETH** | 151-0563 Dynamic Programming | Sampling-based planning, motion primitives | Lattice planner, PRM construction |

## Reference Textbooks

| Textbook | Authors | Coverage |
|----------|---------|----------|
| A Mathematical Introduction to Robotic Manipulation | Murray, Li, Sastry (1994) | Chapters 6-7: Nonholonomic systems, chained form |
| Robot Motion Planning and Control | Laumond (ed.) (1998) | Complete treatment of nonholonomic planning |
| Nonholonomic Mechanics and Control | Bloch (2003) | Geometric foundations, Frobenius, Chow |
| Geometric Control of Mechanical Systems | Bullo & Lewis (2005) | Distribution theory, controllability |
| Planning Algorithms | LaValle (2006) | RRT, PRM, sampling-based planning |

## Course Coverage Matrix

| Topic | MIT | Stan | Berk | CMU | Prin | Calt | Camb | Oxf | ETH |
|-------|-----|------|------|-----|------|------|------|-----|-----|
| Nonholonomic constraints | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Frobenius theorem | ✓ | | ✓ | ✓ | | ✓ | ✓ | | |
| Chow-Rashevskii | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | | |
| Brockett's condition | | ✓ | ✓ | ✓ | | | ✓ | | |
| Lie bracket motion | ✓ | ✓ | | | | ✓ | | | |
| RRT/PRM planning | ✓ | ✓ | | ✓ | | | | ✓ | ✓ |
| Sinusoidal steering | ✓ | | | | | | | | |
| Nilpotent approximation | | | | | | | | | |
| Dubins/Reeds-Shepp paths | | ✓ | | ✓ | | | | | |
| Trailer systems | | | | | | | | | |
