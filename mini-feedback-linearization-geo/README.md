# Feedback Linearization ¡ª Geometric Control Theory

## Module Status: COMPLETE

- **L1-L6**: Complete
- **L7**: Complete (3 applications: Robotics, Aerospace, Automotive)
- **L8**: Partial (2/10 advanced topics implemented)
- **L9**: Partial (documented, not implemented)

## Core Concept

Transform a nonlinear control-affine system into an equivalent linear system
via nonlinear state feedback and coordinate transformation (diffeomorphism).

```
System:  x' = f(x) + g(x) u,  y = h(x)

Feedback:  u = alpha(x) + beta(x) * v

Result:  y^{(rho)} = v  (chain of rho integrators)
```

## Key Definitions

| Term | Definition |
|------|-----------|
| Lie derivative | L_f h(x) = (grad h)^T . f(x) |
| Lie bracket | [f, g](x) = J_g f - J_f g |
| Relative degree rho | Smallest r such that L_g L_f^{r-1} h(x) != 0 |
| Diffeomorphism | Smooth invertible map with smooth inverse |
| Involutive distribution | [X, Y] in Delta for all X, Y in Delta |
| Zero dynamics | Internal dynamics when output y=0 |
| Normal form | Byrnes-Isidori: xi' = A*xi + B(...), eta' = q(xi, eta) |
| Minimum phase | Zero dynamics asymptotically stable |
| Decoupling matrix | E_{ij} = L_{g_j} L_f^{rho_i-1} h_i (MIMO) |

## Key Theorems

| Theorem | Statement |
|---------|----------|
| Frobenius | Delta completely integrable iff involutive |
| Jakubczyk-Respondek (1980) | IS-linearizable iff controllability + involutivity |
| Isidori (1995) | rho well-defined on open set U |
| Byrnes-Isidori (1991) | Internal stability iff zero dynamics stable |
| Chow (1939) | Controllability via Lie algebra rank condition |

## Key Algorithms

| Algorithm | Implementation |
|-----------|---------------|
| Relative degree computation | compute_relative_degree() |
| IO linearization feedback law | input_output_linearize() |
| IS linearization check | input_state_linearize() |
| Involutivity check | check_involutivity() |
| Normal form construction | compute_normal_form() |
| Zero dynamics analysis | analyze_zero_dynamics() |
| MIMO linearization | mimo_linearize() |
| Decoupling matrix | compute_decoupling_matrix() |
| Dynamic extension | free_output_feedback_linearize() |

## Classic Problems (Examples)

| Problem | File | Description |
|---------|------|-------------|
| Pendulum on cart | example_pendulum.c | 4th order, rho=3 |
| Ball and beam | example_ball_beam.c | 4th order, full rel deg |
| Robot manipulator | example_robot_arm.c | MIMO 2-link, computed torque |

## Applications

- **Robotics**: Computed torque control (ABB, KUKA, Fanuc)
- **Aerospace**: Quadrotor attitude control (DJI, Parrot, Boeing)
- **Automotive**: Engine control (Toyota, Bosch, Siemens)
- **Maglev**: Magnetic levitation (Transrapid)
- **Power**: Smart grid stabilization (ISO)
- **Precision**: Wafer stage control (ASML lithography)

## Nine-School Curriculum Mapping

| School | Course | Topics |
|--------|--------|--------|
| MIT | 6.832, 2.152 | Underactuated robotics, nonlinear control |
| Stanford | AA 203, ENGR 205 | Differential flatness, control design |
| Berkeley | ME 237, EE 222 | Geometric methods, passivity |
| CMU | 24-775, 16-745 | Robotic applications, optimization |
| Princeton | MAE 546 | Differential geometric methods |
| Caltech | CDS 240 | Geometric control theory |
| Cambridge | 4F12 | Nonlinear and predictive control |
| Oxford | Nonlinear Control | Lie algebraic methods |
| ETH Zurich | 227-0216-00L | Nonlinear control design |

## Build and Test

```bash
make          # Build static library
make test     # Run test suite
make examples # Build examples
make demo     # Run all examples
make clean    # Remove build artifacts
```

## File Structure

```
mini-feedback-linearization-geo/
  Makefile
  README.md
  include/
    feedback_linearization.h   (core API, data types)
    nonlinear_system.h         (vector, matrix, system)
    lie_derivative.h           (Lie derivative ops)
    lie_bracket.h              (Lie bracket ops)
    coordinate_transform.h     (diffeomorphisms)
    frobenius.h                (distributions, involutivity)
  src/
    feedback_linearization.c   (relative degree, IO/IS/MIMO linearization)
    nonlinear_system.c         (vector, matrix, system operations)
    lie_derivative.c           (Lie derivative computation)
    lie_bracket.c              (Lie bracket, ad operator)
    coordinate_transform.c     (diffeomorphism, system transform)
    frobenius.c                (involutivity, Frobenius, distributions)
  tests/
    test_feedback_linearization.c
  examples/
    example_pendulum.c
    example_ball_beam.c
    example_robot_arm.c
  demos/
    demo_feedback_linearization.c
  benches/
    bench_core.c
  docs/
    knowledge-graph.md
    coverage-report.md
    gap-report.md
    course-alignment.md
    course-tree.md
```

## References

- Isidori, A. (1995). *Nonlinear Control Systems*, 3rd ed. Springer.
- Nijmeijer, H. & van der Schaft, A. (1990). *Nonlinear Dynamical Control Systems*. Springer.
- Khalil, H.K. (2002). *Nonlinear Systems*, 3rd ed. Prentice Hall.
- Slotine, J-J. & Li, W. (1991). *Applied Nonlinear Control*. Prentice Hall.
- Brockett, R.W. (1978). Feedback invariants for nonlinear systems. *IFAC Congress*.
- Jakubczyk, B. & Respondek, W. (1980). On linearization of control systems. *Bull. Acad. Pol. Sci.*
- Byrnes, C.I. & Isidori, A. (1991). Asymptotic stabilization of minimum phase nonlinear systems. *IEEE TAC*.
