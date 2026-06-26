# Knowledge Graph: Feedback Linearization (Geometric Control)

## L1: Definitions
- Feedback linearization: Transform nonlinear system to linear via state feedback + diffeomorphism
- Relative degree (rho): Smallest r where L_g L_f^{r-1} h(x) != 0
- Lie derivative: L_f h(x) = (grad h)^T . f(x), directional derivative along vector field
- Lie bracket: [f,g](x) = J_g f - J_f g, commutator of vector fields
- Diffeomorphism: Smooth invertible map with smooth inverse
- Involutive distribution: [X,Y] in Delta for all X,Y in Delta
- Zero dynamics: Internal dynamics when output is constrained to zero
- Normal form (Byrnes-Isidori): xi prime = A xi + B(alpha + beta u), eta prime = q(xi, eta)
- Minimum phase: Zero dynamics asymptotically stable
- Vector relative degree: Tuple (rho_1, ..., rho_m) for MIMO systems
- Decoupling matrix: E_{ij}(x) = L_{g_j} L_f^{rho_i-1} h_i(x)
- Ad operator: ad_f^k g = [f, ad_f^{k-1} g], iterated Lie bracket
- Nonlinear control-affine system: x prime = f(x) + sum g_i(x) u_i
- Distribution: Smooth assignment of subspace to each point
- Complete integrability: Existence of foliation tangent to distribution

## L2: Core Concepts
- Exact linearization vs approximate linearization
- Input-output linearization vs input-state linearization
- Internal dynamics and their stability
- Byrnes-Isidori normal form
- Dynamic extension for singular decoupling matrix
- Feedback invariance and Brunovsky normal form
- Controllability rank condition for nonlinear systems
- Involutivity condition for exact linearization
- Frobenius integrability and foliations
- MIMO decoupling and non-interacting control
- Partial feedback linearization
- Global vs local results

## L3: Mathematical Structures
- Smooth manifolds and tangent bundles
- Vector fields as derivations on C-infinity functions
- Lie algebra of vector fields
- Differential forms and exterior derivative
- Distributions and codistributions (Pfaffian systems)
- Frobenius theorem (differential geometric version)
- Coordinate charts and transition maps
- Push-forward and pull-back of vector fields
- Jacobian matrices and their non-singularity
- Brunovsky blocks (chain of integrators)
- State-space diffeomorphism group

## L4: Fundamental Theorems
- Frobenius Theorem: Delta integrable iff involutive
- Jakubczyk-Respondek Theorem: IS-linearizable iff controllability + involutivity
- Isidori relative degree theorem: rho well-defined on open set
- Byrnes-Isidori zero dynamics stability theorem
- Chow Theorem: Controllability from Lie algebra rank condition
- Sussmann-Jurdjevic accessibility theorem
- Brockett necessary condition for stabilization
- Descusse-Moog dynamic extension algorithm
- Singh inversion algorithm for nonlinear systems
- Brunovsky normal form for controllable linear systems

## L5: Algorithms/Methods
- Relative degree computation algorithm
- Input-output linearization feedback law
- Input-state linearization coordinate construction
- Involutivity checking via rank test on augmented matrices
- Decoupling matrix computation
- Normal form transformation
- Zero dynamics computation and stability analysis
- Dynamic extension algorithm (Descusse-Moog)
- Linearizing coordinate construction
- Flatness-based trajectory generation

## L6: Canonical Problems
- Pendulum on a cart (4th order, rho=3)
- Ball and beam system (4th order, full relative degree)
- Two-link robot manipulator (MIMO, computed torque)
- Magnetic levitation system (2nd order, rho=1)
- Quadrotor attitude control (MIMO, 6-DOF)
- Inverted pendulum (2nd order)
- VTOL aircraft (MIMO, non-minimum phase)
- Flexible joint robot (4th order per joint)
- Biochemical reactor (relative degree issues)
- DC motor with nonlinear load

## L7: Applications
- Robotics: Computed torque control (ABB, KUKA, Fanuc industrial robots)
- Aerospace: Quadrotor attitude control (DJI, Parrot, Boeing)
- Automotive: Engine air-fuel ratio control (Toyota, Bosch, Siemens)
- Process control: pH neutralization via nonlinear gain scheduling
- Power systems: Excitation control for smart grid stabilization (ISO)
- Maglev: Magnetic levitation train control (Transrapid)
- Precision motion: Wafer stage control in lithography (ASML)
- Biomedical: Insulin-glucose regulation via nonlinear control

## L8: Advanced Topics
- Time-varying feedback linearization (extended state)
- Adaptive feedback linearization (unknown parameters)
- Robust feedback linearization (uncertainty, ISS framework)
- Discrete-time feedback linearization (sampled-data systems)
- Stochastic feedback linearization (Fokker-Planck methods)
- Orbital feedback linearization (limit cycle stabilization)
- Feedback linearization with input constraints
- Differential flatness and equivalence to linear systems
- Lyapunov-based redesign for robust linearization
- Dynamic feedback linearization via integrator chains

## L9: Research Frontiers
- Contraction theory and incremental stability
- Koopman operator methods for nonlinear control
- Learning-based feedback linearization (Gaussian processes)
- Feedback linearization on Lie groups (SE(3) for aerospace)
- Quantum feedback linearization
- Hybrid/switched feedback linearization
- Distributed feedback linearization for multi-agent systems
- Homogeneity-based approximate linearization
- Meta-complexity of nonlinear control verification
