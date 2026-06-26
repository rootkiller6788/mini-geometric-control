# Knowledge Graph ¡ª mini-differential-flatness

## L1: Definitions
- Differentially flat system (FlatControlAffineSystem)
- Flat output (FlatPolynomial)
- Lie derivative (lie_derivative)
- Lie bracket (lie_bracket)
- Relative degree (flat_compute_relative_degree)
- Defect (flat_compute_defect)
- Distribution (FlatDistribution)
- Endogenous dynamic feedback (FlatDynamicCompensator)

## L2: Core Concepts
- Exact linearization via dynamic feedback
- Trajectory generation without ODE integration
- Flat output parameterization of state and input
- Brunovsky normal form (integrator chains)

## L3: Mathematical Structures
- Multivariate polynomials (FlatPolynomial, FlatMonomial)
- Vector fields (FlatVectorField)
- Distributions and codistributions
- Differential forms (FlatOneForm)
- Lie algebras (structure constants, nilpotency, solvability)
- B-splines (FlatSpline)

## L4: Fundamental Laws
- Fliess equivalence theorem (flat_system_is_flat)
- Brunovsky theorem (flat_to_brunovsky)
- Frobenius theorem (flat_distribution_is_involutive)
- Jacobi identity (flat_verify_jacobi_identity)
- Chow-Rashevsky theorem (flat_check_accessibility_rank)
- Controllability Gramian (flat_controllability_gramian)

## L5: Algorithms/Methods
- B-spline trajectory generation (Cox-de Boor)
- Minimum-snap/jerk trajectory optimization
- Dynamic extension algorithm
- Polynomial symbolic algebra engine
- QR/SVD/Cholesky/LDLT decompositions
- Newton/BFGS/Golden section optimization
- Cubic spline interpolation
- Pole placement and LQR design
- Distribution sequence algorithm (Murray-Rathinam-Sluis)

## L6: Canonical Problems
- n-Trailer vehicle (flat output: last trailer position)
- Quadrotor UAV (flat outputs: x, y, z, yaw)
- Overhead crane (flat output: hook position)
- Space robot (flat output: end-effector pose)
- Inverted pendulum, PVTOL, ball-beam, maglev, AUV, satellite
- CSTR chemical reactor

## L7: Applications
- UAV flight control (minimum snap, perching, flip maneuvers)
- Autonomous vehicle parking (trailer, lane change)
- Crane anti-sway material handling
- NASA space robotics on-orbit servicing
- Chemical process control (CSTR temperature)
- Precision agriculture (field coverage waypoints)
- Rocket landing (SpaceX-style)
- Nuclear debris removal (Fukushima context)

## L8: Advanced Topics
- Orbital flatness (trajectory-dependent)
- Flat singularity detection and resolution
- Monte Carlo flatness verification
- Time-varying and switched flat systems
- Structural flatness from symmetries
- Bisimulation of flat systems
- Gap metric robustness analysis

## L9: Research Frontiers
- Machine learning-based flat output discovery (documented)
- Koopman operator theory for flatness (documented)
- Neural network flat output parameterization (documented)
- Hierarchical and compositional flatness (documented)
