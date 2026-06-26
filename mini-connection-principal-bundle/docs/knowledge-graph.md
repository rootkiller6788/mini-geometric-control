# Knowledge Graph: Connection on Principal Bundles

## L1: Definitions
- U(1) group: e^{i*theta} with theta in [0, 2*pi), abelian Lie group
- SU(2) group: unit quaternions a+bi+cj+dk, non-abelian, dim 3
- SO(3) group: 3x3 orthogonal matrices with det=1, dim 3
- Lie algebra u(1): R with zero bracket (abelian)
- Lie algebra su(2): R^3 with cross product bracket
- Lie algebra so(3): R^3 with cross product bracket
- Principal fiber bundle: total space P with free right G-action, pi: P -> M = P/G
- Local trivialization: pi^{-1}(U) ~ U x G
- Transition functions: g_ij: U_i cap U_j -> G
- Section: sigma: U -> P with pi(sigma(x)) = x
- Associated vector bundle: E = P x_{rho} V
- Connection (Ehresmann): horizontal distribution H subset TP
- Connection 1-form: omega: TP -> g
- Curvature 2-form: Omega = d omega + (1/2)[omega, omega]
- Parallel transport: unique horizontal lift of a curve
- Holonomy group: all parallel transports around loops
- Gauge transformation: fiber-preserving bundle automorphism

## L2: Core Concepts
- Polynomial time reduction
- Gauge invariance: physical observables are gauge-invariant
- Wilson loop: gauge-invariant observable W(C) = Tr P exp(i oint A)
- Flat connection: Omega = 0, holonomy depends only on homotopy class
- Horizontal lift of a curve
- Covariant derivative: D = d + [A, .]
- Maurer-Cartan form: theta = g^{-1} dg
- Right-invariance of the connection
- Adjoint representation of G on g
- Gauge group: infinite-dimensional group of local transformations
- Lattice regularization: discretize spacetime as a lattice
- Link variable U_mu(x) in G on each lattice link
- Plaquette: elementary Wilson loop around a square

## L3: Mathematical Structures
- Lattice geometry: L[0..d-1] sizes, periodic boundary conditions
- Site indexing: linear index from multi-dimensional coordinates
- Link indexing: (site, mu) pairs
- Transition function cocycle: g_ij * g_jk * g_ki = id
- Cech cohomology H^1(M, G): classifies principal G-bundles
- Bundle reconstruction from transition functions
- Gauge potential A_mu(x) in g: continuum field
- Link variable U_mu(x) = P exp(-a * A_mu(x))
- Plaquette variable U_{mu,nu}(x)
- Wilson action S = beta * sum (1 - ReTr U_{mu,nu} / N)
- Gauge transformation on lattice: U_mu(x) -> g(x) * U_mu(x) * g(x+mu)^{-1}
- Path on lattice: sequence of (site, direction) steps
- Holonomy along path: ordered product of link variables

## L4: Fundamental Laws
- Cartan structure equation: Omega = d omega + (1/2)[omega, omega]
- Bianchi identity: D Omega = 0
- Ambrose-Singer theorem: holonomy algebra = curvature algebra
- Chern-Weil theory: characteristic classes from invariant polynomials
- First Chern class c_1: integer for U(1) bundles
- Second Chern class c_2: integer for SU(2) bundles (instanton number)
- Elitzur theorem: local gauge symmetry cannot break spontaneously
- Wilson criterion: confinement if area law for Wilson loops
- Gauss-Bonnet: holonomy = curvature integral for 2D surfaces
- Dirac quantization: e*g = n/2 for magnetic monopoles
- Bogomolny bound: M_monopole >= (4*pi/g) * v
- Atiyah-Singer index theorem: relates zero modes to topology

## L5: Algorithms/Methods
- Parallel transport computation (path-ordered product)
- Plaquette computation (elementary Wilson loop)
- Wilson action computation (sum over all plaquettes)
- Creutz ratio computation (force estimator)
- Polyakov loop computation (deconfinement order parameter)
- Chern number c_1 computation (U(1) 2D)
- Instanton number c_2 computation (SU(2) 4D)
- Bianchi identity verification (3-cube product)
- Dirac monopole configuration setup
- BPST instanton configuration setup
- Random gauge transformation generation
- Gauge fixing (not implemented, documented)
- Monte Carlo update (Metropolis, not implemented, documented)

## L6: Canonical Problems
- U(1) constant magnetic field: homogeneous F_{12} = B
- U(1) Dirac monopole: Wu-Yang construction on two patches
- SU(2) BPST instanton: topological charge Q=1 solution
- SO(3) parallel transport on S^2: holonomy = 2*pi*sin(phi)
- 2D compact QED: confinement by monopole condensation
- 4D SU(2) lattice gauge theory: confinement, glueball spectrum
- Wilson loop area law: static quark potential ~ sigma * R

## L7: Applications
- Electromagnetism as U(1) gauge theory on lattice
- Yang-Mills theory as SU(2) gauge theory on lattice
- Geometric phase in quantum mechanics (Berry phase)
- Parallel transport for orientation control (SO(3))
- Magnetic monopoles in grand unified theories
- Confinement in quantum chromodynamics (SU(3))
- Lattice QCD: numerical simulation of strong interactions

## L8: Advanced Topics
- Topological charge and axial anomaly
- Instanton physics and theta vacua
- Confinement mechanism: monopole condensation, center vortices
- Large N expansion ('t Hooft)
- Supersymmetric lattice gauge theory
- Chiral fermions on the lattice (Ginsparg-Wilson relation)
- AdS/CFT correspondence and Wilson loops

## L9: Research Frontiers
- Quantum link models: finite-dimensional Hilbert spaces for gauge fields
- Tensor network methods for lattice gauge theory
- Quantum simulation of lattice gauge theories
- Machine learning for lattice field theory (flow-based sampling)
- Deconfinement transition and heavy-ion collisions
- Lattice gravity and causal dynamical triangulations
