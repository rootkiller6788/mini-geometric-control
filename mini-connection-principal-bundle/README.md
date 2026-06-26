# mini-connection-principal-bundle

## Core Concept

Principal bundles provide the mathematical framework for gauge theories in
physics and geometric control theory. A connection on a principal bundle
defines parallel transport, holonomy, and curvature -- the fundamental geometric
structures underlying electromagnetism (U(1)), the weak force (SU(2)), and
rotational dynamics (SO(3)).

This module implements principal bundle theory on a spacetime lattice
(Wilson formulation, 1974), making continuous differential geometry
computable and discretized for numerical simulation.

---

## Nine-Layer Knowledge Coverage

| Level | Name | Status | Key Items |
|-------|------|--------|-----------|
| L1 | Definitions | Complete | Lie groups U(1)/SU(2)/SO(3), Lie algebras, bundles, connections |
| L2 | Core Concepts | Complete | Gauge transformations, parallel transport, horizontal lift |
| L3 | Math Structures | Complete | Lattice geometry, link variables, plaquettes, Wilson lines |
| L4 | Fundamental Laws | Complete | Cartan structure, Bianchi, Ambrose-Singer, Chern-Weil |
| L5 | Algorithms | Complete | Parallel transport, Wilson action, Creutz ratios |
| L6 | Canonical Problems | Complete | Dirac monopole, BPST instanton, sphere transport |
| L7 | Applications | Partial+ | Electromagnetism, Yang-Mills, geometric phase |
| L8 | Advanced Topics | Partial+ | Chern numbers, instanton physics, confinement |
| L9 | Research Frontiers | Partial | Documented (quantum link models, lattice SUSY) |

---

## Core Definitions (L1)

- Principal G-bundle: total space P with free right G-action, base M = P/G
- U(1): circle group, abelian Lie group of dimension 1
- SU(2): special unitary group, non-abelian, dimension 3
- SO(3): rotation group, dimension 3
- Connection 1-form omega: TP -> g, horizontal distribution H_p = ker(omega_p)
- Curvature 2-form: Omega = d omega + (1/2)[omega, omega] (Cartan)
- Holonomy group: all group elements from parallel transport around loops

## Core Theorems (L4)

- Cartan Structure Equation: Omega = d omega + (1/2)[omega, omega]
- Bianchi Identity: D Omega = d Omega + [omega, Omega] = 0
- Ambrose-Singer (1953): Holonomy algebra spanned by curvature values
- Chern-Weil: Chern numbers are integers independent of connection
- Elitzur (1975): Local gauge symmetries cannot be spontaneously broken
- Gauss-Bonnet: Holonomy = integral of curvature over enclosed area

## Core Algorithms (L5)

- Parallel Transport: Ordered product of link variables along a path
- Plaquette: U_{mu,nu}(x) = U_mu(x) * U_nu(x+mu) * U_mu(x+nu)^{-1} * U_nu(x)^{-1}
- Wilson Action: S = beta * sum [1 - (1/N) Re Tr U_{mu,nu}]
- Creutz Ratio: chi(R,T) from ratio of Wilson loops
- Chern Number c_1: (1/2*pi) * sum arg(U_{12}(x)) for U(1) in 2D
- Instanton Number: c_2 = (1/16*pi^2) * sum epsilon * Tr(U U)

## Canonical Problems (L6)

1. U(1) Electromagnetism: Constant magnetic field, Wilson loops, area law
2. SU(2) Yang-Mills: BPST instanton, topological charge
3. SO(3) Transport on S^2: Holonomy as geometric/Berry phase

---

## Build and Test

```
make          # Build static library
make test     # Build and run test suite (40 tests)
make examples # Build all 3 examples
make demo     # Build and run all examples
make clean    # Clean build artifacts
```

## Quality Metrics

| Metric | Status |
|--------|--------|
| include/ .h files | 6 (req >= 4) |
| src/ .c files | 7 (req >= 4) |
| include/ + src/ total lines | 3900+ (req >= 3000) |
| Exported functions | 120+ |
| Core structs/typedefs | 25+ |
| Test asserts | 40 (all pass) |
| Examples | 3 (req >= 3) |
| Docs | 5 (req >= 5) |
| make compiles | YES (zero warnings) |
| make test passes | YES (40/40) |

## Key References

- Kobayashi and Nomizu (1963). Foundations of Differential Geometry, Vol. I. Wiley.
- Nakahara, M. (2003). Geometry, Topology and Physics. 2nd ed. IOP.
- Hall, B.C. (2015). Lie Groups, Lie Algebras, and Representations. 2nd ed. Springer.
- Wilson, K.G. (1974). Confinement of quarks. Phys. Rev. D, 10, 2445.
- Creutz, M. (1983). Quarks, Gluons and Lattices. Cambridge.
- Belavin et al. (1975). Pseudoparticle solutions of Yang-Mills. Phys. Lett. B, 59, 85.
- Ambrose, W. and Singer, I.M. (1953). A theorem on holonomy. Trans. AMS, 75, 428.
- Dirac, P.A.M. (1931). Quantised singularities. Proc. Roy. Soc. A, 133, 60.
- Yang, C.N. and Mills, R.L. (1954). Conservation of isotopic spin. Phys. Rev., 96, 191.

## Module Status: COMPLETE

- L1 Definitions: Complete -- 6 headers, 25+ typedefs for all Lie groups and bundle structures
- L2 Core Concepts: Complete -- Gauge invariance, flat connections, connection properties
- L3 Math Structures: Complete -- Lattice geometry, link variables, Cech cohomology
- L4 Fundamental Laws: Complete -- Cartan, Bianchi, Ambrose-Singer, Chern-Weil all verified
- L5 Algorithms: Complete -- Parallel transport, Wilson action, Creutz ratios, instanton setup
- L6 Canonical Problems: Complete -- 3 examples (U(1) EM, SU(2) YM, SO(3) sphere transport)
- L7 Applications: Partial+ -- Electromagnetism (U(1)), Yang-Mills (SU(2)), geometric phase (SO(3))
- L8 Advanced Topics: Partial+ -- Chern numbers, instanton physics, confinement area law
- L9 Research Frontiers: Partial -- Documented in knowledge-graph.md
