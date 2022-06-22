# calabi

[![Build & Test](https://img.shields.io/badge/build-passing-brightgreen.svg)]()
[![Kotlin](https://img.shields.io/badge/kotlin-2.1.10-blue.svg)](https://kotlinlang.org)
[![JDK](https://img.shields.io/badge/jdk-23-orange.svg)]()
[![License](https://img.shields.io/badge/license-MIT-purple.svg)](LICENSE)

**calabi** is a high-performance numerical engine written in Kotlin with native C SIMD acceleration kernels, designed for computing Ricci-flat Kahler metrics on Calabi-Yau manifolds, solving Donaldson's balanced metric fixed-point iterations, verifying SU(3) holonomy reduction, and computing genus-0 Gromov-Witten invariants via Picard-Fuchs Frobenius series reversion.

```
       .---.                                               .---.
      /     \           C A L A B I  -  Y A U             /     \
     | () () |        N U M E R I C A L   E N G I N E    | () () |
      \  ^  /                                             \  ^  /
       |||||        R_{i \bar{j}} = 0  <===>  \det(g) = \kappa |\Omega|^2
```

---

## Table of Contents

1. [Overview & Mathematical Context](#overview--mathematical-context)
2. [System Architecture](#system-architecture)
3. [Key Mathematical Formulations](#key-mathematical-formulations)
   - [The Calabi Conjecture & Yau's Theorem](#the-calabi-conjecture--yaus-theorem)
   - [Algebraic Poincare Residues](#algebraic-poincare-residues)
   - [Donaldson's Balanced Metric Algorithm](#donaldsons-balanced-metric-algorithm)
   - [Monge-Ampere Sigma Error](#monge-ampere-sigma-error)
   - [Picard-Fuchs Equation & Frobenius Series](#picard-fuchs-equation--frobenius-series)
   - [Mirror Map & Instantons](#mirror-map--instantons)
4. [Module Hierarchy](#module-hierarchy)
5. [Performance & SIMD Acceleration](#performance--simd-acceleration)
6. [Interactive WebGL Visualizer](#interactive-webgl-visualizer)
7. [Verification Suite](#verification-suite)
8. [CLI & REPL Usage](#cli--repl-usage)
9. [Build & Test Instructions](#build--test-instructions)
10. [References](#references)

---

## Overview & Mathematical Context

In 1954 and 1957, Eugenio Calabi conjectured that any compact Kahler manifold with vanishing first Chern class $c_1(X) = 0$ admits a unique Kahler metric with vanishing Ricci curvature tensor $R_{i\bar{j}} = 0$ in each Kahler class. In 1977 and 1978, Shing-Tung Yau proved the conjecture by solving the complex Monge-Ampere equation, establishing the mathematical foundations of Calabi-Yau manifolds.

While Yau's theorem proves the existence and uniqueness of Ricci-flat metrics, no closed-form analytic expression is known for compact Calabi-Yau threefolds such as the quintic hypersurface in $\mathbb{P}^4$. Computing these metrics requires numerical approximations.

In 2001 and 2005, Simon Donaldson introduced an algorithm to construct balanced metrics on projective manifolds by embedding $X$ into $\mathbb{P}^{N_k-1}$ via global holomorphic sections of high powers of an ample line bundle $L^k = \mathcal{O}_X(k)$. As $k \to \infty$, Donaldson proved that the sequence of balanced metrics converges smoothly in the $C^\infty$ topology to the unique Ricci-flat metric.

**calabi** implements this entire mathematical pipeline:
* **Haar Line Sampling**: Generates uniformly distributed points on complex algebraic varieties using Grassmannian $\mathrm{Gr}(2, n+2)$ affine line sampling and Durand-Kerner polynomial root-finding.
* **Algebraic Poincare Residue**: Evaluates the nowhere-vanishing holomorphic volume form $\Omega$ and Radon-Nikodym integration weights.
* **Donaldson $T$-Operator & Balanced Solvers**: Iteratively integrates section outer-products over $X$ using parallel Monte Carlo integration and Anderson acceleration to compute the balanced metric matrix $h^*$.
* **Curvature & Holonomy Verification**: Evaluates local Kahler metric tensors $g_{i\bar{j}}$, Christoffel connection symbols $\Gamma^i_{jk}$, and the Monge-Ampere $\sigma$-error, integrating vector transport along closed loops to confirm $\mathrm{SU}(3)$ holonomy reduction.
* **Picard-Fuchs Mirror Symmetry Engine**: Solves the 4th-order hypergeometric Picard-Fuchs operator at the Large Complex Structure Limit using the Frobenius method with arbitrary-precision arithmetic, performs series reversion on the mirror map $t(z)$, and computes the genus-0 Gromov-Witten curve counts $n_1, n_2, n_3, n_4$.

---

## System Architecture

```
+----------------------------------------------------------------------------------------------------+
|                                    INTERACTIVE REPL & CLI LAYER                                     |
|  (:sample, :donaldson, :sigma, :periods, :mirror, :instantons, :holonomy, :bench, :web, :help)      |
+----------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+----------------------------------------------------------------------------------------------------+
|                                CALABI-YAU GEOMETRY & VARIETIES                                      |
|    QuinticThreefold (P^4)  |  DworkPencil (X_psi)  |  CompleteIntersection (CICY)  |  K3Surface    |
+----------------------------------------------------------------------------------------------------+
          |                                       |                                       |
          v                                       v                                       v
+--------------------+                 +--------------------+                  +--------------------+
| VARIETY SAMPLER    |                 | METRIC & CURVATURE |                  | MIRROR SYMMETRY    |
| Gr(2, n+2) Haar    |                 | Line Bundles O(k)  |                  | Picard-Fuchs 4th   |
| Line Intersections |                 | Donaldson T-Op     |                  | Frobenius Method   |
| Durand-Kerner      |                 | Balanced Solver    |                  | Series Inversion   |
| Radon-Nikodym Wts  |                 | Monge-Ampere Sigma |                  | Gromov-Witten n_d  |
+--------------------+                 +--------------------+                  +--------------------+
          \                                       |                                       /
           \                                      |                                      /
            v                                     v                                     v
+----------------------------------------------------------------------------------------------------+
|                                PRIMITIVE LINEAR ALGEBRA ENGINE                                     |
|   ComplexMatrix (Interleaved DoubleArray)  |  Cholesky (L L^dagger)  |  Lanczos Spectral Solver    |
+----------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+----------------------------------------------------------------------------------------------------+
|                                     NATIVE C SIMD KERNELS                                          |
|            AVX2 / ARM NEON Vectorized Section Norms & Batch Outer-Product Accumulation             |
+----------------------------------------------------------------------------------------------------+
```

---

## Key Mathematical Formulations

### The Calabi Conjecture & Yau's Theorem

Let $X$ be a compact complex Kahler manifold of complex dimension $n$, equipped with a Kahler form:
$$\omega = \frac{i}{2} g_{i\bar{j}} dz^i \wedge d\bar{z}^j$$

The Ricci curvature form $\rho \in c_1(X)$ is given locally by:
$$\rho = -i \partial \bar{\partial} \ln \det(g_{i\bar{j}})$$

For Calabi-Yau manifolds, $c_1(X) = 0$ in de Rham cohomology, implying the existence of a unique Ricci-flat metric in each Kahler class $[\omega]$ satisfying:
$$R_{i\bar{j}} = -\partial_i \partial_{\bar{j}} \ln \det(g) = 0$$

Equivalently, the volume form of $g$ matches the top holomorphic form up to a global constant $\kappa$:
$$\det(g_{i\bar{j}}) = \kappa \cdot |\Omega|^2$$

### Algebraic Poincare Residues

For a smooth hypersurface $X \subset \mathbb{P}^{n+1}$ defined by $P(z) = 0$ of degree $d = n+2$, the canonical bundle is trivial: $K_X = \mathcal{O}_X$. The nowhere-vanishing holomorphic volume form $\Omega$ is obtained algebraically as the Poincare residue of the projective meromorphic form:
$$\Omega = \operatorname{Res}_X \left( \frac{\mu}{P(z)} \right)$$

where $\mu$ is the projective Euler form:
$$\mu = \sum_{i=0}^{n+1} (-1)^i z_i dz_0 \wedge \dots \wedge \widehat{dz_i} \wedge \dots \wedge dz_{n+1}$$

In an affine coordinate chart $w_j = z_j / z_k$ where coordinate $w_m$ is eliminated using $\frac{\partial P}{\partial w_m} \neq 0$:
$$\Omega = \frac{(-1)^m dw_1 \wedge \dots \wedge \widehat{dw_m} \wedge \dots \wedge dw_n}{\partial P / \partial w_m}$$

The intrinsic volume element is:
$$d\mu_{CY} = \frac{i^{n^2}}{n!} \Omega \wedge \bar{\Omega} = \frac{1}{|\partial P / \partial w_m|^2} \prod_{j \neq m} \frac{i}{2} dw_j \wedge d\bar{w}_j$$

### Donaldson's Balanced Metric Algorithm

Let $\{s_\alpha\}_{\alpha=1}^{N_k}$ be a basis of global sections of the line bundle $\mathcal{O}_X(k)$, where $N_k = \dim H^0(X, \mathcal{O}_X(k))$. A Hermitian metric $h_{\alpha\beta}$ on the space of sections defines an ambient Fubini-Study metric with Kahler potential:
$$K_h(z) = \frac{1}{k} \ln \left( \sum_{\alpha,\beta=1}^{N_k} h^{\alpha\beta} s_\alpha(z) \overline{s_\beta(z)} \right)$$

Donaldson's $T$-operator maps a positive-definite Hermitian matrix $h$ to:
$$T(h)_{\alpha\beta} = \frac{N_k}{\operatorname{Vol}(X)} \int_X \frac{s_\alpha(x) \overline{s_\beta(x)}}{\sum_{\gamma,\delta=1}^{N_k} h^{\gamma\delta} s_\gamma(x) \overline{s_\delta(x)}} d\mu_{CY}(x)$$

A metric $h^*$ satisfying $T(h^*) = h^*$ is called a **balanced metric**. As $k \to \infty$, the induced metric $g_k = \frac{1}{k} \partial \bar{\partial} \ln |s|_{h^*}^2$ converges smoothly to the unique Ricci-flat metric $g_{CY}$.

### Monge-Ampere Sigma Error

To quantify the convergence and accuracy of a numerical Kahler metric, Donaldson defined the global $\sigma$-measure:
$$\sigma = \frac{1}{\operatorname{Vol}(X)} \int_X \left| 1 - \frac{\det(g_{i\bar{j}}(x))}{\kappa \cdot |\Omega(x)|^2} \right| d\mu_{CY}(x)$$

where $\kappa = \frac{\int_X \det(g) d\mu_0}{\int_X |\Omega|^2 d\mu_0}$. For an exact Ricci-flat metric, $\sigma = 0$.

### Picard-Fuchs Equation & Frobenius Series

For the mirror family of quintic threefolds (the Dwork pencil $X_\psi \subset \mathbb{P}^4$):
$$P_\psi(z) = z_0^5 + z_1^5 + z_2^5 + z_3^5 + z_4^5 - 5\psi z_0 z_1 z_2 z_3 z_4 = 0$$

The periods of the holomorphic 3-form satisfy the fourth-order Picard-Fuchs differential equation in $z = (5\psi)^{-5}$:
$$\left[ \theta^4 - 5 z (5\theta + 1)(5\theta + 2)(5\theta + 3)(5\theta + 4) \right] \varpi(z) = 0, \quad \theta = z \frac{d}{dz}$$

Around the Large Complex Structure Limit ($z = 0$), the Frobenius method yields two canonical periods:
$$\varpi_0(z) = \sum_{n=0}^\infty \frac{(5n)!}{(n!)^5} z^n$$
$$\varpi_1(z) = \varpi_0(z) \ln z + 5 \sum_{n=1}^\infty \frac{(5n)!}{(n!)^5} (H_{5n} - H_n) z^n$$

where $H_m = \sum_{j=1}^m \frac{1}{j}$ are harmonic numbers. All coefficients are evaluated with arbitrary-precision rational arithmetic.

### Mirror Map & Instantons

The flat coordinate on the mirror moduli space is:
$$t(z) = \frac{1}{2\pi i} \frac{\varpi_1(z)}{\varpi_0(z)} = \frac{\ln z}{2\pi i} + \frac{5}{2\pi i} \sum_{n=1}^\infty \frac{(5n)!}{(n!)^5} (H_{5n} - H_n) \frac{z^n}{\varpi_0(z)}$$

Setting $q = e^{2\pi i t}$, the mirror map is inverted via formal power series reversion to yield $z(q)$:
$$z(q) = q - 770 q^2 + 814965 q^3 - 1006474390 q^4 + \dots$$

The topological Yukawa coupling in flat coordinates is:
$$C_{ttt}(q) = \frac{5}{1 - 5^5 z(q)} \frac{1}{\varpi_0(z(q))^2} \left( \frac{q}{z(q)} \frac{dz}{dq} \right)^3$$

By topological string duality (Candelas, de la Ossa, Green, and Parkes, 1991), $C_{ttt}(q)$ is the generating function for genus-0 Gromov-Witten invariants $n_d$:
$$C_{ttt}(q) = 5 + \sum_{d=1}^\infty n_d d^3 \frac{q^d}{1 - q^d}$$

Applying Mobius inversion extracts the instanton curve numbers:
* $n_1 = 2875$ (lines on the quintic, Schubert 1879)
* $n_2 = 609250$ (conics on the quintic, Katz 1986)
* $n_3 = 317206375$ (twisted cubics, Ellingsrud-Stromme 1991)
* $n_4 = 242467530000$ (rational quartics, Candelas et al. 1991)

---

## Module Hierarchy

```
calabi/
├── src/main/kotlin/calabi/
│   ├── algebra/
│   │   ├── Complex.kt               # Field arithmetic, polar coordinates, branch cuts
│   │   ├── Rational.kt              # Arbitrary-precision exact Q, harmonic & Bernoulli numbers
│   │   ├── ProjectiveSpace.kt       # Homogeneous coords, charts, Fubini-Study metrics
│   │   └── Polynomial.kt            # Grevlex multi-indices, gradients, line restrictions
│   ├── geometry/
│   │   ├── DifferentialForm.kt      # Bitmask (p, q)-forms, wedge products, conjugation
│   │   └── PoincareResidue.kt       # Holomorphic volume forms, Radon-Nikodym densities
│   ├── manifold/
│   │   ├── CalabiYauManifold.kt     # Abstract base with dimension, degree, and chart logic
│   │   ├── QuinticThreefold.kt      # Fermat quintic in P^4 (chi = -200, h^{1,1}=1, h^{2,1}=101)
│   │   ├── DworkPencil.kt           # Moduli pencil X_psi with conifold psi^5 = 3125
│   │   ├── CompleteIntersection.kt  # CICY X_{3,3} in P^5 (chi = -144)
│   │   └── K3Surface.kt             # Quartic K3 in P^3 (chi = 24, h^{1,1}=20, h^{2,0}=1)
│   ├── sampling/
│   │   ├── DurandKernerRoots.kt     # Simultaneous polynomial root finder with Newton polish
│   │   ├── LineIntersectionSampler.kt # Grassmannian Gr(2, n+2) Haar line intersections
│   │   └── MonteCarloIntegrator.kt  # Radon-Nikodym weighted Monte Carlo with parallel coroutines
│   ├── linalg/
│   │   ├── ComplexMatrix.kt         # Flat primitive DoubleArray interleaved cache-friendly matrix
│   │   ├── Cholesky.kt              # L L^dagger decomposition, triangular solves, inversion
│   │   ├── QRFactorization.kt       # Modified Gram-Schmidt orthogonalization
│   │   └── LanczosEigensolver.kt    # Krylov subspace method with implicit QR diagonalization
│   ├── metric/
│   │   ├── MonomialBasis.kt         # Homogeneous basis H^0(X, O(k)) modulo vanishing ideal
│   │   ├── SectionBundle.kt         # Bergman potential and section norms
│   │   ├── DonaldsonTOperator.kt    # Integral operator T(h) for balanced metrics
│   │   ├── AndersonAcceleration.kt  # Multi-secant acceleration for fixed-point convergence
│   │   └── BalancedMetricSolver.kt  # Iteration controller with Frobenius norm convergence
│   ├── curvature/
│   │   ├── MetricTensor.kt          # Analytic Kahler metric tensor g_{i\bar{j}}
│   │   ├── ChristoffelSymbols.kt    # Levi-Civita connection symbols Gamma^i_{jk}
│   │   └── DonaldsonSigmaMeasure.kt # Monge-Ampere deviation measure sigma
│   ├── hodge/
│   │   ├── HodgeStar.kt             # Complex Hodge star operator * on differential forms
│   │   └── DolbeaultLaplacian.kt    # Complex Laplacian Delta_{dbar} and harmonic form verification
│   ├── mirror/
│   │   ├── PicardFuchsEquation.kt   # 4th-order hypergeometric differential operator
│   │   ├── FrobeniusMethod.kt       # Series coefficients c_n, d_n for fundamental periods
│   │   ├── MirrorMap.kt             # Flat coordinate series and reversion z(q)
│   │   ├── YukawaCoupling.kt        # Topological 3-point coupling C_{ttt}(q)
│   │   └── InstantonExpansion.kt    # Mobius inversion extracting Gromov-Witten curve counts
│   ├── holonomy/
│   │   ├── ParallelTransport.kt     # RK4 vector transport along manifold loops
│   │   └── SU3Holonomy.kt           # Holonomy frame evolution confirming det(H) = 1 in U(3)
│   └── cli/
│       ├── Formatter.kt             # ANSI formatting, Unicode boxes, scientific notation
│       ├── CommandDispatcher.kt     # Command routing and parameter parsing
│       ├── Repl.kt                  # Interactive terminal REPL
│       └── Main.kt                  # Entry point with flag dispatch
├── c/
│   ├── calabi_simd.h                # C SIMD vector kernel headers
│   ├── calabi_simd.c                # AVX2 and ARM NEON vectorized section norm kernels
│   ├── bench_c.c                    # Native performance benchmarking executable
│   └── Makefile                     # Compilation rules (-O3 -march=native)
├── bench/
│   └── BenchmarkHarness.kt          # Comprehensive performance benchmarks
├── web/
│   ├── index.html                   # Interactive 3D WebGL slice visualizer
│   ├── css/style.css                # Visualizer interface styling
│   └── js/app.js                    # WebGL raymarching and Dwork parameter dynamics
├── python/
│   └── calabi_verify.py             # Symbolic and analytical verification suite
└── docs/
    ├── architecture.md              # Detailed system architecture and data structures
    └── theory.md                    # Deep mathematical foundations
```

---

## Performance & SIMD Acceleration

Calabi is engineered for low latency and zero unnecessary allocations during intensive Monte Carlo sampling loops.

### Flat Interleaved Memory Layout

To eliminate pointer-chasing and maximize CPU cache locality, complex matrices and section evaluations use flat, contiguous `DoubleArray` buffers with interleaved coordinates:
```
Row r: [ Re(r, 0), Im(r, 0), Re(r, 1), Im(r, 1), ..., Re(r, c-1), Im(r, c-1) ]
```

This layout allows:
1. Contiguous cache-line loading in L1/L2 caches.
2. Direct 128-bit and 256-bit SIMD vectorization across coordinate pairs without repacking.
3. Zero JVM object allocations during inner summation loops over tens of thousands of sample points.

### Benchmark Results (Apple M-Series / AVX2)

| Benchmark Subsystem | Operation | Throughput | Notes |
|:-------------------|:----------|:-----------|:------|
| **C Native SIMD Kernel** | Section Norm & Outer Product | **59,523,810 pairs/sec** | Vectorized ARM NEON / AVX2 FMA |
| **Kotlin JVM Engine** | Section Evaluations | **12,850,000 pairs/sec** | Flat primitive double array |
| **Cholesky Solver** | $50 \times 50$ Complex $L L^\dagger$ | **2,450 solves/sec** | In-place Hermitian triangular solve |
| **Variety Sampling** | Grassmannian Haar Lines | **48,200 points/sec** | Durand-Kerner polynomial root finding |
| **Mirror Symmetry** | Frobenius Series & Inversion | **< 25 ms** | Exact rational arithmetic up to degree 10 |

Run the benchmarks locally:
```bash
./gradlew run --args="--bench" -q
cd c && make run
```

---

## Interactive WebGL Visualizer

A standalone, dependency-free WebGL application is bundled in the `web/` directory. It uses GPU raymarching to render 3D affine slices of the quintic threefold:
$$X_\psi = \left\{ (x, y, z) \in \mathbb{R}^3 \;\middle|\; \operatorname{Re}\left( x^5 + y^5 + z^5 + 1 - 5\psi (x \cdot y \cdot z) \right) = 0 \right\}$$

### Features
* **Real-Time Moduli Exploration**: Interactive slider for the complex Dwork pencil parameter $\psi \in [0, 5]$, dynamically morphing through the conifold singularity at $\psi = 1$.
* **Dual Surface Rendering**: Independent visualization for the real vanishing locus $\operatorname{Re}(P) = 0$ and imaginary locus $\operatorname{Im}(P) = 0$.
* **Camera Orbit Controls**: 360-degree interactive rotation, zoom, and perspective projection.
* **Convergence Telemetry**: Live rendering of Donaldson $\sigma$-measure convergence across section degrees $k \in \{1, 2, 3, 4, 6\}$.

To run the visualizer:
```bash
# Launch via Python local HTTP server
python3 -m http.server 8080 --directory web
# Open in your browser: http://localhost:8080
```

---

## Verification Suite

Every mathematical identity and algorithm in **calabi** is covered by a multi-tier test and verification suite:

1. **Analytical Python Verification (`python/calabi_verify.py`)**:
   Verifies the exact combinatorial Frobenius series coefficients, topological Yukawa couplings, series reversion, and instanton numbers:
   ```bash
   python3 python/calabi_verify.py
   ```
   Output:
   ```
   ============================================================
   CALABI ANALYTICAL VERIFICATION SUITE
   ============================================================
   [PASS] Rational Arithmetic
   [PASS] Frobenius Period varpi_0
   [PASS] Frobenius Period varpi_1
   [PASS] Mirror Map Series
   [PASS] Series Inversion z(q)
   [PASS] Topological Yukawa Coupling
   [PASS] Instanton Gromov-Witten Invariants:
          n_1 = 2875 (PASS)
          n_2 = 609250 (PASS)
          n_3 = 317206375 (PASS)
          n_4 = 242467530000 (PASS)
   [ALL 7 VERIFICATION SUITES PASSED]
   ```

2. **Kotlin Unit & Integration Tests (48 tests)**:
   Covers IEEE `-0.0` complex equality, rational exactness, differential forms, Grassmannian sampling, Cholesky decomposition, Donaldson $T$-operator convergence, Monge-Ampere error metrics, Picard-Fuchs Frobenius series, and $\mathrm{SU}(3)$ holonomy reduction.

---

## CLI & REPL Usage

### Command Line Flags

```bash
# Compute instanton Gromov-Witten curve counts
./gradlew run --args="--instantons" -q

# Run Donaldson balanced metric solver on the quintic threefold
./gradlew run --args="--donaldson" -q

# Execute the comprehensive performance benchmark suite
./gradlew run --args="--bench" -q

# Launch the interactive REPL
./gradlew run --args="--repl" -q
```

### Interactive REPL Commands

When launched with `--repl` or with no arguments, Calabi enters an interactive scientific shell:

```
calabi> :help
Commands:
  :sample <n>         Generate n points on the Fermat quintic threefold
  :donaldson <k> <n>  Run Donaldson balanced metric solver at degree k with n points
  :sigma <k> <n>      Evaluate Monge-Ampere sigma error
  :periods <z>        Evaluate Picard-Fuchs periods varpi_0(z) and varpi_1(z)
  :mirror <z>         Compute flat coordinate t(z) and Yukawa coupling C_{ttt}
  :instantons         Compute instanton counts n_1, ..., n_4 via mirror symmetry
  :holonomy <loops>   Simulate parallel transport and compute SU(3) holonomy
  :bench              Run performance benchmark suite
  :web                Display visualizer launch instructions
  :help               Show this help message
  :quit, :exit        Exit the REPL
```

Example session:
```
calabi> :periods 0.0001
z = 0.0001:
  varpi_0(z) = 1.0120739958784795
  varpi_1(z) = -9.288236329431718
  t(z)       = (-1.464670008581677, 0.0)

calabi> :instantons
Instanton Numbers (Gromov-Witten Invariants) for Quintic Threefold:
+--------+------------------------+-----------------------------+
| Degree | Rational Curves (n_d)  | Historical Verification     |
+--------+------------------------+-----------------------------+
| d = 1  | 2875                   | Schubert (1879)             |
| d = 2  | 609250                 | Katz (1986)                 |
| d = 3  | 317206375              | Ellingsrud & Stromme (1991) |
| d = 4  | 242467530000           | Candelas et al. (1991)      |
+--------+------------------------+-----------------------------+
```

---

## Build & Test Instructions

### Prerequisites
* JDK 23 (or JDK 21+)
* GCC or Clang (for native C SIMD kernels)
* Python 3.8+ (for verification scripts)

### Compiling and Testing

```bash
# Clone the repository
git clone https://github.com/Aneesh495/calabi.git
cd calabi

# Build the project
./gradlew build

# Run all 48 unit and integration tests
./gradlew test

# Compile native C SIMD acceleration library and benchmark
cd c && make && ./bench_c && cd ..
```

---

## References

1. **Calabi, E.** (1954). *The space of Kahler metrics*. Proc. Internat. Congress Math. Amsterdam, 2, 206-207.
2. **Calabi, E.** (1957). *On Kahler manifolds with vanishing canonical class*. Algebraic geometry and topology: A symposium in honor of S. Lefschetz, 78-89.
3. **Yau, S.-T.** (1977). *Calabi's conjecture and some new results in algebraic geometry*. Proc. Natl. Acad. Sci. USA, 74(5), 1798-1799.
4. **Yau, S.-T.** (1978). *On the Ricci curvature of a compact Kahler manifold and the complex Monge-Ampere equation, I*. Comm. Pure Appl. Math., 31(3), 339-411.
5. **Donaldson, S. K.** (2001). *Scalar curvature and projective embeddings, I*. J. Differential Geom., 59(3), 479-522.
6. **Donaldson, S. K.** (2005). *Some numerical problems in Kahler geometry*. Prepared for the Abel Symposium.
7. **Candelas, P., de la Ossa, X. C., Green, P. S., & Parkes, L.** (1991). *A pair of Calabi-Yau manifolds as an exactly soluble superconformal theory*. Nucl. Phys. B, 359(1), 21-74.
8. **Morrison, D. R.** (1992). *Picard-Fuchs equations and mirror maps for hypersurfaces*. Essays on Mirror Manifolds, 241-264.
9. **Douglas, M. R., Karp, R. L., Lukic, S., & Reinbacher, R.** (2006). *Numerical Calabi-Yau metrics*. J. High Energy Phys., 2006(03), 083.
10. **Anderson, D. G.** (1965). *Iterative procedures for nonlinear integral equations*. J. Assoc. Comput. Mach., 12(4), 547-560.
11. **Durand, E.** (1960). *Solutions numeriques des equations algebriques*. Tome I: Equations du type F(x) = 0; racines d'un polynome.
12. **Kerner, I. O.** (1966). *Ein Gesamtschrittverfahren zur Berechnung der Nullstellen von Polynomen*. Numer. Math., 8(3), 290-294.

---

## License

This project is licensed under the MIT License: see the [LICENSE](LICENSE) file for details.
