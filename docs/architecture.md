# Calabi System Architecture

Calabi is a high-performance numerical engine for computing Ricci-flat Kähler metrics on Calabi-Yau manifolds, solving Donaldson's balanced metric fixed-point iteration, and extracting genus-0 Gromov-Witten invariants via Picard-Fuchs Frobenius series reversion.

```
+----------------------------------------------------------------------------------------------------+
|                                     CALABI REPL & CLI LAYER                                        |
|  (:sample, :donaldson, :sigma, :periods, :mirror, :instantons, :holonomy, :bench, :web)             |
+----------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+----------------------------------------------------------------------------------------------------+
|                                 CALABI-YAU GEOMETRY & MANIFOLDS                                    |
|   QuinticThreefold (P^4)  |  DworkPencil (X_psi)  |  CompleteIntersection (CICY)  |  K3Surface    |
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
|                                 PRIMITIVE LINEAR ALGEBRA ENGINE                                    |
|   ComplexMatrix (Interleaved DoubleArray)  |  Cholesky (L L^dagger)  |  Lanczos Spectral Solver    |
+----------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+----------------------------------------------------------------------------------------------------+
|                                      C SIMD ACCELERATION KERNELS                                   |
|             AVX2 / ARM NEON Vectorized Section Norms & Batch Outer-Product Accumulation             |
+----------------------------------------------------------------------------------------------------+
```

---

## 1. Modular Subsystems

### 1.1 calabi.algebra
Fundamental field arithmetic and projective geometry foundations:
* `Complex`: High-precision complex number arithmetic, polar coordinates, complex powers, and analytic branch cut functions.
* `Rational`: Exact rational arithmetic over Q with arbitrary-precision BigInteger numerators and denominators, canonical Euclidean GCD reduction, and Bernoulli / harmonic numbers.
* `ProjectiveSpace`: Homogeneous coordinates [z_0 : ... : z_n] in P^n(C), canonical Euclidean normalization, dominant affine chart selection maximizing coordinate magnitude, chart transition functions, and Fubini-Study metric calculations.
* `Polynomial`: Multivariate complex polynomials over C[z_0, ..., z_n], graded reverse lexicographic (grevlex) multi-index ordering, exact analytic gradients, Hessian matrices, and restriction to complex projective lines.

### 1.2 calabi.geometry
Exterior calculus and holomorphic volume forms:
* `DifferentialForm`: Differential (p, q)-forms on complex manifolds represented compactly via 32-bit bitmasks for holomorphic and anti-holomorphic differentials. Implements wedge products, signs of coordinate permutations, and complex conjugation.
* `PoincareResidue`: Evaluates the nowhere-vanishing holomorphic (n, 0)-volume form Omega via the algebraic Poincare residue theorem:
  $$\Omega = \operatorname{Res}_X \left( \frac{\mu}{P(z)} \right)$$
  Calculates intrinsic volume densities and Radon-Nikodym weights relative to the ambient Fubini-Study measure.

### 1.3 calabi.sampling
Numerical generation of points on complex algebraic varieties:
* `DurandKernerRoots`: Simultaneous polynomial root-finding algorithm with Aberth-Ehrlich initial circular estimates, quadratic convergence, and Newton-Raphson polishing steps.
* `LineIntersectionSampler`: Samples points uniformly on projective hypersurfaces by picking complex lines in the Grassmannian Gr(2, n+2) according to the Haar measure, restricting polynomials to the line, and extracting intersection points.
* `MonteCarloIntegrator`: Importance sampling numerical integrator equipped with parallel evaluation using Kotlin structured concurrency across multicore CPUs.

### 1.4 calabi.linalg
Zero-allocation high-throughput linear algebra:
* `ComplexMatrix`: Dense matrix over C backed by a contiguous flat primitive `DoubleArray` with interleaved real and imaginary coordinates `[re00, im00, re01, im01, ...]`. Avoids object overhead and enables direct cache locality.
* `Cholesky`: Complex Cholesky decomposition A = L * L^dagger for Hermitian positive-definite matrices with forward and backward triangular solvers and matrix inversion.
* `QRFactorization`: Modified Gram-Schmidt orthogonalization producing unitary factor Q (Q^dagger * Q = I) and upper-triangular factor R.
* `LanczosEigensolver`: Iterative Krylov subspace solver with implicit symmetric QR tridiagonal diagonalization for computing extreme eigenvalues of large Hermitian operators.

### 1.5 calabi.metric & calabi.curvature
Donaldson's balanced metric machinery:
* `MonomialBasis`: Basis of sections H^0(X, O(k)) comprising all monomials of degree k modulo the ideal.
* `SectionBundle`: Equipped with Hermitian metric h on H^0(X, O(k)), evaluating the Bergman / Fubini-Study Kahler potential:
  $$K_h(z) = \frac{1}{k} \ln \sum_{\alpha,\beta} h^{\alpha\beta} s_\alpha(z) \overline{s_\beta(z)}$$
* `DonaldsonTOperator`: Fixed-point operator mapping h to balanced metric T(h).
* `BalancedMetricSolver`: Iterative solver with Anderson acceleration and relative Frobenius convergence monitoring.
* `MetricTensor`: Exact analytic evaluation of the induced metric tensor g_{i bar{j}} = del_i dbar_j K on the manifold tangent bundle.
* `ChristoffelSymbols`: Complex connection symbols Gamma^i_{jk} = g^{i bar{l}} del_j g_{k bar{l}}.
* `DonaldsonSigmaMeasure`: Monge-Ampere equation error metric quantifying deviation from Ricci-flatness:
  $$\sigma = \frac{1}{\operatorname{Vol}(X)} \int_X \left| 1 - \frac{\det(g)}{\kappa \cdot |\Omega|^2} \right| d\mu_{CY}$$

### 1.6 calabi.mirror
Picard-Fuchs differential operators and instanton counting:
* `PicardFuchsEquation`: Fourth-order Picard-Fuchs differential operator:
  $$\mathcal{L} = \theta^4 - 5 z (5\theta + 1)(5\theta + 2)(5\theta + 3)(5\theta + 4)$$
* `FrobeniusMethod`: Computes exact series coefficients c_n and d_n for fundamental periods varpi_0(z) and varpi_1(z).
* `MirrorMap`: Flat coordinate t(z) and series reversion for z(q).
* `YukawaCoupling`: Topological 3-point coupling C_{ttt}(q).
* `InstantonExpansion`: Mobius inversion extracting exact Gromov-Witten curve counts:
  n_1 = 2875, n_2 = 609250, n_3 = 317206375, n_4 = 242467530000.

### 1.7 calabi.holonomy
Holonomy group reduction verification:
* `ParallelTransport`: Integrates tangent vector transport along closed loops on the manifold using 4th-order Runge-Kutta.
* `SU3Holonomy`: Frame transport over closed loops extracting the holonomy matrix H in U(3) and verifying det(H) = 1 (special unitary property).

---

## 2. Memory Layout & Zero Allocation Design

To achieve peak computational throughput during Monte Carlo loops with millions of section evaluations, Calabi stores matrices in flat primitive buffers:

```
Row 0: [ Re(0,0), Im(0,0), Re(0,1), Im(0,1), ..., Re(0,c-1), Im(0,c-1) ]
Row 1: [ Re(1,0), Im(1,0), Re(1,1), Im(1,1), ..., Re(1,c-1), Im(1,c-1) ]
...
```

This ensures:
1. Cache locality: Sequential elements along rows reside contiguously in L1/L2 data cache.
2. SIMD alignment: 128-bit and 256-bit vector units can load pairs of complex coordinates directly with zero memory packing overhead.
3. Garbage collection relief: No short-lived heap allocations occur inside hot integration inner loops.
