package calabi.sampling

import calabi.algebra.Complex
import kotlin.math.PI
import kotlin.math.max

/**
 * DurandKernerRoots implements the Weierstrass (Durand-Kerner) simultaneous root-finding algorithm
 * for arbitrary complex polynomials:
 *
 *   Q(t) = sum_{k=0}^d a_k * t^k
 *
 * All d roots are iteratively refined simultaneously with quadratic convergence.
 */
object DurandKernerRoots {

    /**
     * Solves Q(t) = 0 given coefficients in ascending degree order: [a_0, a_1, ..., a_d].
     */
    fun findRoots(
        coefficients: List<Complex>,
        maxIterations: Int = 100,
        tolerance: Double = 1e-12
    ): List<Complex> {
        val d = coefficients.size - 1
        require(d >= 1) { "Polynomial degree must be at least 1" }
        val lead = coefficients.last()
        require(!lead.isZero(1e-15)) { "Leading coefficient cannot be zero" }

        // Normalize polynomial so leading coefficient is 1: P(t) = t^d + sum_{k=0}^{d-1} c_k t^k
        val monic = coefficients.map { it / lead }

        if (d == 1) {
            // Linear case: t + c_0 = 0 => t = -c_0
            return listOf(-monic[0])
        }

        if (d == 2) {
            // Quadratic formula: t = (-c_1 +/- sqrt(c_1^2 - 4 c_0)) / 2
            val b = monic[1]
            val c = monic[0]
            val disc = (b * b - 4.0 * c).sqrt()
            val r1 = (-b + disc) / 2.0
            val r2 = (-b - disc) / 2.0
            return listOf(r1, r2)
        }

        // Aberth-Ehrlich initial estimate circle radius R = 1 + max |c_k|
        var maxCoeff = 0.0
        for (k in 0 until d) {
            maxCoeff = max(maxCoeff, monic[k].abs())
        }
        val radius = max(1.0, maxCoeff)
        // Center offset to break rotational symmetry
        val z0 = Complex(0.4, 0.9)

        // Initialize roots on a circle of radius R with angle offset
        val roots = Array(d) { k ->
            val angle = 2.0 * PI * k / d + PI / (2.0 * d)
            z0 + Complex.fromPolar(radius, angle)
        }

        // Durand-Kerner iteration
        for (iter in 0 until maxIterations) {
            var maxDelta = 0.0

            for (i in 0 until d) {
                val xi = roots[i]
                // Evaluate P(xi) using Horner's method
                var pVal = Complex.ONE
                for (k in (d - 1) downTo 0) {
                    pVal = pVal * xi + monic[k]
                }

                // Evaluate product prod_{j != i} (xi - xj)
                var prod = Complex.ONE
                for (j in 0 until d) {
                    if (j != i) {
                        prod *= (xi - roots[j])
                    }
                }

                val delta = if (prod.isZero(1e-16)) Complex(1e-5, 1e-5) else (pVal / prod)
                roots[i] = roots[i] - delta
                maxDelta = max(maxDelta, delta.abs())
            }

            if (maxDelta < tolerance) break
        }

        // Newton-Raphson polishing step on each root
        for (i in 0 until d) {
            var r = roots[i]
            for (step in 0 until 5) {
                var pVal = Complex.ONE
                var pDeriv = Complex.ofReal(d.toDouble())
                for (k in (d - 1) downTo 0) {
                    pVal = pVal * r + monic[k]
                    if (k > 0) {
                        pDeriv = pDeriv * r + (monic[k] * k.toDouble())
                    }
                }
                if (pDeriv.isZero(1e-15)) break
                val corr = pVal / pDeriv
                r -= corr
                if (corr.abs() < 1e-14) break
            }
            roots[i] = r
        }

        return roots.toList()
    }
}
