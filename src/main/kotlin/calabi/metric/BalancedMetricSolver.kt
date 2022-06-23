package calabi.metric

import calabi.linalg.ComplexMatrix
import calabi.sampling.SampledPoint
import kotlin.math.sqrt

/**
 * MetricConvergenceReport summarizes the fixed-point progress of the balanced metric solver:
 * - iteration count
 * - relative change ||T(h) - h||_F / ||h||_F
 * - trace normalization error |Tr(h^{-1} T(h)) - N_k|
 */
data class MetricConvergenceReport(
    val iterations: Int,
    val relativeChange: Double,
    val traceError: Double,
    val converged: Boolean,
    val finalMetric: ComplexMatrix
)

/**
 * BalancedMetricSolver orchestrates Donaldson's iteration:
 *
 *   h^{(m+1)} = T(h^{(m)})
 *
 * with optional Anderson acceleration and convergence diagnostics.
 */
class BalancedMetricSolver(
    val bundle: SectionBundle,
    val maxIterations: Int = 50,
    val tolerance: Double = 1e-4,
    val useAcceleration: Boolean = true
) {
    val tOperator = DonaldsonTOperator(bundle)
    val dim: Int = bundle.dim

    /**
     * Solves for the balanced Hermitian metric h^* satisfying T(h^*) = h^*.
     */
    fun solve(
        points: List<SampledPoint>,
        initialH: ComplexMatrix = ComplexMatrix.identity(dim),
        onIteration: ((Int, Double) -> Unit)? = null
    ): MetricConvergenceReport {
        var hCurrent = initialH.copy()
        val accelerator = if (useAcceleration) AndersonAcceleration(dim * dim * 2) else null

        var relChange = 1.0
        var iter = 0

        while (iter < maxIterations && relChange > tolerance) {
            val tVal = tOperator.evaluate(hCurrent, points)

            // Compute Frobenius norm difference ||T(h) - h||_F
            var diffSq = 0.0
            var hNormSq = 0.0
            for (i in hCurrent.data.indices) {
                val diff = tVal.data[i] - hCurrent.data[i]
                diffSq += diff * diff
                hNormSq += hCurrent.data[i] * hCurrent.data[i]
            }

            relChange = sqrt(diffSq) / sqrt(hNormSq)
            onIteration?.invoke(iter + 1, relChange)

            if (relChange <= tolerance) break

            // Next iterate via accelerator or damping
            var nextH: ComplexMatrix? = null
            if (accelerator != null) {
                val cand = accelerator.step(hCurrent, tVal).hermitianPart()
                if (calabi.linalg.Cholesky.factorizeOrNull(cand) != null) {
                    nextH = cand
                }
            }

            if (nextH == null) {
                // Damped fixed-point step: guaranteed positive-definite convex combination
                nextH = (hCurrent * 0.5 + tVal * 0.5).hermitianPart()
            }

            // Enforce trace normalization Tr(h) = N_k
            val trRe = nextH.trace().re
            if (trRe > 1e-12) {
                nextH = nextH * (dim.toDouble() / trRe)
            }

            hCurrent = nextH
            iter++
        }

        // Trace error check
        val cholesky = calabi.linalg.Cholesky.factorizeOrNull(hCurrent)
            ?: calabi.linalg.Cholesky.factorize(hCurrent.hermitianPart() + ComplexMatrix.identity(dim) * 1e-10)
        val hInv = cholesky.inverse()
        val hInvT = hInv * tOperator.evaluate(hCurrent, points)
        val tr = hInvT.trace()
        val traceError = kotlin.math.abs(tr.re - dim.toDouble())

        return MetricConvergenceReport(
            iterations = iter,
            relativeChange = relChange,
            traceError = traceError,
            converged = relChange <= tolerance,
            finalMetric = hCurrent
        )
    }
}

