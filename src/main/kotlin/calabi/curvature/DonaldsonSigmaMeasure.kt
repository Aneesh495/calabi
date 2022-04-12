package calabi.curvature

import calabi.algebra.Complex
import calabi.geometry.PoincareResidue
import calabi.linalg.Cholesky
import calabi.linalg.ComplexMatrix
import calabi.sampling.SampledPoint
import kotlin.math.abs
import kotlin.math.max

/**
 * SigmaMeasureReport encapsulates the Monge-Ampere error metrics:
 * - sigma: L1 error measure int_X |1 - det(g) / (kappa * |Omega|^2)| dmu_CY / Vol(X)
 * - sigmaMax: maximum pointwise deviation
 * - kappa: volume normalization constant
 */
data class SigmaMeasureReport(
    val sigma: Double,
    val sigmaMax: Double,
    val kappa: Double,
    val sampledPointCount: Int
)

/**
 * DonaldsonSigmaMeasure evaluates the Donaldson sigma-measure error:
 *
 *   sigma = (1 / Vol(X)) * int_X | 1 - (det(g_{i bar{j}}) / (kappa * ||Omega||^2)) | dmu_{CY}
 *
 * This measures the exact deviation of the numerical Kahler metric from Ricci-flatness.
 * By Donaldson's theorem, sigma -> 0 as degree k -> infinity.
 */
class DonaldsonSigmaMeasure(
    val metricTensor: MetricTensor,
    val poincareResidue: PoincareResidue
) {
    /**
     * Evaluates the sigma-measure error on a set of Monte Carlo sampled points.
     */
    fun evaluate(
        h: ComplexMatrix,
        points: List<SampledPoint>
    ): SigmaMeasureReport {
        require(points.isNotEmpty()) { "Point set cannot be empty" }
        val cholesky = Cholesky.factorize(h)
        val hInv = cholesky.inverse()

        val detGVals = DoubleArray(points.size)
        val omegaSqVals = DoubleArray(points.size)
        val weights = DoubleArray(points.size)

        var sumDetG = 0.0
        var sumOmegaSq = 0.0
        var totalWeight = 0.0

        for (idx in points.indices) {
            val sp = points[idx]
            val pt = sp.point
            val w = sp.weight
            val affine = pt.toAffineCoordinates()

            val (elimIdx, _) = poincareResidue.selectEliminatedCoordinate(affine)
            val gMat = metricTensor.evaluate(affine, hInv, elimIdx)
            val detG = metricTensor.determinant(gMat)
            val omegaSq = poincareResidue.volumeDensity(affine)

            detGVals[idx] = detG
            omegaSqVals[idx] = omegaSq
            weights[idx] = w

            sumDetG += detG * w
            sumOmegaSq += omegaSq * w
            totalWeight += w
        }

        require(sumOmegaSq > 0.0 && totalWeight > 0.0) { "Zero volume encountered in sigma measure" }
        val kappa = sumDetG / sumOmegaSq

        var sumDeviation = 0.0
        var maxDev = 0.0

        for (idx in points.indices) {
            val detG = detGVals[idx]
            val omegaSq = omegaSqVals[idx]
            val w = weights[idx]

            val ratio = detG / (kappa * omegaSq)
            val dev = abs(1.0 - ratio)

            sumDeviation += dev * w
            maxDev = max(maxDev, dev)
        }

        val sigma = sumDeviation / totalWeight
        return SigmaMeasureReport(
            sigma = sigma,
            sigmaMax = maxDev,
            kappa = kappa,
            sampledPointCount = points.size
        )
    }
}
