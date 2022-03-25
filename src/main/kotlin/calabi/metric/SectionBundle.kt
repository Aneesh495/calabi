package calabi.metric

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import calabi.linalg.Cholesky
import calabi.linalg.ComplexMatrix

/**
 * SectionBundle represents the line bundle L^k = O_X(k) equipped with a basis of sections
 * and a Hermitian metric h on H^0(X, O_X(k)).
 *
 * Given h, the metric defines a local Kahler potential on the manifold:
 *   K_h(z) = (1/k) * ln ( sum_{alpha, beta} h^{alpha beta} s_alpha(z) bar{s_beta(z)} )
 *
 * where h^{alpha beta} is the inverse matrix (h^{-1})_{alpha beta}.
 */
class SectionBundle(
    val basis: MonomialBasis,
    var metricH: ComplexMatrix
) {
    val k: Int = basis.degree
    val dim: Int = basis.size

    init {
        require(metricH.rows == dim && metricH.cols == dim) {
            "Metric matrix dimensions (${metricH.rows} x ${metricH.cols}) must match section basis size ($dim x $dim)"
        }
    }

    /**
     * Computes the section norm squared:
     *   ||s(z)||_h^2 = sum_{alpha, beta} h^{alpha beta} s_alpha(z) bar{s_beta(z)}
     * using the inverse metric hInv.
     */
    fun sectionNormSquared(
        z: HomogeneousPoint,
        hInv: ComplexMatrix,
        tempRe: DoubleArray = DoubleArray(dim),
        tempIm: DoubleArray = DoubleArray(dim)
    ): Double {
        basis.evaluateAllPrimitive(z, tempRe, tempIm)

        var sum = 0.0
        for (alpha in 0 until dim) {
            val saRe = tempRe[alpha]
            val saIm = tempIm[alpha]
            if (saRe == 0.0 && saIm == 0.0) continue

            for (beta in 0 until dim) {
                val sbRe = tempRe[beta]
                val sbIm = -tempIm[beta] // bar{s_beta}

                val hInvRe = hInv.getRe(alpha, beta)
                val hInvIm = hInv.getIm(alpha, beta)

                val saSbRe = saRe * sbRe - saIm * sbIm
                val saSbIm = saRe * sbIm + saIm * sbRe

                sum += (saSbRe * hInvRe - saSbIm * hInvIm)
            }
        }
        return sum
    }

    /**
     * Evaluates the Fubini-Study / Bergman Kahler potential:
     *   K_h(z) = (1/k) * ln ||s(z)||_h^2
     */
    fun kahlerPotential(z: HomogeneousPoint, hInv: ComplexMatrix): Double {
        val normSq = sectionNormSquared(z, hInv)
        require(normSq > 0.0) { "Section norm squared must be positive: $normSq" }
        return (1.0 / k) * kotlin.math.ln(normSq)
    }

    companion object {
        fun identity(basis: MonomialBasis): SectionBundle =
            SectionBundle(basis, ComplexMatrix.identity(basis.size))
    }
}
