package calabi.metric

import calabi.algebra.Complex
import calabi.linalg.Cholesky
import calabi.linalg.ComplexMatrix
import calabi.sampling.SampledPoint

/**
 * DonaldsonTOperator evaluates the Donaldson T-operator T(h) mapping a Hermitian metric
 * h on H^0(X, O(k)) to a new balanced Hermitian metric:
 *
 *   T(h)_{alpha beta} = (N_k / Vol(X)) * int_X (s_alpha bar{s_beta} / ||s||_h^2) dmu_{CY}
 *
 * approximated via Monte Carlo summation over weighted points on the variety:
 *
 *   T(h)_{alpha beta} = (N_k / sum w_i) * sum_{i=1}^M ( (s_alpha(x_i) bar{s_beta(x_i)} / ||s(x_i)||_h^2) * w_i )
 */
class DonaldsonTOperator(
    val bundle: SectionBundle
) {
    val dim: Int = bundle.dim // N_k
    val k: Int = bundle.k

    /**
     * Computes T(h) for a given Hermitian metric matrix h.
     */
    fun evaluate(h: ComplexMatrix, points: List<SampledPoint>): ComplexMatrix {
        require(points.isNotEmpty()) { "Point set cannot be empty for Donaldson T-operator" }
        val cholesky = Cholesky.factorize(h)
        val hInv = cholesky.inverse()

        val tMat = ComplexMatrix.zeros(dim, dim)
        var totalWeight = 0.0

        val tempRe = DoubleArray(dim)
        val tempIm = DoubleArray(dim)

        for (sp in points) {
            val pt = sp.point
            val w = sp.weight
            totalWeight += w

            // Evaluate all sections at this point
            bundle.basis.evaluateAllPrimitive(pt, tempRe, tempIm)

            // Compute ||s(pt)||_h^2 = sum_{gamma, delta} hInv^{gamma, delta} s_gamma bar{s_delta}
            var sNormSq = 0.0
            for (gamma in 0 until dim) {
                val sgRe = tempRe[gamma]
                val sgIm = tempIm[gamma]
                if (sgRe == 0.0 && sgIm == 0.0) continue

                for (delta in 0 until dim) {
                    val sdRe = tempRe[delta]
                    val sdIm = -tempIm[delta]

                    val invRe = hInv.getRe(gamma, delta)
                    val invIm = hInv.getIm(gamma, delta)

                    val prodRe = sgRe * sdRe - sgIm * sdIm
                    val prodIm = sgRe * sdIm + sgIm * sdRe

                    sNormSq += (prodRe * invRe - prodIm * invIm)
                }
            }

            if (sNormSq <= 1e-15) continue
            val factor = w / sNormSq

            // Accumulate into T: s_alpha bar{s_beta} * factor
            for (alpha in 0 until dim) {
                val saRe = tempRe[alpha]
                val saIm = tempIm[alpha]
                if (saRe == 0.0 && saIm == 0.0) continue

                for (beta in 0 until dim) {
                    val sbRe = tempRe[beta]
                    val sbIm = -tempIm[beta] // bar{s_beta}

                    val pRe = (saRe * sbRe - saIm * sbIm) * factor
                    val pIm = (saRe * sbIm + saIm * sbRe) * factor

                    val curRe = tMat.getRe(alpha, beta)
                    val curIm = tMat.getIm(alpha, beta)
                    tMat.set(alpha, beta, curRe + pRe, curIm + pIm)
                }
            }
        }

        require(totalWeight > 0.0) { "Total integration weight must be positive: $totalWeight" }

        // Normalization scale: N_k / totalWeight
        val scale = dim.toDouble() / totalWeight
        val result = tMat * scale

        // Enforce exact Hermiticity
        for (i in 0 until dim) {
            result.set(i, i, result.getRe(i, i), 0.0)
            for (j in (i + 1) until dim) {
                val avgRe = (result.getRe(i, j) + result.getRe(j, i)) * 0.5
                val avgIm = (result.getIm(i, j) - result.getIm(j, i)) * 0.5
                result.set(i, j, avgRe, avgIm)
                result.set(j, i, avgRe, -avgIm)
            }
        }

        return result
    }
}
