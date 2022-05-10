package calabi.mirror

import kotlin.math.roundToLong

/**
 * InstantonReport summarizes the extracted genus-0 Gromov-Witten invariants (instanton numbers):
 * - degree d
 * - count n_d: exact number of rational curves of degree d on the general quintic threefold
 */
data class InstantonReport(
    val degree: Int,
    val instantonCount: Long,
    val classicalName: String
)

/**
 * InstantonExpansion performs Mobius inversion on the topological Yukawa coupling C_{ttt}(q)
 * to extract the genus-0 Gromov-Witten invariants n_d of the quintic threefold:
 *
 *   C_{ttt}(q) = 5 + sum_{d=1}^infty n_d * d^3 * (q^d / (1 - q^d)) = 5 + sum_{k=1}^infty A_k * q^k
 *
 * where A_k = sum_{d | k} n_d * d^3.
 *
 * Reversing gives:
 *   n_k = (1 / k^3) * [ A_k - sum_{d | k, d < k} n_d * d^3 ]
 */
class InstantonExpansion(
    val yukawa: YukawaCoupling = YukawaCoupling()
) {
    /**
     * Computes the first maxDegree instanton numbers n_1, ..., n_{maxDegree}.
     */
    fun computeInstantons(maxDegree: Int = 4): List<InstantonReport> {
        val aCoeffs = yukawa.computeQExpansion(maxDegree)
        val nCounts = LongArray(maxDegree + 1)
        val reports = ArrayList<InstantonReport>(maxDegree)

        val classicalNames = listOf(
            "",
            "Lines on quintic threefold (Schubert 1879)",
            "Conics on quintic threefold (Katz 1986)",
            "Twisted cubics on quintic threefold (Ellingsrud-Stromme 1991)",
            "Rational quartics on quintic threefold (Candelas et al. 1991)",
            "Rational quintics on quintic threefold"
        )

        for (k in 1..maxDegree) {
            val aK = aCoeffs[k]
            var sumDivisors = 0.0

            for (d in 1 until k) {
                if (k % d == 0) {
                    val d3 = (d * d * d).toDouble()
                    sumDivisors += nCounts[d].toDouble() * d3
                }
            }

            val k3 = (k * k * k).toDouble()
            val nRaw = (aK - sumDivisors) / k3
            val nCount = nRaw.roundToLong()
            nCounts[k] = nCount

            val name = if (k < classicalNames.size) classicalNames[k] else "Degree $k rational curves"
            reports.add(InstantonReport(k, nCount, name))
        }

        return reports
    }

    companion object {
        /**
         * Known exact ground truth values from algebraic geometry:
         * n_1 = 2875, n_2 = 609250, n_3 = 317206375, n_4 = 242467530000.
         */
        val GROUND_TRUTH: Map<Int, Long> = mapOf(
            1 to 2875L,
            2 to 609250L,
            3 to 317206375L,
            4 to 242467530000L
        )
    }
}
