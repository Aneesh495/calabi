package calabi.metric

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import calabi.algebra.MultiIndex

/**
 * MonomialBasis represents a basis {s_alpha} of global sections of the line bundle
 * O_X(k) over a projective variety in P^n.
 *
 * For degrees k < deg(P), the sections are given by all monomials of total degree k
 * in (n + 1) homogeneous variables: dimension N_k = binom(n + 1 + k - 1, k).
 */
class MonomialBasis(
    val ambientDimension: Int, // n in P^n
    val degree: Int // k
) {
    val numVariables: Int = ambientDimension + 1
    val multiIndices: List<MultiIndex> = MultiIndex.generateDegreeBasis(numVariables, degree)
    val size: Int = multiIndices.size // N_k

    init {
        require(ambientDimension >= 1) { "Ambient dimension must be at least 1" }
        require(degree >= 1) { "Section degree k must be at least 1" }
    }

    /**
     * Evaluates all N_k basis sections at a homogeneous point z:
     * s_alpha(z) = prod_{i=0}^n z_i^{alpha_i}.
     * Returns a flat complex vector of length N_k.
     */
    fun evaluateAll(z: HomogeneousPoint): List<Complex> {
        val coords = z.coordinates
        val values = ArrayList<Complex>(size)
        for (idx in multiIndices) {
            var prod = Complex.ONE
            for (i in 0 until numVariables) {
                val exp = idx[i]
                if (exp > 0) {
                    prod *= coords[i].pow(exp)
                }
            }
            values.add(prod)
        }
        return values
    }

    /**
     * Evaluates all sections into preallocated primitive double arrays (re and im)
     * for zero-allocation performance in Monte Carlo integration loops.
     */
    fun evaluateAllPrimitive(
        z: HomogeneousPoint,
        outRe: DoubleArray,
        outIm: DoubleArray
    ) {
        val coords = z.coordinates
        for (alpha in 0 until size) {
            val idx = multiIndices[alpha]
            var curRe = 1.0
            var curIm = 0.0

            for (i in 0 until numVariables) {
                val exp = idx[i]
                if (exp == 0) continue

                val zPow = coords[i].pow(exp)
                val newRe = curRe * zPow.re - curIm * zPow.im
                val newIm = curRe * zPow.im + curIm * zPow.re
                curRe = newRe
                curIm = newIm
            }

            outRe[alpha] = curRe
            outIm[alpha] = curIm
        }
    }
}
