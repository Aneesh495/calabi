package calabi.mirror

import calabi.algebra.Complex
import calabi.algebra.Rational
import java.math.BigInteger

/**
 * FrobeniusMethod computes the exact series expansions of the fundamental period solutions
 * to the Picard-Fuchs equation around the Large Complex Structure Limit (z = 0):
 *
 *   varpi_0(z) = sum_{n=0}^infty c_n * z^n
 *   varpi_1(z) = varpi_0(z) * ln(z) + sum_{n=1}^infty d_n * z^n
 *
 * where:
 *   c_n = (5n)! / (n!)^5
 *   d_n = 5 * c_n * (H_{5n} - H_n)
 */
class FrobeniusMethod(
    val maxOrder: Int = 12
) {
    val cCoeffs = ArrayList<BigInteger>(maxOrder + 1)
    val dCoeffs = ArrayList<Rational>(maxOrder + 1)

    init {
        // Precompute exact combinatorial coefficients
        var cCur = BigInteger.ONE
        cCoeffs.add(cCur)
        dCoeffs.add(Rational.ZERO)

        for (n in 1..maxOrder) {
            // c_n = c_{n-1} * prod_{k=0}^4 (5n - k) / n^5
            var num = BigInteger.ONE
            for (k in 0..4) {
                num = num.multiply(BigInteger.valueOf((5 * n - k).toLong()))
            }
            val den = BigInteger.valueOf(n.toLong()).pow(5)
            cCur = cCur.multiply(num).divide(den)
            cCoeffs.add(cCur)

            val h5n = Rational.harmonic(5 * n)
            val hn = Rational.harmonic(n)
            val hDiff = h5n - hn
            val dCur = Rational.of(cCur, BigInteger.ONE) * hDiff * 5
            dCoeffs.add(dCur)
        }
    }

    /**
     * Evaluates the fundamental period varpi_0(z) up to maxOrder.
     */
    fun evaluateVarpi0(z: Complex): Complex {
        var sum = Complex.ZERO
        var zPow = Complex.ONE
        for (n in 0..maxOrder) {
            val cn = cCoeffs[n].toDouble()
            sum += zPow * cn
            zPow *= z
        }
        return sum
    }

    /**
     * Evaluates the regular logarithmic part S_1(z) = sum_{n=1}^infty d_n z^n.
     */
    fun evaluateS1(z: Complex): Complex {
        var sum = Complex.ZERO
        var zPow = z
        for (n in 1..maxOrder) {
            val dn = dCoeffs[n].toDouble()
            sum += zPow * dn
            zPow *= z
        }
        return sum
    }

    /**
     * Evaluates the full second period varpi_1(z) = varpi_0(z) * ln(z) + S_1(z).
     */
    fun evaluateVarpi1(z: Complex): Complex {
        val varpi0 = evaluateVarpi0(z)
        val s1 = evaluateS1(z)
        return varpi0 * z.ln() + s1
    }
}
