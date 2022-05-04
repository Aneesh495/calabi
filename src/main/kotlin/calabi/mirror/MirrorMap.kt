package calabi.mirror

import calabi.algebra.Complex
import kotlin.math.PI

/**
 * MirrorMap computes the flat coordinate t(z) on the complexified Kahler moduli space
 * and the exponentiated mirror map q(z) = exp(2 * pi * i * t(z)), along with its series inversion z(q).
 */
class MirrorMap(
    val frobenius: FrobeniusMethod = FrobeniusMethod(maxOrder = 10)
) {
    /**
     * Evaluates t(z) = (1 / 2*pi*i) * (varpi_1(z) / varpi_0(z)).
     */
    fun evaluateT(z: Complex): Complex {
        val varpi0 = frobenius.evaluateVarpi0(z)
        val varpi1 = frobenius.evaluateVarpi1(z)
        val twoPiI = Complex(0.0, 2.0 * PI)
        return (varpi1 / varpi0) / twoPiI
    }

    /**
     * Evaluates q(z) = exp(2 * pi * i * t(z)) = z * exp(S_1(z) / varpi_0(z)).
     */
    fun evaluateQ(z: Complex): Complex {
        val varpi0 = frobenius.evaluateVarpi0(z)
        val s1 = frobenius.evaluateS1(z)
        val exponent = s1 / varpi0
        return z * exponent.exp()
    }

    /**
     * Computes the Taylor series coefficients of q(z) = z + sum_{k=2}^N a_k z^k.
     */
    fun computeQSeries(order: Int = 8): DoubleArray {
        // q(z) = z * exp( S_1(z) / varpi_0(z) )
        // Compute S_1(z) / varpi_0(z) series
        val c = DoubleArray(order + 1) { i -> frobenius.cCoeffs[i].toDouble() }
        val d = DoubleArray(order + 1) { i -> frobenius.dCoeffs[i].toDouble() }

        // Division of formal power series: R(z) = S_1(z) / varpi_0(z) = sum_{k=1} r_k z^k
        val r = DoubleArray(order + 1)
        for (n in 1..order) {
            var sum = d[n]
            for (k in 1 until n) {
                sum -= r[k] * c[n - k]
            }
            r[n] = sum / c[0]
        }

        // exp(R(z)): E(z) = 1 + sum_{k=1} e_k z^k
        val e = DoubleArray(order + 1)
        e[0] = 1.0
        for (n in 1..order) {
            var sum = 0.0
            for (k in 1..n) {
                sum += k.toDouble() * r[k] * e[n - k]
            }
            e[n] = sum / n.toDouble()
        }

        // q(z) = z * E(z) => q_k = e_{k-1} for k >= 1
        val qSeries = DoubleArray(order + 1)
        qSeries[0] = 0.0
        for (k in 1..order) {
            qSeries[k] = e[k - 1]
        }
        return qSeries
    }

    /**
     * Inverts the power series q(z) = z + sum_{k=2}^N q_k z^k to obtain
     * z(q) = q + sum_{k=2}^N z_k q^k via Lagrange inversion / formal series reversion.
     */
    fun computeInverseZSeries(order: Int = 8): DoubleArray {
        val qCoeffs = computeQSeries(order)
        val zCoeffs = DoubleArray(order + 1)
        zCoeffs[0] = 0.0
        zCoeffs[1] = 1.0 // Coefficient of q is 1

        for (n in 2..order) {
            // Reversion formula
            var sum = 0.0
            // Substitute z(q) into q(z(q)) = q
            // q(z) = z + q_2 z^2 + q_3 z^3 + ...
            val zPowers = Array(n + 1) { DoubleArray(order + 1) }
            zPowers[0][0] = 1.0
            for (p in 1..n) {
                for (i in p..order) {
                    var s = 0.0
                    for (j in 1..(i - p + 1)) {
                        s += zCoeffs[j] * zPowers[p - 1][i - j]
                    }
                    zPowers[p][i] = s
                }
            }

            for (k in 2..n) {
                sum += qCoeffs[k] * zPowers[k][n]
            }
            zCoeffs[n] = -sum
        }

        return zCoeffs
    }
}
