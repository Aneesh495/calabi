package calabi.mirror

import calabi.algebra.Complex

/**
 * YukawaCoupling computes the topological 3-point Yukawa coupling C_{ttt} on the mirror quintic
 * as a function of the flat coordinate q = exp(2 * pi * i * t).
 *
 *   C_{ttt}(q) = (5 / (1 - 5^5 * z(q))) * ( (1 / varpi_0(z(q))) * (q / z(q)) * (dz / dq) )^3
 *
 * The genus-0 topological A-model expansion gives:
 *   C_{ttt}(q) = 5 + sum_{d=1}^infty n_d * d^3 * (q^d / (1 - q^d))
 */
class YukawaCoupling(
    val mirrorMap: MirrorMap = MirrorMap()
) {
    /**
     * Evaluates C_{ttt}(q) for a specific complex value of q.
     */
    fun evaluate(q: Complex, zInvSeries: DoubleArray): Complex {
        // Evaluate z(q)
        var z = Complex.ZERO
        var qPow = Complex.ONE
        for (k in zInvSeries.indices) {
            z += qPow * zInvSeries[k]
            qPow *= q
        }

        // dz/dq series
        var dzDq = Complex.ZERO
        var qPowDer = Complex.ONE
        for (k in 1 until zInvSeries.size) {
            dzDq += qPowDer * (zInvSeries[k] * k.toDouble())
            qPowDer *= q
        }

        val varpi0 = mirrorMap.frobenius.evaluateVarpi0(z)
        val denom1 = Complex.ONE - (z * 3125.0)

        val ratio = (Complex.ONE / varpi0) * (q / z) * dzDq
        val ratio3 = ratio.pow(3)

        return (Complex.ofReal(5.0) / denom1) * ratio3
    }

    /**
     * Computes the Taylor expansion of C_{ttt}(q) = 5 + sum_{k=1}^M A_k * q^k.
     */
    fun computeQExpansion(order: Int = 6): DoubleArray {
        val zCoeffs = mirrorMap.computeInverseZSeries(order + 4)
        val cCoeffs = DoubleArray(order + 4) { i -> mirrorMap.frobenius.cCoeffs[i].toDouble() }

        // Expand z(q), dz/dq, varpi_0(z(q)), and (1 - 3125 z(q))^{-1}
        // dz/dq = sum_{k=0} (k+1) z_{k+1} q^k
        val dzDq = DoubleArray(order + 1)
        for (k in 0..order) {
            dzDq[k] = (k + 1).toDouble() * zCoeffs[k + 1]
        }

        // q / z(q): since z(q) = q (1 + z_2 q + z_3 q^2 + ...), q / z(q) = 1 / (1 + sum_{k=1} z_{k+1} q^k)
        val zNorm = DoubleArray(order + 1)
        for (k in 0..order) {
            zNorm[k] = zCoeffs[k + 1] // zNorm[0] = z_1 = 1.0
        }
        val invZNorm = invertSeries(zNorm, order)

        // varpi_0(z(q)) series
        val varpi0Z = substituteSeries(cCoeffs, zCoeffs, order)
        val invVarpi0 = invertSeries(varpi0Z, order)

        // Factor F = (1 / varpi_0) * (q / z) * (dz / dq)
        val prod1 = multiplySeries(invVarpi0, invZNorm, order)
        val fSeries = multiplySeries(prod1, dzDq, order)

        // F^3
        val f2 = multiplySeries(fSeries, fSeries, order)
        val f3 = multiplySeries(f2, fSeries, order)

        // 1 / (1 - 3125 * z(q))
        val denomP = DoubleArray(order + 1)
        denomP[0] = 1.0
        for (k in 1..order) {
            denomP[k] = -3125.0 * zCoeffs[k]
        }
        val invDenom = invertSeries(denomP, order)

        // C_{ttt} = 5 * invDenom * F^3
        val cttt = multiplySeries(invDenom, f3, order)
        for (k in 0..order) {
            cttt[k] *= 5.0
        }

        return cttt
    }

    private fun invertSeries(a: DoubleArray, order: Int): DoubleArray {
        val b = DoubleArray(order + 1)
        b[0] = 1.0 / a[0]
        for (n in 1..order) {
            var sum = 0.0
            for (k in 1..n) {
                sum += a[k] * b[n - k]
            }
            b[n] = -sum / a[0]
        }
        return b
    }

    private fun multiplySeries(a: DoubleArray, b: DoubleArray, order: Int): DoubleArray {
        val c = DoubleArray(order + 1)
        for (i in 0..order) {
            for (j in 0..(order - i)) {
                c[i + j] += a[i] * b[j]
            }
        }
        return c
    }

    private fun substituteSeries(f: DoubleArray, g: DoubleArray, order: Int): DoubleArray {
        // Computes f(g(q)) where g(0) = 0
        val res = DoubleArray(order + 1)
        res[0] = f[0]

        var gPow = DoubleArray(order + 1)
        gPow[0] = 1.0

        for (n in 1..order) {
            gPow = multiplySeries(gPow, g, order)
            if (n < f.size) {
                val fn = f[n]
                for (k in 0..order) {
                    res[k] += fn * gPow[k]
                }
            }
        }
        return res
    }
}
