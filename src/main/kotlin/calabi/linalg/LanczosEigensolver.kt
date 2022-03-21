package calabi.linalg

import calabi.algebra.Complex
import kotlin.math.abs
import kotlin.math.sqrt
import kotlin.random.Random

/**
 * LanczosEigensolver computes the extreme eigenvalues and eigenvectors of an n x n
 * Hermitian matrix A using the iterative Lanczos algorithm.
 */
class LanczosEigensolver(
    val matrix: ComplexMatrix
) {
    val n: Int = matrix.rows

    init {
        require(matrix.rows == matrix.cols) { "Eigensolver requires a square matrix: ${matrix.rows} x ${matrix.cols}" }
    }

    /**
     * Executes m steps of the Lanczos tridiagonalization algorithm.
     * Returns the diagonal elements alpha and off-diagonal elements beta of the tridiagonal matrix T.
     */
    fun tridiagonalize(maxSteps: Int, random: Random = Random(42)): Pair<DoubleArray, DoubleArray> {
        val m = maxSteps.coerceAtMost(n)
        val alpha = DoubleArray(m)
        val beta = DoubleArray(m - 1)

        // Random initial normalized complex vector v_1
        var vPrev = DoubleArray(n * 2)
        val vCurrent = DoubleArray(n * 2)

        var normSq = 0.0
        for (i in 0 until (n * 2)) {
            val r = random.nextDouble() - 0.5
            vCurrent[i] = r
            normSq += r * r
        }
        val norm = sqrt(normSq)
        for (i in 0 until (n * 2)) {
            vCurrent[i] /= norm
        }

        val vNext = DoubleArray(n * 2)

        for (j in 0 until m) {
            // Compute w = A * vCurrent
            for (r in 0 until n) {
                var sumRe = 0.0
                var sumIm = 0.0
                for (c in 0 until n) {
                    val aRe = matrix.getRe(r, c)
                    val aIm = matrix.getIm(r, c)
                    val vRe = vCurrent[c * 2]
                    val vIm = vCurrent[c * 2 + 1]

                    sumRe += (aRe * vRe - aIm * vIm)
                    sumIm += (aRe * vIm + aIm * vRe)
                }
                vNext[r * 2] = sumRe
                vNext[r * 2 + 1] = sumIm
            }

            // alpha_j = <vCurrent, w> = sum_i bar(vCurrent_i) * w_i
            var alph = 0.0
            for (i in 0 until n) {
                val vRe = vCurrent[i * 2]
                val vIm = -vCurrent[i * 2 + 1]
                val wRe = vNext[i * 2]
                val wIm = vNext[i * 2 + 1]
                alph += (vRe * wRe - vIm * wIm)
            }
            alpha[j] = alph

            // w = w - alpha_j * vCurrent - beta_{j-1} * vPrev
            val bPrev = if (j > 0) beta[j - 1] else 0.0
            for (i in 0 until (n * 2)) {
                vNext[i] -= (alph * vCurrent[i] + bPrev * vPrev[i])
            }

            if (j < m - 1) {
                var bNormSq = 0.0
                for (i in 0 until (n * 2)) {
                    bNormSq += vNext[i] * vNext[i]
                }
                val bet = sqrt(bNormSq)
                beta[j] = bet

                if (bet < 1e-13) {
                    return Pair(alpha.copyOf(j + 1), beta.copyOf(j))
                }

                for (i in 0 until (n * 2)) {
                    vPrev[i] = vCurrent[i]
                    vCurrent[i] = vNext[i] / bet
                }
            }
        }

        return Pair(alpha, beta)
    }

    /**
     * Computes the eigenvalues of the symmetric tridiagonal matrix T(alpha, beta)
     * using the implicit symmetric QR algorithm with Wilkinson shifts.
     */
    fun tridiagonalEigenvalues(alpha: DoubleArray, beta: DoubleArray): DoubleArray {
        val k = alpha.size
        val d = alpha.clone()
        val e = DoubleArray(k)
        for (i in 0 until (k - 1)) {
            e[i] = beta[i]
        }

        for (l in 0 until k) {
            var iter = 0
            while (true) {
                var m = l
                while (m < k - 1) {
                    val dd = abs(d[m]) + abs(d[m + 1])
                    if (abs(e[m]) <= 1e-15 * dd) break
                    m++
                }
                if (m == l) break
                require(iter++ < 100) { "Too many iterations in tridiagonal QR algorithm" }

                var g = (d[l + 1] - d[l]) / (2.0 * e[l])
                var r = kotlin.math.hypot(g, 1.0)
                g = d[m] - d[l] + e[l] / (g + kotlin.math.sign(if (g == 0.0) 1.0 else g) * r)
                var s = 1.0
                var c = 1.0
                var p = 0.0

                var brokenEarly = false
                for (i in (m - 1) downTo l) {
                    var f = s * e[i]
                    val b = c * e[i]
                    r = kotlin.math.hypot(f, g)
                    e[i + 1] = r
                    if (r == 0.0) {
                        d[i + 1] -= p
                        e[m] = 0.0
                        brokenEarly = true
                        break
                    }
                    s = f / r
                    c = g / r
                    g = d[i + 1] - p
                    r = (d[i] - g) * s + 2.0 * c * b
                    p = s * r
                    d[i + 1] = g + p
                    g = c * r - b
                }
                if (brokenEarly) continue
                d[l] -= p
                e[l] = g
                e[m] = 0.0
            }
        }

        d.sort()
        return d
    }

    /**
     * Computes the eigenvalues of the Hermitian matrix using Lanczos.
     */
    fun computeEigenvalues(steps: Int = 30): DoubleArray {
        val (alpha, beta) = tridiagonalize(steps)
        return tridiagonalEigenvalues(alpha, beta)
    }
}
