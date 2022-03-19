package calabi.linalg

import calabi.algebra.Complex
import kotlin.math.sqrt

/**
 * QRFactorization computes A = Q * R for an m x n complex matrix (m >= n)
 * via Modified Gram-Schmidt orthogonalization.
 *
 * Q is an m x n matrix with orthonormal columns: Q^dagger * Q = I_n.
 * R is an n x n upper-triangular matrix with non-negative real diagonal entries.
 */
class QRFactorization(
    val q: ComplexMatrix,
    val r: ComplexMatrix
) {
    val m: Int = q.rows
    val n: Int = q.cols

    companion object {
        fun factorize(a: ComplexMatrix): QRFactorization {
            val m = a.rows
            val n = a.cols
            require(m >= n) { "QR factorization requires rows >= cols: $m < $n" }

            val q = a.copy()
            val r = ComplexMatrix.zeros(n, n)

            for (k in 0 until n) {
                // Compute 2-norm of column k
                var normSq = 0.0
                for (i in 0 until m) {
                    val re = q.getRe(i, k)
                    val im = q.getIm(i, k)
                    normSq += (re * re + im * im)
                }

                val rkk = sqrt(normSq)
                require(rkk > 1e-15) { "Matrix does not have full column rank at column $k" }
                r.set(k, k, rkk, 0.0)

                // Normalize column k
                for (i in 0 until m) {
                    val re = q.getRe(i, k) / rkk
                    val im = q.getIm(i, k) / rkk
                    q.set(i, k, re, im)
                }

                // Orthogonalize remaining columns against normalized column k
                for (j in (k + 1) until n) {
                    // Compute inner product <q_k, q_j> = sum_i bar(q_{i,k}) * q_{i,j}
                    var innerRe = 0.0
                    var innerIm = 0.0
                    for (i in 0 until m) {
                        val qkRe = q.getRe(i, k)
                        val qkIm = -q.getIm(i, k) // conjugate
                        val qjRe = q.getRe(i, j)
                        val qjIm = q.getIm(i, j)

                        innerRe += (qkRe * qjRe - qkIm * qjIm)
                        innerIm += (qkRe * qjIm + qkIm * qjRe)
                    }

                    r.set(k, j, innerRe, innerIm)

                    // q_j = q_j - r_{k,j} * q_k
                    for (i in 0 until m) {
                        val qkRe = q.getRe(i, k)
                        val qkIm = q.getIm(i, k)

                        val subRe = innerRe * qkRe - innerIm * qkIm
                        val subIm = innerRe * qkIm + innerIm * qkRe

                        val curRe = q.getRe(i, j) - subRe
                        val curIm = q.getIm(i, j) - subIm
                        q.set(i, j, curRe, curIm)
                    }
                }
            }

            return QRFactorization(q, r)
        }
    }
}
