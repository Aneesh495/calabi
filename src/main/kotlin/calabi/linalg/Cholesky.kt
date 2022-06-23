package calabi.linalg

import calabi.algebra.Complex
import kotlin.math.sqrt

/**
 * Cholesky performs the complex Cholesky factorization A = L * L^dagger
 * for an n x n Hermitian positive-definite matrix A.
 *
 * L is lower-triangular with strictly positive real diagonal entries.
 */
class Cholesky(
    val lMatrix: ComplexMatrix
) {
    val n: Int = lMatrix.rows

    init {
        require(lMatrix.rows == lMatrix.cols) { "Cholesky matrix must be square" }
    }

    /**
     * Solves L * y = b via forward substitution where L is lower-triangular.
     * b is an n x p matrix (or column vector).
     */
    fun solveLower(b: ComplexMatrix): ComplexMatrix {
        require(b.rows == n) { "Row dimension mismatch: ${b.rows} != $n" }
        val p = b.cols
        val y = ComplexMatrix(n, p)

        for (i in 0 until n) {
            val lDiag = lMatrix.getRe(i, i)
            require(lDiag > 0.0) { "Cholesky diagonal element must be positive: $lDiag at ($i, $i)" }

            for (col in 0 until p) {
                var sumRe = b.getRe(i, col)
                var sumIm = b.getIm(i, col)

                for (k in 0 until i) {
                    val lRe = lMatrix.getRe(i, k)
                    val lIm = lMatrix.getIm(i, k)
                    val yRe = y.getRe(k, col)
                    val yIm = y.getIm(k, col)

                    sumRe -= (lRe * yRe - lIm * yIm)
                    sumIm -= (lRe * yIm + lIm * yRe)
                }

                y.set(i, col, sumRe / lDiag, sumIm / lDiag)
            }
        }
        return y
    }

    /**
     * Solves L^dagger * x = y via back substitution where L^dagger is upper-triangular.
     */
    fun solveUpper(y: ComplexMatrix): ComplexMatrix {
        require(y.rows == n) { "Row dimension mismatch" }
        val p = y.cols
        val x = ComplexMatrix(n, p)

        for (i in (n - 1) downTo 0) {
            val lDiag = lMatrix.getRe(i, i) // L^dagger_{i,i} = bar(L_{i,i}) = L_{i,i} (since real)

            for (col in 0 until p) {
                var sumRe = y.getRe(i, col)
                var sumIm = y.getIm(i, col)

                for (k in (i + 1) until n) {
                    // L^dagger_{i,k} = bar(L_{k,i})
                    val lRe = lMatrix.getRe(k, i)
                    val lIm = -lMatrix.getIm(k, i)
                    val xRe = x.getRe(k, col)
                    val xIm = x.getIm(k, col)

                    sumRe -= (lRe * xRe - lIm * xIm)
                    sumIm -= (lRe * xIm + lIm * xRe)
                }

                x.set(i, col, sumRe / lDiag, sumIm / lDiag)
            }
        }
        return x
    }

    /**
     * Solves A * x = b for x by combining forward and backward substitution:
     * L * y = b, then L^dagger * x = y.
     */
    fun solve(b: ComplexMatrix): ComplexMatrix {
        val y = solveLower(b)
        return solveUpper(y)
    }

    /**
     * Inverts the original Hermitian positive-definite matrix A by solving A * A^{-1} = I.
     */
    fun inverse(): ComplexMatrix =
        solve(ComplexMatrix.identity(n))

    /**
     * Computes the log-determinant ln(det(A)) = 2 * sum_{i=0}^{n-1} ln(L_{i,i}).
     */
    fun logDeterminant(): Double {
        var sum = 0.0
        for (i in 0 until n) {
            val diag = lMatrix.getRe(i, i)
            require(diag > 0.0) { "Non-positive diagonal element in Cholesky factor: $diag" }
            sum += kotlin.math.ln(diag)
        }
        return 2.0 * sum
    }

    companion object {
        /**
         * Computes the Cholesky factor L of an n x n Hermitian positive-definite matrix A.
         */
        fun factorize(a: ComplexMatrix): Cholesky {
            require(a.rows == a.cols) { "Cholesky factorization requires a square matrix: ${a.rows} x ${a.cols}" }
            val n = a.rows
            val l = ComplexMatrix(n, n)

            for (j in 0 until n) {
                var d = a.getRe(j, j)
                for (k in 0 until j) {
                    val lRe = l.getRe(j, k)
                    val lIm = l.getIm(j, k)
                    d -= (lRe * lRe + lIm * lIm)
                }

                require(d > 1e-15) {
                    "Matrix is not Hermitian positive-definite: pivot at index $j is non-positive ($d)"
                }

                val ljj = sqrt(d)
                l.set(j, j, ljj, 0.0)

                for (i in (j + 1) until n) {
                    var sumRe = a.getRe(i, j)
                    var sumIm = a.getIm(i, j)

                    for (k in 0 until j) {
                        val likRe = l.getRe(i, k)
                        val likIm = l.getIm(i, k)
                        val ljkRe = l.getRe(j, k)
                        val ljkIm = -l.getIm(j, k) // bar(L_{j,k})

                        sumRe -= (likRe * ljkRe - likIm * ljkIm)
                        sumIm -= (likRe * ljkIm + likIm * ljkRe)
                    }

                    l.set(i, j, sumRe / ljj, sumIm / ljj)
                }
            }

            return Cholesky(l)
        }

        /**
         * Safely attempts Cholesky factorization, returning null if the matrix
         * is not Hermitian positive-definite or has non-positive pivots.
         */
        fun factorizeOrNull(a: ComplexMatrix): Cholesky? {
            return try {
                factorize(a)
            } catch (e: IllegalArgumentException) {
                null
            }
        }
    }
}

