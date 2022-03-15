package calabi.linalg

import calabi.algebra.Complex
import kotlin.math.sqrt

/**
 * ComplexMatrix represents a dense m x n matrix over C.
 *
 * Internally backed by a flat primitive DoubleArray storing interleaved real and imaginary parts:
 *   [re_{0,0}, im_{0,0}, re_{0,1}, im_{0,1}, ..., re_{m-1,n-1}, im_{m-1,n-1}]
 *
 * This provides zero-allocation access, contiguous memory cache locality, and
 * direct compatibility with native BLAS and SIMD kernels.
 */
class ComplexMatrix(
    val rows: Int,
    val cols: Int,
    val data: DoubleArray = DoubleArray(rows * cols * 2)
) {
    init {
        require(rows > 0 && cols > 0) { "Matrix dimensions must be strictly positive: $rows x $cols" }
        require(data.size == rows * cols * 2) { "Data array length must equal rows * cols * 2" }
    }

    @Suppress("NOTHING_TO_INLINE")
    private inline fun offset(r: Int, c: Int): Int = (r * cols + c) shl 1

    operator fun get(r: Int, c: Int): Complex {
        val off = offset(r, c)
        return Complex(data[off], data[off + 1])
    }

    operator fun set(r: Int, c: Int, v: Complex) {
        val off = offset(r, c)
        data[off] = v.re
        data[off + 1] = v.im
    }

    fun set(r: Int, c: Int, re: Double, im: Double) {
        val off = offset(r, c)
        data[off] = re
        data[off + 1] = im
    }

    fun getRe(r: Int, c: Int): Double = data[offset(r, c)]
    fun getIm(r: Int, c: Int): Double = data[offset(r, c) + 1]

    operator fun plus(other: ComplexMatrix): ComplexMatrix {
        require(rows == other.rows && cols == other.cols) { "Dimension mismatch in addition" }
        val res = ComplexMatrix(rows, cols)
        for (i in data.indices) {
            res.data[i] = data[i] + other.data[i]
        }
        return res
    }

    operator fun minus(other: ComplexMatrix): ComplexMatrix {
        require(rows == other.rows && cols == other.cols) { "Dimension mismatch in subtraction" }
        val res = ComplexMatrix(rows, cols)
        for (i in data.indices) {
            res.data[i] = data[i] - other.data[i]
        }
        return res
    }

    operator fun times(scalar: Complex): ComplexMatrix {
        val res = ComplexMatrix(rows, cols)
        val sRe = scalar.re
        val sIm = scalar.im
        var i = 0
        while (i < data.size) {
            val r = data[i]
            val im = data[i + 1]
            res.data[i] = r * sRe - im * sIm
            res.data[i + 1] = r * sIm + im * sRe
            i += 2
        }
        return res
    }

    operator fun times(scalar: Double): ComplexMatrix {
        val res = ComplexMatrix(rows, cols)
        for (i in data.indices) {
            res.data[i] = data[i] * scalar
        }
        return res
    }

    /**
     * Matrix multiplication C = A * B.
     * Uses cache-friendly loop order (i, k, j).
     */
    operator fun times(other: ComplexMatrix): ComplexMatrix {
        require(cols == other.rows) { "Inner dimension mismatch: ($rows x $cols) * (${other.rows} x ${other.cols})" }
        val p = other.cols
        val res = ComplexMatrix(rows, p)

        for (i in 0 until rows) {
            for (k in 0 until cols) {
                val aOff = offset(i, k)
                val aRe = data[aOff]
                val aIm = data[aOff + 1]
                if (aRe == 0.0 && aIm == 0.0) continue

                for (j in 0 until p) {
                    val bOff = other.offset(k, j)
                    val bRe = other.data[bOff]
                    val bIm = other.data[bOff + 1]

                    val cOff = res.offset(i, j)
                    res.data[cOff] += aRe * bRe - aIm * bIm
                    res.data[cOff + 1] += aRe * bIm + aIm * bRe
                }
            }
        }
        return res
    }

    /**
     * Conjugate transpose (Hermitian adjoint) A^dagger = bar(A)^T.
     */
    fun conjugateTranspose(): ComplexMatrix {
        val res = ComplexMatrix(cols, rows)
        for (r in 0 until rows) {
            for (c in 0 until cols) {
                val srcOff = offset(r, c)
                val dstOff = res.offset(c, r)
                res.data[dstOff] = data[srcOff]
                res.data[dstOff + 1] = -data[srcOff + 1]
            }
        }
        return res
    }

    /**
     * Matrix trace sum_{i=0}^{n-1} A_{i,i}.
     */
    fun trace(): Complex {
        require(rows == cols) { "Trace requires a square matrix: $rows x $cols" }
        var sumRe = 0.0
        var sumIm = 0.0
        for (i in 0 until rows) {
            val off = offset(i, i)
            sumRe += data[off]
            sumIm += data[off + 1]
        }
        return Complex(sumRe, sumIm)
    }

    /**
     * Checks if this matrix is Hermitian: A == A^dagger within tolerance.
     */
    fun isHermitian(tolerance: Double = 1e-12): Boolean {
        if (rows != cols) return false
        for (i in 0 until rows) {
            if (kotlin.math.abs(getIm(i, i)) > tolerance) return false
            for (j in (i + 1) until cols) {
                val offIJ = offset(i, j)
                val offJI = offset(j, i)
                val diffRe = kotlin.math.abs(data[offIJ] - data[offJI])
                val diffIm = kotlin.math.abs(data[offIJ + 1] + data[offJI + 1])
                if (diffRe > tolerance || diffIm > tolerance) return false
            }
        }
        return true
    }

    fun copy(): ComplexMatrix =
        ComplexMatrix(rows, cols, data.clone())

    companion object {
        fun identity(n: Int): ComplexMatrix {
            val mat = ComplexMatrix(n, n)
            for (i in 0 until n) {
                mat.set(i, i, 1.0, 0.0)
            }
            return mat
        }

        fun zeros(rows: Int, cols: Int): ComplexMatrix =
            ComplexMatrix(rows, cols)
    }
}
