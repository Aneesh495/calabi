package calabi.metric

import calabi.linalg.ComplexMatrix
import kotlin.math.sqrt

/**
 * AndersonAcceleration implements multi-secant acceleration for the fixed-point iteration
 * x_{k+1} = T(x_k), mixing m previous iterates to minimize the residual norm.
 */
class AndersonAcceleration(
    val size: Int,
    val memoryDepth: Int = 3,
    val damping: Double = 0.8
) {
    private val xHistory = ArrayList<DoubleArray>()
    private val rHistory = ArrayList<DoubleArray>()

    /**
     * Updates Anderson history with current iterate x and evaluated target T(x).
     * Returns the accelerated next iterate.
     */
    fun step(x: ComplexMatrix, tVal: ComplexMatrix): ComplexMatrix {
        val n = x.data.size
        val xArr = x.data.clone()
        val rArr = DoubleArray(n) { i -> tVal.data[i] - xArr[i] }

        xHistory.add(xArr)
        rHistory.add(rArr)

        if (xHistory.size > memoryDepth + 1) {
            xHistory.removeAt(0)
            rHistory.removeAt(0)
        }

        val k = xHistory.size - 1
        if (k == 0) {
            // First step: standard damped fixed point
            val nextData = DoubleArray(n) { i -> xArr[i] + damping * rArr[i] }
            return ComplexMatrix(x.rows, x.cols, nextData)
        }

        // Form difference matrices Delta R and Delta X
        // Solve min_{gamma} || r_k - sum_{j=0}^{k-1} gamma_j (r_k - r_j) ||_2
        // Using unconstrained least squares / normal equations
        val m = k
        val gMat = Array(m) { DoubleArray(m) }
        val rhs = DoubleArray(m)

        for (i in 0 until m) {
            // diffR_i = r_k - r_i
            var dotRhs = 0.0
            for (p in 0 until n) {
                val diffI = rHistory[k][p] - rHistory[i][p]
                dotRhs += diffI * rHistory[k][p]
            }
            rhs[i] = dotRhs

            for (j in 0 until m) {
                var dotG = 0.0
                for (p in 0 until n) {
                    val diffI = rHistory[k][p] - rHistory[i][p]
                    val diffJ = rHistory[k][p] - rHistory[j][p]
                    dotG += diffI * diffJ
                }
                gMat[i][j] = dotG
                if (i == j) gMat[i][j] += 1e-8 // Regularization
            }
        }

        // Solve small m x m system G * gamma = rhs via Gaussian elimination
        val gamma = solveLinearSystem(gMat, rhs)

        // Compute accelerated point: xNext = (1 - damping) * sum + damping * ...
        var sumGamma = 0.0
        for (j in 0 until m) sumGamma += gamma[j]
        val alphaK = 1.0 - sumGamma

        val nextData = DoubleArray(n)
        for (p in 0 until n) {
            var xMix = alphaK * xHistory[k][p]
            var rMix = alphaK * rHistory[k][p]
            for (j in 0 until m) {
                xMix += gamma[j] * xHistory[j][p]
                rMix += gamma[j] * rHistory[j][p]
            }
            nextData[p] = xMix + damping * rMix
        }

        return ComplexMatrix(x.rows, x.cols, nextData)
    }

    private fun solveLinearSystem(a: Array<DoubleArray>, b: DoubleArray): DoubleArray {
        val n = b.size
        val aug = Array(n) { i -> DoubleArray(n + 1) { j -> if (j < n) a[i][j] else b[i] } }

        for (i in 0 until n) {
            var maxRow = i
            for (k in (i + 1) until n) {
                if (kotlin.math.abs(aug[k][i]) > kotlin.math.abs(aug[maxRow][i])) {
                    maxRow = k
                }
            }
            val temp = aug[i]
            aug[i] = aug[maxRow]
            aug[maxRow] = temp

            val pivot = aug[i][i]
            if (kotlin.math.abs(pivot) < 1e-14) continue
            for (j in i until (n + 1)) {
                aug[i][j] /= pivot
            }

            for (k in 0 until n) {
                if (k != i) {
                    val factor = aug[k][i]
                    for (j in i until (n + 1)) {
                        aug[k][j] -= factor * aug[i][j]
                    }
                }
            }
        }

        return DoubleArray(n) { i -> aug[i][n] }
    }
}
