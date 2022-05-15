package calabi.holonomy

import calabi.algebra.AffineCoordinates
import calabi.algebra.Complex
import calabi.curvature.ChristoffelSymbols
import calabi.curvature.MetricTensor
import calabi.linalg.ComplexMatrix
import kotlin.math.PI

/**
 * ParallelTransport integrates the parallel transport differential equation for complex tangent vectors:
 *
 *   d v^i / dt + Gamma^i_{jk}(gamma(t)) * (d gamma^j / dt) * v^k(t) = 0
 *
 * along parameterized curves gamma(t) on the Calabi-Yau manifold using 4th-order Runge-Kutta.
 */
class ParallelTransport(
    val christoffelSymbols: ChristoffelSymbols
) {
    val manifoldDim: Int = christoffelSymbols.manifoldDim

    /**
     * Transports a complex tangent vector v along a closed loop gamma(t) discretized into N points.
     * Returns the transported vector at the end of the loop.
     */
    fun transportAlongLoop(
        loopPoints: List<AffineCoordinates>,
        eliminatedIndices: List<Int>,
        hInv: ComplexMatrix,
        initialVector: List<Complex>
    ): List<Complex> {
        val numSteps = loopPoints.size - 1
        require(numSteps >= 3) { "Loop must contain at least 3 segments" }

        var currentV = initialVector.toMutableList()

        for (step in 0 until numSteps) {
            val affCurrent = loopPoints[step]
            val affNext = loopPoints[step + 1]
            val elimIdx = eliminatedIndices[step]

            // Tangent velocity vector d gamma / dt
            val dt = 1.0 / numSteps.toDouble()
            val dGamma = Array(manifoldDim) { a ->
                (affNext[a] - affCurrent[a]) / dt
            }

            // Connection at current point: Gamma^i_{jk}
            val gammaMat = christoffelSymbols.evaluate(affCurrent, hInv, elimIdx)

            // RK4 step for dv/dt = - Gamma(gamma, dgamma, v)
            fun evalDeriv(v: List<Complex>): List<Complex> {
                val dv = ArrayList<Complex>(manifoldDim)
                for (i in 0 until manifoldDim) {
                    var sum = Complex.ZERO
                    for (j in 0 until manifoldDim) {
                        for (k in 0 until manifoldDim) {
                            sum += gammaMat[i][j][k] * dGamma[j] * v[k]
                        }
                    }
                    dv.add(-sum)
                }
                return dv
            }

            val k1 = evalDeriv(currentV)
            val vMid1 = (0 until manifoldDim).map { i -> currentV[i] + k1[i] * (0.5 * dt) }

            val k2 = evalDeriv(vMid1)
            val vMid2 = (0 until manifoldDim).map { i -> currentV[i] + k2[i] * (0.5 * dt) }

            val k3 = evalDeriv(vMid2)
            val vEnd = (0 until manifoldDim).map { i -> currentV[i] + k3[i] * dt }

            val k4 = evalDeriv(vEnd)

            for (i in 0 until manifoldDim) {
                currentV[i] = currentV[i] + (k1[i] + k2[i] * 2.0 + k3[i] * 2.0 + k4[i]) * (dt / 6.0)
            }
        }

        return currentV
    }
}
