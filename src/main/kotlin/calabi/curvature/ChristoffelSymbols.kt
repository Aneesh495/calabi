package calabi.curvature

import calabi.algebra.AffineCoordinates
import calabi.algebra.Complex
import calabi.linalg.Cholesky
import calabi.linalg.ComplexMatrix

/**
 * ChristoffelSymbols computes the complex Christoffel connection symbols of the second kind:
 *
 *   Gamma^i_{jk} = g^{i bar{l}} * del_j g_{k bar{l}}
 *
 * on the intrinsic tangent bundle of the Kahler manifold.
 *
 * For a Kahler manifold, the metric connection is torsion-free and symmetric in lower indices:
 *   Gamma^i_{jk} = Gamma^i_{kj}.
 */
class ChristoffelSymbols(
    val metricTensor: MetricTensor
) {
    val manifoldDim: Int = metricTensor.manifoldDim

    /**
     * Evaluates Gamma^i_{jk} at an affine point on the manifold using numerical directional derivatives
     * of the exact analytic metric tensor g_{a bar{b}}.
     *
     * Returns an array of size manifoldDim x manifoldDim x manifoldDim.
     */
    fun evaluate(
        affine: AffineCoordinates,
        hInv: ComplexMatrix,
        eliminatedAffineIndex: Int,
        stepSize: Double = 1e-5
    ): Array<Array<Array<Complex>>> {
        val gBase = metricTensor.evaluate(affine, hInv, eliminatedAffineIndex)
        val chol = Cholesky.factorize(gBase)
        val gInv = chol.inverse()

        // Numerical tangent derivative: del_j g_{k bar{l}}
        val dg = Array(manifoldDim) { Array(manifoldDim) { Array(manifoldDim) { Complex.ZERO } } }

        val tangentIndices = ArrayList<Int>(manifoldDim)
        for (i in 0 until affine.dimension) {
            if (i != eliminatedAffineIndex) tangentIndices.add(i)
        }

        for (j in 0 until manifoldDim) {
            val affJ = tangentIndices[j]

            // Forward point: affine + h in coordinate affJ
            val coordsFwd = affine.coords.toMutableList()
            coordsFwd[affJ] = coordsFwd[affJ] + Complex.ofReal(stepSize)
            val affFwd = AffineCoordinates(affine.chartIndex, affine.dimension, coordsFwd)

            val coordsBwd = affine.coords.toMutableList()
            coordsBwd[affJ] = coordsBwd[affJ] - Complex.ofReal(stepSize)
            val affBwd = AffineCoordinates(affine.chartIndex, affine.dimension, coordsBwd)

            val gFwd = metricTensor.evaluate(affFwd, hInv, eliminatedAffineIndex)
            val gBwd = metricTensor.evaluate(affBwd, hInv, eliminatedAffineIndex)

            val inv2H = 1.0 / (2.0 * stepSize)
            for (k in 0 until manifoldDim) {
                for (l in 0 until manifoldDim) {
                    val diff = (gFwd[k, l] - gBwd[k, l]) * inv2H
                    dg[j][k][l] = diff
                }
            }
        }

        // Gamma^i_{jk} = sum_l gInv^{i l} * dg_j(k, l)
        val gamma = Array(manifoldDim) { Array(manifoldDim) { Array(manifoldDim) { Complex.ZERO } } }

        for (i in 0 until manifoldDim) {
            for (j in 0 until manifoldDim) {
                for (k in 0 until manifoldDim) {
                    var sum = Complex.ZERO
                    for (l in 0 until manifoldDim) {
                        val gInvIL = gInv[i, l]
                        val dgJkl = dg[j][k][l]
                        sum += gInvIL * dgJkl
                    }
                    gamma[i][j][k] = sum
                }
            }
        }

        // Symmetrize in lower indices (torsion-free Kahler condition)
        for (i in 0 until manifoldDim) {
            for (j in 0 until manifoldDim) {
                for (k in (j + 1) until manifoldDim) {
                    val avg = (gamma[i][j][k] + gamma[i][k][j]) * 0.5
                    gamma[i][j][k] = avg
                    gamma[i][k][j] = avg
                }
            }
        }

        return gamma
    }
}
