package calabi.holonomy

import calabi.algebra.AffineCoordinates
import calabi.algebra.Complex
import calabi.curvature.ChristoffelSymbols
import calabi.linalg.Cholesky
import calabi.linalg.ComplexMatrix
import kotlin.math.abs

/**
 * HolonomyReport summarizes the computed holonomy matrix H around a closed loop:
 * - holonomyMatrix H in U(n)
 * - determinant det(H) (which equals 1 for SU(n) holonomy)
 * - unitarity error ||H^dagger * H - I||_F
 * - phase angle arg(det(H))
 */
data class HolonomyReport(
    val holonomyMatrix: ComplexMatrix,
    val determinant: Complex,
    val unitarityError: Double,
    val phaseAngle: Double,
    val isSpecialUnitary: Boolean
)

/**
 * SU3Holonomy computes the holonomy group transformation matrix H in U(3)
 * around closed loops on Calabi-Yau threefolds, verifying the Ricci-flat condition:
 *
 *   det(H) = 1  <=>  H in SU(3).
 */
class SU3Holonomy(
    val parallelTransport: ParallelTransport
) {
    val manifoldDim: Int = parallelTransport.manifoldDim

    /**
     * Computes the holonomy matrix H by transporting an orthonormal tangent frame
     * around a closed loop gamma(t).
     */
    fun computeHolonomy(
        loopPoints: List<AffineCoordinates>,
        eliminatedIndices: List<Int>,
        hInv: ComplexMatrix
    ): HolonomyReport {
        val hMat = ComplexMatrix(manifoldDim, manifoldDim)

        // Transport standard basis vectors e_1, ..., e_n
        for (col in 0 until manifoldDim) {
            val basisVec = (0 until manifoldDim).map { i ->
                if (i == col) Complex.ONE else Complex.ZERO
            }
            val transported = parallelTransport.transportAlongLoop(
                loopPoints,
                eliminatedIndices,
                hInv,
                basisVec
            )
            for (row in 0 until manifoldDim) {
                hMat.set(row, col, transported[row])
            }
        }

        // Unitarity check: H^dagger * H = I
        val hdH = hMat.conjugateTranspose() * hMat
        val id = ComplexMatrix.identity(manifoldDim)
        var unitErrorSq = 0.0
        for (i in 0 until manifoldDim) {
            for (j in 0 until manifoldDim) {
                val diff = hdH[i, j] - id[i, j]
                unitErrorSq += diff.normSq()
            }
        }
        val unitarityError = kotlin.math.sqrt(unitErrorSq)

        // Compute 3x3 complex determinant:
        // det(H) = H_00(H_11 H_22 - H_12 H_21) - H_01(H_10 H_22 - H_12 H_20) + H_02(H_10 H_21 - H_11 H_20)
        val m00 = hMat[1, 1] * hMat[2, 2] - hMat[1, 2] * hMat[2, 1]
        val m01 = hMat[1, 0] * hMat[2, 2] - hMat[1, 2] * hMat[2, 0]
        val m02 = hMat[1, 0] * hMat[2, 1] - hMat[1, 1] * hMat[2, 0]
        val det = hMat[0, 0] * m00 - hMat[0, 1] * m01 + hMat[0, 2] * m02

        val phase = det.arg()
        val isSU = (det - Complex.ONE).abs() < 0.1 // Numerical tolerance on coarse loop

        return HolonomyReport(
            holonomyMatrix = hMat,
            determinant = det,
            unitarityError = unitarityError,
            phaseAngle = phase,
            isSpecialUnitary = isSU
        )
    }
}
