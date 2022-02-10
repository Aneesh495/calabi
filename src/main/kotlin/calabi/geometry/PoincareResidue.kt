package calabi.geometry

import calabi.algebra.AffineCoordinates
import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import calabi.algebra.Polynomial
import kotlin.math.pow

/**
 * PoincareResidue computes the nowhere-vanishing holomorphic (n, 0)-volume form Omega
 * on a Calabi-Yau hypersurface X = {P = 0} in P^{n+1} via the algebraic Poincare residue theorem:
 *
 *   Omega = Res_X ( mu / P(z) )
 *
 * where mu = sum_{i=0}^{n+1} (-1)^i z_i dz_0 ^ ... ^ \widehat{dz_i} ^ ... ^ dz_{n+1}.
 */
class PoincareResidue(
    val definingPolynomial: Polynomial
) {
    val ambientDimension: Int = definingPolynomial.numVariables - 1 // P^{n+1} -> ambientDimension = n + 1
    val manifoldDimension: Int = ambientDimension - 1 // n-fold

    init {
        require(definingPolynomial.isHomogeneous) { "Defining polynomial must be homogeneous" }
        require(definingPolynomial.degree == ambientDimension + 1) {
            "For vanishing first Chern class c_1(X) = 0, degree (${definingPolynomial.degree}) must equal ambient dimension + 1 (${ambientDimension + 1})"
        }
    }

    /**
     * Identifies the optimal coordinate to eliminate on the affine chart U_k,
     * chosen to maximize the magnitude of the partial derivative |del P / del w_m|.
     */
    fun selectEliminatedCoordinate(affine: AffineCoordinates): Pair<Int, Complex> {
        val homo = affine.toHomogeneous()
        val grad = definingPolynomial.gradient(homo.coordinates)

        var bestIdx = -1
        var maxDerivativeMod = 0.0

        var affineIdx = 0
        for (i in homo.coordinates.indices) {
            if (i != affine.chartIndex) {
                val mod = grad[i].normSq()
                if (mod > maxDerivativeMod) {
                    maxDerivativeMod = mod
                    bestIdx = affineIdx
                }
                affineIdx++
            }
        }

        require(bestIdx >= 0 && maxDerivativeMod > 1e-14) {
            "Hypersurface is singular at coordinate: partial derivatives vanish"
        }

        val homoElimIdx = if (bestIdx >= affine.chartIndex) bestIdx + 1 else bestIdx
        val derivativeVal = grad[homoElimIdx]
        return Pair(bestIdx, derivativeVal)
    }

    /**
     * Evaluates the holomorphic volume form Omega at an affine point on the variety.
     * Returns the complex coefficient of the volume form on the independent affine tangent directions:
     * Omega = coeff * dw_{i_1} ^ ... ^ dw_{i_n} (omitting dw_m).
     */
    fun evaluateHolomorphicFormCoeff(affine: AffineCoordinates): Pair<Int, Complex> {
        val (elimIdx, deriv) = selectEliminatedCoordinate(affine)
        // Sign is (-1)^elimIdx according to the interior product definition of Poincare residue
        val sign = if ((elimIdx and 1) == 0) 1.0 else -1.0
        val coeff = Complex.ofReal(sign) / deriv
        return Pair(elimIdx, coeff)
    }

    /**
     * Evaluates the invariant Calabi-Yau volume density |Omega|^2 relative to the Lebesgue measure
     * on the remaining n affine coordinates:
     *
     *   dmu_CY = (1 / |del P / del w_m|^2) * prod_{j != m} (i/2 dw_j ^ dbar{w}_j)
     */
    fun volumeDensity(affine: AffineCoordinates): Double {
        val (_, deriv) = selectEliminatedCoordinate(affine)
        val denom = deriv.normSq()
        return 1.0 / denom
    }

    /**
     * Evaluates the invariant ratio between the Calabi-Yau volume form and the Fubini-Study volume form
     * restricted to the tangent space of the hypersurface at this point:
     *
     *   weight(x) = dmu_CY / dmu_FS
     */
    fun localMeasureWeight(affine: AffineCoordinates): Double {
        val cyDensity = volumeDensity(affine)
        val fsDensity = affine.fubiniStudyVolumeDensity()
        // Account for projective scaling factor: sum |z_i|^2
        val homo = affine.toHomogeneous()
        val normSq = homo.coordinates.sumOf { it.normSq() }
        val scaling = normSq.pow(ambientDimension.toDouble())
        return cyDensity * scaling / fsDensity
    }
}
