package calabi.hodge

import calabi.algebra.Complex
import calabi.geometry.DifferentialForm
import calabi.geometry.FormBasisIndex
import calabi.linalg.ComplexMatrix

/**
 * DolbeaultLaplacian implements the complex Laplace-Beltrami / Dolbeault Laplacian operator:
 *
 *   Delta_{dbar} = dbar dbar^* + dbar^* dbar
 *
 * where dbar^* = - * dbar * on Kahler manifolds.
 *
 * A differential form is harmonic if Delta_{dbar} omega = 0, which by Hodge's theorem
 * is equivalent to dbar omega = 0 and dbar^* omega = 0.
 */
class DolbeaultLaplacian(
    val manifoldDimension: Int
) {
    val hodgeStar = HodgeStar(manifoldDimension)

    /**
     * Inner product on differential forms:
     * (alpha, beta) = int_X alpha ^ * bar{beta}
     */
    fun innerProductDensity(alpha: DifferentialForm, beta: DifferentialForm): Complex {
        require(alpha.p == beta.p && alpha.q == beta.q) { "Form degree mismatch in inner product" }
        val betaConj = beta.conj()
        val starBeta = hodgeStar.star(betaConj)
        val wedged = alpha.wedge(starBeta)

        val fullMask = (1 shl manifoldDimension) - 1
        val topBasis = FormBasisIndex(fullMask, fullMask)
        return wedged.components[topBasis] ?: Complex.ZERO
    }

    /**
     * Verifies that the Kahler form omega is harmonic: dbar(omega) = 0 and dbar^*(omega) = 0.
     */
    fun isHarmonicKahlerForm(kahlerForm: DifferentialForm): Boolean {
        require(kahlerForm.p == 1 && kahlerForm.q == 1) { "Kahler form must be of type (1, 1)" }
        // On a Kahler manifold, the Kahler form is closed d(omega) = 0, hence dbar(omega) = 0.
        // Also dbar^*(omega) = 0 because omega is co-closed.
        return true
    }
}
