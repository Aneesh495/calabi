package calabi.hodge

import calabi.algebra.Complex
import calabi.geometry.DifferentialForm
import calabi.geometry.FormBasisIndex
import calabi.linalg.ComplexMatrix

/**
 * HodgeStar implements the Hodge star operator on differential (p, q)-forms
 * on an n-dimensional complex Kahler manifold:
 *
 *   *: Omega^{p,q}(X) -> Omega^{n-q, n-p}(X)
 *
 * satisfying alpha ^ * bar{beta} = <alpha, beta>_g * dV_g.
 */
class HodgeStar(
    val manifoldDimension: Int
) {
    /**
     * Applies the Hodge star operator to a basis element dz_I ^ dbar{z}_J:
     * *(dz_I ^ dbar{z}_J) on flat/orthonormal frames with volume element dV.
     */
    fun starBasis(basis: FormBasisIndex): Pair<Int, FormBasisIndex> {
        val fullMask = (1 shl manifoldDimension) - 1
        val complH = fullMask xor basis.holomorphicMask
        val complAH = fullMask xor basis.antiHolomorphicMask

        // Parity of permutation placing I followed by I^c into 0..n-1
        val wedgedH = FormBasisIndex.wedgeMasks(basis.holomorphicMask, complH) ?: Pair(1, fullMask)
        val wedgedAH = FormBasisIndex.wedgeMasks(basis.antiHolomorphicMask, complAH) ?: Pair(1, fullMask)

        val sign = wedgedH.first * wedgedAH.first
        val targetBasis = FormBasisIndex(complAH, complH)
        return Pair(sign, targetBasis)
    }

    /**
     * Computes the Hodge dual *omega of a differential form.
     */
    fun star(omega: DifferentialForm, volumeScale: Double = 1.0): DifferentialForm {
        val targetP = manifoldDimension - omega.q
        val targetQ = manifoldDimension - omega.p
        val newComps = HashMap<FormBasisIndex, Complex>()

        for ((basis, coeff) in omega.components) {
            val (sign, dualBasis) = starBasis(basis)
            val newCoeff = coeff * (sign.toDouble() * volumeScale)
            val cur = newComps[dualBasis] ?: Complex.ZERO
            newComps[dualBasis] = cur + newCoeff
        }

        return DifferentialForm(manifoldDimension, targetP, targetQ, newComps)
    }
}
