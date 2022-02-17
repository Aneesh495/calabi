package calabi.manifold

import calabi.algebra.Complex
import calabi.algebra.Polynomial
import calabi.algebra.ProjectiveSpace
import calabi.geometry.PoincareResidue

/**
 * QuinticThreefold represents the smooth hypersurface X_5 subset P^4 defined by
 * P(z) = sum_{i=0}^4 z_i^5 = 0 (the Fermat quintic), or with deformations.
 *
 * Hodge numbers: h^{1,1} = 1, h^{2,1} = 101.
 * Euler characteristic: chi = 2 * (1 - 101) = -200.
 */
class QuinticThreefold(
    val deformationParameter: Complex = Complex.ZERO
) : CalabiYauManifold {

    override val name: String =
        if (deformationParameter.isZero(1e-15)) "Fermat Quintic Threefold" else "Deformed Quintic (psi=$deformationParameter)"

    override val complexDimension: Int = 3
    override val ambientSpace: ProjectiveSpace = ProjectiveSpace(4)

    val polynomial: Polynomial = Polynomial.dworkQuintic(deformationParameter)
    override val definingPolynomials: List<Polynomial> = listOf(polynomial)

    override val hodgeNumbers: Pair<Int, Int> = Pair(1, 101)
    override val eulerCharacteristic: Int = 2 * (hodgeNumbers.first - hodgeNumbers.second)

    val poincareResidue: PoincareResidue = PoincareResidue(polynomial)

    /**
     * Discriminant Delta(psi) = (5^5 - psi^5)^5.
     * Conifold singular locus occurs at psi^5 = 5^5 = 3125.
     */
    fun isSingular(tolerance: Double = 1e-6): Boolean {
        val psi5 = deformationParameter.pow(5)
        val conifoldDist = (psi5 - Complex.ofReal(3125.0)).abs()
        return conifoldDist <= tolerance
    }
}
