package calabi.manifold

import calabi.algebra.Complex
import calabi.algebra.Polynomial
import calabi.algebra.ProjectiveSpace
import calabi.geometry.PoincareResidue

/**
 * K3Surface represents the 2-dimensional Calabi-Yau manifold (K3 surface).
 *
 * Implemented here as the quartic hypersurface in P^3:
 *   P(z) = z_0^4 + z_1^4 + z_2^4 + z_3^4 = 0
 *
 * Ambient space: P^3 (N=3), degree = 4 = N + 1.
 * Complex dimension: dim_C = 2.
 * Hodge diamond:
 *         1
 *       0   0
 *     1  20   1
 *       0   0
 *         1
 *
 * Hodge numbers: h^{1,1} = 20, h^{2,0} = 1.
 * Euler characteristic: chi = 2 - 0 + 20 - 0 + 2 = 24.
 */
class K3Surface(
    val quarticPolynomial: Polynomial = Polynomial.fermat(4, 4)
) : CalabiYauManifold {

    override val name: String = "Quartic K3 Surface in P^3"
    override val complexDimension: Int = 2
    override val ambientSpace: ProjectiveSpace = ProjectiveSpace(3)
    override val definingPolynomials: List<Polynomial> = listOf(quarticPolynomial)

    override val hodgeNumbers: Pair<Int, Int> = Pair(20, 1) // (h^{1,1}, h^{2,0})
    override val eulerCharacteristic: Int = 24

    val poincareResidue: PoincareResidue = PoincareResidue(quarticPolynomial)

    companion object {
        fun fermat(): K3Surface = K3Surface(Polynomial.fermat(4, 4))
    }
}
