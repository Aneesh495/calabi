package calabi.geometry

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import calabi.algebra.Polynomial
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class PoincareResidueTest {

    @Test
    fun testQuinticPoincareResidueEvaluation() {
        val quintic = Polynomial.fermat(5, 5)
        val residue = PoincareResidue(quintic)

        assertEquals(4, residue.ambientDimension) // P^4
        assertEquals(3, residue.manifoldDimension) // Calabi-Yau 3-fold

        // A known point on the Fermat quintic:
        // [1 : exp(i * pi / 5) : 0 : 0 : 0]
        // 1^5 + (e^{i pi/5})^5 = 1 + (-1) = 0
        val z1 = Complex.fromPolar(1.0, kotlin.math.PI / 5.0)
        val pt = HomogeneousPoint.of(Complex.ONE, z1, Complex.ZERO, Complex.ZERO, Complex.ZERO)

        assertEquals(0.0, quintic.evaluate(pt).abs(), 1e-12)

        val affine = pt.toAffineCoordinates(0)
        val (elimIdx, deriv) = residue.selectEliminatedCoordinate(affine)

        // Deriv is 5 * z1^4, which has norm 5.0
        assertEquals(5.0, deriv.abs(), 1e-12)
        assertEquals(0, elimIdx) // first affine coord (z_1) is eliminated

        val density = residue.volumeDensity(affine)
        assertEquals(1.0 / 25.0, density, 1e-12)

        val weight = residue.localMeasureWeight(affine)
        assertTrue(weight > 0.0, "Local measure weight must be strictly positive: $weight")
    }
}
