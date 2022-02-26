package calabi.manifold

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import org.junit.jupiter.api.Test
import kotlin.math.PI
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class QuinticThreefoldTest {

    @Test
    fun testQuinticTopologicalInvariants() {
        val quintic = QuinticThreefold()
        assertEquals(3, quintic.complexDimension)
        assertEquals(4, quintic.ambientSpace.dimension)
        assertEquals(Pair(1, 101), quintic.hodgeNumbers)
        assertEquals(-200, quintic.eulerCharacteristic)
    }

    @Test
    fun testPointsOnFermatQuintic() {
        val quintic = QuinticThreefold()
        val zPhase = Complex.fromPolar(1.0, PI / 5.0)

        // Point [1 : zPhase : 0 : 0 : 0]
        val pt1 = HomogeneousPoint.of(Complex.ONE, zPhase, Complex.ZERO, Complex.ZERO, Complex.ZERO)
        assertTrue(quintic.contains(pt1, 1e-12))

        // Point not on the quintic: [1 : 1 : 1 : 1 : 1] -> sum is 5 != 0
        val ptNot = HomogeneousPoint.of(Complex.ONE, Complex.ONE, Complex.ONE, Complex.ONE, Complex.ONE)
        assertFalse(quintic.contains(ptNot, 1e-6))
    }

    @Test
    fun testDworkPencilConifoldSingularities() {
        // Conifold at psi = 5
        val pencilConifold = DworkPencil(Complex.ofReal(5.0))
        assertTrue(pencilConifold.manifold.isSingular(1e-4))
        assertTrue(pencilConifold.isNearConifold(0.01))

        // Smooth point at psi = 0.5
        val pencilSmooth = DworkPencil(Complex.ofReal(0.5))
        assertFalse(pencilSmooth.manifold.isSingular(1e-4))
        assertFalse(pencilSmooth.isNearConifold(0.01))

        // Moduli coordinate z = (5 * 0.5)^(-5) / 3125 = (2.5)^(-5) / 3125
        val z = pencilSmooth.moduliCoordinate()
        assertTrue(z.re > 0.0)
    }
}
