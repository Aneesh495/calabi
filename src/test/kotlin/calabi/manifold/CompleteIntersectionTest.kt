package calabi.manifold

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import org.junit.jupiter.api.Test
import kotlin.math.PI
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class CompleteIntersectionTest {

    @Test
    fun testK3SurfaceInvariants() {
        val k3 = K3Surface.fermat()
        assertEquals(2, k3.complexDimension)
        assertEquals(3, k3.ambientSpace.dimension)
        assertEquals(Pair(20, 1), k3.hodgeNumbers)
        assertEquals(24, k3.eulerCharacteristic)

        // Point [1 : e^{i pi / 4} : 0 : 0]
        val zPhase = Complex.fromPolar(1.0, PI / 4.0)
        val pt = HomogeneousPoint.of(Complex.ONE, zPhase, Complex.ZERO, Complex.ZERO)
        assertTrue(k3.contains(pt, 1e-12))
    }

    @Test
    fun testBicubicThreefoldInP5() {
        val bicubic = CompleteIntersectionCalabiYau.bicubicInP5()
        assertEquals(3, bicubic.complexDimension)
        assertEquals(5, bicubic.ambientSpace.dimension)
        assertEquals(2, bicubic.definingPolynomials.size)
        assertEquals(-144, bicubic.eulerCharacteristic)
        assertEquals(Pair(1, 73), bicubic.hodgeNumbers)
    }
}
