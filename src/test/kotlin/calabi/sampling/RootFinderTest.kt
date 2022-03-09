package calabi.sampling

import calabi.algebra.Complex
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class RootFinderTest {

    @Test
    fun testKnownRootsQuartic() {
        // (t - 1)(t + 2)(t - 3i)(t + 3i) = t^4 + t^3 + 7t^2 + 9t - 18
        val coeffs = listOf(
            Complex.ofReal(-18.0),
            Complex.ofReal(9.0),
            Complex.ofReal(7.0),
            Complex.ofReal(1.0),
            Complex.ofReal(1.0)
        )

        val roots = DurandKernerRoots.findRoots(coeffs)
        assertEquals(4, roots.size)

        val expectedRoots = listOf(
            Complex.ofReal(1.0),
            Complex.ofReal(-2.0),
            Complex(0.0, 3.0),
            Complex(0.0, -3.0)
        )

        for (expected in expectedRoots) {
            val matched = roots.any { it.isCloseTo(expected, 1e-8) }
            assertTrue(matched, "Expected root $expected not found in computed roots: $roots")
        }
    }

    @Test
    fun testRootsOfUnityFifthDegree() {
        // t^5 - 1 = 0 => coeffs = [-1, 0, 0, 0, 0, 1]
        val coeffs = listOf(
            Complex.ofReal(-1.0),
            Complex.ZERO,
            Complex.ZERO,
            Complex.ZERO,
            Complex.ZERO,
            Complex.ONE
        )

        val roots = DurandKernerRoots.findRoots(coeffs)
        assertEquals(5, roots.size)

        for (r in roots) {
            assertEquals(1.0, r.abs(), 1e-10)
            val r5 = r.pow(5)
            assertTrue(r5.isCloseTo(Complex.ONE, 1e-10))
        }
    }
}
