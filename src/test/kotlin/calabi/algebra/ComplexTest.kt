package calabi.algebra

import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import kotlin.math.PI
import kotlin.math.sqrt
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class ComplexTest {

    @Test
    fun testBasicArithmetic() {
        val a = Complex(3.0, 4.0)
        val b = Complex(1.0, -2.0)

        val sum = a + b
        assertEquals(4.0, sum.re, 1e-12)
        assertEquals(2.0, sum.im, 1e-12)

        val diff = a - b
        assertEquals(2.0, diff.re, 1e-12)
        assertEquals(6.0, diff.im, 1e-12)

        val prod = a * b
        // (3+4i)(1-2i) = 3 - 6i + 4i - 8i^2 = 11 - 2i
        assertEquals(11.0, prod.re, 1e-12)
        assertEquals(-2.0, prod.im, 1e-12)

        val div = a / b
        // (3+4i)/(1-2i) = (3+4i)(1+2i)/5 = (3 + 6i + 4i - 8)/5 = -1 + 2i
        assertEquals(-1.0, div.re, 1e-12)
        assertEquals(2.0, div.im, 1e-12)
    }

    @Test
    fun testModulusAndConjugate() {
        val z = Complex(3.0, 4.0)
        assertEquals(5.0, z.abs(), 1e-12)
        assertEquals(25.0, z.normSq(), 1e-12)
        assertEquals(Complex(3.0, -4.0), z.conj())

        // z * conj(z) = |z|^2
        val prod = z * z.conj()
        assertEquals(25.0, prod.re, 1e-12)
        assertEquals(0.0, prod.im, 1e-12)
    }

    @Test
    fun testEulerIdentity() {
        // e^(i * pi) + 1 = 0
        val eIpi = (Complex.I * PI).exp()
        val result = eIpi + Complex.ONE
        assertTrue(result.abs() < 1e-12, "Euler identity error: $result")
    }

    @Test
    fun testPowersAndRoots() {
        val z = Complex(1.0, 1.0)
        val z8 = z.pow(8)
        // (1+i)^2 = 2i, (1+i)^4 = -4, (1+i)^8 = 16
        assertEquals(16.0, z8.re, 1e-12)
        assertEquals(0.0, z8.im, 1e-12)

        val w = Complex(0.0, 4.0)
        val sqrtW = w.sqrt()
        // sqrt(4i) = sqrt(2) + sqrt(2)i
        val expected = sqrt(2.0)
        assertEquals(expected, sqrtW.re, 1e-12)
        assertEquals(expected, sqrtW.im, 1e-12)
    }

    @Test
    fun testRootsOfUnity() {
        val roots5 = Complex.rootsOfUnity(5)
        assertEquals(5, roots5.size)
        // All have norm 1 and 5th power is 1
        for (r in roots5) {
            assertEquals(1.0, r.abs(), 1e-12)
            val r5 = r.pow(5)
            assertTrue(r5.isCloseTo(Complex.ONE, 1e-10))
        }
    }

    @Test
    fun testDivisionByZeroThrows() {
        val a = Complex(1.0, 2.0)
        assertThrows<IllegalArgumentException> {
            a / Complex.ZERO
        }
    }
}
