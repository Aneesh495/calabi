package calabi.algebra

import org.junit.jupiter.api.Test
import org.junit.jupiter.api.assertThrows
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class RationalTest {

    @Test
    fun testCanonicalReduction() {
        val r1 = Rational.of(6, 8)
        assertEquals("3/4", r1.toString())

        val r2 = Rational.of(-10, -20)
        assertEquals("1/2", r2.toString())

        val r3 = Rational.of(15, -5)
        assertEquals("-3", r3.toString())
    }

    @Test
    fun testArithmeticOperations() {
        val a = Rational.of(1, 3)
        val b = Rational.of(1, 6)

        assertEquals(Rational.of(1, 2), a + b)
        assertEquals(Rational.of(1, 6), a - b)
        assertEquals(Rational.of(1, 18), a * b)
        assertEquals(Rational.of(2, 1), a / b)
    }

    @Test
    fun testHarmonicNumbers() {
        // H_1 = 1, H_2 = 3/2, H_3 = 11/6, H_4 = 25/12
        assertEquals(Rational.of(1, 1), Rational.harmonic(1))
        assertEquals(Rational.of(3, 2), Rational.harmonic(2))
        assertEquals(Rational.of(11, 6), Rational.harmonic(3))
        assertEquals(Rational.of(25, 12), Rational.harmonic(4))
    }

    @Test
    fun testBernoulliNumbers() {
        // B_0 = 1, B_1 = -1/2, B_2 = 1/6, B_3 = 0, B_4 = -1/30
        assertEquals(Rational.of(1, 1), Rational.bernoulli(0))
        assertEquals(Rational.of(-1, 2), Rational.bernoulli(1))
        assertEquals(Rational.of(1, 6), Rational.bernoulli(2))
        assertEquals(Rational.ZERO, Rational.bernoulli(3))
        assertEquals(Rational.of(-1, 30), Rational.bernoulli(4))
    }

    @Test
    fun testDivisionByZero() {
        assertThrows<IllegalArgumentException> {
            Rational.of(1, 2) / Rational.ZERO
        }
    }
}
