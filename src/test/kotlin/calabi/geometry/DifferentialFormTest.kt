package calabi.geometry

import calabi.algebra.Complex
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class DifferentialFormTest {

    @Test
    fun testWedgeAntisymmetry() {
        val dz0 = DifferentialForm.dz(4, 0)
        val dz1 = DifferentialForm.dz(4, 1)

        val w01 = dz0.wedge(dz1)
        val w10 = dz1.wedge(dz0)

        // dz0 ^ dz1 = - dz1 ^ dz0
        val sum = w01 + w10
        assertTrue(sum.components.isEmpty(), "Wedge product of 1-forms must be antisymmetric: sum=$sum")

        // dz0 ^ dz0 = 0
        val w00 = dz0.wedge(dz0)
        assertTrue(w00.components.isEmpty(), "Self-wedge of 1-form must be zero")
    }

    @Test
    fun testAssociativityOfWedgeProduct() {
        val dz0 = DifferentialForm.dz(4, 0)
        val dz1 = DifferentialForm.dz(4, 1)
        val dz2 = DifferentialForm.dz(4, 2)

        val wA = (dz0.wedge(dz1)).wedge(dz2)
        val wB = dz0.wedge(dz1.wedge(dz2))

        val diff = wA - wB
        assertTrue(diff.components.isEmpty(), "Wedge product must be associative")
    }

    @Test
    fun testStandardKahlerForm() {
        val kahler = DifferentialForm.standardKahlerForm(3)
        assertEquals(1, kahler.p)
        assertEquals(1, kahler.q)
        assertEquals(3, kahler.components.size)

        // Kahler form is real: bar(omega) = omega
        val conj = kahler.conj()
        for ((basis, coeff) in kahler.components) {
            val conjCoeff = conj.components[basis] ?: Complex.ZERO
            assertTrue(coeff.isCloseTo(conjCoeff, 1e-14), "Kahler form must be real")
        }
    }
}
