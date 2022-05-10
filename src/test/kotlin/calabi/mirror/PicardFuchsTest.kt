package calabi.mirror

import calabi.algebra.Complex
import java.math.BigInteger
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class PicardFuchsTest {

    @Test
    fun testCombinatorialCoefficients() {
        val frob = FrobeniusMethod(maxOrder = 4)

        // c_0 = 1, c_1 = 120, c_2 = 113400, c_3 = 168168000
        assertEquals(BigInteger.valueOf(1), frob.cCoeffs[0])
        assertEquals(BigInteger.valueOf(120), frob.cCoeffs[1])
        assertEquals(BigInteger.valueOf(113400), frob.cCoeffs[2])
        assertEquals(BigInteger.valueOf(168168000), frob.cCoeffs[3])

        // d_1 = 770, d_2 = 810225
        assertEquals(770.0, frob.dCoeffs[1].toDouble(), 1e-10)
        assertEquals(810225.0, frob.dCoeffs[2].toDouble(), 1e-10)
    }

    @Test
    fun testConifoldSingularity() {
        assertEquals(1.0 / 3125.0, PicardFuchsEquation.CONIFOLD_SINGULARITY.re, 1e-15)
    }
}
