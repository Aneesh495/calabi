package calabi.mirror

import org.junit.jupiter.api.Test
import kotlin.test.assertEquals

class InstantonCountingTest {

    @Test
    fun testGromovWittenInstantonInvariants() {
        val instanton = InstantonExpansion()
        val reports = instanton.computeInstantons(maxDegree = 3)

        // n_1 = 2875 (number of lines on the quintic threefold, Schubert 1879)
        assertEquals(2875L, reports[0].instantonCount)

        // n_2 = 609250 (number of conics on the quintic threefold, Katz 1986)
        assertEquals(609250L, reports[1].instantonCount)

        // n_3 = 317206375 (twisted cubics on the quintic threefold, Ellingsrud-Stromme 1991)
        assertEquals(317206375L, reports[2].instantonCount)
    }
}
