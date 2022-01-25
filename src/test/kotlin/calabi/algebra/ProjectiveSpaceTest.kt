package calabi.algebra

import org.junit.jupiter.api.Test
import kotlin.math.PI
import kotlin.random.Random
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class ProjectiveSpaceTest {

    @Test
    fun testHomogeneousEquivalenceAndNormalization() {
        val p1 = HomogeneousPoint.ofDoubles(1.0, 2.0, 3.0)
        val p2 = HomogeneousPoint.ofDoubles(2.0, 4.0, 6.0)

        val norm1 = p1.normalized()
        val norm2 = p2.normalized()

        // After normalization, coordinates should match within tolerance
        for (i in 0..2) {
            assertTrue(norm1[i].isCloseTo(norm2[i], 1e-12), "Coordinate mismatch at $i")
        }
        assertEquals(0.0, p1.fubiniStudyDistance(p2), 1e-12)
    }

    @Test
    fun testAffineChartTransitions() {
        val p = HomogeneousPoint.of(
            Complex(2.0, 1.0),
            Complex(0.0, 4.0),
            Complex(1.0, -1.0)
        )

        val chart0 = p.toAffineCoordinates(0)
        val chart1 = p.toAffineCoordinates(1)

        val pRecovered0 = chart0.toHomogeneous()
        val pRecovered1 = chart1.toHomogeneous()

        assertEquals(0.0, p.fubiniStudyDistance(pRecovered0), 1e-12)
        assertEquals(0.0, p.fubiniStudyDistance(pRecovered1), 1e-12)

        // Direct transition
        val transitioned = chart0.transitionTo(1)
        for (i in 0 until chart1.dimension) {
            assertTrue(chart1[i].isCloseTo(transitioned[i], 1e-12))
        }
    }

    @Test
    fun testFubiniStudyMetricHermiticity() {
        val p = HomogeneousPoint.ofDoubles(1.0, 2.0, 3.0, 4.0)
        val affine = p.toAffineCoordinates(0)
        val metric = affine.fubiniStudyMetric()

        // Check g_{i bar{j}} = conj(g_{j bar{i}})
        for (i in 0 until affine.dimension) {
            for (j in 0 until affine.dimension) {
                val gij = metric[i][j]
                val gjiConj = metric[j][i].conj()
                assertTrue(gij.isCloseTo(gjiConj, 1e-14), "Hermiticity failed at ($i, $j)")
            }
            // Diagonal elements must be strictly positive real
            assertTrue(metric[i][i].re > 0.0)
            assertTrue(metric[i][i].isReal(1e-14))
        }
    }

    @Test
    fun testProjectiveLineOrthonormality() {
        val p4 = ProjectiveSpace(4)
        val rng = Random(42)
        val (p, q) = p4.sampleProjectiveLine(rng)

        assertEquals(4, p.dimension)
        assertEquals(4, q.dimension)

        // Check <p, q> = 0
        var inner = Complex.ZERO
        for (i in 0..4) {
            inner += p[i].conj() * q[i]
        }
        assertTrue(inner.abs() < 1e-12, "Spanning points p and q must be orthonormal: inner=$inner")
        assertEquals(PI / 2.0, p.fubiniStudyDistance(q), 1e-12)
    }
}
