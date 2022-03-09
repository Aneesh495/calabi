package calabi.sampling

import calabi.algebra.Complex
import calabi.manifold.QuinticThreefold
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class MonteCarloTest {

    @Test
    fun testSampledPointsLieOnHypersurface() {
        val quintic = QuinticThreefold()
        val integrator = MonteCarloIntegrator(quintic, seed = 12345L)

        val points = integrator.sampleWeightedPoints(30)
        assertEquals(30, points.size)

        for (sp in points) {
            val residual = quintic.polynomial.evaluate(sp.point).abs()
            assertTrue(residual < 1e-6, "Point does not lie on hypersurface: residual=$residual")
            assertTrue(sp.weight > 0.0, "Weight must be positive: weight=${sp.weight}")
            assertTrue(sp.weight.isFinite(), "Weight must be finite")
        }

        // Integral of constant 1 should equal 1
        val integralOfOne = integrator.integrate(points) { Complex.ONE }
        assertEquals(1.0, integralOfOne.re, 1e-10)
        assertEquals(0.0, integralOfOne.im, 1e-10)
    }

    @Test
    fun testParallelSampling() {
        val quintic = QuinticThreefold()
        val integrator = MonteCarloIntegrator(quintic, seed = 999L)

        val parallelPoints = integrator.sampleWeightedPointsParallel(25)
        assertEquals(25, parallelPoints.size)

        for (sp in parallelPoints) {
            val residual = quintic.polynomial.evaluate(sp.point).abs()
            assertTrue(residual < 1e-6, "Parallel point residual too high: $residual")
        }
    }
}
