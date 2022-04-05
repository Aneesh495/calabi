package calabi.metric

import calabi.linalg.Cholesky
import calabi.linalg.ComplexMatrix
import calabi.manifold.QuinticThreefold
import calabi.sampling.MonteCarloIntegrator
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class DonaldsonConvergenceTest {

    @Test
    fun testDonaldsonTOperatorProperties() {
        val quintic = QuinticThreefold()
        val basis = MonomialBasis(ambientDimension = 4, degree = 1)
        val bundle = SectionBundle.identity(basis)
        val tOp = DonaldsonTOperator(bundle)

        val integrator = MonteCarloIntegrator(quintic, seed = 777L)
        val points = integrator.sampleWeightedPoints(40)

        val id = ComplexMatrix.identity(5)
        val tId = tOp.evaluate(id, points)

        // Verify T(I) is Hermitian
        assertTrue(tId.isHermitian(1e-10), "T(h) must be Hermitian")

        // Verify T(I) is positive-definite via Cholesky
        val cholesky = Cholesky.factorize(tId)
        for (i in 0 until 5) {
            assertTrue(cholesky.lMatrix[i, i].re > 0.0, "Cholesky pivots of T(I) must be positive")
        }

        // Trace Tr(h^{-1} T(h)) = Tr(I^{-1} T(I)) = Tr(T(I)) should approximate N_k = 5
        val tr = tId.trace()
        assertEquals(5.0, tr.re, 0.2) // Monte Carlo variance on 40 points
    }

    @Test
    fun testBalancedMetricSolverStep() {
        val quintic = QuinticThreefold()
        val basis = MonomialBasis(ambientDimension = 4, degree = 1)
        val bundle = SectionBundle.identity(basis)
        val solver = BalancedMetricSolver(bundle, maxIterations = 3, tolerance = 1e-3, useAcceleration = false)

        val integrator = MonteCarloIntegrator(quintic, seed = 888L)
        val points = integrator.sampleWeightedPoints(35)

        val report = solver.solve(points)
        assertEquals(3, report.iterations)
        assertTrue(report.finalMetric.isHermitian(1e-8))
    }
}
