package calabi.curvature

import calabi.linalg.ComplexMatrix
import calabi.manifold.QuinticThreefold
import calabi.metric.DonaldsonTOperator
import calabi.metric.MonomialBasis
import calabi.metric.SectionBundle
import calabi.sampling.MonteCarloIntegrator
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class MongeAmpereTest {

    @Test
    fun testDonaldsonSigmaMeasureEvaluation() {
        val quintic = QuinticThreefold()
        val basis = MonomialBasis(ambientDimension = 4, degree = 1)
        val bundle = SectionBundle.identity(basis)
        val metricTensor = MetricTensor(bundle, quintic.polynomial)
        val sigmaMeasure = DonaldsonSigmaMeasure(metricTensor, quintic.poincareResidue)

        val integrator = MonteCarloIntegrator(quintic, seed = 54321L)
        val points = integrator.sampleWeightedPoints(30)

        // Evaluate on initial identity metric
        val h0 = ComplexMatrix.identity(5)
        val report0 = sigmaMeasure.evaluate(h0, points)

        assertTrue(report0.sigma > 0.0, "Sigma error must be positive: ${report0.sigma}")
        assertTrue(report0.sigma.isFinite(), "Sigma error must be finite")
        assertTrue(report0.kappa > 0.0, "Volume constant kappa must be positive: ${report0.kappa}")
        assertEquals(30, report0.sampledPointCount)

        // Perform one step of Donaldson T-operator: h_1 = T(h_0)
        val tOp = DonaldsonTOperator(bundle)
        val h1 = tOp.evaluate(h0, points)
        val report1 = sigmaMeasure.evaluate(h1, points)

        assertTrue(report1.sigma > 0.0)
        assertTrue(report1.kappa > 0.0)
    }
}
