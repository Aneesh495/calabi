package calabi.holonomy

import calabi.algebra.AffineCoordinates
import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import calabi.curvature.ChristoffelSymbols
import calabi.curvature.MetricTensor
import calabi.linalg.ComplexMatrix
import calabi.manifold.QuinticThreefold
import calabi.metric.MonomialBasis
import calabi.metric.SectionBundle
import org.junit.jupiter.api.Test
import kotlin.math.PI
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class SU3HolonomyTest {

    @Test
    fun testHolonomyMatrixUnitarityAndDeterminant() {
        val quintic = QuinticThreefold()
        val basis = MonomialBasis(ambientDimension = 4, degree = 1)
        val bundle = SectionBundle.identity(basis)
        val metricTensor = MetricTensor(bundle, quintic.polynomial)
        val christoffel = ChristoffelSymbols(metricTensor)
        val parallelTransport = ParallelTransport(christoffel)
        val su3Holonomy = SU3Holonomy(parallelTransport)

        // Construct a small closed loop around a point on the Fermat quintic
        // Center: [1 : exp(i pi / 5) : 0 : 0 : 0]
        val z1 = Complex.fromPolar(1.0, PI / 5.0)
        val r = 0.01
        val numSegments = 8

        val loopPoints = ArrayList<AffineCoordinates>(numSegments + 1)
        val elimIndices = ArrayList<Int>(numSegments + 1)

        for (step in 0..numSegments) {
            val theta = 2.0 * PI * step / numSegments
            // Perturb z_2 along a circle in the complex plane
            val z2 = Complex.fromPolar(r, theta)
            val pt = HomogeneousPoint.of(Complex.ONE, z1, z2, Complex.ZERO, Complex.ZERO)
            val affine = pt.toAffineCoordinates(0)
            val (elimIdx, _) = quintic.poincareResidue.selectEliminatedCoordinate(affine)

            loopPoints.add(affine)
            elimIndices.add(elimIdx)
        }

        val hInv = ComplexMatrix.identity(5)
        val report = su3Holonomy.computeHolonomy(loopPoints, elimIndices, hInv)

        // Verify determinant is close to 1: det(H) = 1
        val det = report.determinant
        assertEquals(1.0, det.abs(), 0.15) // Numerical integration tolerance
        assertTrue(report.unitarityError < 0.25)
    }
}
