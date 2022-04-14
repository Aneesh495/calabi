package calabi.curvature

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import calabi.linalg.Cholesky
import calabi.linalg.ComplexMatrix
import calabi.manifold.QuinticThreefold
import calabi.metric.MonomialBasis
import calabi.metric.SectionBundle
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class MetricTensorTest {

    @Test
    fun testMetricTensorHermiticityAndPositiveDefiniteness() {
        val quintic = QuinticThreefold()
        val basis = MonomialBasis(ambientDimension = 4, degree = 1)
        val bundle = SectionBundle.identity(basis)
        val metricTensor = MetricTensor(bundle, quintic.polynomial)

        // Point on Fermat quintic: [1 : e^{i pi / 5} : 0 : 0 : 0]
        val z1 = Complex.fromPolar(1.0, kotlin.math.PI / 5.0)
        val pt = HomogeneousPoint.of(Complex.ONE, z1, Complex.ZERO, Complex.ZERO, Complex.ZERO)
        val affine = pt.toAffineCoordinates(0)

        val (elimIdx, _) = quintic.poincareResidue.selectEliminatedCoordinate(affine)
        val hInv = ComplexMatrix.identity(5)

        val gMat = metricTensor.evaluate(affine, hInv, elimIdx)

        // Dimension of Calabi-Yau 3-fold is 3
        assertEquals(3, gMat.rows)
        assertEquals(3, gMat.cols)

        // Check Hermiticity
        assertTrue(gMat.isHermitian(1e-12), "Induced metric tensor must be Hermitian")

        // Check positive definiteness via Cholesky
        val chol = Cholesky.factorize(gMat)
        for (i in 0 until 3) {
            assertTrue(chol.lMatrix[i, i].re > 0.0, "Diagonal pivot of metric tensor must be positive: index $i")
        }

        // Check determinant
        val detG = metricTensor.determinant(gMat)
        assertTrue(detG > 0.0, "Metric determinant must be strictly positive: $detG")
    }

    @Test
    fun testChristoffelSymbolsTorsionFree() {
        val quintic = QuinticThreefold()
        val basis = MonomialBasis(ambientDimension = 4, degree = 1)
        val bundle = SectionBundle.identity(basis)
        val metricTensor = MetricTensor(bundle, quintic.polynomial)
        val christoffel = ChristoffelSymbols(metricTensor)

        val z1 = Complex.fromPolar(1.0, kotlin.math.PI / 5.0)
        val pt = HomogeneousPoint.of(Complex.ONE, z1, Complex.ZERO, Complex.ZERO, Complex.ZERO)
        val affine = pt.toAffineCoordinates(0)
        val (elimIdx, _) = quintic.poincareResidue.selectEliminatedCoordinate(affine)
        val hInv = ComplexMatrix.identity(5)

        val gamma = christoffel.evaluate(affine, hInv, elimIdx)

        // Verify lower index symmetry: Gamma^i_{j k} == Gamma^i_{k j}
        for (i in 0 until 3) {
            for (j in 0 until 3) {
                for (k in 0 until 3) {
                    val gJK = gamma[i][j][k]
                    val gKJ = gamma[i][k][j]
                    assertTrue(gJK.isCloseTo(gKJ, 1e-12), "Christoffel symmetry failed at ($i, $j, $k)")
                }
            }
        }
    }
}
