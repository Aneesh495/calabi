package calabi.metric

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class SectionBundleTest {

    @Test
    fun testBasisDimensions() {
        val basis1 = MonomialBasis(ambientDimension = 4, degree = 1)
        assertEquals(5, basis1.size)

        val basis2 = MonomialBasis(ambientDimension = 4, degree = 2)
        assertEquals(15, basis2.size)

        val basis3 = MonomialBasis(ambientDimension = 4, degree = 3)
        assertEquals(35, basis3.size)
    }

    @Test
    fun testSectionEvaluationAndNorm() {
        val basis = MonomialBasis(ambientDimension = 4, degree = 1)
        val bundle = SectionBundle.identity(basis)

        val pt = HomogeneousPoint.ofDoubles(1.0, 2.0, 0.5, -1.0, 3.0)
        val sections = basis.evaluateAll(pt)
        assertEquals(5, sections.size)
        val vals = sections.map { it.re }.toSet()
        assertEquals(setOf(1.0, 2.0, 0.5, -1.0, 3.0), vals)

        val normSq = bundle.sectionNormSquared(pt, bundle.metricH)
        val expectedNormSq = 1.0*1.0 + 2.0*2.0 + 0.5*0.5 + 1.0*1.0 + 3.0*3.0 // 1 + 4 + 0.25 + 1 + 9 = 15.25
        assertEquals(expectedNormSq, normSq, 1e-10)

        val potential = bundle.kahlerPotential(pt, bundle.metricH)
        assertEquals(kotlin.math.ln(expectedNormSq), potential, 1e-10)
    }
}
