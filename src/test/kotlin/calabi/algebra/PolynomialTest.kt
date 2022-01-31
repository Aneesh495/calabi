package calabi.algebra

import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class PolynomialTest {

    @Test
    fun testMultiIndexGenerationAndOrdering() {
        // binom(5 + 2 - 1, 2) = binom(6, 2) = 15 quadratic monomials in 5 variables
        val basis2 = MultiIndex.generateDegreeBasis(5, 2)
        assertEquals(15, basis2.size)

        // binom(5 + 3 - 1, 3) = binom(7, 3) = 35 cubic monomials in 5 variables
        val basis3 = MultiIndex.generateDegreeBasis(5, 3)
        assertEquals(35, basis3.size)

        // Verify all have correct total degree
        for (m in basis2) assertEquals(2, m.totalDegree)
        for (m in basis3) assertEquals(3, m.totalDegree)
    }

    @Test
    fun testEulerHomogeneousFunctionTheorem() {
        // For any homogeneous polynomial of degree d: sum z_i * (del P / del z_i) = d * P(z)
        val p = Polynomial.dworkQuintic(Complex(0.5, 0.2))
        assertTrue(p.isHomogeneous)
        assertEquals(5, p.degree)

        val pt = listOf(
            Complex(1.0, 0.5),
            Complex(-0.5, 1.2),
            Complex(0.3, -0.7),
            Complex(1.1, 0.2),
            Complex(-0.9, -0.4)
        )

        val pVal = p.evaluate(pt)
        val grad = p.gradient(pt)

        var eulerSum = Complex.ZERO
        for (i in 0..4) {
            eulerSum += pt[i] * grad[i]
        }

        val expected = pVal * 5.0
        assertTrue(eulerSum.isCloseTo(expected, 1e-12), "Euler theorem failed: $eulerSum vs $expected")
    }

    @Test
    fun testHessianSymmetry() {
        val p = Polynomial.dworkQuintic(Complex(1.5, -0.8))
        val pt = listOf(
            Complex(0.7, 0.1),
            Complex(0.2, -0.4),
            Complex(1.3, 0.9),
            Complex(-0.6, 0.3),
            Complex(0.5, -0.5)
        )

        val hessian = p.hessian(pt)
        for (i in 0..4) {
            for (j in 0..4) {
                assertTrue(
                    hessian[i][j].isCloseTo(hessian[j][i], 1e-12),
                    "Hessian symmetry failed at ($i, $j): ${hessian[i][j]} != ${hessian[j][i]}"
                )
            }
        }
    }

    @Test
    fun testProjectiveLineRestriction() {
        val p = Polynomial.dworkQuintic(Complex(0.8, -0.3))
        val ptA = HomogeneousPoint.ofDoubles(1.0, 0.0, 2.0, -1.0, 0.5)
        val ptB = HomogeneousPoint.ofDoubles(0.2, 1.5, -0.5, 0.8, 1.0)

        val coeffs = p.restrictToProjectiveLine(ptA, ptB)
        assertEquals(6, coeffs.size) // degree 5 has 6 coefficients (a_0 to a_5)

        // Evaluate univariate polynomial at t = 1.7 + 0.4i and compare with P(ptA + t * ptB)
        val t = Complex(1.7, 0.4)
        var polyVal = Complex.ZERO
        var tPow = Complex.ONE
        for (k in coeffs.indices) {
            polyVal += coeffs[k] * tPow
            tPow *= t
        }

        val linePoint = (0..4).map { i -> ptA[i] + t * ptB[i] }
        val directVal = p.evaluate(linePoint)

        assertTrue(
            polyVal.isCloseTo(directVal, 1e-10),
            "Line restriction mismatch: polyVal=$polyVal vs directVal=$directVal"
        )
    }
}
