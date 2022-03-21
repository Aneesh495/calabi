package calabi.linalg

import calabi.algebra.Complex
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class EigensolverTest {

    @Test
    fun testLanczosDiagonalMatrixEigenvalues() {
        // Known diagonal matrix with eigenvalues 1.0, 3.0, 7.0, 15.0
        val n = 4
        val diag = doubleArrayOf(1.0, 3.0, 7.0, 15.0)
        val mat = ComplexMatrix(n, n)
        for (i in 0 until n) {
            mat.set(i, i, diag[i], 0.0)
        }

        val solver = LanczosEigensolver(mat)
        val computed = solver.computeEigenvalues(steps = 4)

        assertEquals(n, computed.size)
        // Verify minimum and maximum match
        assertEquals(1.0, computed.first(), 1e-6)
        assertEquals(15.0, computed.last(), 1e-6)
    }

    @Test
    fun testHermitianMatrixSpectralBounds() {
        val mat = ComplexMatrix(3, 3)
        mat.set(0, 0, 5.0, 0.0)
        mat.set(0, 1, 1.0, 1.0)
        mat.set(0, 2, 0.0, 2.0)
        mat.set(1, 0, 1.0, -1.0)
        mat.set(1, 1, 8.0, 0.0)
        mat.set(1, 2, -1.0, 0.0)
        mat.set(2, 0, 0.0, -2.0)
        mat.set(2, 1, -1.0, 0.0)
        mat.set(2, 2, 3.0, 0.0)

        assertTrue(mat.isHermitian())

        val solver = LanczosEigensolver(mat)
        val evs = solver.computeEigenvalues(steps = 3)

        // All eigenvalues of a Hermitian matrix are real numbers
        // Trace = sum of eigenvalues = 5 + 8 + 3 = 16
        val trace = evs.sum()
        assertEquals(16.0, trace, 1e-4)
    }
}
