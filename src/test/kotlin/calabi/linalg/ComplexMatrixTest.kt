package calabi.linalg

import calabi.algebra.Complex
import org.junit.jupiter.api.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class ComplexMatrixTest {

    @Test
    fun testBasicMultiplicationAndAdjoint() {
        val a = ComplexMatrix(2, 2)
        a.set(0, 0, Complex(1.0, 2.0))
        a.set(0, 1, Complex(3.0, 0.0))
        a.set(1, 0, Complex(0.0, -1.0))
        a.set(1, 1, Complex(2.0, 1.0))

        val adj = a.conjugateTranspose()
        assertEquals(Complex(1.0, -2.0), adj[0, 0])
        assertEquals(Complex(0.0, 1.0), adj[0, 1])
        assertEquals(Complex(3.0, 0.0), adj[1, 0])
        assertEquals(Complex(2.0, -1.0), adj[1, 1])

        // Verify trace
        assertEquals(Complex(3.0, 3.0), a.trace())

        // Product A * I = A
        val id = ComplexMatrix.identity(2)
        val prod = a * id
        for (i in 0..1) {
            for (j in 0..1) {
                assertTrue(a[i, j].isCloseTo(prod[i, j], 1e-12))
            }
        }
    }

    @Test
    fun testCholeskyFactorizationAndInverse() {
        // Construct a strictly positive-definite 3x3 Hermitian matrix: A = B * B^dagger + 2 * I
        val b = ComplexMatrix(3, 3)
        b.set(0, 0, Complex(2.0, 1.0))
        b.set(0, 1, Complex(1.0, -1.0))
        b.set(0, 2, Complex(0.0, 3.0))
        b.set(1, 0, Complex(-1.0, 2.0))
        b.set(1, 1, Complex(3.0, 0.0))
        b.set(1, 2, Complex(1.0, 1.0))
        b.set(2, 0, Complex(0.0, 1.0))
        b.set(2, 1, Complex(2.0, -1.0))
        b.set(2, 2, Complex(4.0, 0.0))

        val a = (b * b.conjugateTranspose()) + (ComplexMatrix.identity(3) * 2.0)
        assertTrue(a.isHermitian())

        val cholesky = Cholesky.factorize(a)
        val l = cholesky.lMatrix

        // Verify L is lower triangular
        for (i in 0..2) {
            for (j in (i + 1)..2) {
                assertTrue(l[i, j].isZero(1e-14))
            }
            assertTrue(l[i, i].re > 0.0)
            assertTrue(l[i, i].isReal(1e-14))
        }

        // Verify L * L^dagger = A
        val reconstructed = l * l.conjugateTranspose()
        for (i in 0..2) {
            for (j in 0..2) {
                assertTrue(
                    a[i, j].isCloseTo(reconstructed[i, j], 1e-10),
                    "Cholesky reconstruction mismatch at ($i, $j): ${a[i, j]} != ${reconstructed[i, j]}"
                )
            }
        }

        // Verify A * A^{-1} = I
        val inv = cholesky.inverse()
        val prodId = a * inv
        val id3 = ComplexMatrix.identity(3)
        for (i in 0..2) {
            for (j in 0..2) {
                assertTrue(
                    id3[i, j].isCloseTo(prodId[i, j], 1e-10),
                    "Matrix inversion check failed at ($i, $j): ${prodId[i, j]}"
                )
            }
        }
    }

    @Test
    fun testQRFactorization() {
        val a = ComplexMatrix(3, 2)
        a.set(0, 0, Complex(1.0, 1.0))
        a.set(0, 1, Complex(0.0, 2.0))
        a.set(1, 0, Complex(2.0, 0.0))
        a.set(1, 1, Complex(1.0, -1.0))
        a.set(2, 0, Complex(-1.0, 1.0))
        a.set(2, 1, Complex(3.0, 0.0))

        val qr = QRFactorization.factorize(a)
        val q = qr.q
        val r = qr.r

        // Verify Q^dagger * Q = I_2
        val qdq = q.conjugateTranspose() * q
        val id2 = ComplexMatrix.identity(2)
        for (i in 0..1) {
            for (j in 0..1) {
                assertTrue(id2[i, j].isCloseTo(qdq[i, j], 1e-12))
            }
        }

        // Verify Q * R = A
        val reconstructed = q * r
        for (i in 0..2) {
            for (j in 0..1) {
                assertTrue(a[i, j].isCloseTo(reconstructed[i, j], 1e-12))
            }
        }
    }
}
