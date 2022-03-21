package calabi.algebra

import java.math.BigDecimal
import java.math.MathContext
import kotlin.math.*

/**
 * Complex represents an element of the field of complex numbers C.
 *
 * Implements standard field operations, complex conjugation, polar decompositions,
 * and complex analytic functions with standard branch cut conventions along the
 * negative real axis.
 *
 * All operations are pure and immutable.
 */
data class Complex(
    val re: Double,
    val im: Double
) : Comparable<Complex> {

    operator fun plus(other: Complex): Complex =
        Complex(re + other.re, im + other.im)

    operator fun plus(other: Double): Complex =
        Complex(re + other, im)

    operator fun plus(other: Number): Complex =
        this + other.toDouble()

    operator fun minus(other: Complex): Complex =
        Complex(re - other.re, im - other.im)

    operator fun minus(other: Double): Complex =
        Complex(re - other, im)

    operator fun minus(other: Number): Complex =
        this - other.toDouble()

    operator fun unaryMinus(): Complex =
        Complex(-re, -im)

    operator fun times(other: Complex): Complex =
        Complex(re * other.re - im * other.im, re * other.im + im * other.re)

    operator fun times(scalar: Double): Complex =
        Complex(re * scalar, im * scalar)

    operator fun times(scalar: Number): Complex =
        this * scalar.toDouble()

    operator fun div(other: Complex): Complex {
        val denom = other.re * other.re + other.im * other.im
        require(denom != 0.0) { "Division by complex zero" }
        return Complex(
            (re * other.re + im * other.im) / denom,
            (im * other.re - re * other.im) / denom
        )
    }

    operator fun div(scalar: Double): Complex {
        require(scalar != 0.0) { "Division by zero scalar" }
        return Complex(re / scalar, im / scalar)
    }

    operator fun div(scalar: Number): Complex =
        this / scalar.toDouble()

    /**
     * Squared Euclidean modulus |z|^2 = re^2 + im^2.
     */
    fun normSq(): Double = re * re + im * im

    /**
     * Euclidean modulus |z| = sqrt(re^2 + im^2).
     */
    fun abs(): Double = hypot(re, im)

    /**
     * Argument / phase angle in radians in the interval (-pi, pi].
     */
    fun arg(): Double = atan2(im, re)

    /**
     * Complex conjugate bar(z) = re - i * im.
     */
    fun conj(): Complex = Complex(re, -im)

    /**
     * Multiplicative inverse 1 / z.
     */
    fun inv(): Complex {
        val d = normSq()
        require(d != 0.0) { "Inverse of zero complex number" }
        return Complex(re / d, -im / d)
    }

    /**
     * Complex exponential exp(z) = exp(re) * (cos(im) + i * sin(im)).
     */
    fun exp(): Complex {
        val r = kotlin.math.exp(re)
        return Complex(r * kotlin.math.cos(im), r * kotlin.math.sin(im))
    }

    /**
     * Principal branch of the complex logarithm ln(z) = ln|z| + i * arg(z).
     * Branch cut on the negative real axis (-infinity, 0].
     */
    fun ln(): Complex {
        val r = abs()
        require(r > 0.0) { "Logarithm of zero" }
        return Complex(kotlin.math.ln(r), arg())
    }

    /**
     * Principal square root sqrt(z) = sqrt|z| * exp(i * arg(z) / 2).
     */
    fun sqrt(): Complex {
        if (re == 0.0 && im == 0.0) return ZERO
        val r = abs()
        val u = kotlin.math.sqrt((r + re) / 2.0)
        val v = kotlin.math.sign(if (im == 0.0) 1.0 else im) * kotlin.math.sqrt(max(0.0, (r - re) / 2.0))
        return Complex(u, v)
    }

    /**
     * Integer power z^n via binary exponentiation.
     */
    fun pow(n: Int): Complex {
        if (n == 0) return ONE
        if (n < 0) return ONE / pow(-n)
        var base = this
        var exp = n
        var result = ONE
        while (exp > 0) {
            if ((exp and 1) == 1) {
                result *= base
            }
            base *= base
            exp = exp shr 1
        }
        return result
    }

    /**
     * Arbitrary complex power z^w = exp(w * ln(z)).
     */
    fun pow(w: Complex): Complex {
        if (this == ZERO) {
            if (w.re > 0.0) return ZERO
            throw IllegalArgumentException("0^w is undefined for Re(w) <= 0")
        }
        return (w * this.ln()).exp()
    }

    fun pow(p: Double): Complex = pow(Complex(p, 0.0))

    fun sin(): Complex =
        Complex(kotlin.math.sin(re) * kotlin.math.cosh(im), kotlin.math.cos(re) * kotlin.math.sinh(im))

    fun cos(): Complex =
        Complex(kotlin.math.cos(re) * kotlin.math.cosh(im), -kotlin.math.sin(re) * kotlin.math.sinh(im))

    fun tan(): Complex = sin() / cos()

    fun sinh(): Complex =
        Complex(kotlin.math.sinh(re) * kotlin.math.cos(im), kotlin.math.cosh(re) * kotlin.math.sin(im))

    fun cosh(): Complex =
        Complex(kotlin.math.cosh(re) * kotlin.math.cos(im), kotlin.math.sinh(re) * kotlin.math.sin(im))

    /**
     * Approximate equality check with absolute and relative tolerance.
     */
    fun isCloseTo(other: Complex, tolerance: Double = 1e-12): Boolean {
        val diffSq = (re - other.re) * (re - other.re) + (im - other.im) * (im - other.im)
        if (diffSq <= tolerance * tolerance) return true
        val scale = max(normSq(), other.normSq())
        return diffSq <= tolerance * tolerance * scale
    }

    fun isZero(tolerance: Double = 1e-14): Boolean =
        abs() <= tolerance

    fun isReal(tolerance: Double = 1e-14): Boolean =
        abs(im) <= tolerance

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is Complex) return false
        val rA = if (re == 0.0) 0.0 else re
        val rB = if (other.re == 0.0) 0.0 else other.re
        val iA = if (im == 0.0) 0.0 else im
        val iB = if (other.im == 0.0) 0.0 else other.im
        return rA == rB && iA == iB
    }

    override fun hashCode(): Int {
        val r = if (re == 0.0) 0.0 else re
        val i = if (im == 0.0) 0.0 else im
        return 31 * r.hashCode() + i.hashCode()
    }

    override fun compareTo(other: Complex): Int {
        val c = re.compareTo(other.re)
        return if (c != 0) c else im.compareTo(other.im)
    }

    override fun toString(): String {
        return when {
            abs(im) < 1e-15 -> String.format("%.8g", re)
            abs(re) < 1e-15 -> String.format("%.8gi", im)
            im < 0.0 -> String.format("%.8g - %.8gi", re, -im)
            else -> String.format("%.8g + %.8gi", re, im)
        }
    }

    companion object {
        val ZERO = Complex(0.0, 0.0)
        val ONE = Complex(1.0, 0.0)
        val I = Complex(0.0, 1.0)
        val PI = kotlin.math.PI
        val TWO_PI = 2.0 * kotlin.math.PI

        fun ofReal(re: Double): Complex = Complex(re, 0.0)
        fun ofImag(im: Double): Complex = Complex(0.0, im)

        fun fromPolar(r: Double, theta: Double): Complex {
            require(r >= 0.0) { "Radial distance must be non-negative" }
            return Complex(r * kotlin.math.cos(theta), r * kotlin.math.sin(theta))
        }

        fun rootsOfUnity(n: Int): List<Complex> {
            require(n >= 1) { "Root order must be at least 1" }
            return (0 until n).map { k ->
                val angle = TWO_PI * k / n
                fromPolar(1.0, angle)
            }
        }
    }
}

operator fun Double.plus(c: Complex): Complex = Complex(this + c.re, c.im)
operator fun Double.minus(c: Complex): Complex = Complex(this - c.re, -c.im)
operator fun Double.times(c: Complex): Complex = Complex(this * c.re, this * c.im)
operator fun Double.div(c: Complex): Complex = Complex(this, 0.0) / c

operator fun Number.plus(c: Complex): Complex = this.toDouble() + c
operator fun Number.minus(c: Complex): Complex = this.toDouble() - c
operator fun Number.times(c: Complex): Complex = this.toDouble() * c
operator fun Number.div(c: Complex): Complex = this.toDouble() / c
