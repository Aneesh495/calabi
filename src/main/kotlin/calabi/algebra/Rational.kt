package calabi.algebra

import java.math.BigDecimal
import java.math.BigInteger
import java.math.MathContext

/**
 * Rational represents an exact element of the field of rational numbers Q.
 *
 * Always maintained in canonical reduced form:
 * - denominator is strictly positive
 * - gcd(numerator, denominator) == 1
 * - zero is represented uniquely as 0/1
 */
data class Rational private constructor(
    val num: BigInteger,
    val den: BigInteger
) : Number(), Comparable<Rational> {

    operator fun plus(other: Rational): Rational {
        val n = num.multiply(other.den).add(other.num.multiply(den))
        val d = den.multiply(other.den)
        return of(n, d)
    }

    operator fun plus(other: Long): Rational = this + of(other)
    operator fun plus(other: Int): Rational = this + of(other.toLong())

    operator fun minus(other: Rational): Rational {
        val n = num.multiply(other.den).subtract(other.num.multiply(den))
        val d = den.multiply(other.den)
        return of(n, d)
    }

    operator fun minus(other: Long): Rational = this - of(other)
    operator fun minus(other: Int): Rational = this - of(other.toLong())

    operator fun unaryMinus(): Rational =
        Rational(num.negate(), den)

    operator fun times(other: Rational): Rational {
        val g1 = num.gcd(other.den)
        val g2 = den.gcd(other.num)
        val n = (num.divide(g1)).multiply(other.num.divide(g2))
        val d = (den.divide(g2)).multiply(other.den.divide(g1))
        return of(n, d)
    }

    operator fun times(other: Long): Rational = this * of(other)
    operator fun times(other: Int): Rational = this * of(other.toLong())

    operator fun div(other: Rational): Rational {
        require(other.num != BigInteger.ZERO) { "Division by rational zero" }
        return this * other.inv()
    }

    operator fun div(other: Long): Rational = this / of(other)
    operator fun div(other: Int): Rational = this / of(other.toLong())

    fun inv(): Rational {
        require(num != BigInteger.ZERO) { "Inverse of zero rational" }
        return if (num.signum() < 0) {
            Rational(den.negate(), num.negate())
        } else {
            Rational(den, num)
        }
    }

    fun abs(): Rational = if (num.signum() < 0) -this else this

    fun pow(exp: Int): Rational {
        if (exp == 0) return ONE
        if (exp < 0) return inv().pow(-exp)
        return Rational(num.pow(exp), den.pow(exp))
    }

    fun toComplex(): Complex = Complex(toDouble(), 0.0)

    fun toBigDecimal(mc: MathContext = MathContext.DECIMAL128): BigDecimal =
        BigDecimal(num).divide(BigDecimal(den), mc)

    override fun toDouble(): Double =
        num.toDouble() / den.toDouble()

    override fun toFloat(): Float = toDouble().toFloat()
    override fun toLong(): Long = num.divide(den).toLong()
    override fun toInt(): Int = toLong().toInt()
    override fun toShort(): Short = toLong().toShort()
    override fun toByte(): Byte = toLong().toByte()

    override fun compareTo(other: Rational): Int =
        num.multiply(other.den).compareTo(other.num.multiply(den))

    override fun toString(): String =
        if (den == BigInteger.ONE) num.toString() else "$num/$den"

    companion object {
        val ZERO = Rational(BigInteger.ZERO, BigInteger.ONE)
        val ONE = Rational(BigInteger.ONE, BigInteger.ONE)
        val TWO = Rational(BigInteger.valueOf(2), BigInteger.ONE)
        val HALF = Rational(BigInteger.ONE, BigInteger.valueOf(2))

        fun of(numerator: BigInteger, denominator: BigInteger): Rational {
            require(denominator != BigInteger.ZERO) { "Denominator cannot be zero" }
            if (numerator == BigInteger.ZERO) return ZERO

            var n = numerator
            var d = denominator
            if (d.signum() < 0) {
                n = n.negate()
                d = d.negate()
            }
            val gcd = n.gcd(d)
            return Rational(n.divide(gcd), d.divide(gcd))
        }

        fun of(numerator: Long, denominator: Long = 1L): Rational =
            of(BigInteger.valueOf(numerator), BigInteger.valueOf(denominator))

        fun of(numerator: Int, denominator: Int = 1): Rational =
            of(BigInteger.valueOf(numerator.toLong()), BigInteger.valueOf(denominator.toLong()))

        /**
         * Harmonic number H_n = sum_{k=1}^n 1/k.
         */
        fun harmonic(n: Int): Rational {
            require(n >= 0) { "n must be non-negative" }
            var sum = ZERO
            for (k in 1..n) {
                sum += of(1, k)
            }
            return sum
        }

        /**
         * Computes the Bernoulli numbers B_n for small to medium n.
         * B_0 = 1, B_1 = -1/2, B_2 = 1/6, B_4 = -1/30, B_6 = 1/42, ...
         */
        fun bernoulli(n: Int): Rational {
            require(n >= 0) { "Index must be non-negative" }
            val a = Array(n + 1) { ZERO }
            for (m in 0..n) {
                a[m] = of(1, m + 1)
                for (j in m downTo 1) {
                    a[j - 1] = (a[j - 1] - a[j]) * of(j)
                }
            }
            return a[0]
        }
    }
}
