package calabi.algebra

import java.util.*

/**
 * MultiIndex represents an exponent vector alpha = (alpha_0, alpha_1, ..., alpha_n) for monomials z^alpha.
 */
data class MultiIndex(
    val exponents: IntArray
) : Comparable<MultiIndex> {

    val numVariables: Int = exponents.size
    val totalDegree: Int = exponents.sum()

    init {
        for (e in exponents) {
            require(e >= 0) { "Exponents must be non-negative integers: found $e" }
        }
    }

    operator fun get(index: Int): Int = exponents[index]

    /**
     * Partial derivative of the exponent vector with respect to variable at index j.
     * Returns null if exponent at j is zero.
     */
    fun derive(j: Int): Pair<Int, MultiIndex>? {
        if (exponents[j] == 0) return null
        val newExps = exponents.clone()
        val coeff = newExps[j]
        newExps[j]--
        return Pair(coeff, MultiIndex(newExps))
    }

    /**
     * Graded reverse lexicographic (grevlex) monomial ordering:
     * 1. Compare total degrees
     * 2. If equal, compare exponents from right to left, reversing the outcome
     */
    override fun compareTo(other: MultiIndex): Int {
        val degComp = totalDegree.compareTo(other.totalDegree)
        if (degComp != 0) return degComp
        for (i in (numVariables - 1) downTo 0) {
            val c = exponents[i].compareTo(other.exponents[i])
            if (c != 0) return -c
        }
        return 0
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is MultiIndex) return false
        return exponents.contentEquals(other.exponents)
    }

    override fun hashCode(): Int = exponents.contentHashCode()

    override fun toString(): String =
        exponents.joinToString(prefix = "(", postfix = ")", separator = ", ")

    companion object {
        fun of(vararg exps: Int): MultiIndex = MultiIndex(exps)

        /**
         * Generates all multi-indices of total degree d in numVars variables.
         * The total number of such multi-indices is binom(numVars + d - 1, d).
         */
        fun generateDegreeBasis(numVars: Int, degree: Int): List<MultiIndex> {
            val result = ArrayList<MultiIndex>()
            val current = IntArray(numVars)

            fun backtrack(varIndex: Int, remainingDegree: Int) {
                if (varIndex == numVars - 1) {
                    current[varIndex] = remainingDegree
                    result.add(MultiIndex(current.clone()))
                    return
                }
                for (deg in 0..remainingDegree) {
                    current[varIndex] = deg
                    backtrack(varIndex + 1, remainingDegree - deg)
                }
            }

            backtrack(0, degree)
            result.sort()
            return result
        }
    }
}

/**
 * Monomial represents a term c * z^alpha where c is a non-zero complex coefficient.
 */
data class Monomial(
    val coefficient: Complex,
    val multiIndex: MultiIndex
) {
    init {
        require(!coefficient.isZero(1e-15)) { "Monomial coefficient cannot be zero" }
    }

    val degree: Int = multiIndex.totalDegree

    fun evaluate(z: List<Complex>): Complex {
        var product = coefficient
        for (i in 0 until multiIndex.numVariables) {
            val exp = multiIndex[i]
            if (exp > 0) {
                product *= z[i].pow(exp)
            }
        }
        return product
    }

    fun derive(variableIndex: Int): Monomial? {
        val derived = multiIndex.derive(variableIndex) ?: return null
        val newCoeff = coefficient * derived.first.toDouble()
        if (newCoeff.isZero(1e-15)) return null
        return Monomial(newCoeff, derived.second)
    }
}

/**
 * Polynomial represents a multivariate polynomial over C[z_0, ..., z_n].
 */
data class Polynomial(
    val numVariables: Int,
    val terms: Map<MultiIndex, Complex>
) {
    val degree: Int = terms.keys.maxOfOrNull { it.totalDegree } ?: 0
    val isHomogeneous: Boolean = terms.isEmpty() || terms.keys.map { it.totalDegree }.distinct().size == 1

    operator fun plus(other: Polynomial): Polynomial {
        require(numVariables == other.numVariables) { "Variable count mismatch" }
        val newTerms = HashMap(terms)
        for ((idx, coeff) in other.terms) {
            val current = newTerms[idx] ?: Complex.ZERO
            val updated = current + coeff
            if (updated.isZero(1e-15)) {
                newTerms.remove(idx)
            } else {
                newTerms[idx] = updated
            }
        }
        return Polynomial(numVariables, newTerms)
    }

    operator fun minus(other: Polynomial): Polynomial =
        this + (-other)

    operator fun unaryMinus(): Polynomial =
        Polynomial(numVariables, terms.mapValues { -it.value })

    operator fun times(scalar: Complex): Polynomial {
        if (scalar.isZero(1e-15)) return zero(numVariables)
        val newTerms = terms.mapValues { it.value * scalar }.filterValues { !it.isZero(1e-15) }
        return Polynomial(numVariables, newTerms)
    }

    operator fun times(scalar: Double): Polynomial =
        this * Complex.ofReal(scalar)

    operator fun times(other: Polynomial): Polynomial {
        require(numVariables == other.numVariables) { "Variable count mismatch" }
        val newTerms = HashMap<MultiIndex, Complex>()
        for ((idxA, coeffA) in terms) {
            for ((idxB, coeffB) in other.terms) {
                val newExps = IntArray(numVariables) { i -> idxA[i] + idxB[i] }
                val newIdx = MultiIndex(newExps)
                val prod = coeffA * coeffB
                val current = newTerms[newIdx] ?: Complex.ZERO
                newTerms[newIdx] = current + prod
            }
        }
        val clean = newTerms.filterValues { !it.isZero(1e-15) }
        return Polynomial(numVariables, clean)
    }

    /**
     * Evaluates the polynomial at complex coordinate vector z.
     */
    fun evaluate(z: List<Complex>): Complex {
        require(z.size == numVariables) { "Coordinate size ${z.size} != numVariables $numVariables" }
        var sum = Complex.ZERO
        for ((multiIndex, coeff) in terms) {
            var termVal = coeff
            for (i in 0 until numVariables) {
                val exp = multiIndex[i]
                if (exp > 0) {
                    termVal *= z[i].pow(exp)
                }
            }
            sum += termVal
        }
        return sum
    }

    fun evaluate(pt: HomogeneousPoint): Complex =
        evaluate(pt.coordinates)

    /**
     * Evaluates the complex gradient nabla P(z) = (del P / del z_0, ..., del P / del z_n).
     */
    fun gradient(z: List<Complex>): List<Complex> {
        val grad = ArrayList<Complex>(numVariables)
        for (i in 0 until numVariables) {
            val partial = partialDerivative(i)
            grad.add(partial.evaluate(z))
        }
        return grad
    }

    /**
     * Partial derivative polynomial del P / del z_j.
     */
    fun partialDerivative(j: Int): Polynomial {
        val newTerms = HashMap<MultiIndex, Complex>()
        for ((idx, coeff) in terms) {
            val derived = idx.derive(j) ?: continue
            val newCoeff = coeff * derived.first.toDouble()
            if (!newCoeff.isZero(1e-15)) {
                val current = newTerms[derived.second] ?: Complex.ZERO
                newTerms[derived.second] = current + newCoeff
            }
        }
        return Polynomial(numVariables, newTerms)
    }

    /**
     * Hessian matrix H_{ij} = del^2 P / del z_i del z_j evaluated at z.
     */
    fun hessian(z: List<Complex>): Array<Array<Complex>> {
        val matrix = Array(numVariables) { Array(numVariables) { Complex.ZERO } }
        for (i in 0 until numVariables) {
            val pi = partialDerivative(i)
            for (j in 0 until numVariables) {
                val pij = pi.partialDerivative(j)
                matrix[i][j] = pij.evaluate(z)
            }
        }
        return matrix
    }

    /**
     * Restricts the homogeneous polynomial P to a projective line parametrized by
     * z(s, t) = s * p + t * q for spanning points p, q in P^{numVars-1}.
     * Returns coefficients [a_0, a_1, ..., a_d] of the univariate polynomial in t with s=1:
     * Q(t) = P(p + t * q) = sum_{k=0}^d a_k * t^k.
     */
    fun restrictToProjectiveLine(p: HomogeneousPoint, q: HomogeneousPoint): List<Complex> {
        require(isHomogeneous) { "Line restriction requires a homogeneous polynomial" }
        val d = degree
        val coeffs = Array(d + 1) { Complex.ZERO }

        // Expand each monomial term: coeff * prod_{i=0}^{numVars-1} (p_i + t * q_i)^{alpha_i}
        for ((multiIndex, coeff) in terms) {
            // Compute univariate expansion of product of linear factors
            var univariatePoly = Array(1) { Complex.ONE } // degree 0: 1.0

            for (i in 0 until numVariables) {
                val alpha = multiIndex[i]
                if (alpha == 0) continue

                // Expand (p_i + t * q_i)^alpha using binomial theorem
                val factorPoly = Array(alpha + 1) { Complex.ZERO }
                for (k in 0..alpha) {
                    val binom = Rational.binomial(alpha, k).toDouble()
                    val pPow = p[i].pow(alpha - k)
                    val qPow = q[i].pow(k)
                    factorPoly[k] = binom * pPow * qPow
                }

                // Multiply univariatePoly by factorPoly
                val newDeg = univariatePoly.size - 1 + alpha
                val nextPoly = Array(newDeg + 1) { Complex.ZERO }
                for (degA in univariatePoly.indices) {
                    for (degB in factorPoly.indices) {
                        nextPoly[degA + degB] += univariatePoly[degA] * factorPoly[degB]
                    }
                }
                univariatePoly = nextPoly
            }

            for (k in univariatePoly.indices) {
                if (k <= d) {
                    coeffs[k] += coeff * univariatePoly[k]
                }
            }
        }
        return coeffs.toList()
    }

    companion object {
        fun zero(numVariables: Int): Polynomial =
            Polynomial(numVariables, emptyMap())

        fun monomial(coeff: Complex, multiIndex: MultiIndex): Polynomial =
            Polynomial(multiIndex.numVariables, mapOf(multiIndex to coeff))

        /**
         * Fermat hypersurface polynomial of degree d in n variables: sum_{i=0}^{n-1} z_i^d = 0.
         */
        fun fermat(numVariables: Int, degree: Int): Polynomial {
            val terms = HashMap<MultiIndex, Complex>()
            for (i in 0 until numVariables) {
                val exps = IntArray(numVariables)
                exps[i] = degree
                terms[MultiIndex(exps)] = Complex.ONE
            }
            return Polynomial(numVariables, terms)
        }

        /**
         * Dwork pencil quintic threefold: sum_{i=0}^4 z_i^5 + psi * prod_{i=0}^4 z_i = 0.
         */
        fun dworkQuintic(psi: Complex): Polynomial {
            val terms = HashMap<MultiIndex, Complex>()
            // Fermat sum: z_0^5 + z_1^5 + z_2^5 + z_3^5 + z_4^5
            for (i in 0..4) {
                val exps = IntArray(5)
                exps[i] = 5
                terms[MultiIndex(exps)] = Complex.ONE
            }
            // Product term: psi * z_0 * z_1 * z_2 * z_3 * z_4
            if (!psi.isZero(1e-15)) {
                val prodExps = IntArray(5) { 1 }
                terms[MultiIndex(prodExps)] = psi
            }
            return Polynomial(5, terms)
        }
    }
}
