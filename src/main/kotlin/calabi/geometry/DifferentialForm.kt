package calabi.geometry

import calabi.algebra.Complex

/**
 * FormBasisIndex represents a basis element dz_I wedge dbar{z}_J.
 * Represented compactly via bitmasks:
 * - holomorphicMask: bits set for indices in I
 * - antiHolomorphicMask: bits set for indices in J
 */
data class FormBasisIndex(
    val holomorphicMask: Int,
    val antiHolomorphicMask: Int
) : Comparable<FormBasisIndex> {

    val p: Int = Integer.bitCount(holomorphicMask)
    val q: Int = Integer.bitCount(antiHolomorphicMask)

    override fun compareTo(other: FormBasisIndex): Int {
        val pComp = p.compareTo(other.p)
        if (pComp != 0) return pComp
        val qComp = q.compareTo(other.q)
        if (qComp != 0) return qComp
        val hComp = holomorphicMask.compareTo(other.holomorphicMask)
        if (hComp != 0) return hComp
        return antiHolomorphicMask.compareTo(other.antiHolomorphicMask)
    }

    override fun toString(): String {
        val dzList = maskToIndices(holomorphicMask).map { "dz_$it" }
        val dbarZList = maskToIndices(antiHolomorphicMask).map { "dbar{z}_$it" }
        val all = dzList + dbarZList
        return if (all.isEmpty()) "1" else all.joinToString(" ^ ")
    }

    companion object {
        val SCALAR = FormBasisIndex(0, 0)

        fun maskToIndices(mask: Int): List<Int> {
            val list = ArrayList<Int>()
            var m = mask
            var idx = 0
            while (m > 0) {
                if ((m and 1) == 1) list.add(idx)
                m = m ushr 1
                idx++
            }
            return list
        }

        fun indicesToMask(indices: Collection<Int>): Int {
            var mask = 0
            for (i in indices) {
                require(i in 0..30) { "Index out of range for 32-bit mask: $i" }
                val bit = 1 shl i
                require((mask and bit) == 0) { "Duplicate index in basis form: $i" }
                mask = mask or bit
            }
            return mask
        }

        /**
         * Computes the wedge product of two holomorphic or anti-holomorphic index masks:
         * maskA wedge maskB.
         * Returns Pair(sign, unionMask) or null if they share any index (dz_i wedge dz_i = 0).
         */
        fun wedgeMasks(maskA: Int, maskB: Int): Pair<Int, Int>? {
            if ((maskA and maskB) != 0) return null // Shared coordinate vanishes
            val union = maskA or maskB

            var inversions = 0
            var tempA = maskA
            var idxA = 0
            while (tempA > 0) {
                if ((tempA and 1) == 1) {
                    val smallerInB = Integer.bitCount(maskB and ((1 shl idxA) - 1))
                    val largerInB = Integer.bitCount(maskB) - smallerInB
                    inversions += largerInB
                }
                tempA = tempA ushr 1
                idxA++
            }
            val sign = if ((inversions and 1) == 0) 1 else -1
            return Pair(sign, union)
        }
    }
}

/**
 * DifferentialForm represents a differential (p, q)-form on a complex manifold.
 * omega = sum_{I, J} c_{I, J} dz_I wedge dbar{z}_J
 */
data class DifferentialForm(
    val manifoldDimension: Int,
    val p: Int,
    val q: Int,
    val components: Map<FormBasisIndex, Complex>
) {
    init {
        for (k in components.keys) {
            require(k.p == p && k.q == q) {
                "Heterogeneous degree in differential form: expected ($p, $q) but found (${k.p}, ${k.q})"
            }
        }
    }

    operator fun plus(other: DifferentialForm): DifferentialForm {
        require(manifoldDimension == other.manifoldDimension) { "Manifold dimension mismatch" }
        if (components.isEmpty()) return other
        if (other.components.isEmpty()) return this
        require(p == other.p && q == other.q) { "Form degree mismatch: ($p, $q) vs (${other.p}, ${other.q})" }
        val newComps = HashMap(components)
        for ((basis, coeff) in other.components) {
            val cur = newComps[basis] ?: Complex.ZERO
            val updated = cur + coeff
            if (updated.isZero(1e-15)) {
                newComps.remove(basis)
            } else {
                newComps[basis] = updated
            }
        }
        return DifferentialForm(manifoldDimension, p, q, newComps)
    }

    operator fun minus(other: DifferentialForm): DifferentialForm =
        this + (-other)

    operator fun unaryMinus(): DifferentialForm =
        DifferentialForm(manifoldDimension, p, q, components.mapValues { -it.value })

    operator fun times(scalar: Complex): DifferentialForm {
        if (scalar.isZero(1e-15)) return zero(manifoldDimension, p, q)
        return DifferentialForm(manifoldDimension, p, q, components.mapValues { it.value * scalar })
    }

    operator fun times(scalar: Double): DifferentialForm =
        this * Complex.ofReal(scalar)

    /**
     * Exterior wedge product: omega_A wedge omega_B.
     * Yields a form of degree (p_A + p_B, q_A + q_B).
     */
    fun wedge(other: DifferentialForm): DifferentialForm {
        require(manifoldDimension == other.manifoldDimension) { "Manifold dimension mismatch" }
        val newP = p + other.p
        val newQ = q + other.q
        if (newP > manifoldDimension || newQ > manifoldDimension) {
            return zero(manifoldDimension, newP.coerceAtMost(manifoldDimension), newQ.coerceAtMost(manifoldDimension))
        }

        val newComps = HashMap<FormBasisIndex, Complex>()

        for ((basisA, coeffA) in components) {
            for ((basisB, coeffB) in other.components) {
                val wedgedH = FormBasisIndex.wedgeMasks(basisA.holomorphicMask, basisB.holomorphicMask) ?: continue
                val wedgedAH = FormBasisIndex.wedgeMasks(basisA.antiHolomorphicMask, basisB.antiHolomorphicMask) ?: continue

                val commuteParity = (basisA.q * basisB.p) and 1
                val commuteSign = if (commuteParity == 0) 1 else -1

                val totalSign = wedgedH.first * wedgedAH.first * commuteSign
                val targetBasis = FormBasisIndex(wedgedH.second, wedgedAH.second)
                val targetCoeff = coeffA * coeffB * totalSign.toDouble()

                val cur = newComps[targetBasis] ?: Complex.ZERO
                val updated = cur + targetCoeff
                if (updated.isZero(1e-15)) {
                    newComps.remove(targetBasis)
                } else {
                    newComps[targetBasis] = updated
                }
            }
        }

        return DifferentialForm(manifoldDimension, newP, newQ, newComps)
    }

    /**
     * Complex conjugation of differential form:
     * bar(c * dz_I wedge dbar{z}_J) = bar(c) * (-1)^(p*q) * dz_J wedge dbar{z}_I
     */
    fun conj(): DifferentialForm {
        val newComps = HashMap<FormBasisIndex, Complex>()
        val sign = if (((p * q) and 1) == 0) 1.0 else -1.0
        for ((basis, coeff) in components) {
            val swapped = FormBasisIndex(basis.antiHolomorphicMask, basis.holomorphicMask)
            newComps[swapped] = coeff.conj() * sign
        }
        return DifferentialForm(manifoldDimension, q, p, newComps)
    }

    companion object {
        fun zero(dim: Int, p: Int = 0, q: Int = 0): DifferentialForm =
            DifferentialForm(dim, p, q, emptyMap())

        fun scalar(dim: Int, value: Complex): DifferentialForm =
            DifferentialForm(dim, 0, 0, mapOf(FormBasisIndex.SCALAR to value))

        fun dz(dim: Int, index: Int): DifferentialForm {
            require(index in 0 until dim) { "Index $index out of bounds for dim $dim" }
            val basis = FormBasisIndex(1 shl index, 0)
            return DifferentialForm(dim, 1, 0, mapOf(basis to Complex.ONE))
        }

        fun dbarZ(dim: Int, index: Int): DifferentialForm {
            require(index in 0 until dim) { "Index $index out of bounds for dim $dim" }
            val basis = FormBasisIndex(0, 1 shl index)
            return DifferentialForm(dim, 0, 1, mapOf(basis to Complex.ONE))
        }

        /**
         * Kähler form omega_0 = i/2 * sum_{k=0}^{n-1} dz_k wedge dbar{z}_k
         */
        fun standardKahlerForm(dim: Int): DifferentialForm {
            var form = zero(dim, 1, 1)
            val coeff = Complex.I * 0.5
            for (k in 0 until dim) {
                val term = DifferentialForm(dim, 1, 1, mapOf(FormBasisIndex(1 shl k, 1 shl k) to coeff))
                form += term
            }
            return form
        }
    }
}
