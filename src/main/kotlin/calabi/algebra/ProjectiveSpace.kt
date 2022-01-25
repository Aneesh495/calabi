package calabi.algebra

import kotlin.math.ln
import kotlin.math.max
import kotlin.random.Random

/**
 * HomogeneousPoint represents an equivalence class [z_0 : z_1 : ... : z_n] in P^n(C).
 *
 * Two points are equivalent under [z] ~ [lambda * z] for any non-zero complex scalar lambda.
 */
data class HomogeneousPoint(
    val coordinates: List<Complex>
) {
    val dimension: Int = coordinates.size - 1

    init {
        require(coordinates.size >= 2) { "Projective space dimension must be at least 1 (at least 2 coordinates)" }
        require(!coordinates.all { it.isZero(1e-15) }) { "Homogeneous coordinates cannot be identically zero" }
    }

    operator fun get(index: Int): Complex = coordinates[index]

    /**
     * Normalizes homogeneous coordinates such that the Euclidean norm is 1: sum |z_i|^2 = 1.
     * Also fixes the global phase by making the first non-zero coordinate strictly positive real.
     */
    fun normalized(): HomogeneousPoint {
        val sumSq = coordinates.sumOf { it.normSq() }
        val norm = kotlin.math.sqrt(sumSq)
        val firstNonZero = coordinates.firstOrNull { !it.isZero(1e-12) } ?: Complex.ONE
        val phase = firstNonZero.arg()
        val phaseRot = Complex.fromPolar(1.0, -phase)

        val scaled = coordinates.map { (it * phaseRot) / norm }
        return HomogeneousPoint(scaled)
    }

    /**
     * Determines the optimal affine chart U_k = { [z] in P^n | z_k != 0 },
     * selected by maximizing |z_k| to ensure numerical conditioning.
     */
    fun dominantChartIndex(): Int {
        var maxIdx = 0
        var maxMod = coordinates[0].normSq()
        for (i in 1 until coordinates.size) {
            val mod = coordinates[i].normSq()
            if (mod > maxMod) {
                maxMod = mod
                maxIdx = i
            }
        }
        return maxIdx
    }

    /**
     * Converts to local affine coordinates on chart U_k:
     * w_j = z_j / z_k for j != k.
     */
    fun toAffineCoordinates(chartIndex: Int = dominantChartIndex()): AffineCoordinates {
        val zk = coordinates[chartIndex]
        require(!zk.isZero(1e-14)) { "Cannot construct affine chart at zero coordinate: index $chartIndex" }
        val affineList = ArrayList<Complex>(dimension)
        for (i in coordinates.indices) {
            if (i != chartIndex) {
                affineList.add(coordinates[i] / zk)
            }
        }
        return AffineCoordinates(chartIndex, dimension, affineList)
    }

    /**
     * Computes the Fubini-Study potential K_FS = ln(sum_{i=0}^n |z_i|^2)
     */
    fun fubiniStudyPotential(): Double {
        val sumSq = coordinates.sumOf { it.normSq() }
        return ln(sumSq)
    }

    /**
     * Distance in the Fubini-Study metric:
     * d_FS([p], [q]) = arccos( |<p, q>| / (||p|| * ||q||) )
     */
    fun fubiniStudyDistance(other: HomogeneousPoint): Double {
        require(dimension == other.dimension) { "Dimension mismatch: $dimension vs ${other.dimension}" }
        var innerProduct = Complex.ZERO
        var normP = 0.0
        var normQ = 0.0
        for (i in 0..dimension) {
            innerProduct += coordinates[i].conj() * other.coordinates[i]
            normP += coordinates[i].normSq()
            normQ += other.coordinates[i].normSq()
        }
        val cosTheta = (innerProduct.abs() / (kotlin.math.sqrt(normP) * kotlin.math.sqrt(normQ))).coerceIn(0.0, 1.0)
        return kotlin.math.acos(cosTheta)
    }

    companion object {
        fun of(vararg coords: Complex): HomogeneousPoint = HomogeneousPoint(coords.toList())
        fun ofDoubles(vararg coords: Double): HomogeneousPoint =
            HomogeneousPoint(coords.map { Complex.ofReal(it) })
    }
}

/**
 * AffineCoordinates represents a point on the standard affine chart U_k = { z_k != 0 }.
 *
 * Coordinates are w = (w_1, ..., w_n) where w_j = z_j / z_k for j != k.
 */
data class AffineCoordinates(
    val chartIndex: Int,
    val dimension: Int,
    val coords: List<Complex>
) {
    init {
        require(coords.size == dimension) {
            "Coordinate vector size (${coords.size}) must match affine dimension ($dimension)"
        }
    }

    operator fun get(index: Int): Complex = coords[index]

    /**
     * Converts back to homogeneous coordinates [w_0 : ... : 1 : ... : w_n] with 1 at chartIndex.
     */
    fun toHomogeneous(): HomogeneousPoint {
        val list = ArrayList<Complex>(dimension + 1)
        var affineCursor = 0
        for (i in 0..dimension) {
            if (i == chartIndex) {
                list.add(Complex.ONE)
            } else {
                list.add(coords[affineCursor++])
            }
        }
        return HomogeneousPoint(list)
    }

    /**
     * Transition function from chart U_k to chart U_m.
     */
    fun transitionTo(targetChartIndex: Int): AffineCoordinates {
        if (targetChartIndex == chartIndex) return this
        val homo = toHomogeneous()
        return homo.toAffineCoordinates(targetChartIndex)
    }

    /**
     * Evaluates the local Fubini-Study metric tensor g_{i bar{j}} on this affine chart:
     * g_{i bar{j}} = ((1 + ||w||^2) * delta_{ij} - bar{w}_i * w_j) / (1 + ||w||^2)^2
     */
    fun fubiniStudyMetric(): Array<Array<Complex>> {
        val normSq = coords.sumOf { it.normSq() }
        val denom = (1.0 + normSq) * (1.0 + normSq)
        val matrix = Array(dimension) { Array(dimension) { Complex.ZERO } }

        for (i in 0 until dimension) {
            for (j in 0 until dimension) {
                val delta = if (i == j) (1.0 + normSq) else 0.0
                val wiWjBar = coords[i] * coords[j].conj()
                val elem = (Complex(delta, 0.0) - wiWjBar) / denom
                matrix[i][j] = elem
            }
        }
        return matrix
    }

    /**
     * Fubini-Study volume determinant:
     * det(g_FS) = 1 / (1 + ||w||^2)^(n + 1)
     */
    fun fubiniStudyVolumeDensity(): Double {
        val normSq = coords.sumOf { it.normSq() }
        return 1.0 / (1.0 + normSq).let { v ->
            var p = 1.0
            repeat(dimension + 1) { p *= v }
            p
        }
    }
}

/**
 * ProjectiveSpace represents P^n(C) and provides geometric sampling and coordinate operations.
 */
class ProjectiveSpace(
    val dimension: Int
) {
    init {
        require(dimension >= 1) { "Projective dimension must be at least 1" }
    }

    /**
     * Samples a point uniformly at random from P^n(C) according to the Haar/Fubini-Study measure.
     * Accomplished by sampling standard Gaussian complex coordinates and normalizing.
     */
    fun sampleUniformPoint(random: Random = Random.Default): HomogeneousPoint {
        val coords = ArrayList<Complex>(dimension + 1)
        for (i in 0..dimension) {
            val u1 = max(1e-15, random.nextDouble())
            val u2 = random.nextDouble()
            val r = kotlin.math.sqrt(-2.0 * kotlin.math.ln(u1))
            val theta = 2.0 * kotlin.math.PI * u2
            val re = r * kotlin.math.cos(theta)

            val v1 = max(1e-15, random.nextDouble())
            val v2 = random.nextDouble()
            val r2 = kotlin.math.sqrt(-2.0 * kotlin.math.ln(v1))
            val theta2 = 2.0 * kotlin.math.PI * v2
            val im = r2 * kotlin.math.cos(theta2)

            coords.add(Complex(re, im))
        }
        return HomogeneousPoint(coords).normalized()
    }

    /**
     * Samples a projective line L = span(p, q) in P^n(C), representing a point in Gr(2, n+1).
     * The two spanning points p and q are chosen orthonormal with respect to the standard Hermitian inner product.
     */
    fun sampleProjectiveLine(random: Random = Random.Default): Pair<HomogeneousPoint, HomogeneousPoint> {
        val p = sampleUniformPoint(random)
        val qRaw = sampleUniformPoint(random)

        // Gram-Schmidt orthogonalization: q = qRaw - <p, qRaw> * p
        var inner = Complex.ZERO
        for (i in 0..dimension) {
            inner += p[i].conj() * qRaw[i]
        }
        val qOrth = (0..dimension).map { i ->
            qRaw[i] - inner * p[i]
        }
        val q = HomogeneousPoint(qOrth).normalized()
        return Pair(p, q)
    }
}
