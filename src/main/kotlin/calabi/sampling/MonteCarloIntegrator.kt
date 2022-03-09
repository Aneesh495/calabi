package calabi.sampling

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import calabi.geometry.PoincareResidue
import calabi.manifold.CalabiYauManifold
import calabi.manifold.QuinticThreefold
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.runBlocking
import kotlin.random.Random

/**
 * MonteCarloIntegrator performs numerical integration over a Calabi-Yau manifold
 * with respect to the unique Calabi-Yau volume form dmu_{CY} = Omega ^ bar{Omega}.
 *
 * Employs importance sampling via the Fubini-Study measure on Gr(2, n+2) with
 * Radon-Nikodym weight corrections:
 *
 *   int_X f(x) dmu_{CY} = (1 / N) * sum_{i=1}^N f(x_i) * w(x_i)
 */
class MonteCarloIntegrator(
    val manifold: QuinticThreefold,
    val seed: Long = 42L
) {
    val sampler = LineIntersectionSampler(manifold.polynomial, manifold.ambientSpace)
    val poincareResidue = manifold.poincareResidue

    /**
     * Generates a collection of points equipped with local Calabi-Yau Radon-Nikodym weights.
     */
    fun sampleWeightedPoints(count: Int): List<SampledPoint> {
        val rawPoints = sampler.samplePoints(count, Random(seed))
        val weighted = ArrayList<SampledPoint>(rawPoints.size)
        for (pt in rawPoints) {
            val affine = pt.toAffineCoordinates()
            val w = poincareResidue.localMeasureWeight(affine)
            if (w.isFinite() && w > 0.0) {
                weighted.add(SampledPoint(pt, w))
            }
        }
        return weighted
    }

    /**
     * Parallel sampling of weighted points using Kotlin Coroutines across all available CPU cores.
     */
    fun sampleWeightedPointsParallel(count: Int, numChunks: Int = Runtime.getRuntime().availableProcessors()): List<SampledPoint> = runBlocking {
        val perChunk = (count + numChunks - 1) / numChunks
        val deferreds = (0 until numChunks).map { chunkIdx ->
            async(Dispatchers.Default) {
                val chunkSampler = LineIntersectionSampler(manifold.polynomial, manifold.ambientSpace)
                val rng = Random(seed + chunkIdx * 10007L)
                val pts = chunkSampler.samplePoints(perChunk, rng)
                val list = ArrayList<SampledPoint>(pts.size)
                for (pt in pts) {
                    val affine = pt.toAffineCoordinates()
                    val w = poincareResidue.localMeasureWeight(affine)
                    if (w.isFinite() && w > 0.0) {
                        list.add(SampledPoint(pt, w))
                    }
                }
                list
            }
        }
        deferreds.awaitAll().flatten().take(count)
    }

    /**
     * Integrates a scalar complex function f: X -> C over the Calabi-Yau manifold.
     */
    fun integrate(
        points: List<SampledPoint>,
        f: (HomogeneousPoint) -> Complex
    ): Complex {
        require(points.isNotEmpty()) { "Point set cannot be empty for integration" }
        var sumF = Complex.ZERO
        var totalWeight = 0.0
        for (sp in points) {
            val fVal = f(sp.point)
            sumF += fVal * sp.weight
            totalWeight += sp.weight
        }
        return sumF / totalWeight
    }

    /**
     * Total volume of the manifold estimated via Monte Carlo:
     * Vol(X) = (1/N) * sum_{i=1}^N w(x_i).
     */
    fun estimateVolume(points: List<SampledPoint>): Double {
        if (points.isEmpty()) return 0.0
        return points.sumOf { it.weight } / points.size
    }
}
