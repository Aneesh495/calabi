package calabi.sampling

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import calabi.algebra.Polynomial
import calabi.algebra.ProjectiveSpace
import kotlin.random.Random

/**
 * SampledPoint represents a point lying on the algebraic hypersurface together with
 * its local invariant measure integration weight:
 *
 *   weight = dmu_CY / dmu_FS
 */
data class SampledPoint(
    val point: HomogeneousPoint,
    val weight: Double
)

/**
 * LineIntersectionSampler generates points on a Calabi-Yau hypersurface X = {P = 0} in P^{n+1}
 * by sampling complex lines uniformly according to the Haar measure on Gr(2, n+2) and finding
 * all intersection points via Durand-Kerner polynomial root finding.
 */
class LineIntersectionSampler(
    val polynomial: Polynomial,
    val ambientSpace: ProjectiveSpace
) {
    init {
        require(polynomial.isHomogeneous) { "Polynomial must be homogeneous" }
        require(polynomial.numVariables == ambientSpace.dimension + 1) {
            "Variable count must match projective coordinates"
        }
    }

    /**
     * Samples intersection points from a single random complex projective line.
     * Yields deg(P) intersection points on the hypersurface.
     */
    fun sampleLineIntersections(random: Random = Random.Default): List<HomogeneousPoint> {
        val (p, q) = ambientSpace.sampleProjectiveLine(random)
        val coeffs = polynomial.restrictToProjectiveLine(p, q)
        val roots = DurandKernerRoots.findRoots(coeffs)

        val points = ArrayList<HomogeneousPoint>(roots.size)
        val numVars = polynomial.numVariables
        for (root in roots) {
            val coords = (0 until numVars).map { i ->
                p[i] + root * q[i]
            }
            val hp = HomogeneousPoint(coords).normalized()
            // Verify point is on hypersurface within tolerance
            val residual = polynomial.evaluate(hp).abs()
            if (residual < 1e-6) {
                points.add(hp)
            }
        }
        return points
    }

    /**
     * Generates a batch of at least targetCount valid points on the hypersurface.
     */
    fun samplePoints(targetCount: Int, random: Random = Random.Default): List<HomogeneousPoint> {
        val result = ArrayList<HomogeneousPoint>(targetCount)
        while (result.size < targetCount) {
            val linePoints = sampleLineIntersections(random)
            for (pt in linePoints) {
                result.add(pt)
                if (result.size >= targetCount) break
            }
        }
        return result
    }
}
