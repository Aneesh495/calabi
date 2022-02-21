package calabi.manifold

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import calabi.algebra.Polynomial
import calabi.algebra.ProjectiveSpace

/**
 * CompleteIntersectionCalabiYau represents a Calabi-Yau manifold defined by the simultaneous vanishing
 * of a system of transverse homogeneous polynomials in projective space:
 *
 *   X = { z in P^N | P_1(z) = 0, ..., P_k(z) = 0 }
 *
 * The Calabi-Yau condition c_1(X) = 0 requires sum_{j=1}^k deg(P_j) = N + 1.
 * Complex dimension: dim_C(X) = N - k.
 */
class CompleteIntersectionCalabiYau(
    override val name: String,
    val ambientN: Int,
    override val definingPolynomials: List<Polynomial>,
    override val hodgeNumbers: Pair<Int, Int>,
    override val eulerCharacteristic: Int
) : CalabiYauManifold {

    override val ambientSpace: ProjectiveSpace = ProjectiveSpace(ambientN)
    override val complexDimension: Int = ambientN - definingPolynomials.size

    init {
        require(definingPolynomials.isNotEmpty()) { "Complete intersection requires at least one polynomial" }
        require(complexDimension >= 1) { "Manifold dimension must be at least 1" }
        for (p in definingPolynomials) {
            require(p.isHomogeneous) { "All defining polynomials must be homogeneous" }
            require(p.numVariables == ambientN + 1) {
                "Polynomial variable count (${p.numVariables}) must equal ambient dimension + 1 (${ambientN + 1})"
            }
        }
        val totalDegree = definingPolynomials.sumOf { it.degree }
        require(totalDegree == ambientN + 1) {
            "Calabi-Yau condition c_1(X) = 0 requires total degree ($totalDegree) == ambient N + 1 (${ambientN + 1})"
        }
    }

    /**
     * Evaluates the Jacobian matrix J_{ij} = del P_i / del z_j at a point.
     * Transversality requires that J has maximal rank k at all points on X.
     */
    fun jacobian(point: HomogeneousPoint): Array<Array<Complex>> {
        val k = definingPolynomials.size
        val nVars = ambientN + 1
        val jac = Array(k) { Array(nVars) { Complex.ZERO } }
        for (i in 0 until k) {
            val grad = definingPolynomials[i].gradient(point.coordinates)
            for (j in 0 until nVars) {
                jac[i][j] = grad[j]
            }
        }
        return jac
    }

    companion object {
        /**
         * The standard codimension-2 Calabi-Yau threefold X_{3,3} in P^5:
         * intersection of two cubic hypersurfaces in P^5.
         * N = 5, k = 2, total degree 3 + 3 = 6 = N + 1.
         * dim_C = 5 - 2 = 3.
         * Hodge numbers: h^{1,1} = 1, h^{2,1} = 73, chi = 2 * (1 - 73) = -144.
         */
        fun bicubicInP5(): CompleteIntersectionCalabiYau {
            val p1 = Polynomial.fermat(6, 3)
            // Second cubic with shifted variables
            val terms = HashMap<calabi.algebra.MultiIndex, Complex>()
            for (i in 0..5) {
                val exps = IntArray(6)
                exps[(i + 1) % 6] = 2
                exps[i] = 1
                terms[calabi.algebra.MultiIndex(exps)] = Complex.ONE
            }
            val p2 = Polynomial(6, terms)
            return CompleteIntersectionCalabiYau(
                name = "Bicubic Threefold X_{3,3} in P^5",
                ambientN = 5,
                definingPolynomials = listOf(p1, p2),
                hodgeNumbers = Pair(1, 73),
                eulerCharacteristic = -144
            )
        }
    }
}
