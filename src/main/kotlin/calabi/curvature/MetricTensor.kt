package calabi.curvature

import calabi.algebra.AffineCoordinates
import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import calabi.algebra.Polynomial
import calabi.linalg.Cholesky
import calabi.linalg.ComplexMatrix
import calabi.metric.SectionBundle

/**
 * MetricTensor computes the induced local Kahler metric tensor g_{i bar{j}} = del_i dbar_j K
 * and its determinant det(g) on the tangent space of the Calabi-Yau hypersurface.
 */
class MetricTensor(
    val bundle: SectionBundle,
    val definingPolynomial: Polynomial
) {
    val k: Int = bundle.k
    val dim: Int = bundle.dim
    val ambientDim: Int = bundle.basis.ambientDimension // e.g. 4 for P^4
    val manifoldDim: Int = ambientDim - 1 // e.g. 3 for threefold

    /**
     * Evaluates the n x n complex Hermitian metric tensor g_{i bar{j}} on the intrinsic
     * tangent space of the Calabi-Yau manifold at the given affine point.
     *
     * Uses analytic differentiation of sections:
     *   del_i s_alpha = (del s_alpha / del u_i) + (del s_alpha / del w_m) * (del w_m / del u_i)
     * where del w_m / del u_i = - (del P / del u_i) / (del P / del w_m) by the implicit function theorem.
     */
    fun evaluate(
        affine: AffineCoordinates,
        hInv: ComplexMatrix,
        eliminatedAffineIndex: Int
    ): ComplexMatrix {
        val homo = affine.toHomogeneous()
        val numVars = definingPolynomial.numVariables

        // Implicit function theorem derivatives: dw_m / du_a
        val gradP = definingPolynomial.gradient(homo.coordinates)
        val homoElimIdx = if (eliminatedAffineIndex >= affine.chartIndex) eliminatedAffineIndex + 1 else eliminatedAffineIndex
        val dPdwM = gradP[homoElimIdx]
        require(!dPdwM.isZero(1e-15)) { "Singular Jacobian in metric computation" }

        // Tangent coordinate indices in affine coords (excluding eliminatedAffineIndex)
        val tangentIndices = ArrayList<Int>(manifoldDim)
        for (i in 0 until affine.dimension) {
            if (i != eliminatedAffineIndex) tangentIndices.add(i)
        }

        // Projective gradient factors dw_m / du_a = - (del P / del u_a) / (del P / del w_m)
        val dWmDu = Array(manifoldDim) { a ->
            val affIdx = tangentIndices[a]
            val homoIdx = if (affIdx >= affine.chartIndex) affIdx + 1 else affIdx
            -(gradP[homoIdx] / dPdwM)
        }

        // Evaluate basis sections s_alpha and their intrinsic tangent derivatives del_a s_alpha
        val sVals = bundle.basis.evaluateAll(homo)

        // del_a s_alpha for a = 0..manifoldDim-1 and alpha = 0..dim-1
        val dsVals = Array(manifoldDim) { Array(dim) { Complex.ZERO } }

        for (alpha in 0 until dim) {
            val multiIdx = bundle.basis.multiIndices[alpha]
            for (a in 0 until manifoldDim) {
                val affIdx = tangentIndices[a]
                val homoIdx = if (affIdx >= affine.chartIndex) affIdx + 1 else affIdx

                // Explicit derivative wrt homoIdx
                val derDirect = multiIdx.derive(homoIdx)
                var explicitDer = Complex.ZERO
                if (derDirect != null) {
                    var term = Complex.ofReal(derDirect.first.toDouble())
                    for (v in 0 until numVars) {
                        val exp = derDirect.second[v]
                        if (exp > 0) term *= homo[v].pow(exp)
                    }
                    explicitDer = term
                }

                // Chain rule derivative through eliminated coordinate w_m
                val derElim = multiIdx.derive(homoElimIdx)
                var chainDer = Complex.ZERO
                if (derElim != null) {
                    var term = Complex.ofReal(derElim.first.toDouble())
                    for (v in 0 until numVars) {
                        val exp = derElim.second[v]
                        if (exp > 0) term *= homo[v].pow(exp)
                    }
                    chainDer = term * dWmDu[a]
                }

                dsVals[a][alpha] = explicitDer + chainDer
            }
        }

        // S = sum h^{alpha beta} s_alpha bar{s_beta}
        var sSum = 0.0
        for (alpha in 0 until dim) {
            val sa = sVals[alpha]
            if (sa.isZero(1e-15)) continue
            for (beta in 0 until dim) {
                val sbConj = sVals[beta].conj()
                val hInvElem = hInv[alpha, beta]
                sSum += (sa * sbConj * hInvElem).re
            }
        }

        require(sSum > 0.0) { "Section norm sum S must be strictly positive" }
        val invS = 1.0 / sSum
        val invS2 = invS * invS

        // V_a = sum h^{alpha beta} (del_a s_alpha) bar{s_beta}
        val vA = Array(manifoldDim) { Complex.ZERO }
        for (a in 0 until manifoldDim) {
            var sum = Complex.ZERO
            for (alpha in 0 until dim) {
                val dsa = dsVals[a][alpha]
                if (dsa.isZero(1e-15)) continue
                for (beta in 0 until dim) {
                    val sbConj = sVals[beta].conj()
                    val hInvElem = hInv[alpha, beta]
                    sum += dsa * sbConj * hInvElem
                }
            }
            vA[a] = sum
        }

        // g_{a bar{b}} = (1/k) * [ (1/S) sum h^{alpha beta} (del_a s_alpha) bar{del_b s_beta} - (1/S^2) V_a bar{V_b} ]
        val gMat = ComplexMatrix(manifoldDim, manifoldDim)
        val invK = 1.0 / k.toDouble()

        for (a in 0 until manifoldDim) {
            for (b in 0 until manifoldDim) {
                var term1 = Complex.ZERO
                for (alpha in 0 until dim) {
                    val dsa = dsVals[a][alpha]
                    if (dsa.isZero(1e-15)) continue
                    for (beta in 0 until dim) {
                        val dsbConj = dsVals[b][beta].conj()
                        val hInvElem = hInv[alpha, beta]
                        term1 += dsa * dsbConj * hInvElem
                    }
                }

                val term2 = vA[a] * vA[b].conj()
                val valG = (term1 * invS - term2 * invS2) * invK

                gMat.set(a, b, valG)
            }
        }

        // Enforce exact Hermiticity
        for (i in 0 until manifoldDim) {
            gMat.set(i, i, gMat.getRe(i, i), 0.0)
            for (j in (i + 1) until manifoldDim) {
                val avgRe = (gMat.getRe(i, j) + gMat.getRe(j, i)) * 0.5
                val avgIm = (gMat.getIm(i, j) - gMat.getIm(j, i)) * 0.5
                gMat.set(i, j, avgRe, avgIm)
                gMat.set(j, i, avgRe, -avgIm)
            }
        }

        return gMat
    }

    /**
     * Computes the metric determinant det(g_{i bar{j}}) via Cholesky factorization.
     */
    fun determinant(gMat: ComplexMatrix): Double {
        val chol = Cholesky.factorize(gMat)
        var det = 1.0
        for (i in 0 until gMat.rows) {
            val d = chol.lMatrix.getRe(i, i)
            det *= (d * d)
        }
        return det
    }
}
