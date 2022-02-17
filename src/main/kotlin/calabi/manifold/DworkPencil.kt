package calabi.manifold

import calabi.algebra.Complex
import kotlin.math.PI

/**
 * DworkPencil represents the 1-parameter family of Calabi-Yau quintic threefolds:
 *
 *   X_psi: z_0^5 + z_1^5 + z_2^5 + z_3^5 + z_4^5 + psi * z_0 z_1 z_2 z_3 z_4 = 0
 *
 * This pencil admits a large discrete symmetry group G = (Z_5)^3.
 * The quotient W = X_psi / G is the mirror quintic manifold, which has inverted
 * Hodge numbers: h^{1,1}(W) = 101, h^{2,1}(W) = 1, chi(W) = +200.
 */
class DworkPencil(
    val psi: Complex
) {
    val manifold: QuinticThreefold = QuinticThreefold(psi)

    /**
     * Large Complex Structure Limit (LCSL) coordinate:
     * z = (5 * psi)^{-5} = 1 / (3125 * psi^5)
     */
    fun moduliCoordinate(): Complex {
        require(!psi.isZero(1e-15)) { "LCSL moduli coordinate z undefined at Fermat point psi=0" }
        val denom = psi.pow(5) * 3125.0
        return denom.inv()
    }

    /**
     * Checks if the parameter lies near the conifold singular point (psi^5 = 3125).
     */
    fun isNearConifold(threshold: Double = 0.05): Boolean {
        val psi5 = psi.pow(5)
        return (psi5 - Complex.ofReal(3125.0)).abs() < threshold * 3125.0
    }

    /**
     * Five conifold singular values of psi: psi_k = 5 * exp(2 * pi * i * k / 5).
     */
    companion object {
        fun conifoldPoints(): List<Complex> =
            (0 until 5).map { k ->
                Complex.fromPolar(5.0, 2.0 * PI * k / 5.0)
            }
    }
}
