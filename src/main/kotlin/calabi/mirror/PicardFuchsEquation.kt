package calabi.mirror

import calabi.algebra.Complex

/**
 * PicardFuchsEquation represents the 4th-order Picard-Fuchs differential operator
 * governing the period integrals of the Dwork pencil of quintic threefolds:
 *
 *   L = theta^4 - 5 * z * (5*theta + 1)(5*theta + 2)(5*theta + 3)(5*theta + 4)
 *
 * where theta = z * (d / dz) is the Euler logarithmic derivative operator,
 * and z = (5 * psi)^{-5} is the moduli coordinate around the Large Complex Structure Limit (LCSL).
 */
class PicardFuchsEquation {

    /**
     * Evaluates the action of L on a function f(z) with known Euler derivatives:
     * theta_k = theta^k(f)(z) for k = 0, 1, 2, 3, 4.
     */
    fun evaluateOperator(
        z: Complex,
        theta0: Complex,
        theta1: Complex,
        theta2: Complex,
        theta3: Complex,
        theta4: Complex
    ): Complex {
        // Expand (5*theta + 1)(5*theta + 2)(5*theta + 3)(5*theta + 4):
        // Product of (5*th + 1)(5*th + 4) = 25*th^2 + 25*th + 4
        // Product of (5*th + 2)(5*th + 3) = 25*th^2 + 25*th + 6
        // Product = (25*th^2 + 25*th)^2 + 10*(25*th^2 + 25*th) + 24
        //         = 625*th^4 + 1250*th^3 + 875*th^2 + 250*th + 24
        // Multiplied by 5*z:
        // 5*z * (625*th^4 + 1250*th^3 + 875*th^2 + 250*th + 24)
        // = z * (3125*th^4 + 6250*th^3 + 4375*th^2 + 1250*th + 120)

        val polyTheta = theta4 * 3125.0 + theta3 * 6250.0 + theta2 * 4375.0 + theta1 * 1250.0 + theta0 * 120.0
        val term2 = z * polyTheta
        return theta4 - term2
    }

    companion object {
        /**
         * Returns the conifold singular points of the differential equation:
         * 1 - 5^5 * z = 0 => z_c = 1 / 3125.
         */
        val CONIFOLD_SINGULARITY: Complex = Complex.ofReal(1.0 / 3125.0)
    }
}
