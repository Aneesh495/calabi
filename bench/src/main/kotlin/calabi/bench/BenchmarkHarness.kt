package calabi.bench

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint
import calabi.linalg.Cholesky
import calabi.linalg.ComplexMatrix
import calabi.manifold.QuinticThreefold
import calabi.metric.DonaldsonTOperator
import calabi.metric.MonomialBasis
import calabi.metric.SectionBundle
import calabi.mirror.FrobeniusMethod
import calabi.mirror.InstantonExpansion
import calabi.mirror.MirrorMap
import calabi.sampling.MonteCarloIntegrator
import kotlin.system.measureNanoTime

object BenchmarkHarness {

    fun runAll() {
        println("=================================================================")
        println("           Calabi-Yau Metric & Mirror Symmetry Benchmark         ")
        println("=================================================================\n")

        benchMonteCarloSampling()
        benchMonomialSectionEvaluation()
        benchCholeskyAndInversion()
        benchDonaldsonIteration()
        benchMirrorSymmetryInversion()

        println("\n=================================================================")
        println("             All Benchmark Suites Completed Successfully         ")
        println("=================================================================")
    }

    private fun benchMonteCarloSampling() {
        print("[1] Benchmarking Monte Carlo Grassmannian Sampling... ")
        val quintic = QuinticThreefold()
        val integrator = MonteCarloIntegrator(quintic, seed = 123L)

        // Warmup
        integrator.sampleWeightedPoints(20)

        val targetCount = 100
        val nano = measureNanoTime {
            integrator.sampleWeightedPoints(targetCount)
        }
        val sec = nano / 1e9
        val rate = targetCount / sec
        println("DONE")
        println("    Sampled $targetCount points in ${String.format("%.4f", sec)} sec (${String.format("%.1f", rate)} points/sec)")
    }

    private fun benchMonomialSectionEvaluation() {
        print("[2] Benchmarking Monomial Section Evaluation (k=1, 2, 3)... ")
        val pt = HomogeneousPoint.ofDoubles(1.0, 0.5, -0.8, 1.2, -0.3)
        val trials = 100000

        for (k in 1..3) {
            val basis = MonomialBasis(ambientDimension = 4, degree = k)
            val tempRe = DoubleArray(basis.size)
            val tempIm = DoubleArray(basis.size)

            // Warmup
            repeat(1000) { basis.evaluateAllPrimitive(pt, tempRe, tempIm) }

            val nano = measureNanoTime {
                for (t in 0 until trials) {
                    basis.evaluateAllPrimitive(pt, tempRe, tempIm)
                }
            }
            val sec = nano / 1e9
            val mps = (trials * basis.size) / (sec * 1e6)
            println("k=$k (N_$k = ${basis.size}): ${String.format("%.2f", mps)} million evaluations/sec")
            print("    ")
        }
        println("DONE")
    }

    private fun benchCholeskyAndInversion() {
        print("[3] Benchmarking Complex Cholesky & Matrix Inversion... ")
        val n = 35 // N_3 = 35 sections
        val b = ComplexMatrix(n, n)
        for (i in 0 until n) {
            for (j in 0 until n) {
                b.set(i, j, Complex(1.0 / (i + j + 1.0), 0.5 / (i + j + 2.0)))
            }
        }
        val a = (b * b.conjugateTranspose()) + (ComplexMatrix.identity(n) * 2.0)

        // Warmup
        repeat(50) { Cholesky.factorize(a).inverse() }

        val trials = 500
        val nano = measureNanoTime {
            for (t in 0 until trials) {
                val chol = Cholesky.factorize(a)
                chol.inverse()
            }
        }
        val sec = nano / 1e9
        val ops = trials / sec
        println("DONE")
        println("    Matrix dimension $n x $n: ${String.format("%.1f", ops)} factorize+invert ops/sec")
    }

    private fun benchDonaldsonIteration() {
        print("[4] Benchmarking Donaldson T-Operator Fixed-Point Step... ")
        val quintic = QuinticThreefold()
        val basis = MonomialBasis(ambientDimension = 4, degree = 2) // N_2 = 15
        val bundle = SectionBundle.identity(basis)
        val tOp = DonaldsonTOperator(bundle)

        val integrator = MonteCarloIntegrator(quintic, seed = 555L)
        val points = integrator.sampleWeightedPoints(50)
        val h = ComplexMatrix.identity(basis.size)

        // Warmup
        repeat(2) { tOp.evaluate(h, points) }

        val trials = 10
        val nano = measureNanoTime {
            for (t in 0 until trials) {
                tOp.evaluate(h, points)
            }
        }
        val msPerIter = (nano / 1e6) / trials
        println("DONE")
        println("    k=2 (N_2=15), 50 points: ${String.format("%.2f", msPerIter)} ms/iteration")
    }

    private fun benchMirrorSymmetryInversion() {
        print("[5] Benchmarking Picard-Fuchs Frobenius Series & Reversion... ")
        val instanton = InstantonExpansion()

        // Warmup
        repeat(5) { instanton.computeInstantons(3) }

        val trials = 50
        val nano = measureNanoTime {
            for (t in 0 until trials) {
                instanton.computeInstantons(4)
            }
        }
        val msPerRun = (nano / 1e6) / trials
        println("DONE")
        println("    Frobenius reversion up to degree 4: ${String.format("%.3f", msPerRun)} ms/run")
    }
}
