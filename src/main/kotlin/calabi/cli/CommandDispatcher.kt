package calabi.cli

import calabi.algebra.Complex
import calabi.curvature.DonaldsonSigmaMeasure
import calabi.curvature.MetricTensor
import calabi.holonomy.ParallelTransport
import calabi.holonomy.SU3Holonomy
import calabi.linalg.ComplexMatrix
import calabi.manifold.DworkPencil
import calabi.manifold.QuinticThreefold
import calabi.metric.BalancedMetricSolver
import calabi.metric.MonomialBasis
import calabi.metric.SectionBundle
import calabi.mirror.FrobeniusMethod
import calabi.mirror.InstantonExpansion
import calabi.mirror.MirrorMap
import calabi.sampling.MonteCarloIntegrator
import kotlin.math.PI

class CommandDispatcher {

    private val quintic = QuinticThreefold()
    private val integrator = MonteCarloIntegrator(quintic, seed = 42L)
    private val instanton = InstantonExpansion()

    fun handleCommand(input: String): Boolean {
        val trimmed = input.trim()
        if (trimmed.isEmpty()) return true

        val parts = trimmed.split("\\s+".toRegex())
        val cmd = parts[0].lowercase()

        when (cmd) {
            ":help", "help", ":h" -> printHelp()
            ":sample", "sample" -> {
                val count = parts.getOrNull(1)?.toIntOrNull() ?: 5
                samplePoints(count)
            }
            ":instantons", "instantons" -> {
                val maxDeg = parts.getOrNull(1)?.toIntOrNull() ?: 4
                displayInstantons(maxDeg)
            }
            ":donaldson", "donaldson" -> {
                val k = parts.getOrNull(1)?.toIntOrNull() ?: 1
                val iters = parts.getOrNull(2)?.toIntOrNull() ?: 5
                runDonaldson(k, iters)
            }
            ":sigma", "sigma" -> {
                val k = parts.getOrNull(1)?.toIntOrNull() ?: 1
                computeSigma(k)
            }
            ":periods", "periods" -> {
                val psiVal = parts.getOrNull(1)?.toDoubleOrNull() ?: 0.5
                displayPeriods(psiVal)
            }
            ":mirror", "mirror" -> {
                val psiVal = parts.getOrNull(1)?.toDoubleOrNull() ?: 0.5
                displayMirror(psiVal)
            }
            ":holonomy", "holonomy" -> {
                val segs = parts.getOrNull(1)?.toIntOrNull() ?: 8
                computeHolonomy(segs)
            }
            ":web", "web" -> {
                println("${Formatter.ANSI_GREEN}WebGL Laboratory available in ./web directory.${Formatter.ANSI_RESET}")
                println("Run 'make web-serve' or 'open web/index.html' in your browser.")
            }
            ":quit", ":exit", "quit", "exit" -> return false
            else -> {
                println("${Formatter.ANSI_YELLOW}Unknown command: $cmd. Type :help for commands.${Formatter.ANSI_RESET}")
            }
        }
        return true
    }

    private fun printHelp() {
        val headers = listOf("Command", "Arguments", "Description")
        val rows = listOf(
            listOf(":sample", "[count]", "Samples points on quintic threefold via Grassmannian line intersection"),
            listOf(":donaldson", "[k] [iters]", "Solves for Donaldson balanced metric at degree k"),
            listOf(":sigma", "[k]", "Computes Monge-Ampere Ricci flatness error sigma(k)"),
            listOf(":instantons", "[maxDegree]", "Displays genus-0 Gromov-Witten curve counts n_1, ..., n_d"),
            listOf(":periods", "[psi]", "Evaluates Picard-Fuchs periods varpi_0, varpi_1 for Dwork pencil"),
            listOf(":mirror", "[psi]", "Computes flat coordinate t(z) and exponentiated mirror map q"),
            listOf(":holonomy", "[segments]", "Integrates parallel transport around closed loop and verifies SU(3)"),
            listOf(":web", "", "Displays WebGL interactive visualization laboratory instructions"),
            listOf(":quit", "", "Exits the Calabi interactive session")
        )
        println(Formatter.formatTable(headers, rows))
    }

    private fun samplePoints(count: Int) {
        println("${Formatter.ANSI_CYAN}Sampling $count points on ${quintic.name}...${Formatter.ANSI_RESET}")
        val points = integrator.sampleWeightedPoints(count)
        val headers = listOf("#", "Coordinates [z_0 : z_1 : z_2 : z_3 : z_4]", "Weight w(x)")
        val rows = points.mapIndexed { idx, sp ->
            val coordStr = sp.point.coordinates.joinToString(", ") { c -> String.format("%.3g", c.re) }
            listOf((idx + 1).toString(), "[$coordStr]", String.format("%.5e", sp.weight))
        }
        println(Formatter.formatTable(headers, rows))
    }

    private fun displayInstantons(maxDegree: Int) {
        println("${Formatter.ANSI_CYAN}Gromov-Witten Invariants (Instantons) on Quintic Threefold:${Formatter.ANSI_RESET}")
        val reports = instanton.computeInstantons(maxDegree.coerceIn(1, 4))
        val headers = listOf("Degree (d)", "Instanton Count (n_d)", "Classification")
        val rows = reports.map { r ->
            listOf("d = ${r.degree}", String.format("%,d", r.instantonCount), r.classicalName)
        }
        println(Formatter.formatTable(headers, rows))
    }

    private fun runDonaldson(k: Int, iters: Int) {
        val basis = MonomialBasis(ambientDimension = 4, degree = k)
        val bundle = SectionBundle.identity(basis)
        val solver = BalancedMetricSolver(bundle, maxIterations = iters, tolerance = 1e-4)

        println("${Formatter.ANSI_CYAN}Solving Donaldson balanced metric for O($k) (dim N_$k = ${basis.size})...${Formatter.ANSI_RESET}")
        val points = integrator.sampleWeightedPoints(40)

        val report = solver.solve(points) { step, change ->
            println("  Iteration $step: relative change = ${String.format("%.4e", change)}")
        }

        println("${Formatter.ANSI_GREEN}Completed ${report.iterations} iterations.${Formatter.ANSI_RESET}")
        println("Final relative change : ${String.format("%.5e", report.relativeChange)}")
        println("Trace normalization   : ${String.format("%.4f", basis.size.toDouble() + report.traceError)} / ${basis.size}")
    }

    private fun computeSigma(k: Int) {
        val basis = MonomialBasis(ambientDimension = 4, degree = k)
        val bundle = SectionBundle.identity(basis)
        val metricTensor = MetricTensor(bundle, quintic.polynomial)
        val sigmaMeasure = DonaldsonSigmaMeasure(metricTensor, quintic.poincareResidue)

        println("${Formatter.ANSI_CYAN}Computing Monge-Ampere sigma error for O($k)...${Formatter.ANSI_RESET}")
        val points = integrator.sampleWeightedPoints(30)
        val report = sigmaMeasure.evaluate(ComplexMatrix.identity(basis.size), points)

        println("Mean sigma error : ${Formatter.ANSI_BOLD}${String.format("%.5f", report.sigma)}${Formatter.ANSI_RESET}")
        println("Max local error  : ${String.format("%.5f", report.sigmaMax)}")
        println("Volume factor    : ${String.format("%.5f", report.kappa)}")
    }

    private fun displayPeriods(psiVal: Double) {
        val psi = Complex.ofReal(psiVal)
        val pencil = DworkPencil(psi)
        val z = pencil.moduliCoordinate()
        val frob = FrobeniusMethod(maxOrder = 8)

        val varpi0 = frob.evaluateVarpi0(z)
        val varpi1 = frob.evaluateVarpi1(z)

        println("${Formatter.ANSI_CYAN}Picard-Fuchs Periods for psi = $psiVal:${Formatter.ANSI_RESET}")
        println("  Moduli coordinate z = (5 psi)^(-5) : $z")
        println("  Period varpi_0(z)                  : $varpi0")
        println("  Period varpi_1(z)                  : $varpi1")
    }

    private fun displayMirror(psiVal: Double) {
        val psi = Complex.ofReal(psiVal)
        val pencil = DworkPencil(psi)
        val z = pencil.moduliCoordinate()
        val mirror = MirrorMap()

        val t = mirror.evaluateT(z)
        val q = mirror.evaluateQ(z)

        println("${Formatter.ANSI_CYAN}Mirror Map Coordinates for psi = $psiVal:${Formatter.ANSI_RESET}")
        println("  Flat coordinate t(z)    : $t")
        println("  Exponentiated coord q(z): $q")
    }

    private fun computeHolonomy(numSegments: Int) {
        val basis = MonomialBasis(ambientDimension = 4, degree = 1)
        val bundle = SectionBundle.identity(basis)
        val metricTensor = MetricTensor(bundle, quintic.polynomial)
        val christoffel = calabi.curvature.ChristoffelSymbols(metricTensor)
        val transport = ParallelTransport(christoffel)
        val su3 = SU3Holonomy(transport)

        println("${Formatter.ANSI_CYAN}Computing SU(3) holonomy matrix over $numSegments closed loop segments...${Formatter.ANSI_RESET}")
        val z1 = Complex.fromPolar(1.0, PI / 5.0)
        val loop = (0..numSegments).map { s ->
            val theta = 2.0 * PI * s / numSegments
            val z2 = Complex.fromPolar(0.01, theta)
            val pt = calabi.algebra.HomogeneousPoint.of(Complex.ONE, z1, z2, Complex.ZERO, Complex.ZERO)
            pt.toAffineCoordinates(0)
        }
        val elims = loop.map { quintic.poincareResidue.selectEliminatedCoordinate(it).first }

        val report = su3.computeHolonomy(loop, elims, ComplexMatrix.identity(5))
        println("Holonomy Determinant det(H) : ${Formatter.ANSI_BOLD}${report.determinant}${Formatter.ANSI_RESET}")
        println("Unitarity Error ||H^† H - I|| : ${String.format("%.4e", report.unitarityError)}")
        println("Phase Angle arg(det(H))     : ${String.format("%.4e rad", report.phaseAngle)}")
    }
}
