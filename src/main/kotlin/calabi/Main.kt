package calabi

import calabi.cli.CommandDispatcher
import calabi.cli.Repl

fun main(args: Array<String>) {
    if (args.isEmpty()) {
        Repl().run()
        return
    }

    val dispatcher = CommandDispatcher()
    when (args[0]) {
        "--bench", "-b" -> {
            calabi.bench.BenchmarkHarness.runAll()
        }
        "--instantons", "-i" -> {
            dispatcher.handleCommand(":instantons 4")
        }
        "--donaldson", "-d" -> {
            val k = args.getOrNull(1)?.toIntOrNull() ?: 1
            dispatcher.handleCommand(":donaldson $k 5")
        }
        "--help", "-h" -> {
            dispatcher.handleCommand(":help")
        }
        else -> {
            val cmd = args.joinToString(" ")
            dispatcher.handleCommand(cmd)
        }
    }
}
