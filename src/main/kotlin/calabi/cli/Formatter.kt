package calabi.cli

import calabi.algebra.Complex
import calabi.algebra.HomogeneousPoint

object Formatter {

    const val ANSI_RESET = "\u001B[0m"
    const val ANSI_BOLD = "\u001B[1m"
    const val ANSI_CYAN = "\u001B[36m"
    const val ANSI_GREEN = "\u001B[32m"
    const val ANSI_YELLOW = "\u001B[33m"
    const val ANSI_BLUE = "\u001B[34m"
    const val ANSI_MAGENTA = "\u001B[35m"

    fun banner(): String = """
$ANSI_CYAN$ANSI_BOLD
   ██████╗ █████╗ ██╗      █████╗ ██████╗ ██╗
  ██╔════╝██╔══██╗██║     ██╔══██╗██╔══██╗██║
  ██║     ███████║██║     ███████║██████╔╝██║
  ██║     ██╔══██║██║     ██╔══██║██╔══██╗██║
  ╚██████╗██║  ██║███████╗██║  ██║██████╔╝██║
   ╚═════╝╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═════╝ ╚═╝
$ANSI_RESET
  Ricci-Flat Kähler Metrics, Donaldson Balanced Iterations & Mirror Symmetry
  Type $ANSI_YELLOW:help$ANSI_RESET for commands, or $ANSI_YELLOW:quit$ANSI_RESET to exit.
"""

    fun formatTable(
        headers: List<String>,
        rows: List<List<String>>
    ): String {
        val colWidths = headers.mapIndexed { idx, h ->
            val maxRow = rows.maxOfOrNull { it.getOrNull(idx)?.length ?: 0 } ?: 0
            kotlin.math.max(h.length, maxRow) + 2
        }

        val sb = StringBuilder()

        // Header
        val hLine = colWidths.joinToString("+", "+", "+") { "-".repeat(it) }
        sb.appendLine(hLine)
        sb.append("|")
        for (i in headers.indices) {
            val pad = colWidths[i] - headers[i].length - 1
            sb.append(" ").append(ANSI_BOLD).append(headers[i]).append(ANSI_RESET).append(" ".repeat(kotlin.math.max(0, pad))).append("|")
        }
        sb.appendLine()
        sb.appendLine(hLine)

        // Rows
        for (r in rows) {
            sb.append("|")
            for (i in headers.indices) {
                val cell = r.getOrNull(i) ?: ""
                val pad = colWidths[i] - cell.length - 1
                sb.append(" ").append(cell).append(" ".repeat(kotlin.math.max(0, pad))).append("|")
            }
            sb.appendLine()
        }
        sb.append(hLine)

        return sb.toString()
    }
}
