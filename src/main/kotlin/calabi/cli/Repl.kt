package calabi.cli

import java.util.Scanner

class Repl {
    private val dispatcher = CommandDispatcher()

    fun run() {
        print(Formatter.banner())
        val scanner = Scanner(System.`in`)

        while (true) {
            print("${Formatter.ANSI_BOLD}calabi> ${Formatter.ANSI_RESET}")
            if (!scanner.hasNextLine()) break
            val line = scanner.nextLine()
            val continueRunning = dispatcher.handleCommand(line)
            if (!continueRunning) {
                println("${Formatter.ANSI_CYAN}Exiting Calabi. Goodbye!${Formatter.ANSI_RESET}")
                break
            }
        }
    }
}
