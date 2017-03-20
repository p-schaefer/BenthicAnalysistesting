#' User Interface
#'
#' Launches the user interface.
#' @export
#' @examples
#' runUI()

runUI <- function() {
  appDir <- system.file("shiny", "myapp", package = "BenthicAnalysis")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing BenthicAnalysis.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}
