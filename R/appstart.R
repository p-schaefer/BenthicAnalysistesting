#' User Interface
#'
#' Launches the user interface.
#' @export
#' @examples
#' runUI()

runUI <- function() {
  appDir <- system.file("shiny", "myapp", package = "BenthicAnalysistesting")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing BenthicAnalysistesting.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}
