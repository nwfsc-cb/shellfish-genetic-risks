#' load shellfish model
#'
#' @return
#' @export
#'
#' @examples
#'
#' \dontrun{
#' load_shellfish()
#' }
#'
#'
load_shellfish <- function(){

  reticulate::source_python(system.file("spam.py", package = "shellfishrisks"), envir = .GlobalEnv)

}
