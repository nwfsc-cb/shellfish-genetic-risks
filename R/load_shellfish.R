#' load shellfish model
#'
#' @param condaenv  the Conda environment used
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
load_shellfish <- function(condaenv = "r-reticulate"){

  reticulate::use_condaenv(condaenv)

  message(paste0("shellfishrisks is using Conda environment ", condaenv))

  reticulate::source_python(system.file("spam.py", package = "shellfishrisks"), envir = .GlobalEnv)

}
