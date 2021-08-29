#' Generate a control file from the default template
#'
#' @param ctrl_file_name the desired name of the control file
#' @param destdir the directory to store the control file in
#'
#' @return the file path of the control file
#' @export
#'
#' @examples
#' \dontrun{
#'
#' generate_ctrl_file()
#' }
#'
generate_ctrl_file <-
  function(ctrl_file_name = "ctrl",
           destdir = ".") {
    file.copy(
      system.file("ctrl.csv", package = "shellfishrisks"),
      to = file.path(destdir, paste0(ctrl_file_name, ".csv")),
      overwrite = TRUE,
      recursive = FALSE
    )

    out <- file.path(destdir, paste0(ctrl_file_name, ".csv"))
  }
