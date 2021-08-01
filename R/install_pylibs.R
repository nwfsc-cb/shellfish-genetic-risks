#' Install Required Python Libraries
#'
#' @param method python install method
#' @param conda path to the conda executable
#'
#' @export
#'
#' @examples
#' \dontrun{
#' install_pylibs()
#' }
#'
install_pylibs <- function(method = "auto", conda = "auto") {

  # vector of required python packages
  pypkgs <-
    c(
      "argparse",
      "pandas",
      "numpy",
      "itertools",
      "random",
      "time",
      "datetime",
      "os",
      "subprocess",
      "simuPOP"
    )


  #installation wrapper
  instfoo <- function(pkg) {
    has <- reticulate::py_module_available(pkg)
    if (!has) {
      reticulate::py_install(pkg, method = method, conda = conda)
    }
  }

  for (i in pypkgs) {
    instfoo(i)
  }


}
