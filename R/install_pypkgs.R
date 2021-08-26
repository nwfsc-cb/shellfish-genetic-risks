#' Install Required Python Packages
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
install_pypkgs <- function(method = "auto", conda = "auto") {

  reticulate::use_condaenv("r-reticulate")

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
      "simuPOP",
      "pickleshare"
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
