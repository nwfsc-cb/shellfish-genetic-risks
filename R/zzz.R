.onLoad <- function(libname, pkgname) {

  library(reticulate)

  reticulate::use_condaenv("r-reticulate")

  message("shellfishrisks is using Conda environment 'r-reticulate'")

  reticulate::source_python(system.file("spam.py", package = "shellfishrisks"), envir = .GlobalEnv)


  # # global reference to scipy (will be initialized in .onLoad)
  # scipy <- NULL
  #
  # .onLoad <- function(libname, pkgname) {
  #   # use superassignment to update global reference to scipy
  #   scipy <<- reticulate::import("scipy", delay_load = TRUE)
  # }
}
