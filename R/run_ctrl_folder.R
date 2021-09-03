#' Run all control files in a folder
#'
#' @param ctrl_folder the location of the folder containing the control files
#' @param reps the number of reps per control file
#' @param cores the number of cores. Set to greater than 1 to run in parallel
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' run_ctrl_folder("my_control_files", cores = 6)
#'
#' }
#'
run_ctrl_folder <- function(ctrl_folder, reps = 1, cores = 1){


  `%dopar%` <- foreach::`%dopar%`

  ctrl_files <- list.files(ctrl_folder)

  ctrl_files <- ctrl_files[grepl("\\.csv", ctrl_files)]

  batches <- gsub("\\.csv", "",ctrl_files)

  doParallel::registerDoParallel(cores = cores)

  on.exit(doParallel::registerDoParallel(cores = 1))

  foreach::foreach(i = seq_along(ctrl_files)) %dopar% {


    run_shellfishrisk_ctrlfile(batch = batches[i], reps = 1, ctrl_file_path = here::here(ctrl_folder,ctrl_files[i]))


  }


}
