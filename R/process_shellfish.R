#' Process Shellfish Genetics Results
#'
#' @param batch the name of the batch of results
#' @param results_dir the name of the directory where results are stored
#'
#' @return a list containing each of the results objects
#' @export
#'
#' @examples
#' \dontrun{
#'
#' results <- process_shellfish(batch = "dev")
#'
#' }
process_shellfish <- function(batch, results_dir = "."){

  # batch <- "dev"

  # slowly and painfully separate out results since my regex skills are poor

  fls <- list.files(results_dir)

  fls <- fls[grepl(batch,fls) & grepl(".txt",fls)] # get just files that match that results naming conventions.


  results <- gsub("(.*\\d_)","", fls)

  coreid <- as.integer(gsub("\\D","", fls))

  coreid <- unique(coreid)

  loader <- function(result, fls, results_dir, coreid) {

    storage <-  lapply(1:coreid, function(r) read.delim(file.path(results_dir, fls[grepl(result, fls) & grepl(r, fls)])))

    storage <- dplyr::bind_rows(storage, .id = "coreid")

    storage$coreid <- as.integer(storage$coreid)

    return(storage)

  }
  out <- lapply(results, loader, fls = fls, results_dir = results_dir, coreid = coreid)

  out <- setNames(out, gsub("\\.txt","",results))

  return(out)

}
