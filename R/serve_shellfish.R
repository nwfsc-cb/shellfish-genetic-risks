#' Load and Serve Shellfish Genetics Results
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
#' results <- serve_shellfish(batch = "dev")
#'
#' }
serve_shellfish <- function(batches, results_dir = NA){

  # batch <- "dev"

  # slowly and painfully separate out results since my regex skills are poor


  tmp <- vector(mode = "list", length = length(batches))

  if (is.na(results_dir)) {

    results_dir <- "results"

  }

  for (b in seq_along(batches)) {

      tmp_results_dir <- file.path(results_dir, batches[b])


  fls <- list.files(tmp_results_dir)

  if (length(fls) == 0) {
    stop("No results found in batch folder, check back name or results folder file path")
  }

  fls <- fls[grepl(batches[b],fls) & grepl(".txt",fls)] # get just files that match that results naming conventions.

  tmp_fls <- gsub(paste0("batch_",batches[b],"_"),"", fls) # trim out batch stuff

  results <- gsub(paste0("batch_",batches[b],"_","(.*\\d_)"),"", fls) # all results files are preceeded by a underscore followed by an integer

  coreid <- unique(as.integer(gsub("\\D","", tmp_fls)))

  loader <- function(result, fls, tmp_results_dir, coreid) {

    storage <-  lapply(1:coreid, function(r) read.delim(file.path(tmp_results_dir, fls[grepl(result, fls) & grepl(r, fls)])))

    storage <- dplyr::bind_rows(storage, .id = "coreid")

    storage$coreid <- as.integer(storage$coreid)

    storage$batch <- batches[b]

    return(storage)

  }
  out <- lapply(results, loader, fls = fls, tmp_results_dir = tmp_results_dir, coreid = coreid)

  out <- setNames(out, gsub("\\.txt","",results))

  tmp[[b]] <- out

  } # close batches loop

  obj_names <- names(tmp[[1]])

  out <- purrr::map(obj_names, ~ purrr::map_df(tmp, .x))

  out <- setNames(out, gsub("\\.txt","",results))

  # calculate survivorship table

  out$survival <- out$life_hist_report %>%
    dplyr::mutate(age = floor(Age)) %>%
    dplyr::group_by(IndID, age, batch) %>%
    dplyr::summarise(alive = 1, .groups = "drop") %>%
    dplyr::group_by(IndID, batch) %>%
    dplyr::mutate(died = alive - dplyr::lead(alive, default = 0)) %>%
    dplyr::group_by(age,batch) %>%
    dplyr::summarise(survival = mean(died == 0), .groups = "drop") %>%
    dplyr::group_by(batch) %>%
    tidyr::nest() %>%
    dplyr::mutate(survivorship = purrr::map(data, ~ c(1,cumprod(.x$survival[1:(length(.x$survival) - 1)])))) %>%
    tidyr::unnest(cols = c(data,survivorship))

  return(out)

}
