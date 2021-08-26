#' Clean up results of shellfishrisks
#'
#' @param batch
#'
#' @export
#'
#' @examples
#' \dontrun{
#' clean_shellfish(batch = "dev")
#' }
clean_shellfish <- function(batch){

newdir <- file.path("results",batch)

if (!dir.exists(newdir)){

  dir.create(newdir, recursive = TRUE)

}

fls <- list.files(".")

fls <- fls[grepl(batch,fls) & (grepl(".txt",fls) | grepl(".pop",fls)) | grepl("progress.txt", fls)] # get just files that match that results naming conventions.

file.copy(fls, to = newdir, overwrite = TRUE)

file.remove(fls)

message(paste("results moved to",newdir))

}
