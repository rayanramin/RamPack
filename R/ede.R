#' ede
#' @description ensure directory exists
#' for example ede("/path/to/directory/")
#' @export
#'

ede <- function(dir){
  if(!dir.exists(paths = dir)){dir.create(path = dir,recursive = T,showWarnings = F)}
}
