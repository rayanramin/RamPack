#' clear the environment
#' @example clear()
#'
#' @export


clear <- function(){
  command=paste0('rm(list=ls())')
  cat("removing everything with: ",command,"\n")
  eval(parse(text = command), envir=.GlobalEnv)
}
