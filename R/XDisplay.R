#' Set DISPLAY settings for X11 forwarding
#' @example XDisplay(11)
#' sets display environment to localhost:11.0
#'
#' @param version int (default=10)
#' @export


XDisplay <- function(version=10, from_system=F){
  eval(parse(text = paste0('Sys.setenv("DISPLAY"="localhots:', version , '.0")')))
  if(from_system){
    try(system("echo $DISPLAY")) -> version
    eval(parse(text = paste0('Sys.setenv("DISPLAY"=', version , '")')))
  }
}
