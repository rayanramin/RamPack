#' Set DISPLAY settings for X11 forwarding
#' @example XDisplay(11)
#' sets display environment to localhost:11.0
#'
#' @param version int (default=10)
#' @export


XDisplay <- function(version=10, from_system=F){
  if(!from_system){
    command=paste0('Sys.setenv("DISPLAY"="localhost:', version , '.0")')
    cat(command,"\n")
  eval(parse(text = command))
  }  else if(from_system){
    print("setting system environment")
    system("echo $DISPLAY", intern = TRUE) -> version
    command=paste0('Sys.setenv("DISPLAY"="localhost:', version , '")')
    cat(command,"\n")
    eval(parse(text = paste0('Sys.setenv("DISPLAY"="', version , '")')))
  }
}
