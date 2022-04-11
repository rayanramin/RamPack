#' duplicated_all
#' @description gets all duplicated values
#' for example in duplicated_all(c("a","b","b","c")), the output will be c(F,T,T,F)
#' @export
#'


duplicated_all <- function(x){
  tmp <- x[duplicated(x)]
  x %in% tmp
}
