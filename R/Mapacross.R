#' Mapacross, similar to Mike Lawrence's Matlab function
#' @description map an element to a vector A and return matching value from another vector B
#' @param x input vector
#' @param A the vector to map x to
#' @param B the vector the extract ouput from
#' output will be of from B or m=NA if missing
#' @example mapacross(c(1,4,3),c(1:10), letters[1:10]) -> out
#' @export

# write a function to map an element to a vector A and return matching value from another vector B
mapacross <- function(x, A, B, m=NA) {
  ifelse(x %in% A, B[match(x, A)], m)
}

