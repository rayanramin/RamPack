#' melt
#' @description melt from reshape package with a small fix
#'
#'
melt <- function (data, varnames = names(dimnames(data)), ...)
{
  values <- as.vector(data)
  dn <- dimnames(data)
  if (is.null(dn))
    dn <- vector("list", length(dim(data)))
  dn_missing <- sapply(dn, is.null)
  dn[dn_missing] <- lapply(dim(data), function(x) 1:x)[dn_missing]
  char <- sapply(dn, is.character)
  dn[char] <- lapply(dn[char], type.convert, as.is = TRUE)
  indices <- do.call(expand.grid, dn)
  names(indices) <- varnames
  data.frame(indices, value = values)
}
