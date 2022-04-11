  #' Calculating Uracilation index
  #'
  #' input's column names should be:("Chr.","POS","REF","COV","A","C","G","T")
  #' output is a vector with Uracilation fraction values at C:G positions, the length of the vector is the same as the nrow of the input
  #' this is ~300 times faster than UX_TC for getting only CtoT_GtoA vector
  #' @param Input is readcount data.frame
  #' @param context the nucleotide context of cytosines; "NC" , "TC" , "WRCY"; Defaults to "NC" for all C:G positions.
  #' @keywords Uracilation Fraction Index
  #' @export
  #' @examples
  #' UX(BRC) -> ui
  #'@aliases ux

UX <- function(Input,context = "NC") {
  X<-Input
  if(missing(X)) { stop("ERROR: input X is required") }
  if(any(! c("POS","REF","COV","A","C","G","T") %in% names(X))) {stop("ERROR : input format is incorrect.")}
  if(! context %in% c("NC","TC","WRCY")) stop("Undefined context")
  # defining conditions
  cond1 <- logical(nrow(X))
  cond2 <- logical(nrow(X))
  switch (context,
"NC"= {
  cond1 <- (X$REF == "C")
  cond2 <- (X$REF == "G")},
"TC"={
   for (i in 2:nrow(X)) { cond1[i] <- (( X$REF[i] == "C") & (X$REF[i-1] == "T"))}
  for (i in 1:(nrow(X)-1)) { cond2[i] <- (( X$REF[i] == "G") & (X$REF[i+1] == "A")) }},
  "WRCY"={
    for (i in 3:nrow(X)) { cond1[i] <- (( X$REF[i] == "C") & (X$REF[i-2] == "T" | X$REF[i-2] == "A") & (X$REF[i-1] == "A" | X$REF[i-1] == "G") & (X$REF[i+1] == "C" | X$REF[i+1] == "T")) }
    for (i in 2:(nrow(X)-2)){ cond2[i] <- (( X$REF[i] == "G") & (X$REF[i-1] == "A" | X$REF[i-1] == "G") & (X$REF[i+1] == "C" | X$REF[i+1] == "T") & (X$REF[i+2] == "A" | X$REF[i+2] == "T"))}
    })
# legacy code for NC :
  # cond1 <- logical(nrow(X)) ; cond1 <- X$REF == "C"
  # cond2 <- logical(nrow(X)) ; cond2 <- X$REF == "G"

  # C to T changes only
  t <- ifelse(cond1, as.integer(X$T), 0) / X$COV
  # G to A changes only
  a <- ifelse(cond2, as.integer(X$A), 0) / X$COV
  # All changes
  c(a+t)
}

