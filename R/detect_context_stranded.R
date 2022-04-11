#' Detect sequence context (stranded)
#' @description  This function detects NC, TC, WRC, WRCY contexts and returns a vector of -1, 1 and 0 showing the presense and absense of the selected context. -1 and +1 shows the strand that cytosine is on.
#' This function is not optimized for speed.
#' @param X data.frame with the following columns: POS (numeric) & REF ("A", "C", "G", "T")
#' @param context context to be detected, corrently only "NC","TC","WRC" and "WRCY" are covered.
#'
#' @export
#' @examples
#' detect_context_stranded(Data,"TC") -> out
#'


detect_context_stranded <- function(X, context = "NC"){
  if (missing(X)) {
    stop("ERROR: input X is required")}
  if (any(!c("POS", "REF", "COV") %in% names(X))) {
    stop("ERROR : input format is incorrect.")}
  if (!context %in% c("NC", "TC","WRC", "WRCY")){
    stop("Undefined context")}
  output <- rep(0, nrow(X))

  switch(context, NC = {
    ifelse(X$REF == "C", 1, ifelse(X$REF == "G",-1,0)) -> output
  }, TC = {
    for (i in 1:nrow(X)) {
      pos <- X$POS[i]
      posp1 <- X$REF[X$POS == pos + 1]
      posm1 <- X$REF[X$POS == pos - 1]
      if((X$REF[i] == "C") && length(posm1)==1  && (posm1 =="T"))1 -> output[i]
      if((X$REF[i] == "G") && length(posp1)==1 && (posp1 =="A"))-1 -> output[i]
    }
  }, WRC = {
    for (i in 1:nrow(X)) {
      pos <- X$POS[i]
      posp2 <- X$REF[X$POS == pos + 2]
      posp1 <- X$REF[X$POS == pos + 1]
      posm1 <- X$REF[X$POS == pos - 1]
      posm2 <- X$REF[X$POS == pos - 2]
      if((X$REF[i] == "C") &&
         ((length(posm1)+length(posm2))==2) &&
         (posm2 %in% c("T","A")) &&
         (posm1 %in% c("A","G")))  1 -> output[i]
      if((X$REF[i] == "G") &&
         ((length(posp2)+length(posp1))==2) &&
         (posp1 %in% c("C","T")) &&
         (posp2 %in% c("A","T"))) -1 -> output[i]
    }
  }, WRCY = {
    for (i in 1:nrow(X)) {
      pos <- X$POS[i]
      posp2 <- X$REF[X$POS == pos + 2]
      posp1 <- X$REF[X$POS == pos + 1]
      posm1 <- X$REF[X$POS == pos - 1]
      posm2 <- X$REF[X$POS == pos - 2]
      if((X$REF[i] == "C") &&
         ((length(posp1)+length(posm1)+length(posm2))==3) &&
         (posm2 %in% c("T","A")) &&
         (posm1 %in% c("A","G")) &&
         (posp1 %in% c("C","T")))  1 -> output[i]
      if((X$REF[i] == "G") &&
         ((length(posp2)+length(posp1)+length(posm1))==3) &&
         (posm1 %in% c("A","G")) &&
         (posp1 %in% c("C","T")) &&
         (posp2 %in% c("A","T"))) -1 -> output[i]
    }
  })
  output
}
