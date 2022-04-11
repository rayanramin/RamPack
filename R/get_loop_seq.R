#' Get hairpin loop sequence
#' @description gets the sequence of hairpin loops from a data.frame
#' @param X input must be a data.frame of Hairpin_survey output with columns looplen,looppos, minus2,minus1,minus0,plus1,plus2,plus3
#' output will be of character strings like "(NATC.C.C)" with the targeted C shown as .C.
#' @example get_loop_seq(Data) -> Data$loopseq
#' @export



get_loop_seq <- function(X){
  ACGT <- c("A","C","G","T")
  tmp <- paste0("NNNNNNNNNN",ACGT[match(X$minus2,1:4)],ACGT[match(X$minus1,1:4)],ACGT[match(X$minus0,1:4)],".C.",ACGT[match(X$plus1,1:4)],ACGT[match(X$plus2,1:4)],ACGT[match(X$plus3,1:4)],"NNNNNNNNNN")
  out <- rep(NA,nrow(X))
  ifelse(X$looppos==0,paste0(substr(tmp,14,16),"(",paste0(substr(tmp,17,17+X$looplen-1),")")),
         ifelse(X$looppos==X$looplen+1,paste0("(",substr(tmp,14-X$looplen,13),")",substr(tmp,14,16)),
                paste0("(",substr(tmp,15-X$looppos,16+X$looplen-X$looppos),")")))
}
duplicated_all <- function(x){
  tmp <- x[duplicated(x)]
  x %in% tmp
}
