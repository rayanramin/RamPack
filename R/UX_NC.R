#' Calculating Uracilation index
#'
#' input's column names should be:("Chr.","POS","REF","COV","A","C","G","T")
#' output is a list:  CtoT_indx , GtoA_indx , U_indx , Number_of_Cs , Number_of_Gs , Number_of_C_G ,  Out = UDF
#' this runs very slow with big inputs
#' @param Input is readcount data.frame
#' @keywords Uracilation Fraction Index
#' @export
#' @examples
#' UX_NC(BRC) -> ui
#'@aliases ux_nc

UX_NC <- function(Y) {
  X <- Y
  # defining conditions    
  cond1 <- logical(nrow(X))
  cond2 <- logical(nrow(X))
  for (i in 1:nrow(X)) { cond1[i] <- ( X$REF[i] == "C")  } 
  for (i in 1:nrow(X)) { cond2[i] <- ( X$REF[i] == "G")  }
  # C to T changes only
  t <- ifelse(cond1, as.integer(X$T), 0) / X$COV
  nc <- sum(cond1)
  AllC <- sum(t) / nc * 1e3
  # G to A changes only
  a <- ifelse(cond2, as.integer(X$A), 0) / X$COV
  ng <- sum(cond2)
  AllG <- sum(a) / ng * 1e3
  # All changes
  u <- (a+t) 
  z <- sum(cond2 | cond1)
  AllU <- sum(u) / z * 1e3
  
  UDF <-  data.frame("POS" = X$POS , "CtoT" = t ,  "GtoA" = a , "CtoT_GtoA" = u)
  
  return(list(CtoT_indx = AllC , GtoA_indx = AllG ,  U_indx = AllU , Number_of_Cs = nc, Number_of_Gs = ng , Number_of_C_G = z,  Out = UDF  ))
}
