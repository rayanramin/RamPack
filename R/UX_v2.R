#' UX version2
#' @description combines UX function with detect_context_stranded
#' new function to calculate uracilation index at the selected context.
#' @param X BRC-style data.frame with columns (POS,REF,COV,A,T,C,G)
#' @param context context to be detected, corrently only "NC","TC","WRC" and "WRCY" are covered.
#' @example UX_v2(A_BRC, "WRC") -> A_BRC$UI
#' @aliases UX_V2 ux_v2 UX_v2
#' @export


UX_v2 <- function(X, context="NC"){
  cntxt <- RamPack::detect_context_stranded(X,context)
  t <- ifelse(cntxt==1, as.integer(X$T), 0)/X$COV
  nc <- sum(cntxt==1)
  AllC <- sum(t)/nc * 1000
  a <- ifelse(cntxt==-1, as.integer(X$A), 0)/X$COV
  ng <- sum(cntxt==-1)
  AllG <- sum(a)/ng * 1000
  u <- (a + t)
  z <- sum(cntxt%in%c(-1,1))
  AllU <- sum(u)/z * 1000
  eps <- 0.001
  bias <- log((AllC+eps)/(AllG+eps),2)
  return(list(CtoT_indx = AllC, GtoA_indx = AllG, U_indx = AllU, strand_bias=bias ,
              Number_of_context_top_strand= nc, Number_of_context_bottom_strand = ng, Total_number_of_detected_context = z ))
}


#' UX version3
#' @description combines UX function with detect_context_stranded, it filters for position with coverage > 0
#' new function to calculate uracilation index at the selected context.
#' @param X BRC-style data.frame with columns (POS,REF,COV,A,T,C,G)
#' @param context context to be detected, corrently only "NC","TC","WRC" and "WRCY" are covered.
#' @example UX_v3(A_BRC, "WRC") -> A_BRC$UI
#' @aliases UX_V3 ux_v3 UX_v3
#' @export


UX_v3 <- function(X, context="NC"){
  ### new function to calculate uracilation index at the selected context.
  cntxt <- RamPack::detect_context_stranded(X,context)
  t <- ifelse(cntxt==1 & X$COV >0, as.integer(X$T)/X$COV, 0)
  nc <- sum(cntxt==1)
  AllC <- sum(t)/nc * 1000
  a <- ifelse(cntxt==-1 & X$COV >0, as.integer(X$A)/X$COV, 0)
  ng <- sum(cntxt==-1)
  AllG <- sum(a)/ng * 1000
  u <- (a + t)
  z <- sum(cntxt%in%c(-1,1))
  AllU <- sum(u)/z * 1000
  eps <- 0.001
  bias <- log((AllC+eps)/(AllG+eps),2)
  return(list(CtoT_indx = AllC, GtoA_indx = AllG, U_indx = AllU, strand_bias=bias ,
              Number_of_context_top_strand= nc, Number_of_context_bottom_strand = ng, Total_number_of_detected_context = z ))
}




#' UX version4
#' @description combines UX function with detect_context_stranded, it filters for position with coverage > 0
#' new function to calculate uracilation index at the selected context.
#' @param X BRC-style data.frame with columns (POS,REF,COV,A,T,C,G)
#' @param test_vector numeric vector with 1,-1 and 0 values indicating which rows to include, similar to the output of Detect_Context
#' @example UX_v3(A_BRC, "WRC") -> A_BRC$UI
#' @aliases UX_V3 ux_v3 UX_v3
#' @export


UX_v4 <- function(X, test_vector){
  cntxt <- test_vector
  t <- ifelse(cntxt==1, as.integer(X$T), 0)/X$COV
  nc <- sum(cntxt==1)
  AllC <- sum(t)/nc * 1000
  a <- ifelse(cntxt==-1, as.integer(X$A), 0)/X$COV
  ng <- sum(cntxt==-1)
  AllG <- sum(a)/ng * 1000
  u <- (a + t)
  z <- sum(cntxt%in%c(-1,1))
  AllU <- sum(u)/z * 1000
  eps <- 0.001
  bias <- log((AllC+eps)/(AllG+eps),2)
  return(list(CtoT_indx = AllC, GtoA_indx = AllG, U_indx = AllU, strand_bias=bias ,
              Number_of_context_top_strand= nc, Number_of_context_bottom_strand = ng, Total_number_of_detected_context = z ))
}


