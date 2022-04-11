#'Normalized Differntial Coverage
#' @description NDC1 is the original version that uses the whole genome for normalization ; NDC2 divide the genome into smaller parts for normalization ; NDC3 is similar to NDC2 but better for when coverage is 0 at many positions.
#'
#' default values are optimized based on e.coli experiments
#'
#' @param x coverage vector from the treatment sample
#' @param y coverage vector from the control sample
#' @param w width of normalization window (for NDC2|3); w=1e5
#' @param b step size of w; b =1e4
#' @param mavn size of moving average window; mavn = 120
#' @param NDC_Version is the version of NDC; NDC_Version = 2
#' @import dplyr
#' @keywords Normalized differential coverage NDC NDC1 NDC2 NDC3
#' @export
#' @examples
#' NDC(treatment$cov , control$cov , w=5e5 , b = 2e4 , mavn = 150 ) -> NDC
#'@aliases NDC ndc



NDC <- function(x, y, w=1e5 , b =1e4 , mavn = 120 , NDC_Version = 2){
#Functions (default values are set based on optimized signal to noise ratio)
#NDC1 : old version of NDC good for genomes with more uniform coverage
mav <- function(x,n=100){stats::filter(x,rep(1/n,n),sides=2, circular=T)}
NDC1 <- function(x, y){
    mx <- mean(x)
    my <- mean(y)
    ((mav(x) /mx ) - (mav(y)/ my ))
}
#NDC2 : Newer version of NDC , more sensitive and higher signal to noise ratio
NDC2 <- function(x, y, w , b  , mavn ){
    mx <- x %>% zoo::rollapply( width = w , by = b , FUN = "mean" , fill = NA) %>% zoo::na.fill("extend")
    my <- y %>% zoo::rollapply( width = w , by = b , FUN = "mean" , fill = NA) %>% zoo::na.fill("extend")
    ((mav(x, mavn) /mx ) - (mav(y ,mavn )/ my ))
}
#NDC3 : Similar to NDC2 , good for genomes with lots of zero coverage areas.
NDC3 <- function(x, y, w , b , mavn  ){
    mxp1 <- c(x +1) %>% zoo::rollapply( width = w , by = b , FUN = "mean" , fill = NA) %>% zoo::na.fill("extend")
    myp1 <- c(y +1) %>% zoo::rollapply( width = w , by = b , FUN = "mean" , fill = NA) %>% zoo::na.fill("extend")
    ((mav(x, mavn) /mxp1 ) - (mav(y ,mavn )/ myp1 ))
}

if( NDC_Version == 1){
    NDC1(x, y)
} else if (NDC_Version == 2) {
   NDC2(x, y, w , b  , mavn)
} else if (NDC_Version == 3) {
   NDC3(x, y, w , b  , mavn)
}
}
