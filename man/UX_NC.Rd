% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UX_NC.R
\name{UX_NC}
\alias{UX_NC}
\alias{ux_nc}
\title{Calculating Uracilation index}
\usage{
UX_NC(Y)
}
\arguments{
\item{Input}{is readcount data.frame}
}
\description{
input's column names should be:("Chr.","POS","REF","COV","A","C","G","T")
output is a list:  CtoT_indx , GtoA_indx , U_indx , Number_of_Cs , Number_of_Gs , Number_of_C_G ,  Out = UDF
this runs very slow with big inputs
}
\examples{
UX_NC(BRC) -> ui
}
\keyword{Fraction}
\keyword{Index}
\keyword{Uracilation}
