#' discreteHawkes: To Be Continued.
#'
#' The foo package provides three categories of important functions:
#' foo, bar and baz.
#' 
#' @section Foo functions:
#' The foo functions ...
#'
#' @docType package
#' @name discreteHawkes
#' @useDynLib discreteHawkes
NULL

loadModule("MyModule", TRUE)

#' function that calls the class of Rcpp
#'
#' @export
exportRcppClass <- function(class)
{
  res = new(class)
  return(res)
}