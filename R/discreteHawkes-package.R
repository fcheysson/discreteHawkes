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

loadModule("HawkesModule", TRUE)

.onUnload <- function (libpath) {
  library.dynam.unload("discreteHawkes", libpath)
}