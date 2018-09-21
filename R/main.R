loadModule("MyModule", TRUE)

#' function that calls the class of Rcpp
#'
#' @export
exportRcppClass <- function(class)
{
  res = new(class)
  return(res)
}