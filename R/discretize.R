#' discrete
#'
#' @param data 
#' @param T
#' @param length 
#' @param binSize 
#'
#' @return
#' @export
#'
#' @examples
#' @export
discrete <- function(data, T, length=NULL, binSize=NULL) {
  if (is.null(length) & is.null(binSize))
    stop("One of length or binSize must be specified.")
  if (!is.null(binSize) & !is.null(length))
    cat("Warning: Both length and binSize specified. binSize will be ignored.\n")
  if (!is.null(length)) {
    if (length <= 0)
      stop("length is less or equal to zero.")
    ti <- seq(0, T, length.out=length+1)
    bin <- sapply(1:length, function(i) {
      sum(data > ti[i] & data < ti[i+1])
    })
  } else {
    if (T %% binSize != 0)
      cat("Warning: T is not a multiplier of binSize. Last bin will be discarded.\n")
    ti <- seq(0, T, by=binSize)
    bin <- sapply(1:(length(ti)-1), function(i) {
      sum(data > ti[i] & data < ti[i+1])
    })
  }
  return(bin)
}