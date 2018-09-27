#' Title
#'
#' @param hawkes 
#'
#' @return
#' @export
#'
#' @examples
#' @export
whittle <- function(data, binSize, trunc=5) {
  model <- new(ExpHawkes, data, binSize)
  
  n <- length(data)
  
  # Periodogram
  dft <- fft(data - mean(data))   # Since centering, need to ignore omega=0 frequency
  I <- Mod(dft)^2 / n
  omega <- 0:(n-1) * 2 * pi / n
  
  # Whittle pseudo likelihood function
  wlik <- function(param_) {
    param <- param_
    param[2] <- param_[2] * param_[3]
    model$param <- param
    return( -model$wlik(I, trunc) )
    # spec <- sapply(omega, function(w) {
    #   model$gammaf1(w, trunc)
    # })
    # return( sum(log(spec) + I / spec) )
  }
  
  opt <- optim(par=c(1, .5, 2), fn = wlik, lower = rep(.0001, 3), upper = c(2e16, .9999, 2e16), method = "L-BFGS-B")
  
  return( opt )
}

#' whittle_cov
#'
#' @param hawkes 
#'
#' @return
#' @export
#'
#' @examples
#' @export
whittle_cov <- function(data, binSize) {
  n <- length(data)
  
  # Periodogram
  dft <- fft(data - mean(data))
  I <- Mod(dft)^2 / n
  omega <- 0:(n-1) * 2 * pi / n
  
  # Whittle pseudo likelihood function
  # using debiased
  # Da Fonseca, J., & Zaatour, R. (2014). Hawkes process: Fast calibration, application to trade clustering, and diffusive limit. Journal of Futures Markets (Vol. 34). http://doi.org/10.1002/fut.21644
  wlik <- function(param_) {
    param <- param_
    param[2] <- param_[2] * param_[3]
    model$param <- param
    return( -model$wlikCov(I) )
    # hvar <- hawkes::jumpVariance(param[1], param[2], param[3], binSize)
    # hcov <- sapply(1:(n-1), function(tau) {
    #   (1 - tau / n) * hawkes::jumpAutocorrelation(param[1], param[2], param[3], binSize, tau-1)
    # }) * hvar
    # spec <- 2 * Re(sapply(omega, function(w) {
    #   sum( hcov * exp( -1i * w * 1:(n-1) ) )
    # })) + hvar
    # return( sum(log(spec) + I / spec) )
  }
  
  opt <- optim(par=c(1, .5, 2), fn = wlik, lower = rep(.0001, 3), upper = c(2e16, .9999, 2e16), method = "L-BFGS-B")
  
  return( opt )
}
