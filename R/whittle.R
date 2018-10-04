#' Title
#'
#' @param hawkes 
#'
#' @return
#' @export
#'
#' @examples
whittle <- function(data, binsize, init=c(1,.5,2), trunc=5, ...) {
  model <- new(ExpHawkes)
  model$ddata <- data
  model$binsize <- binsize
  
  n <- length(data)
  
  # Periodogram
  dft <- fft(data - mean(data))
  I <- Mod(dft)^2 / n

  # Whittle pseudo likelihood function (for optim)
  wlik <- function(param_) {
    # Maybe reparameterize (lambda, alpha/beta, beta) to put upper conditions on optim c(Inf, .9999, Inf)
    param <- param_
    param[2] <- param[2] * param[3]
    model$param <- param
    return( -model$whittleLik(I, trunc) )
  }
  
  opt <- optim(par=init, fn = wlik, hessian = TRUE, lower = rep(.0001, 3), upper = c(Inf, .9999, Inf), method = "L-BFGS-B", ...)
  
  # param <- c(opt$par[1], opt$par[2] * opt$par[3], opt$par[3])
  derivative <- matrix(byrow = TRUE, nrow = 3, ncol = 3,
                       data = c(1,          0, 0,
                                0, opt$par[3], 0,
                                0, opt$par[2], 1))
  vcov <- t(derivative) %*% solve(opt$hessian) %*% derivative
  
  model$vcov <- vcov
  model$opt <- opt

  return( model )
}

# #' whittle_cov
# #'
# #' @param hawkes 
# #'
# #' @return
# #' @export
# #'
# #' @examples
# whittle_cov <- function(data, binsize, init=c(1,.5,2), ...) {
#   model <- new(ExpHawkes)
#   model$ddata <- data
#   model$binsize <- binsize
# 
#   n <- length(data)
# 
#   # Periodogram
#   dft <- fft(data - mean(data))
#   I <- Mod(dft)^2 / n
# 
#   # Debiased Whittle pseudo likelihood function
#   # Shumway, R. H., & Stoffer, D. S. (2011). Time Series Analysis and Its Applications. http://doi.org/10.1007/978-1-4419-7865-3
#   # Da Fonseca, J., & Zaatour, R. (2014). Hawkes process: Fast calibration, application to trade clustering, and diffusive limit. Journal of Futures Markets (Vol. 34). http://doi.org/10.1002/fut.21644
#   wlik <- function(param_) {
#     param <- param_
#     param[2] <- param[2] * param[3]
#     model$param <- param
#     dcov <- (1-0:(n-1)/n) * c( model$var(), model$cov_(0:(n-2)) )
#     dspectrum <- 2 * Re( fft(dcov) ) - model$var()
#     return( +sum( log(dspectrum) + I / dspectrum ) )
#   }
# 
#   opt <- optim(par=init, fn = wlik, lower = rep(.0001, 3), upper = c(Inf, .9999, Inf), method = "L-BFGS-B", ...)
# 
#   return( opt )
# }
