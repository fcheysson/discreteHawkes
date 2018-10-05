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
    param <- param_
    param[2] <- param[2] * param[3]
    model$param <- param
    return( -model$whittleLik(I, trunc) )
  }
  
  opt <- optim(par=init, fn = wlik, hessian = TRUE, lower = rep(.0001, 3), upper = c(Inf, .9999, Inf), method = "L-BFGS-B", ...)
  
  derivative <- matrix(byrow = TRUE, nrow = 3, ncol = 3,
                       data = c(1,          0, 0,
                                0, opt$par[3], 0,
                                0, opt$par[2], 1))
  vcov <- t(derivative) %*% solve(opt$hessian) %*% derivative
  
  model$vcov <- vcov
  model$opt <- opt

  return( model )
}

#' Title
#'
#' @param trunc 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
whittle_ <- function(model, trunc=5, ...) {
  data <- model$ddata
  n <- length(data)
  
  # Periodogram
  dft <- fft(data - mean(data))
  I <- Mod(dft)^2 / n
  
  # Whittle pseudo likelihood function (for optim)
  wlik <- function(param_) {
    param <- param_
    param[2] <- param[2] * param[3]
    model$param <- param
    return( -model$whittleLik(I, trunc) )
  }
  
  opt <- optim(par=model$param, fn = wlik, hessian = TRUE, lower = rep(.0001, 3), upper = c(Inf, .9999, Inf), method = "L-BFGS-B", ...)
  
  derivative <- matrix(byrow = TRUE, nrow = 3, ncol = 3,
                       data = c(1,          0, 0,
                                0, opt$par[3], 0,
                                0, opt$par[2], 1))
  vcov <- t(derivative) %*% solve(opt$hessian) %*% derivative
  
  model$vcov <- vcov
  model$opt <- opt
}
