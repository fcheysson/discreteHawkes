##############################################################################
#                        SIMULATION OF HAWKES PROCESS                        #
##############################################################################

# Hawkes process by thinning (Ogata's modified thinning algorithm)

#' @export
hawkes <- function(T, lambda, alpha, beta, lambda0=NULL) UseMethod("hawkes")

#' @export
hawkes.default <- function(T, lambda, alpha, beta, lambda0=NULL) {

	# If no lambda0 specified, take lambda
	# NOTE: in the future, find the correct distribution for lambda0
	if (is.null(lambda0)) lambda0 <- lambda
	
	# lambda0 must be higher than lambda
	# (else Ogata's algorithm does not work)
	else if (lambda0 < lambda) 
		stop("For Ogata's algorithm, lambda0 must be higher than lambda.")

	# Object to return	
	sim <- list(p=numeric(0), T=T, 
			lambda=lambda, alpha=alpha, beta=beta, lambda0=lambda0,
			n=0, A=numeric(0), M=lambda, 
			t=numeric(0), e=numeric(0), U=numeric(0), m=0)
	class(sim) <- "hawkes"

	# We do as if there was a point in t = 0 that we will
	# remove at the end of the algorithm
	M <- lambda0
	A <- (lambda0 - lambda) / alpha - 1
	p <- 0
	e <- NA
	t <- 0
	U <- NA

	# Initialize
	i <- 1
	n <- 1
	M[2] <- lambda0
	
	# Ogata's modified thinning algorithm
	while (t[n] < T) {
		n <- n + 1
		e[n] <- rexp(1, rate=M[n])
		t[n] <- t[n-1] + e[n]
		U[n] <- runif(1, min=0, max=M[n])

		if (t[n] < T && U[n] < (lambda + alpha * (A[i]+1) * exp(-beta*(t[n]-p[i])))) {
			A[i+1] <- exp(-beta*(t[n]-p[i]))*(1+A[i])
			p[i+1] <- t[n]
			M[n+1] <- lambda + alpha*(A[i+1]+1)
			i <- i + 1
		} else {
			M[n+1] <- lambda + (M[n] - lambda) * exp(-beta * e[n])
		}
	}
	
	# Remove the first artificial point
	sim$p <- p[-1]
	sim$M <- M[-1]
	sim$A <- A[-1]
	sim$t <- t[-1]
	sim$e <- e[-1]
	sim$U <- U[-1]
	sim$n <- i-1
	sim$m <- n-1

	return( sim )
}

#' @export
plot.hawkes <- function(hawkes, precision=1e3, ...) {
	# Conditional intensity
	matplot(z <- seq(0, hawkes$T, by=hawkes$T / precision), 
		intensity(hawkes, z), 
		type="l", ylim=c(0, max(hawkes$M)),
		xlab=expression(italic(t)), ylab=expression(italic(U)), ...)
	# Upper bound M of Ogata's modified thinning algorithm
	segments(x0=c(0, hawkes$t[-hawkes$m]), x1=c(hawkes$t[-hawkes$m], hawkes$T), 
		y0=hawkes$M[-1-hawkes$m], 
		lwd=4, col="blue")
	# All points considered and the corresponding value for U
	points(x=hawkes$t[-hawkes$m], 
		y=hawkes$U[-hawkes$m], 
		pch=ifelse(hawkes$U[-hawkes$m] > sapply(hawkes$t[-hawkes$m], function(t) {intensity(hawkes, t)}), 3, 1),
		col=ifelse(hawkes$U[-hawkes$m] > sapply(hawkes$t[-hawkes$m], function(t) {intensity(hawkes, t)}), "firebrick1", "green4"))
	# The resulting points of the Hawkes process
	points(x=hawkes$p, y=rep(0, hawkes$n), col="green4", pch=15)
	# Nice lines showing which point considered resulted in which point of the Hawkes process
	segments(x0=hawkes$p, y0=rep(0, hawkes$n), 
		y1=(hawkes$U[-hawkes$m])[which(hawkes$U[-hawkes$m] <= sapply(hawkes$t[-hawkes$m], function(t) {intensity(hawkes, t)}))],
		col="green4", lty=2)
	# Legends
	legend("topleft", legend=c(expression(italic(lambda(t))), expression(italic(M)), expression(Accepted), expression(Rejected)),
		lty=c(1,1,NA,NA), col=c("black", "blue", "green", "red"), lwd=c(1,4,NA,NA), pch=c(NA, NA, 1, 3))
}

#' @export
intensity <- function(object, t) UseMethod("intensity")

#' @export
intensity.hawkes <- function(hawkes, t) {
	# If outside of bounds, return error
	if (any(t < 0) || any(t > hawkes$T)) 
		stop("t is out of bound.")
	# Create vector of past closest points of hawkes
	index <- sapply(t, function(j) {
		Position(function(p) p >= j, hawkes$p, nomatch=hawkes$n+1)-1
	})
	int <- rep(hawkes$lambda, length(t)) 
	pind <- index > 0
	int[!pind] <- int[!pind] + 
		(hawkes$lambda0 - hawkes$lambda) * exp(-hawkes$beta * t[!pind])
	int[pind] <- int[pind] + 
		hawkes$alpha * (hawkes$A[index[pind]]+1) * 
		exp(- hawkes$beta * (t[pind] - hawkes$p[index[pind]]))
	return( int )
}

#' @export
compensator <- function(object, t, ...) UseMethod("compensator")

#' @export
compensator.hawkes<- function(hawkes, t, ...) {
	# If outside of bounds, return error
	if (any(t < 0) || any(t > hawkes$T)) 
		stop("t is out of bound.")
	# Create vector of past closest points of hawkes
	index <- sapply(t, function(j) {
		Position(function(p) p >= j, hawkes$p, nomatch=hawkes$n+1)-1
	})
	int <- hawkes$lambda * t
	pind <- (index > 0)
	int[!pind] <- int[!pind] +
		(hawkes$lambda0 - hawkes$lambda) / hawkes$beta * (1 - exp(-hawkes$beta * t[!pind]))
	int[pind] <- int[pind] + (hawkes$lambda0 - hawkes$lambda) / hawkes$beta +
		hawkes$alpha / hawkes$beta * 
		(index[pind] - 
		exp(-hawkes$beta * (t[pind] - hawkes$p[index[pind]])) * (hawkes$A[index[pind]] + 1))
	return( int )
}

#' @export
residuals.hawkes <- function(hawkes) {
	return( hawkes$lambda * hawkes$p + (hawkes$lambda0 - hawkes$lambda) / hawkes$beta +
		hawkes$alpha / hawkes$beta * (1:hawkes$n - 1 - hawkes$A) )
}
