
# x = a msat matrix
festim <- function(x, verbose = TRUE) {
  if(is.null(x@epsilon)) x <- set.log.emiss(x)
  N <- nrow(x)
  res <- data.frame(f = numeric(N), a = numeric(N), likelihood = numeric(N), convergence = numeric(N))
  for(i in 1:nrow(x)) {
    if(verbose) cat("Individual #",i,"\n")
    logEmission <- x@log.emiss[ c(2*i-1,2*i), ]
    last_theta <- c(0.05, 0.01)
    last_likelihood <- .Call('festim_logLikelihood_gradient', PACKAGE = "FEstim", logEmission, x@delta.dist, last_theta[1], last_theta[2])

    f <- function(theta) { 
      if(theta[1] < 0) theta[1] <- 0
      if(verbose) cat("theta = ", theta)
      last_theta <<- theta
      last_likelihood <<- .Call('festim_logLikelihood_gradient', PACKAGE = "FEstim", logEmission, x@delta.dist, theta[1], theta[2])
      if(verbose) cat(" likelihood + grad = ", last_likelihood, "\n")
      last_likelihood[1]
    }

    gradf <- function(theta) {
      if(theta[1] < 0) theta[1] <- 0
      if(all(theta == last_theta)) return(last_likelihood[-1])
      warn("ououps")
      last_theta <<- theta
      last_likelihood <<- .Call('festim_logLikelihood_gradient', PACKAGE = "FEstim", logEmission, x@delta.dist, theta[1], theta[2])
      last_likelihood[-1]
    }

    xx <- optim( last_theta, f, gradf, method="L-BFGS-B", lower = c(0,0), upper = c(Inf, 1), control = list(fnscale = -1))
    res$a[i] <- xx$par[1]
    res$f[i] <- xx$par[2]
    res$likelihood[i] <- xx$value
    res$convergence[i] <- xx$convergence    
  }
  res
}


