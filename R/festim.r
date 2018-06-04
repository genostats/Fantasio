##################################################################################
#This function allow the computation of a, f, likelihood0, likelihood1 and p.lrt #                                                                                #
#!!! x : a submap matrix                                                         #                                       
#!!! verbose : if you want to have little infos about the computation process    #
#!!! debug : if you want to have precise infos about the computation process     #
#                                                                                #
#*** return the submap matrix with each slots filled                             #
##################################################################################


#probs=TRUE,

festim <- function(x, verbose=TRUE, debug=FALSE) {
  if(is.null(x@epsilon)) x <- set.log.emiss(x)
  N <- nrow(x)
  x@a <- numeric(N)
  x@f <- numeric(N)
  x@likelihood0 <- numeric(N)
  x@likelihood1 <- numeric(N)
  x@p.lrt <- numeric(N)

  # res <- data.frame(f = numeric(N), a = numeric(N), likelihood = numeric(N), convergence = numeric(N), famid = numeric(N), id = numeric(N))
  for(i in 1:nrow(x)) {
    logEmission <- get.log.emiss(x,i)
    likelihood0 <- .Call('festim_logLikelihood', PACKAGE = "FEstim", logEmission, x@delta.dist, 0.01, 0)

    if(verbose) cat("Estimation f and a for individual #",i,"\n")

    last_theta <- c(0.05, 0.01)
    last_likelihood <- .Call('festim_logLikelihood_gradient', PACKAGE = "FEstim", logEmission, x@delta.dist, last_theta[1], last_theta[2])

    f <- function(theta) { 
      if(theta[1] < 0) theta[1] <- 1e-2
      if(theta[2] < 0) { 
        theta[2] <- 0
        if(debug) cat("Negative f set to 0")
      }
      if(theta[2] > 0.999) { 
        theta[2] <- 0.999
        if(debug) cat("f larger than 0.999 set to 0.999")
      }
      if(debug) cat("a = ", theta[1], " f = ", theta[2])
      last_theta <<- theta
      last_likelihood <<- .Call('festim_logLikelihood_gradient', PACKAGE = "FEstim", logEmission, x@delta.dist, theta[1], theta[2])
      if(debug) cat(" likelihood = ", last_likelihood[1], " gradient = ", last_likelihood[-1], "\n")
      last_likelihood[1]
    }

    gradf <- function(theta) {
      if(theta[1] < 0) theta[1] <- 1e-2
      if(theta[2] < 0) { 
        theta[2] <- 0
      }
      if(theta[2] > 0.999) { 
        theta[2] <- 0.999
      }
      if(all(theta == last_theta)) return(last_likelihood[-1])
      warning("Optimization algorithm might been having a glitch")
      last_theta <<- theta
      last_likelihood <<- .Call('festim_logLikelihood_gradient', PACKAGE = "FEstim", logEmission, x@delta.dist, theta[1], theta[2])
      last_likelihood[-1]
    }

    xx <- optim( last_theta, f, gradf, method="L-BFGS-B", lower = c(1e-2,0), upper = c(Inf, 0.999), control = list(fnscale = -1))
    if(xx$convergence != 0) {
      warning("Individual #",i, ", id = ", x@ped$id[i], ", optimization algorithm did not converge");
      xx$par <- c(NA, NA)
    }

    x@a[i] <- xx$par[1]
    x@f[i] <- xx$par[2]
    x@likelihood1[i] <- xx$value
    x@likelihood0[i] <- likelihood0
    x@p.lrt[i] <- pchisq( 2*(xx$value - likelihood0), df = 1, lower.tail = FALSE)
  }
  #if(probs) {
  #  if(verbose) cat("Computing posterior HBD probabilities\n")
  #  #x <- set.HBD.prob(x, all=all.individual)
  #  if(verbose) cat("Computing FLOD and HFLOD\n")
  #  #x <- set.FLOD(x)
  #  #x <- set.HFLOD(x)
  #}
  x
}


