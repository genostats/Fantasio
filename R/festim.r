#' Computation of a and f 
#' 
#' This function is used to compute a and f and other statistical values (p.lrt, likelihood0, likelihood1)
#' 
#' @param x a submap matrix   
#' @param verbose whether you want informations about the computation process (default is TRUE)
#' @param debug  whether you want advanced output about the computation process (default is FALSE)
#' 
#' @details This function is using the optim function with an L-BGFS method to find the best value for a and f
#' 
#' @return returns the submap matrix with the slot a, f, p.lrt, likelihood0 and likelihood1 filled.
#' 
#' 
#' @examples  
#' #Please refer to vignette 
#'
#' 
#' @export
festim <- function(x, n.cores = 1, verbose = FALSE, debug = FALSE) { # should be handled through a method ...
  if(is(x, "fMatrix"))
    festim_fmatrix(x, verbose, debug)
  else if(is(x, "atlas"))
    festim_submapslist(x, n.cores, verbose, debug)
  else
    stop("festim should be applied of objects of class fMatrix or atlas")
}

festim_submapslist <- function(x, n.cores, verbose, debug) {

  if(n.cores != 1 & .Platform$OS.type != "unix") {
    warning("FORK cluster unavailable only one core used")
    n.cores <- 1
  }

  if(n.cores == 1) {
    new_submaps_list <- lapply(x@submaps_list, festim_fmatrix, verbose = verbose, debug = debug)
  } else {
    cl <- makeForkCluster(n.cores) 
    new_submaps_list <- parLapply(cl, x@submaps_list, festim_fmatrix, verbose = FALSE, debug = FALSE)
    stopCluster(cl)
  }
  x@submaps_list <- new_submaps_list
  x
}


# this functions uses slots that are common to 
# it msat matrices, snpsMatrices (submaps by distance), HotspotsMatrices
festim_fmatrix <- function(x, verbose = TRUE, debug = FALSE) {
  N <- nrow(x)
  x@a <- numeric(N)
  x@f <- numeric(N)
  x@likelihood0 <- numeric(N)
  x@likelihood1 <- numeric(N)
  x@p.lrt <- numeric(N)

  # res <- data.frame(f = numeric(N), a = numeric(N), likelihood = numeric(N), convergence = numeric(N), famid = numeric(N), id = numeric(N))
  for(i in seq_len(nrow(x))) {
    logEmission <- getLogEmiss(x,i)
    likelihood0 <- .Call('festim_logLikelihood', PACKAGE = "Fantasio", logEmission, x@delta.dist, 0.01, 0)

    if(verbose) cat("Estimation f and a for individual #",i,"\n")

    last_theta <- c(0.05, 0.05)
    last_likelihood <- .Call('festim_logLikelihood_gradient', PACKAGE = "Fantasio", logEmission, x@delta.dist, last_theta[1], last_theta[2])

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
      last_likelihood <<- .Call('festim_logLikelihood_gradient', PACKAGE = "Fantasio", logEmission, x@delta.dist, theta[1], theta[2])
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
      last_likelihood <<- .Call('festim_logLikelihood_gradient', PACKAGE = "Fantasio", logEmission, x@delta.dist, theta[1], theta[2])
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
    x@p.lrt[i] <- pchisq( 2*(xx$value - likelihood0), df = 2, lower.tail = FALSE)
  }
  x
}


