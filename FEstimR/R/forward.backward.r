forward.backward <- function(logEmission, delta.Distance, a, f) {
  if( a < 0 ) stop("a should be >= 0")
  if( f < 0 | f > 1) stop("f should be between 0 and 1")
  if( f == 0) 
    return( matrix( c(1,0), nrow = 2, ncol = ncol(logEmission) ) )
  if( f == 1) 
    return( matrix( c(0,1), nrow = 2, ncol = ncol(logEmission) ) )
  .Call('festim_forward_backward', PACKAGE = "FEstim", logEmission, delta.Distance, a, f)
}

