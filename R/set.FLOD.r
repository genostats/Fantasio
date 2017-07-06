set.FLOD <- function(x, q = 0.00001)
{
  #Computation of FLOD score with the formula
  FLOD_prob <- matrix(0.0, nrow = x@nrow, ncol = x@ncol)
  for (i in 1:x@nrow)
  {
      FLOD_prob[i,1:x@ncol] <- log10( (x@HBD.prob[i,] + q * ( 1 - x@HBD.prob[i,])) / (x@f[i] + q * ( 1 - x@f[i])) )
  }
  x@FLOD <- FLOD_prob 
  return(x)
  
}

