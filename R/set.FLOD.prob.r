set.FLOD.prob <- function(x)
{
  FLOD_prob <- matrix(0.0, nrow = x@nrow, ncol = x@ncol)
  q = 0.0001
  for (i in 1:x@nrow)
  {
      FLOD_prob[i,1:x@ncol] <- log10( (x@HBD.prob[i,] + q * ( 1 - x@HDB.prob[i,]) / (x@f + q * ( 1 - x@f)) ))
  }
  x@FLOD.prob <- FLOD_prob 
  return(x)
  
}

