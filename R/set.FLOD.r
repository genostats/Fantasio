set.FLOD <- function(x, q = 0.00001)
{
  #Computation of FLOD score with the formula
  FLOD_prob <- matrix(0.0,nrow=nrow(x@HBD.prob), ncol=x@ncol)
  for (i in 1:nrow(FLOD_prob))
  {
      FLOD_prob[i,1:x@ncol] <- log10( (x@HBD.prob[i,] + q * ( 1 - x@HBD.prob[i,])) / (x@f[i] + q * ( 1 - x@f[i])) )
  }
  x@FLOD <- FLOD_prob 
  return(x)
  
}

#require(Matrix)
#?Matrix
#Matrix(nrow=4,ncol=10,sparse=TRUE)
#Matrix(0,nrow=4,ncol=10,sparse=TRUE)
#x <- Matrix(0,nrow=4,ncol=10,sparse=TRUE)
#x[1,] <- rep(0,10)
#xs