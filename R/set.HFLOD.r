set.HFLOD <- function(x)
{
  HFLOD <- matrix(0.0, nrow = x@ncol, ncol = 2)
  for (j in 1:x@ncol)
  {
    # function h(alpha)
    h <- function(alpha)
    {
      return(sum(log10(alpha*exp(x@FLOD[,j] *log(10))+(1-alpha)), na.rm = TRUE))
      #return(sum(log(alpha * (exp(x@FLOD[,j]*log(10))-1) + 1)/log(10), na.rm = TRUE))
    }
  
    #optimization of h(alpha) ; 
    res <- optimize( h, c(0,1), maximum = TRUE, tol = 0.001 )
    
    
    HFLOD[j,1]<-res$objective # HFLOD = h(alpha max)
    HFLOD[j,2]<-res$maximum  #alpha max 
  }
  x@HFLOD <- HFLOD
  return(x)
}

