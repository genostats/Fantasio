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
    
    #optimization of h(alpha)
    res <- optim(par = 1, fn = h, method="L-BFGS-B", lower = 0, upper = 1, control = list(fnscale = -1)) 
    #        optim( last_theta, f, gradf, method="L-BFGS-B", lower = c(0,0), upper = c(Inf, 1), control = list(fnscale = -1))
    HFLOD[j,1]<-res$value
    HFLOD[j,2]<-res$par        
  }
  x@HFLOD <- HFLOD
  return(x)
}

