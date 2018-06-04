set.HFLOD <- function(submaps)
{
  HFLOD <- matrix(0.0, nrow = ncol(submaps@FLOD_recap), ncol = 2)
  
  for (j in 1:nrow(HFLOD))
  {
    # function h(alpha)
    h <- function(alpha)
      return(sum(log10( alpha*exp(submaps@FLOD_recap[,j]*log(10))+(1-alpha) ), na.rm = TRUE))

    #maximisation of h(alpha) ; 
    res <- optimize( h, c(0,1), maximum = TRUE, tol = 0.001 )
    
    
    HFLOD[j,1]<-res$objective # HFLOD = h(alpha max)
    HFLOD[j,2]<-res$maximum  #alpha max 
  }
  return(HFLOD)
}

