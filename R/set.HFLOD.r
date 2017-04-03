set.HFLOD <- function(x, step = 0.01)
{
  #initialiser alpha
  alpha <- matrix(0.0)
  nb <- 0
  i <- 0
  while(nb * step < 1)
  {
    alpha[i] <- nb * step
    nb = nb + 1
    i = i + 1
  }
  alpha[i] <- 1
  
  
  # calculer les HFLOD 
  HFLOD <- matrix(0.0, nrow = x@ncol, ncol = x@nrow) #en ligne les marqueurs en colonnes les individus ??? 
  sum_HFLOD <- 0 
  last <- function(x) { return( length(x) ) } #obtenir le dernier elements de alpha
  
  for (i in 1:x@ncol)
  {
    i <- 0
    while(i != last(alpha))
    {
      x <- log10(alpha[i] * exp(x@FLOD[i] * log(10)) + (1 - alpha[i]))
      sum_HFLOD = sum_HFLOD + x
    }
    HFLOD[i,1:x@nrow] <- sum_HFLOD
  }
  
  
  x@HFLOD <- HFLOD 
  return(x)
  
}

