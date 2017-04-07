set.HFLOD <- function(x, step = 0.01)
{
  #initialiser alpha
  alpha <- matrix(0.0)
  nb <- 0
  i <- 0
  while(nb * step < 1)
  {
    alpha[i] <- nb * step
    nb <- nb + 1
    i <- i + 1
  }
  alpha[i] <- 1
  
  
  # calculer les HFLOD pour chaque marqueurs 
  HFLOD <- matrix(0.0, nrow = x@ncol, ncol = 2) #en ligne les marqueurs en colonnes HFLOD et alpha
  
#######################################################################################################################################  

  temp <- matrix(alpha)
  #sum_HFLOD <- 0 
  
  ####### obtenir les FLOD pour l'individus i 
  #ind <- x@nrow
  #p <- 0
  #all_FLOD <- NULL
  
   

  last <- function(x) { return( length(x) ) } #obtenir le dernier elements de alpha
  FLOD_ind <- 0
  for (i in 1:x@ncol)
  {
      FLOD_ind <- x@FLOD[,i]#FLOD pour le marqueurs i
      j <- 1
      g <- 0
      while(j <= last(alpha))
      {
          g <- temp[j] * (exp(FLOD_ind * log(10))-1)+1
          temp[j] <- max(g, 0.0001, na.rm = TRUE)
          temp[j] <- log(temp[j])/log(10)
          
          for(v in 1:last(temp))
          {
            HFLOD[i] <- HFLOD[i] + temp[i]
          }
         j <- j + 1
          
          
        
    
      }
    }
    
   
  

  
  x@HFLOD <- HFLOD 
  return(x)
  
}

