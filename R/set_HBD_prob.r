#x est une msat_matrix avec a et f calculés

set_HBD_Prob <- function(x)
{
  
  get.log.emiss <- function(x, i) #lignes de log.emiss qui concernent l'individu i
  {
    x@log.emiss[ c(2*i-1, 2*i), ]
  }
  
  HBD_prob <- matrix(0.0, nrow = x@nrow, ncol = x@ncol) # matrice contenant les probabilités pour les individus

  for (i in 1:x@nrow)
  {
    HBD_prob[i,1:x@ncol] <-forward.backward(get.log.emiss(x, i), x@delta.dist, x@a[i], x@f[i] )[2,]
  }
  
  
  return(HBD_prob) 
  
}
