#x est une msat_matrix avec a et f calculés
HBD.prob <- function(x)
{
  HBD_prob <- matrix(0.0, nrow = x@nrow, ncol = x@ncol) # matrice contenant les probabilités pour les individus

  for (i in 1:x@nrow)
  {
    HBD_prob[i,1:x@ncol] <-forward.backward(get.log.emiss(x, i), x@delta.dist, x@a[i], x@f[i] )[2,]
  }
  x@HBD.prob <- HBD_prob 
  return(x)
}
