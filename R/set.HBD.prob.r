#x an msat.matrix with a and f computed
set.HBD.prob <- function(x)
{
  HBD_prob <- matrix(0.0, nrow = x@nrow, ncol = x@ncol) # matrice contenant les probabilités pour les individus
  dimnames(HBD_prob) <- list(rownames(HBD_prob, do.NULL = FALSE, prefix = "individual_"), colnames(HBD_prob) <- c(x@map$id))
  
  for (i in 1:x@nrow)
  {
    if(!is.na(x@a[i]) & !is.na(x@f[i]))
      HBD_prob[i,1:x@ncol] <-forward.backward(get.log.emiss(x, i), x@delta.dist, x@a[i], x@f[i] )[2,]
  }
  x@HBD.prob <- HBD_prob 
  return(x)
}

