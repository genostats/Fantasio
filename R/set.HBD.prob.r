#x an msat.matrix with a and f computed
set.HBD.prob <- function(x, inbred_individuals = which(x@f > 0))
{
  l <- inbred_individuals
  
  HBD_prob <- matrix(NA, nrow= length(l), ncol = x@ncol)#HBD matrix
  dimnames(HBD_prob) <- list(rownames(HBD_prob, do.NULL = FALSE, prefix = "individual_"), colnames(HBD_prob) <- c(x@map$id))
  
  for (i in 1:nrow(HBD_prob))
  {
    if(!is.na(x@a[l[i]]) & !is.na(x@f[l[i]]))
    {
        HBD_prob[i,1:x@ncol] <-forward.backward(get.log.emiss(x, i), x@delta.dist, x@a[l[i]], x@f[l[i]] )[2,]
    }
      
  }
  x@HBD.prob <- HBD_prob 
  return(x)
}


