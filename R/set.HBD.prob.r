#x an matrix with a and f computed 
# inbred_individual à utiliser
set.HBD.prob <- function(x)
{
  l <- which(x@f > 0 & x@a < 1)
  
  id    <- c(x@ped$id[l])
  famid <- c(x@ped$famid[l])
  
  namevector <- c()
  
  for(i in 1:length(id))
  {
    namevector <- c(namevector, paste("individual", id[i], famid[i], sep = "_"))
  }
  
  HBD_prob <- matrix(NA, nrow = length(l), ncol = x@ncol)#HBD matrix
  dimnames(HBD_prob) <- list(rownames(HBD_prob) <- namevector, colnames(HBD_prob) <- c(x@map$id))
  
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




