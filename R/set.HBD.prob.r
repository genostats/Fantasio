#x an matrix with a and f computed 
# a modifier condition
set.HBD.prob <- function(x, condition = which(x@a <= 1))
{
  id    <- c(x@ped$id[condition])
  famid <- c(x@ped$famid[condition])
  
  namevector <- c()
  
  for(i in 1:length(id))
  {
    namevector <- c(namevector, id[i])
  }
  
  HBD_prob <- matrix(NA, nrow = length(condition), ncol = x@ncol)#HBD matrix
  dimnames(HBD_prob) <- list(rownames(HBD_prob) <- namevector, colnames(HBD_prob) <- c(x@map$id))
  
  for (i in 1:nrow(HBD_prob))
  {
    if(!is.na(x@a[condition[i]]) & !is.na(x@f[condition[i]]))
    {
        HBD_prob[i,1:x@ncol] <-forward.backward(get.log.emiss(x, i), x@delta.dist, x@a[condition[i]], x@f[condition[i]] )[2,]
    }
  }
  x@HBD.prob <- HBD_prob 
  return(x)
}




