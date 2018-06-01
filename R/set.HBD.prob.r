# x a matrix with a and f computed 
set.HBD.prob <- function(x, which.inds = which(x@a <= 1)) {
  id    <- c(x@ped$id[which.inds])
  famid <- c(x@ped$famid[which.inds])
  
  HBD_prob <- matrix(NA, nrow = length(which.inds), ncol = x@ncol)#HBD matrix
  rownames(HBD_prob) <- id
  colnames(HBD_prob) <- c(x@map$id)
  
  for (i in 1:nrow(HBD_prob)) {
    i1 <- which.inds[i]
    if(!is.na(x@a[i1]) & !is.na(x@f[i1]))
      HBD_prob[i,1:x@ncol] <-forward.backward(get.log.emiss(x, i1), x@delta.dist, x@a[i1], x@f[i1] )[2,]
  }
  x@HBD.prob <- HBD_prob 
  return(x)
}




