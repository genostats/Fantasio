#x an msat.matrix with a and f computed
set.HBD.prob <- function(x)
{
  l <- list()
  for(i in 1:length(x@f))
  {
    n <- 1
    if(x@f[i] > 0 )
    {
      l[[i]] <- i
      n <- n + 1
    }else{
      l[[i]] <- 0
    }
  }
  
  HBD_prob <- matrix(NA, nrow= sum(l > 0), ncol = x@ncol)#matrice contenant les probabilités pour les individus
  #dimnames(HBD_prob) <- list(rownames(HBD_prob, do.NULL = FALSE, prefix = "individual_"), colnames(HBD_prob) <- c(x@map$id))
  
  for (i in 1:nrow(HBD_prob))
  {
    if(!is.na(x@a[i]) & !is.na(x@f[i]))#a laisser ??? 
    {
      if(l[[i]] == 0)
      {
        next()#if the number at this position is equal to 0 pass the test
      }
        
      HBD_prob[i,1:x@ncol] <-forward.backward(get.log.emiss(x, i), x@delta.dist, x@a[i], x@f[i] )[2,]
    }
      
  }
  x@HBD.prob <- HBD_prob 
  return(x)
}


