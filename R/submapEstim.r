submapEstim <- function(submaps)
{
  estimation <- list()
  dfA <- data.frame("estimation of a.submap" = sapply(submaps, function(x) x@a))
  dfF <- data.frame("estimation of f.submap" = sapply(submaps, function(x) x@f))
  
  estimation[[1]] <- dfF
  estimation[[2]] <- dfA
   
  return(estimation)
}