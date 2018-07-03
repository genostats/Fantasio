submapEstim <- function(submaps)
{
  if(class(submaps[[1]])[1] != "snps.matrix" & class(submaps[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps to eat.") 
  
  estimation <- list()
  dfA <- data.frame("estimation of a.submap" = sapply(submaps, function(x) x@a))
  dfF <- data.frame("estimation of f.submap" = sapply(submaps, function(x) x@f))
  
  estimation[[1]] <- dfF
  estimation[[2]] <- dfA
   
  return(estimation)
}