summaryMap <- function(submaps)
{
  if(class(submaps[[1]])[1] != "snps.matrix" & class(submaps[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps to eat.")  
  
  
  snps <- unlist(sapply(submaps, function(x) x@map$id)) #every markers that has been choosen on submaps
  b <- as.data.frame(table(snps),  stringsAsFactors=FALSE)
  b
}

