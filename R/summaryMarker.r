summaryMarker <- function(submaps, bedmatrix)
{
  if(class(submaps[[1]])[1] != "snps.matrix" & class(submaps[[1]])[1] != "hotspots.matrix")
    stop("need either an hotspots.segments list of submaps or a snps.segments list of submaps to eat.") 
  
  if(class(bedmatrix)[1] != "bed.matrix")
  {
    stop("Need a bed.matrix to eat")
  }
  
  
  b <- summaryMap(submaps)
  res <- c()
  for(i in 1:length(submaps))
  {
    res <- c(res, length(which(b$Freq == i)))
  }
  zero <- length(bedmatrix@snps$chr) - sum(res)
  
  df <- data.frame(
    number_of_time_picked = 0:length(submaps),
    number_of_markers = c(zero, res)
  )
  df
}