# a summary for a particular submap because depending on the submap markers differs
marker.summary <- function(x)
{
  
  res_HBD <- c()
  res_FLOD <- c()
  
  for (i in 1:ncol(x@HBD.prob))
  {
    res_HBD <- c(res_HBD, mean(x@HBD.prob[,i], na.rm = TRUE))
    res_FLOD <- c(res_FLOD, mean(x@FLOD[,i], na.rm = TRUE))
  }
  
  df <- data.frame(
    marker           = x@map$id,
    chromosom        = x@map$chr,
    distance_in_chr  = x@map$distance,
    HBD              = res_HBD,
    FLOD             = res_FLOD
  )
}

