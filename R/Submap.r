# x = bedmatrix, n = the number of submap
submap <- function(x, n = 100, intensity = 10 , hotspot_version = "hg17")
{
  submap <- matrix(list(), nrow=n, ncol=1)
  for ( i in 1:n)
  {
    cat("creating submap number : ", i, "\n" )
    spider <- createSubmap(x, intensity, hotspot_version)
    submap[[i,1]] <- spider
  }
  cat("Finish ! \n")
  return(submap)
}