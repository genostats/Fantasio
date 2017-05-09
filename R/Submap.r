# x = bedmatrix, n = the number of submap
submap <- function(x, n = 100, intensity = 10 , hotspot_version = "hg17")
{
  submap <- array(list(), c(n,1))
  for ( i in 1:n)
  {
    cat("creating submap number : ", i, "\n" )
    spider <- createSubmap(x, intensity, hotspot_version)
    submap[[i,1]] <- spider
  }
  cat("Done ! \n")
  return(submap)
}

#3 min pour 5 sous cartes en moyennes
#60 min pour les 100 sous cartes en moyennes