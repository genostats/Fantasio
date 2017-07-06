
getMarkerChromosom <- function(map)
{
  #find the markers index in the interval
  submap <- c()
  
  #Step 4 : choose a random mkr
  for( i in 1:length(map))  # on parcourt les segments
  {
    if(length(map[[i]]) == 0) next  # segment vide...
    if(length(map[[i]]) == 1) 
    {
      s <- map[[i]]
    } else { 
      s <- sample(map[[i]], 1)
    }
    submap <- c(submap, s)
  }
  
  return(submap)
}

# x = a bed matrix
# mkr_map = list qui donne pour chaque chromosome la carte des segments
#  une carte de segments = une liste d'indices...
# return an" msat matrix without genotypes but with log emiss...

createSubmap <- function(x, mkr_map, epsilon = 1e-3)
{
  submap <- c()
  for(chr in 1:22)
  {
    map <- mkr_map[[chr]]
    v <- getMarkerChromosom(map)
    submap <- c(submap, v)
  }
  
  map <- x@snps[submap , c("id","chr")]
  if(all(x@snps$dist == 0)) {
    map$distance <- x@snps$dist[submap]
  } else {
    map$distance <- x@snps$pos[submap]*1e-6
  }
  
  log.emiss <- bed.logEmiss(x, submap, epsilon)

  new("submap.matrix", length(submap), nrow(x), submap, 
      x@ped[,c("famid", "id", "father", "mother", "sex", "pheno")],
      map, log.emiss, epsilon)

}


