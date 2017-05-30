
getMarkerChromosom <- function(map)
{
  #find the markers index in the interval
  submap <- c()

  
  #Step 4 : choose a random mkr
  for( i in 1:length(map))
  {
    b <- map[[i]]
    if(length(b) == 1) 
    {
      s <- b
    } else { 
        s <- sample(b, 1)
      }
   submap <- c(submap, s)
  }
  
  return(submap)
}





# x = a bed matrix
# return an" msat matrix without genotypes but with log emiss...

createSubmap <- function(x, mkr_map)
{
  submap <- c()
  for(chr in 1:22)
  {
    n <- mkr_map[[chr]]
    v <- getMarkerChromosom(n)
    submap <- c(submap, v)
  }
  
  map <- x@snps[submap , c("id","chr")]
  if(all(x@snps$dist == 0)) {
    map$distance <- x@snps$dist[submap]
  } else {
    map$distance <- x@snps$pos[submap]*1e-6
  }
  
  res <- new("msat.matrix", length(submap), nrow(x), 
             x@ped[,c("famid", "id", "father", "mother", "sex", "pheno")]
             ,matrix(nrow = 0, ncol = 0), map, matrix(0, nrow = 0 , ncol =0))
  
  res@log.emiss <- bed.logEmiss(x, submap, 1e-3)
  res@epsilon <- 1e-3
  #res <- festim(res)
  #res <- HBD.prob(res)
  #res <- FLOD.prob(res)
  #res <- set.HFLOD(res)
  return(res)
}


