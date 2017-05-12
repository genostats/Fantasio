# x = bedmatrix, n = the number of submap
submap <- function(x, n = 100, intensity = 10 , hotspot_version = "hg17")
{
  hotspot <- switch(hotspot_version,
                    hg17 = { data(hotspot_hg17); hotspot_hg17;},
                    hg18 = { data(hotspot_hg18); hotspot_hg18;},
                    hg19 = { data(hotspot_hg19); hotspot_hg19; })
  #Etape 1 : liste des hotspots pour le genome
  VI <- list()
  for ( i in 1:22)
  {
    chr_hotspot <- hotspot[which(hotspot[,1]==i),]
    w <- which(chr_hotspot$IntensitycMMb > intensity)
    segment <- cbind(c(0,chr_hotspot$End[w]),
                     c(chr_hotspot$Start[w],Inf) )
    VI[[i]] <- segment
  }
  
  #Etape 2 : liste des marqueurs pour le genome
  VII <- list()
  for( j in 1:22)
  {
    v <- x@snps$pos[x@snps$chr==j]
    VII[[j]] <- v
  }
  
  
  submap <- array(list(), c(n,1))
  for ( i in 1:n)
  {
    cat("creating submap number : ", i,"/", n, "\n" )
    spider <- createSubmap(x, VI, VII)
    submap[[i,1]] <- spider
  }
  cat("Done ! \n")
  return(submap)
}

