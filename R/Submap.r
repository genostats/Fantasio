# x = bedmatrix, n = the number of submap
submap <- function(x, n = 100, intensity = 10 , hotspot_version = "hg17")
{
  hotspot <- switch(hotspot_version,
                    hg17 = { data(hotspot_hg17); hotspot_hg17;},
                    hg18 = { data(hotspot_hg18); hotspot_hg18;},
                    hg19 = { data(hotspot_hg19); hotspot_hg19; })
  
  #Step 1 : list of all the genome's hotspot
  
  #print("Gathering all informations about genome's hotspot : ")
  VI <- list()
  for ( i in 1:22)
  {
   # print(".")
    chr_hotspot <- hotspot[which(hotspot[,1]==i),]
    w <- which(chr_hotspot$IntensitycMMb > intensity)
    segment <- cbind(c(0,chr_hotspot$End[w]),
                     c(chr_hotspot$Start[w],Inf) )
    VI[[i]] <- segment
  }
  
  #Step 2 : list of all the genome's markers
  
  #print("Gathering all the genome's markers : ")
  VII <- list()
  for( j in 1:22)
  {
    #print(".")
    v <- x@snps$pos[x@snps$chr==j]
    VII[[j]] <- v
  }
  
  #Step 3 : list of all the marker for a segment for the genome
  
  #print("Gathering all the segments for the genome : ")
  VIII <- list()
  for(i in 1:22)
  {
    #print(".")
    chr_segment <- VI[[i]]
    mkr <- VII[[i]]
    chr <- list()
    for( j in 1:nrow(chr_segment))
    {
      b <- which(mkr > chr_segment[j,1] & mkr < chr_segment[j,2])
      if (length(b)== 0) next
      
      chr[[j]] <- b
      
    }
    VIII[[i]] <- chr
  }
  
  
  
  submap <- array(list(), c(n,1))
  for ( i in 1:n)
  {
    cat("creating submap number : ", i,"/", n, "\n" )
    spider <- createSubmap(x, VIII)
    submap[[i,1]] <- spider
  }
  cat("Done ! \n")
  return(submap)
}

