getSegment <- function(chr, intensity, hotspot_version) 
{
  #select the correct hotspot file
  hotspot <- switch(hotspot_version,
                    hg17 = { data(hotspot_hg17); hotspot_hg17;},
                    hg18 = { data(hotspot_hg18); hotspot_hg18;},
                    hg19 = { data(hotspot_hg19); hotspot_hg19; })
  
  #save all  current chromosom's hotspots
  chr_hotspot <- hotspot[ which(hotspot[,1]==chr),]
  #hotspot_hg17[which(hotspot_hg17==chr),]
  
  #save only the hotspot with the right intensity
  w <- which(chr_hotspot$IntensitycMMb > intensity)
  #w <- match(chr_hotspot$IntensitycMMb > intensity, chr_hotspot$IntensitycMMb, nomatch = 0)
  #create a matrix with segments delimitations
  segment_table <- cbind( c(0,chr_hotspot$End[w]),
                          c(chr_hotspot$Start[w],Inf) )
  return(segment_table)
}

getMarkersChromosom <- function(x, chr,intensity, hotspot_version)
{
  submap <- c() #save the marker's index for the current chromosom
  v <- x@snps$pos[x@snps$chr==chr]  
  segment_table <- getSegment(chr, intensity, hotspot_version) #get the hotspot for this chr
  
  #find the markers index in the interval
  for ( i in 1:nrow(segment_table))
  {
    b <- which(v > segment_table[i,1] & v < segment_table[i,2])
    #select a random index in those markers
    if (length(b)== 0) next
    s <- sample(b, 1, replace= TRUE)
    if(length(b) == 1)
      s <- b
    submap <- c(submap, s)
  }
    return(submap)
}
# x = une bed matrix
# renvoie une "pseudo" msat matrix sans genotypes mais avec log emiss...
createSubmap <- function(x, intensity = 10, hotspot_version = "hg17")
{
    final_submap <- c()
    for(chr in 1:22)
    {
      #get marker's index for all the chromosom
      submap <- getMarkersChromosom(x, chr, intensity, hotspot_version) 
      final_submap <- c(final_submap, submap)
    }
    #return(final_submap)#a submap has been created
    
    map <- x@snps[ final_submap , c("id","chr")]
    if(all(x@snps$dist == 0)) {
      map$distance <- x@snps$dist[final_submap]
    } else {
      map$distance <- x@snps$pos[final_submap]*1e-6
    }
    
    res <- new("msat.matrix", length(final_submap), nrow(x), 
              x@ped[,c("famid", "id", "father", "mother", "sex", "pheno")]
              ,matrix(nrow = 0, ncol = 0), map, matrix(0, nrow = 0 , ncol =0))
    res@log.emiss <- bed.logEmiss(x, final_submap, 1e-3)
    res@epsilon <- 1e-3
    return(res)
}


