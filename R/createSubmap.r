#getSegment <- function(chr, intensity, hotspot_version) 
#{
#  #select the correct hotspot file
#  hotspot <- switch(hotspot_version,
#                    hg17 = { data(hotspot_hg17); hotspot_hg17;},
#                    hg18 = { data(hotspot_hg18); hotspot_hg18;},
#                    hg19 = { data(hotspot_hg19); hotspot_hg19; })
#  
#  #save all  current chromosom's hotspots
#  chr_hotspot <- hotspot[which(hotspot[,1]==chr),]
#  #hotspot_hg17[which(hotspot_hg17==chr),]
#  
#  #save only the hotspot with the right intensity
#  w <- which(chr_hotspot$IntensitycMMb > intensity)
#  
#  #create a matrix with segments delimitations
#  segment_table <- cbind(c(0,chr_hotspot$End[w]),
#                         c(chr_hotspot$Start[w],Inf) )
#  return(segment_table)
#}

getMarkersChromosom <- function(chr, segment_map, mkrmap)
{
  #find the markers index in the interval
  submap <- c()
  segment_table <- segment_map[[chr]]
  mkr <- mkrmap[[chr]]
  #trouver les mkr dans le segment puis en choisir un au hasard
  for ( i in 1:nrow(segment_table))
  {
    b <- which(mkr > segment_table[i,1] & mkr < segment_table[i,2])
    if (length(b)== 0) next
    if(length(b) == 1) s <- b
    s <- sample(b, 1)
    submap <- c(submap, s)
  }
  return(submap)
}
# x = une bed matrix
# renvoie une "pseudo" msat matrix sans genotypes mais avec log emiss...
createSubmap <- function(x, segment_map, mkrmap)
{
  submap <- c()
  for(chr in 1:22)
  {
    v <- getMarkersChromosom(chr, segment_map, mkrmap)
    submap <- c(submap, v)
  }
  #return(submap)#a submap has been created
  
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
  res <- festim(res)
  res <- HBD.prob(res)
  res <- FLOD.prob(res)
  res <- set.HFLOD(res)
  return(res)
}


