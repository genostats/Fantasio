select.marker.uphill <- function(x, vector, pas, random, unit)
{
  random_snp <- random
  if(unit=="Bases")                                  #recuperer les distances pour pouvoir trouver le marqueur le plus proche a chaque saut
    dist <- x@snps$pos[vector]
  else
    dist <- x@snps$dist[vector]
  
  random_snp <- dist[random_snp]
  seg <- c()
  marker  <- which.min(abs(dist-(random_snp+pas)))   #trouver la position dans vector du marqueur le plus proche en aval
  
  if((dist[marker] - random_snp) < pas)
  {
    if(marker == length(vector))
    {
      return(seg)
    }
    marker <- marker + 1
  }
    
  seg <- c(seg, vector[marker])
  while(marker < length(vector))
  {
    marker <- which.min(abs(dist-(dist[marker]+pas)))
    #distance inferieur
    if((dist[marker] - dist[which(vector == seg[length(seg)])]) < pas)
    {
      if(marker == length(vector))
      {
        return(seg)
      }
      else{
        marker <- marker + 1
      }
    }
    if(vector[marker] %in% seg)
    {
      return(seg)
    }
    seg <- c(seg, vector[marker]) 

  }
  return(seg)
}