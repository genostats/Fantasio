select.marker.downstream <- function(x, vector, pas, random, unit)
{
  random_snp <- random
  if(unit=="Bases")                                     #get the distance to find which are the nearest markers
    dist <- x@snps$pos[vector]
  else
    dist <- x@snps$dist[vector] 
  
  seg <- c(vector[random_snp])                               
  random_snp <- dist[random_snp]                        #get the dist of the random marker
  
  marker <- which.min(abs(dist-(random_snp-pas)))       #find the nearest marker upstream
  
  if((random_snp - dist[marker]) < pas)
  {
    if(marker == 1)
    {
      return(seg)
    }
    marker <- marker - 1
  }
  seg <- c(vector[marker], seg)                          #find the position in the vector
  
  while(marker > 1)
  {
    marker <- which.min(abs(dist-(dist[marker]-pas)))
    #distance diff
    if((dist[which(vector == seg[1])] - dist[marker]) < pas) #the first element in the vector is the last entered
    {
      if(marker == 1)
      {
        return(seg)
      }
      marker <- marker - 1
    }
    #same distance
    if(vector[marker] %in% seg)
    {
      return(seg)
    }
    seg <- c(vector[marker], seg)
  }
  return(seg)
}