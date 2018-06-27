select.marker.downstream <- function(x, vector, pas, random, unit)
{
  random_snp <- random
  if(unit=="Bases")                                #recuperer les distances pour pouvoir trouver le marqueur le plus proche a chaque saut
    dist <- x@snps$pos[vector]
  else
    dist <- x@snps$dist[vector] 
  
  seg <- c(vector[random_snp])                               
  random_snp <- dist[random_snp]                   #recuperer la distance du marqueur choisit de maniere aleatoire
  
  marker <- which.min(abs(dist-(random_snp-pas)))  #trouver la position dans vector du marqueur le plus proche en amont
  
  if((random_snp - dist[marker]) < pas)
  {
    if(marker == 1)
    {
      return(seg)
    }
    marker <- marker - 1
  }
  seg <- c(vector[marker], seg)                     #touver la position du marqueur dans le vector 
  
  while(marker > 1)
  {
    marker <- which.min(abs(dist-(dist[marker]-pas)))
    #distance diff
    if((dist[which(vector == seg[1])] - dist[marker]) < pas) # the first element in the vector is the last entered
    {
      if(marker == 1)
      {
        return(seg)
      }
      marker <- marker - 1
    }
    #distance non diff
    if(vector[marker] %in% seg)
    {
      return(seg)
    }
    seg <- c(vector[marker], seg)
  }
  return(seg)
}