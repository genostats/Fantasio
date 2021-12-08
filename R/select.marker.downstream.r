# cette fonction est supposée 
# partie de la position 'random' et choisir à gauche de cette position des 
# autres SNP avec une distance >= pas entre tous les SNPs

select.marker.downstream <- function(x, vector, pas, random, unit)
{
  random_snp <- random
  if(unit=="Bases")                                     #get the distance to find which are the nearest markers
    dist <- x@snps$pos[vector]
  else
    dist <- x@snps$dist[vector] 
  
  seg <- c(vector[random_snp])                               
  random_snp <- dist[random_snp]                        #get the dist of the random marker
  
  # 'marker' est le SNP le plus proche de la position (random_snp - pas)
  # il faudrait en fait prendre le SNP le plus proche *en dessous* (à gauche) de cette position
  marker <- which.min(abs(dist-(random_snp-pas)))       #find the nearest marker upstream
  
  # ici on teste si marker est à une distance < pas de random_snp
  # si c'est le cas on descend de 1
  # -> ça marche mais hyper alambiqué, corriger le choix de 'marker' ci-dessus
  if((random_snp - dist[marker]) < pas)
  {
    if(marker == 1)
    {
      return(seg)
    }
    marker <- marker - 1
  }
  seg <- c(vector[marker], seg)                          #find the position in the vector
 
  # La dessous on recommence en boucle : faire une seule boucle ! 
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
