# à partir de la position 'random' la fonction choisit à gauche de cette position les 
# autres SNP avec une distance >= pas entre tous les SNPs

select.marker.downstream <- function(x, vector, pas, random, unit)
{
  marker <- random
  if(unit=="Bases")				#get the distance to find which are the nearest markers
    dist <- x@snps$pos[vector[1]:vector[2]]
  else
    dist <- x@snps$dist[vector[1]:vector[2]]
    
  seg <- c((vector[1]:vector[2])[marker])
  while (marker > 1) {
  	marker_dist <- dist[marker]		#get the dist of the random marker
  	d <- abs( dist - marker_dist )
  	ds <- which ( d[1:marker] > 0.5 )
	if (length(ds) == 0) {
		marker <- 0
	} else {
		min_marker <- which ( d == min(d[ds]) )
		if (length(min_marker) > 1) {
			marker <- sample(min_marker,1)
		} else { 
		marker <- which ( d == min(d[ds]) )
		}
	}
  	seg <- c( (vector[1]:vector[2])[marker], seg )
  }
 return(seg)
}
