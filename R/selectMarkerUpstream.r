select.marker.upstream <- function(x, vector, pas, random, unit)
{
  marker <- random
  if(unit=="Bases")				#get the distance to find which are the nearest markers
    dist <- x@snps$pos[vector[1]:vector[2]]
  else
    dist <- x@snps$dist[vector[1]:vector[2]]
    
  seg <- c()
  while (marker < length(vector[1]:vector[2])) {
  	marker_dist <- dist[marker]		#get the dist of the random marker
  	d <- ( dist - marker_dist )
  	ds <- which ( d > 0.5 )
	if (length(ds) == 0) {
		marker <- Inf
	} else {
		min_marker <- which ( d == min(d[ds]) )
		if (length(min_marker) > 1) {
			marker <- sample(min_marker,1)
		} else { 
		marker <- which ( d == min(d[ds]) )
		}
	seg <- c( seg, (vector[1]:vector[2])[marker] )
	}
  }
 return(seg)
}
