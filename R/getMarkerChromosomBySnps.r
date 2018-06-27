getMarkerChromosomBySnps <- function(x, map, pas, unit="cM")
{
  submap <- c()
  for(i in 1:length(map))#parcourt des segments du chr
  {
    #cat("segment ", i)
    for(j in 1:length(map[[i]]))#parcourt des mini segments du segment du chr
    {
      #cat("mini ", j)
      if(length(map[[i]]) == 0)
        next()
      mini_segments <- map[[i]][[j]]# vecteur d'indice du segment i mini segment j
      if(length(mini_segments) == 0) next                                          #mini_segment vide
      if(length(mini_segments) == 1) 
      {
        res <- mini_segments
      }
      else
      {
        random <- which(mini_segments == sample(mini_segments, 1))                  #indice du marqueur aleatoire dans mini_segment
        while( random == 1 & random == length(mini_segments))
        {
          random <- which(mini_segments == sample(mini_segments, 1))
        }
        start <- select.marker.downstream(x, mini_segments, pas=pas, random, unit)      #parcout en amont
        end   <- select.marker.uphill(x, mini_segments, pas=pas, random, unit)          #parcout en aval
        res <- c(start, end)
      }
      submap <- c(submap, res)
    }
  }
  submap
}