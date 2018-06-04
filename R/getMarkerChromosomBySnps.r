getMarkerChromosomBySnps <- function(x, map, pas, unit="cM")
{
  submap <- c()
  for(i in 1:length(map))
  {
    #cat("segment ", i)
    for(j in 1:length(map[[i]]))
    {
      #cat("mini ", j)
      if(length(map[[i]]) == 0)
        next()
      mini_segments <- map[[i]][[j]]
      if(length(mini_segments) == 0) next                                          #mini_segment vide
      if(length(mini_segments) == 1) 
      {
        res <- mini_segments
      }
      else
      {
        random <- which(mini_segments == sample(mini_segments, 1))                  #marqueur aleatoire dans mini_segment
        while( random == 1 & random == length(mini_segments))
        {
          random <- which(mini_segments == sample(mini_segments, 1))
        }
        start <- select.marker.downstream(x, mini_segments, pas, random, unit)      #parcout en amont
        end   <- select.marker.uphill(x, mini_segments, pas, random, unit)          #parcout en aval
        res <- c(start, end)
      
        
      }
      submap <- c(submap, res)
    }
  }
  submap
}