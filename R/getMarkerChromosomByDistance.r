getMarkerChromosomByDistance <- function(x, map, pas, unit="cM")
{
  submap <- c()
  for(i in seq_along(map))                    #loop over the chromosome segmetns
  {
    isItAList <- is.list(map[[i]])
    
    if(isItAList)
      looping <- length(map[[i]])
    else
      looping <- 1
    
    for(j in seq_len(looping))                       #loop over the mini segments in the segment
    {
      if(length(map[[i]]) == 0)               #test whether the segment is empty
        next()
      
      if(isItAList)
        mini_segments <- map[[i]][[j]]        #vector index of the mini segment j of the segment i
      else
        mini_segments <- map[[i]]
      
      if(length(mini_segments) == 0) next     #empty mini_segment 
      if(length(mini_segments) == 1) 
      {
        res <- mini_segments
      }
      else
      {
        random <- which(mini_segments == sample(mini_segments, 1))
        # random <- sample.int(length(mini_segments), 1)                
        ###
        ### la condition ci-dessous n'est jamais vraie
        ### si il faut un | à la place du & (pour éviter les extrémites)
        ### alors si length(min_segments) est 2 : bcle infinie
        ### random <- 1 + sample.int(length(mini_segments) - 2, 1) ferait l'affaire sinon
        while( random == 1 & random == length(mini_segments))
        {
          random <- which(mini_segments == sample(mini_segments, 1))
        }
        start <- select.marker.downstream(x, mini_segments, pas, random, unit)      #upstream selection
        end   <- select.marker.uphill(x, mini_segments, pas, random, unit)          #downhill selection
        res <- c(start, end)
      
        
      }
      submap <- c(submap, res)
    }
  }
  submap
}
