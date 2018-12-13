getMarkerChromosomByDistance <- function(x, map, pas, unit="cM")
{
  submap <- c()
  for(i in 1:length(map))                    #loop over the chromosome segmetns
  {
    isItAList <- is.list(map[[i]])
    
    if(isItAList)
      looping <- length(map[[i]])
    else
      looping <- 1
    
    for(j in 1:looping)                       #loop over the mini segments in the segment
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