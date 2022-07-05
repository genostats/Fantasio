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
    
    for(j in seq_len(looping))                #loop over the mini segments in the segment
    {
      if(length(map[[i]]) == 0)               #test whether the segment is empty
        next()
      
      if(isItAList)
        mini_segments <- map[[i]][[j]]        #vector index of the mini segment j of the segment i
      else
        mini_segments <- map[[i]]
      
      if(length(mini_segments) == 0) next     #empty mini_segment 
      else if(length(mini_segments) == 1) { 		
      	res <- mini_segments
      } else {
      	 	l <- length(mini_segments[1]:mini_segments[2])  # choose one of them randomly
        	random <- sample(l, 1)				    
   
        	start 	<- selectMarkerDownstream(x, mini_segments, pas, random, unit)	#downstream selection
        	end   	<- selectMarkerUpstream(x, mini_segments, pas, random, unit)       #upstream selection
        	res 	<- c(start, end)
      }
      submap <- c(submap, res)
    }
  }
  submap
}
